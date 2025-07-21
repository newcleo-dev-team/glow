"""
Module containing classes that deals with the extraction of the geometric
information from the lattice built in SALOME. The functionalities in this
module serve for preparing all the data for the output TDT file generation.
"""
import math

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple, Self

from glow.support.types import BoundaryType, CellType, GeometryType, \
    LatticeGeometryType, PropertyType, SymmetryType
from glow.geometry_layouts.lattices import Lattice
from glow.support.utility import build_compound_borders, \
    check_shape_expected_types, get_id_from_name, get_id_from_shape, \
    translate_wrt_reference
from glow.interface.geom_interface import ShapeType, \
    extract_sorted_sub_shapes, extract_sub_shapes, get_in_place, \
    get_kind_of_shape, get_min_distance, get_point_coordinates, \
    get_shape_name, get_shape_type, is_point_inside_shape, make_compound, \
    make_face, make_vertex, make_vertex_inside_face, make_vertex_on_curve, \
    set_shape_name


# Sufficiently small value used to determine face-edge connectivity by
# building a point on an edge and shifting to one side of the epsilon.
EPSILON = 1e-05


@dataclass(order=True)
class Face():
    """
    Providing a data representation for a subface of the lattice, which
    represents a calculation zone containing a list of properties.
    This dataclass can be ordered on the basis of the ``no`` attribute, which
    provides a global index for the faces in the lattice.

    Attributes
    ----------
    no          : int
                  Global index of the face
    face        : Any
                  A GEOM object representing a subface of the lattice
    property    : str
                  The value of the property associated to the subface
    inner_point : Tuple[float, float, float]
                  Coordinates of a point within the subface
    edge_vs_id  : Dict[Any, str]
                  Dictionary associating the GEOM edge objects of the lattice
                  subface to the corresponding ID, based on the edge
                  geometrical characteristics
    """
    face        : Any
    property    : str

    sort_index  : int = field(init=False, repr=False)
    no          : int = field(init=False)
    inner_point : Tuple[float, float, float] = field(init=False)
    edge_vs_id  : Dict[Any, str] = field(default_factory=dict, init=False)

    def __post_init__(self):
        """
        Method that is automatically run after the dataclass initialization
        for setting all the attributes that depends on others. In addition,
        the attribute that allows to order instances of this class on the
        basis of the 'no' attribute is set as well.
        """
        # Build a point inside the face
        self.inner_point = get_point_coordinates(
            make_vertex_inside_face(self.face))
        try:
            # Check the type of the received face is correct
            check_shape_expected_types(self.face, [ShapeType.FACE])
            # Extract the face number from the corresponding GEOM name
            # attribute
            self.no = get_id_from_shape(self.face)
        except RuntimeError as e:
            raise RuntimeError(
                f"Error with 'Face' whose inner point is: {self.inner_point}"
            ) from e
        # Define the attribute for sorting instances of this class
        self.sort_index = self.no
        # Extract the edges and associate an ID to each
        edges = extract_sorted_sub_shapes(self.face, ShapeType.EDGE)
        for edge in edges:
            # Build the ID of the edge, based on the geometrical
            # characteristics  --> two edges share the same entry
            edge_id = build_edge_id(edge)
            self.edge_vs_id[edge] = edge_id
            # print(f"Subface no. {self.no}, edge ID: {id}")

    def __str__(self):
        return f"Region {self.no}, name={get_shape_name(self.face)}, " + \
               f"property={self.property}"


class Edge():
    """
    Class that provides a data representation for an edge of a face in
    the lattice.
    It provides an global index allowing to uniquely identify the edge,
    and two `Face` attributes, allowing to identify the face on the left
    and right of the edge.
    The identification of the left/right faces is performed by building
    a point on the edge's normal so that it is sligtly on the edge left.

    Attributes
    ----------
    no : int
        The edge number extracted from the GEOM edge object name.
    edge : Any
        A GEOM object of type EDGE representing an edge
    data : Any
        Characteristic data of the given GEOM edge object.
    kind : str
        Indicating the type of edge with admitted values being `CIRCLE`,
        `ARC_CIRCLE`, `SEGMENT`.
    right : Face | None
        A `Face` object providing the information for the face to the right
        of the edge, or `None` if the edge does not have any face on its
        right.
    left : Face | None
        A `Face` object providing the information for the face to the left
        of the edge, or `None` if the edge does not have any face on its
        left.
    """
    def __init__(self, edge: Any, *faces: List[Face]):
        # Store all the information of the given GEOM edge object
        self.data: List[Any] = get_kind_of_shape(edge)
        try:
            # Check the type of the received edge is correct
            check_shape_expected_types(edge, [ShapeType.EDGE])
            # Get the number of the edge directly from the name attribute of
            # the corresponding GEOM edge object
            self.no : int = get_id_from_shape(edge)
        except RuntimeError as e:
            raise RuntimeError(
                f"Error with 'Edge' whose data is: {self.data}") from e
        # Store the edge object
        self.edge: Any  = edge
        # Get the type of edge as a string
        self.kind: str = str(self.data[0])
        # Initialize to 'None' both the right and the left 'Face' objects
        self.right: Face | None = None
        self.left: Face | None = None
        # Loop through all the given GEOM face objects associated to the
        # current edge. This allows to define the faces on the right and
        # those on the left wrt the edge.
        for f in faces:
            # Add a face connected to the edge
            self.__add_face(f)
        # Check that both right and left faces have been assigned, if the
        # edge is shared by two faces
        if len(faces) > 1 and (not self.right or not self.left):
            raise RuntimeError(
                f"The edge no. {self.no} have 2 faces (no. "
                f"{[face.no for face in faces]}), but they have not be "
                "correctly assigned to the left and right attributes "
                f"(left = {self.left}, right = {self.right}).")

    def __origin(self) -> Tuple[float, float, float]:
        """
        Method that retrieves the 'origin' point (X-Y-Z coordinates) of
        the edge (i.e. the edge starting point) according to its type.
        This information is retrieved from the instance attribute storing
        the characteristic data of the edge.

        Returns
        -------
        A tuple with the X-Y-Z coordinates of the edge starting point
        """
        # Handle the point retrieval differently if the edge is a circle
        if self.kind == "CIRCLE":
            # Get the X-Y-Z coordinates of the circle center
            x1, y1, z1  = self.data[1:4]
            # Get the circle radius
            r = self.data[-1]
            # Return the point (cx + r, cy, cz) with cx, cy, cz the centre
            # of the circle and r its radius.
            return x1 + r, y1, z1
        # Handle the other edge types by returning the X-Y-Z coordinates of
        # the first point of the edge
        return tuple(self.data[-6:-3])

    def __add_face(self, face: Face) -> None:
        """
        Method that allows to define whether the given GEOM face object,
        connected to the edge, is placed to the right or to left of the
        edge object the instance refers to.

        This analysis is based on the creation of a point at a very small
        distance from the edge middle point so that it belongs to the
        left-oriented normal vector for the edge.
        The identification of the face position (left or right) relative
        to the edge depends on the edge nature:
        - for SEGMENT-type edges, the given face is considered on the left
          of the edge if the point belongs to the face; on the contrary,
          we have a right face;
        - for ARC_CIRCLE and CIRCLE-type edges, the criteria is the opposite,
          i.e. the face is considered on the right of the edge if the point
          belongs to the face; on the contrary, we have a left face.

        Parameters
        ----------
        face    : Face
                  'Face' object whose position (right or left) relative
                  to the edge has to be determined
        epsilon : float
                  Margin small enough to place a point wrt the edge
        """
        # Build the point on the left of the edge and use it to identify the
        # face position relative to the edge
        if not is_point_inside_shape(self.__build_point_on_edge_normal(),
                                     face.face):
            if self.kind == "SEGMENT":
                self.right = face
            else:
                self.left = face
        else:
            if self.kind == "SEGMENT":
                self.left = face
            else:
                self.right = face

    def __build_point_on_edge_normal(self, epsilon: float = EPSILON) -> Any:
        """
        Method that builds a vertex object positioned at an infinitesimal
        distance from the edge object this instance refers to.
        Depending on the edge type, this point is built according to the
        following rules:
        - CIRCLE: the point is positioned slightly on the left wrt the
          X-coordinate of the point laying on the right-most position of
          the circle (identified by the center X-coordinate + the radius).
        - ARC_CIRCLE: given a point positioned at the middle of the arc,
          another point is built on the vector connecting the arc center and
          the first point. This second point has both its X-Y coordinates
          slightly scaled down by a reduction factor so that its distance
          from the arc center is less than the radius.
        - SEGMENT: the normalized left-oriented vector normal to the edge is
          calculated. Its X-Y components are used to determine a point
          slightly to the left of a point positioned at the middle of the
          segment.

        Parameters
        ----------
        epsilon : float
            Indicating a value small enough to determine a point to the left
            of the edge

        Returns
        -------
        A vertex object representing a point slightly to the left of the
        middle point of the edge.
        """
        if self.kind == "CIRCLE":
            # Extract the X-Y-Z coordinates of the circle center
            (xc, yc, zc) = self.data[1:4]
            # Extract the circle radius
            radius = self.data[7]
            # Build a point (as a GEOM object) positioned at an infinitesimal
            # distance from the starting point of the circle
            return make_vertex((xc+radius-epsilon, yc, zc))
        if self.kind == "ARC_CIRCLE":
            # Extract the X-Y-Z coordinates of the circle center the arc
            # belongs to
            (xc, yc, zc) = self.data[1:4]
            # Extract the radius of the arc
            radius = self.data[7]
            # Build a GEOM point on the arc, positioned at its middle
            pm = make_vertex_on_curve(self.edge, 0.5)
            # Get the X-Y-Z coordinates of the point
            coord = get_point_coordinates(pm)
            # Build a 2D vector from the arc center to the arc middle point,
            # which is slightly scaled down by a reduction factor so that
            # its length is less than the radius
            vm = ((coord[0] - xc)*(1.0 - epsilon/radius),
                  (coord[1] - yc)*(1.0 - epsilon/radius))
            # Build a GEOM point along the direction of the just built vector
            return make_vertex((xc + vm[0], yc + vm[1], zc))
        if self.kind == "SEGMENT":
            # Extract the X-Y-Z coordinates of the extremes of the segment
            (x1, y1, z1, x2, y2, _) = self.data[1:7]
            # Build a GEOM point on the segment, positioned at its middle
            pm = make_vertex_on_curve(self.edge, 0.5)
            # Get the X-Y-Z coordinates of the point
            coord = get_point_coordinates(pm)
            # Define the left-oriented normal vector of the segment
            n = ((y1 - y2), (x2 - x1))
            # Get the length of the normal vector
            l = math.sqrt(n[0]*n[0] + n[1]*n[1])
            # Get the coordinates of a point positioned at a distance epsilon
            # from the segment along its left-oriented normalized normal
            # vector
            (xm, ym) = ((coord[0] + epsilon*n[0]/l),
                        (coord[1] + epsilon*n[1]/l))
            # Build a GEOM point positioned at the calculated coordinates
            return make_vertex((xm, ym, z1))
        # If here, raise an exception as the edge has not a valid type
        raise ValueError(f"The kind of shape {self.kind} is not valid!")

    def __str__(self):
        lname = rname = "None"
        ename = get_shape_name(self.edge)
        if self.left:
            rname = str(self.left)
        if self.right:
            lname = str(self.right)
        return (f"Edge {ename} \n\torigin = {self.__origin()},"
                f"\n\tleft = {lname}, \n\tright = {rname}")

    # ------------------------------------------
    # Methods for comparing two 'Edge' instances
    # ------------------------------------------
    def __lt__(self, other: Self):
        return self.no < other.no

    def __le__(self, other: Self):
        return self.no <= other.no

    def __eq__(self, other: Self):
        return self.no == other.no

    def __ne__(self, other: Self):
        return self.no != other.no

    def __gt__(self, other: Self):
        return self.no > other.no

    def __ge__(self, other: Self):
        return self.no >= other.no


class Boundary:
    """
    Class that provides the data structure for a GEOM edge object being a
    border of the lattice.

    Parameters
    ----------
    border : Any
        A GEOM edge object representing a border of the lattice
    type_geo : LatticeGeometryType
        Providing the lattice type of geometry
    lattice_o : Any
        A vertex object representing the lattice center point
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the lattice the border refers to

    Attributes
    ----------
    border      : Any
                  A GEOM edge object representing a border of the lattice
    type        : BoundaryType
                  Providing the type of BCs
    angle       : float
                  Angle (in degrees) of the border defined from its origin to
                  its second point, so that it is always positive
    edge_indxs  : List[int]
                  List of the indices of the edges the border is made of
    tx          : float
                  The X-component of the border axis
    ty          : float
                  The Y-component of the border axis
    """
    def __init__(self,
                 border: Any,
                 type_geo: LatticeGeometryType,
                 lattice_o: Any,
                 dimensions: Tuple[float, float]) -> None:
        # Check the received border is an edge
        try:
            check_shape_expected_types(border, [ShapeType.EDGE])
        except RuntimeError as e:
            raise RuntimeError("Error while initializing the 'Border' "
                               "instance.") from e
        # Initialize instance attributes
        self.type       : BoundaryType
        self.border     : Any = border
        self.angle      : float = 0.0
        self.edge_indxs : List[int] = []
        self.tx         : float = 0.0
        self.ty         : float = 0.0
        # Set the border characteristics in terms of type and border axis
        self.__build_border_characteristics(*dimensions, lattice_o, type_geo)

    def get_bc_type_number(self) -> int:
        """
        Method that returns the index associated to the BC type in the
        corresponding instance attribute dictionary.

        Returns
        -------
        An integer describing the type for the BC this class instance
        refers to.
        """
        return self.type.value

    def __build_border_characteristics(
            self,
            lx: float,
            ly: float,
            lattice_o: Any,
            type_geo: LatticeGeometryType) -> None:
        """
        Method that defines the characteristics of a lattice border that
        represents a boundary for the lattice itself.
        These characteristics are defined in terms of the associated BC type
        (value of the 'BoundaryType' enumeration), the border axes and its
        angle; these depend on the specific lattice type of geometry (as
        value of the 'LatticeGeometryType' enumeration).

        Parameters
        ----------
        lx  : float
              The X-characteristic dimension of the lattice
        ly  : float
              The Y-characteristic dimension of the lattice
        lattice_o : Any
            Vertex object representing the lattice center point
        """
        # Get the X-Y coordinates of the two end points of the edge
        x1, y1, _, x2, y2 = get_kind_of_shape(self.border)[1:6]
        # Calculate the X-Y distances between the start-end points
        dx = x2 - x1
        dy = y2 - y1
        # Calculate the angle of the border in degrees
        self.angle = math.degrees(math.atan2(dy, dx))

        # The border origin must be defined so that the angle between
        # the start-end points is positive. If not, the edge start point
        # is inverted.
        if self.angle < -EPSILON or math.isclose(self.angle, 180.0):
            self.angle = math.degrees(math.atan2(-dy, -dx))
            x1 = x2
            y1 = y2

        # Initialize the border axes to the edge start point
        self.tx = x1
        self.ty = y1

        # Assign the BC type depending on the lattice type of geometry. The
        # geometries identified by 'RECTANGLE_TRAN' and 'HEXAGON_TRAN' need
        # to re-evaluate the border axes.
        match type_geo:
            case (LatticeGeometryType.SYMMETRIES_TWO |
                  LatticeGeometryType.RECTANGLE_SYM |
                  LatticeGeometryType.RECTANGLE_EIGHT |
                  LatticeGeometryType.SA60 |
                  LatticeGeometryType.S30):
                # Identifying a border characterized by the 'AXIAL_SYMMETRY'
                # type of BC, which corresponds to the 'REFL' case in DRAGON5
                self.type = BoundaryType.AXIAL_SYMMETRY
            case (LatticeGeometryType.RA60 |
                  LatticeGeometryType.R120 |
                  LatticeGeometryType.ROTATION):
                # Identifying a border characterized by any of the 'ROTATION'
                # or 'TRANSLATION' types of BC, which correspond to the 'ROTA'
                # or 'TRAN' cases respectively in DRAGON5. The position of the
                # border wrt the lattice center guides the choice.
                if get_min_distance(lattice_o, self.border) > 1e-7:
                    self.type = BoundaryType.TRANSLATION
                else:
                    self.type = BoundaryType.ROTATION
            case LatticeGeometryType.RECTANGLE_TRAN:
                # The BC information for the case of a cartesian geometry
                # with TRAN BCs follows the axes definition below:
                #              M=3 (0,-ly)
                #             ************
                #  M=2 (lx,0) *          * M=4 (-lx,0)
                #             ************
                #              M=1 (0,ly)
                if math.isclose(math.sin(math.radians(self.angle)), 0.0):
                    # The sign of 'dx' discriminates between the M=1 (dx > 0)
                    # and M=3 (dx < 0)
                    self.tx = 0.0
                    self.ty = (dx/abs(dx)) * ly
                elif math.isclose(math.sin(math.radians(self.angle)), 1.0):
                    # The sign of 'dy' discriminates between the M=4 (dy > 0)
                    # and M=2 (dy < 0)
                    self.tx = -(dy/abs(dy)) * lx
                    self.ty = 0.0
                else:
                    raise RuntimeError(
                        f"The border has an angle of {self.angle}° which is "
                        "not one of the admitted values (0°, 90°) for a "
                        "cartesian geometry with TRAN as BC.")
                # Assign the BC type
                self.type = BoundaryType.TRANSLATION
            case LatticeGeometryType.HEXAGON_TRAN:
                # The BC information for the case of an hexagonal geometry
                # with translation on its sides follows the axes definition
                # below:
                #                     M=4 (0,-2ly)
                #                    *****
                #   M=3 (3/2lx,-ly) *     *  M=5 (-3/2lx,-ly)
                #                  *       *
                #   M=2 (3/2lx, ly) *     *  M=6 (-3/2lx, ly)
                #                    *****
                #                     M=1 (0, 2ly)
                if math.isclose(math.sin(math.radians(self.angle)), 1.0):
                    raise RuntimeError(
                        "The border refers to a Y-oriented hexagon which "
                        "is not admitted for tracking.")
                if abs(dy) < 1e-7:
                    # The sign of 'dx' discriminates between the M=1 (dx > 0)
                    # and M=4 (dx < 0)
                    self.tx = 0.0
                    self.ty = (dx/abs(dx)) * 2*ly
                elif dy > 1e-7:
                    # The sign of 'dx' discriminates between the M=6 (dx > 0)
                    # and M=5 (dx < 0)
                    self.tx = -3/2 * lx
                    self.ty = (dx/abs(dx)) * ly
                else:
                    # The sign of 'dx' discriminates between the M=2 (dx > 0)
                    # and M=3 (dx < 0)
                    self.tx = 3/2 * lx
                    self.ty = (dx/abs(dx)) * ly
                # Assign the BC type
                self.type = BoundaryType.TRANSLATION
            case _:
                raise RuntimeError(
                    f"The {type_geo} lattice geometry type is not "
                    "currently handled.")

    def find_edges_on_border(self,
                             boundaries: Any,
                             id_vs_edge: Dict[str, Any]) -> None:
        """
        Method that finds and associates all the GEOM edges related to the
        lattice border this instance refers to.
        These edges are the ones that are connected to a single face only,
        as these are the edges on the borders of the lattice.

        Parameters
        ----------
        boundaries : Any
            A GEOM compound object made from the list of GEOM edge objects
            belonging to the lattice boundaries
        id_vs_edge : Dict[Any, str]
            A dictionary of all the lattice edge IDs VS the corresponding
            GEOM edge objects
        """
        # Loop through all the sub-edges that are part of the lattice border
        # this class instance refers
        for shape in extract_sorted_sub_shapes(
            get_in_place(boundaries, self.border), ShapeType.EDGE):
            # Extract the shape type from all the information about the current
            # one
            shape_type = get_kind_of_shape(shape)[0]
            # Check the retrieved shape is of type 'SEGMENT'
            # TODO add support also for ARC_CIRCLE
            if str(shape_type) == 'SEGMENT':
                # Retrieve the GEOM edge object corresponding to its ID
                edge_ref = id_vs_edge[build_edge_id(shape)]
                # Append the edge's global index number to the list
                self.edge_indxs.append(
                    get_id_from_name(get_shape_name(edge_ref)))
            else:
                # Raise an exception if the found sub-shape type is not
                # 'SEGMENT'
                raise RuntimeError(
                    "Only edges of type 'SEGMENT' can be contained in a "
                    f"lattice border! (found {shape_type})")


class LatticeDataExtractor():
    """
    Class that extracts the needed geometrical data from the given lattice
    to be used further on for generating the output TDT file.
    It determines the association of faces with properties, that of edges
    with connected faces and the edges that represents the lattice borders.

    Parameters
    ----------
    lattice : Lattice
        An instance of the 'Lattice' class storing the lattice information
    geom_type : GeometryType
        The type of geometry of the lattice cells used to extract the regions

    Attributes
    ----------
    lattice             : Lattice
                          The 'Lattice' object storing the lattice information
    borders             : Any
                          The GEOM compound object grouping the edges
                          representing the lattice borders
    boundaries          : List[Boundary]
                          The list of 'Boundary' objects storing the lattice
                          borders information
    subfaces            : List[Any]
                          The list of 'Face' objects storing information about
                          each region of the lattice
    edges               : List[Edge]
                          The list of 'Edge' objects storing information about
                          each edge of the lattice
    id_vs_edge          : Dict[str, Any]
                          A dictionary of the lattice edge ID VS the GEOM edge
                          object reference
    lattice_edges       : List[Any]
                          A list of the GEOM edge objects for the lattice
    """
    # Identifying the symmetry types for each type of lattice cells for which
    # the lattice translation should be evaluated
    CASES_FOR_TRANSLATION = {
        CellType.RECT: [SymmetryType.FULL, SymmetryType.HALF],
        CellType.HEX: [SymmetryType.THIRD, SymmetryType.SIXTH]
    }

    def __init__(self, lattice: Lattice, geom_type: GeometryType) -> None:
        # Initialize the instance attributes
        self.lattice: Lattice = deepcopy(lattice)
        self.borders: List[Any] = []
        self.lattice_edges: List[Any] = []
        self.boundaries: List[Boundary] = []
        self.subfaces: List[Face] = []
        self.edges: List[Edge] = []
        # Extract the information to be stored about borders and edges of the
        # lattice according to the applied symmetry and geometry
        self.__preprocess(geom_type)
        # Associate each edge with an index and build a dictionary
        self.id_vs_edge: Dict[str, Any] = classify_lattice_edges(
            self.lattice_edges)

    def build_boundaries(self) -> None:
        """
        Method that constructs a list of 'Boundary' objects representing the
        lattice boundary edges, i.e. those connected to a single face.
        All the GEOM edges being part of each boundary are associated to the
        same 'Boundary' object.
        """
        # No boundaries to extract if an 'ISOTROPIC' type of geometry, meaning
        # 'VOID' or 'ALBE 1.0' BCs in DRAGON5.
        if self.lattice.type_geo == LatticeGeometryType.ISOTROPIC:
            return
        # Initialize the list of 'Boundary' objects
        boundary_edges = []

        # Loop through all the stored 'Edge' objects
        print("LEN EDGES", len(self.edges))
        for edge in self.edges:
            # Check if the current edge has no faces on the left or on the
            # right, i.e. it represents a boundary edge
            if not edge.left or not edge.right:
                # Append the corresponding GEOM edge object to the list
                boundary_edges.append(edge.edge)

        # Build a compound from all the edges placed on the lattice borders
        boundary_edgs_cmpd = make_compound(boundary_edges)
        # Loop through all the edge objects representing the lattice borders
        print("LEN BORDERS:", len(self.borders))
        for border in self.borders:
            # Build an object of the 'Boundary' class
            boundary = Boundary(border=border,
                                type_geo=self.lattice.type_geo,
                                lattice_o=self.lattice.lattice_center,
                                dimensions=(self.lattice.lx, self.lattice.ly))
            # Store all the indices of the edges belonging to the border
            boundary.find_edges_on_border(boundary_edgs_cmpd, self.id_vs_edge)
            # Append the built 'Boundary' object to the corresponding list
            self.boundaries.append(boundary)

    def build_edges(
            self, edge_names_vs_faces: Dict[str, List[Any | Face]]) -> None:
        """
        Method that builds a list of 'Edge' objects from the given dictionary.
        It associates for each edge name a list containing the corresponding
        GEOM edge and the 'Face' objects; the latter represent the faces
        sharing the same edge.

        Parameters
        ----------
        edge_names_vs_faces : Dict[str, List[Any | Face]]
            Dictionary of edge names VS the list with the corresponding
            GEOM edge and the connected 'Face' objects
        """
        # Loop through all the lists of objects associated to each edge
        for shapes in edge_names_vs_faces.values():
            # Instantiate an object of the 'Edge' class and append to the
            # corresponding list
            self.edges.append(Edge(*shapes))

    def build_edges_and_faces_association(self) -> Dict[str, List[Face]]:
        """
        Method that associates the faces sharing the same edge with the edge
        name. These names contain a global index to identify the edge in the
        lattice.

        Returns
        -------
        A dictionary whose entries associate a list of adjacent faces (as
        'Face' objects) to the name of the corresponding shared edge.
        """
        # Initialize the dictionary storing the edges names VS the list of
        # connected faces
        edges_name_vs_faces : Dict[str, List[Any | Face]] = {}
        # Loop through all the lattice subfaces ('Face' objects)
        print("LEN SUBFACES", len(self.subfaces))
        for subface in self.subfaces:
            # Log the 'Face' characteristics
            print(subface)
            # Loop through all the edges of the current subface
            for subface_edge, edge_id in subface.edge_vs_id.items():
                # Extract the corresponding GEOM edge object(s)
                unique_edges = self.__get_unique_edges(subface_edge, edge_id)
                # Update the dictionary of edge names VS connected faces
                self.__update_edge_face_association(
                    edges_name_vs_faces, subface, unique_edges)
        # Return the dictionary of edge names VS the list of connected faces
        return edges_name_vs_faces

    def build_faces(self, property_type: PropertyType) -> None:
        """
        Method that builds a list of 'Face' objects from the GEOM face
        objects, extracted from the lattice regions, and the property
        values. Each region in the lattice corresponds to a region in a
        cell, according to the type of geometry (either technological or
        sectorized).
        Each region must be associated with a value for the given property
        type. If not the case, an exception is raised.

        Parameters
        ----------
        property_type : PropertyType = PropertyType.MATERIAL
            The type of property associated to the lattice regions

        Raises
        ------
        RuntimeError:
            - If no properties are associated to a lattice region.
            - If no value for the given property type is associated to a
              lattice region.
        """
        # Index identifying the 'Face' object
        subface_indx = 0
        # Loop through all the lattice regions and build the corresponding
        # data structure storing the face object and the value of the given
        # type of property
        print("LEN REGIONS:", len(self.lattice.regions))
        for region in self.lattice.regions:
            # Update the subface index
            subface_indx += 1
            # Set the subface name by providing its index
            set_shape_name(region.face, f"FACE_{subface_indx}")
            # Get the value of the given property type associated to the
            # region, if any
            if not region.properties:
                raise RuntimeError(
                    "The lattice analysis failed: no properties have been "
                    f"assigned for region '{region.name}'.")
            try:
                value = region.properties[property_type]
            except KeyError:
                raise RuntimeError(
                    f"The lattice analysis failed: no {property_type.name} "
                    "property type has been defined for region "
                    f"'{region.name}'")
            if not value:
                raise RuntimeError(
                    "The lattice analysis failed: no value for the property "
                    f"type {property_type.name} has been defined for region "
                    f"'{region.name}'")
            # Build a 'Face' object and append to the corresponding list
            self.subfaces.append(Face(region.face, value))

    def print_log_analysis(self, edge_name_vs_faces: Dict[str, List[Any]]):
        """
        Method that prints on the stdout the log of the data extraction from
        the lattice.

        Parameters
        ----------
        edge_name_vs_faces : Dict[str, List[Any]]
            A dictionary of the name of the lattice edges VS the GEOM faces
            connected to it
        """
        # Displays the number of faces, the number of edges and the edges
        # associated with one face and two faces.
        n0 = 0
        n1 = 0
        n2 = 0
        for f in edge_name_vs_faces.values():
            if len(f) == 2:
                n1 += 1
            elif len(f) == 3:
                n2 += 1
            else:
                n0 += 1
        print("\t# subparts (faces) : ", len(self.subfaces))
        print("\t# edges found : ", len(self.edges))
        print("\t# edges w/one face :", n1)
        print("\t# edges w/two faces :", n2)
        print("\t# edges w errors :", n0)

    def __apply_lattice_elements_translation(
            self,
            lattice_cmpd: Any,
            new_center: Tuple[float, float, float]) -> Any:
        """
        Method that translates the lattice regions and the given lattice
        compound so that they are positioned according to the provided new
        center.

        If the current lattice center differs from the provided one, all
        the regions of the lattice are translated so to keep their relative
        distance from the translated lattice center. The same goes for the
        provided lattice compound and the vertex object representing the
        center of the stored `Lattice` instance.

        Parameters
        ----------
        lattice_cmpd : Any
            The lattice compound object to be translated
        new_center : Tuple[float, float, float]
            The coordinates for the new lattice center

        Returns
        -------
        The translated lattice compound, or the same compound if the center
        has not changed.
        """
        # Procede only if the lattice center has changed
        if all(math.isclose(c, nc) for c, nc in zip(
            get_point_coordinates(self.lattice.lattice_center), new_center)):
            return lattice_cmpd
        pre_center = self.lattice.lattice_center
        # Rebuild the lattice center in the new position
        self.lattice.lattice_center = make_vertex(new_center)
        # Translate the regions of the lattice
        for region in self.lattice.regions:
            region.face = translate_wrt_reference(
                region.face, pre_center, new_center)
        # Translate the given lattice compound and return it
        return translate_wrt_reference(
            lattice_cmpd, pre_center, new_center)

    def __evaluate_lattice_center(self) -> Tuple[float, float, float]:
        """
        Method that evaluates and returns the coordinates of the lattice
        center so that the lower-left corner of the lattice is positioned
        in the XYZ space origin.

        The method constructs a face from the lattice's borders, extracts
        its vertices, and determines the lower-left corner vertex based on
        the minimum X, Y, and Z coordinates.
        If the lower-left corner does not coincide with the XYZ space origin,
        it computes and returns the coordinates the center should have so
        that the lower-left corner is placed in the origin.
        Otherwise, it returns the current lattice center coordinates.

        Returns
        -------
        A tuple providing the coordinates of the lattice center so that the
        lower-left corner coincides with the XYZ space origin.
        """
        # Get the lattice face vertices
        lattice_vertices = extract_sub_shapes(make_face(self.borders),
                                              ShapeType.VERTEX)
        # Get the lower-left corner vertex as the one having the minimum value
        # for the X-Y-Z coordinates
        coords = [get_point_coordinates(p) for p in lattice_vertices]
        lower_left = min(coords, key=lambda x: (x[0], x[1], x[2]))
        # Check if the lower-left corner coincides with the XYZ origin; if
        # not, evaluate the new lattice center to fulfill the condition
        if any(not math.isclose(c, 0.0) for c in lower_left):
            print("!!! Lower-left corner not in O")
            # Return the coordinates of the new center
            return (0.0 - lower_left[0]), (0.0 - lower_left[1]), 0.0
        # Return the current lattice center
        return get_point_coordinates(self.lattice.lattice_center)

    def __get_lattice_compound(self) -> Any:
        """
        Method that returns the lattice compound that corresponds to the
        currently applied symmetry type. It identifies either the full
        lattice of a part of it, if any symmetry is applied.

        Returns
        -------
          The lattice compound that corresponds to the currently applied
          symmetry type.
        """
        lattice_cmpd = self.lattice.lattice_cmpd
        if not self.lattice.symmetry_type == SymmetryType.FULL:
            lattice_cmpd = self.lattice.lattice_symm
        return lattice_cmpd

    def __get_unique_edges(
            self, subface_edge: Any, edge_id: str) -> List[Any]:
        """
        Method that retrieves the GEOM edge object associated to the given
        ID (second argument) from the attribute dictionary of IDs VS GEOM
        edges.
        In case of an exception, due to a missing entry in the `id_vs_edge`
        dictionary, a further analysis is performed. This could happen for
        edges shared by two adjacent faces: the same edge could have a
        different orientation in the two faces, i.e. starting and ending
        points are inverted, or an edge having a single face from one side
        can can be associated to different faces from the other.
        In both cases, it is important to retrieve all the corresponding sub
        edges by exploiting the GEOM function `GetInPlace`; this is used to
        extract the sub-shape(s) of the lattice unique edges, which are
        coincident with, or could be a part of, the GEOM edge provided as
        first argument.
        If any edge is retrieved, the corresponding ID is built and used to
        get the corresponding GEOM edges stored in the `id_vs_edge` attribute.

        Parameters
        ----------
        subface_edge : Any
            The GEOM edge object to retrieve the corresponding sub-edges from
        edge_id : str
            The ID of the edge to look for in the dictionary of IDs VS edges

        Raises
        ------
        RuntimeError
            If no corresponding edge is found in the lattice

        Returns
        -------
        A list of GEOM edge objects that is either directly associated to the
        given ID, or representing the sub-edges the edge can be subdivided
        into.
        """
        try:
            return [self.id_vs_edge[edge_id]]
        except KeyError as exc:
            # Get the sub-edges the given edge can be subdivided into: these
            # are associated to a different face of the lattice
            edge = get_in_place(make_compound(self.lattice_edges),
                                subface_edge)
            # Only edge and compound of edges are treated
            edge_type = get_shape_type(edge)
            error_message = "No corresponding edge in the lattice could " +\
                f"be retrieved for the subface edge whose data is {edge_id}"
            if not edge or edge_type not in [ShapeType.COMPOUND,
                                             ShapeType.EDGE]:
                raise RuntimeError(error_message) from exc
            # EDGE-type case
            if not edge_type == ShapeType.COMPOUND:
                return [self.id_vs_edge[build_edge_id(edge)]]
            # COMPOUND-type case
            edges_in_place = extract_sub_shapes(edge, ShapeType.EDGE)
            if edges_in_place:
                return [
                    self.id_vs_edge[build_edge_id(e)] for e in edges_in_place]
            else:
                # Raise an exception if the compound does not have edges
                raise RuntimeError(error_message) from exc

    def __preprocess(self, geom_type: GeometryType) -> None:
        """
        Method that initializes the information to be stored from the lattice
        according to the applied symmetry and the given geometry.
        It re-builds the regions of the lattice according to the given
        `GeometryType`, if it is needed or the indicated geometry type is
        different from the one used to display them in the SALOME viewer.
        This allows that analysis on the lattice is always performed on the
        up to date regions.
        In case the lattice falls in one of the cases identified by the
        `CASES_FOR_TRANSLATION` attribute, its compound is translated together
        with its regions and the center, so that the lower-left corner
        coincides with the origin of the XYZ space.

        Parameters
        ----------
        geom_type : GeometryType
            The type of geometry of the lattice cells used to extract the
            regions
        """
        # Check if the lattice geometry layout needs to be rebuilt
        if (self.lattice.is_update_needed or
            geom_type != self.lattice.displayed_geom):
            self.lattice.build_regions(geom_type)
        # Get the GEOM compound identifying either the full lattice of a part
        # of it, if a symmetry is applied
        lattice_cmpd = self.__get_lattice_compound()
        # Extract the lattice borders
        self.borders = build_compound_borders(lattice_cmpd)
        # Handle the lattice translation so that the lower-left corner is in
        # the XYZ space origin; this is valid for specific symmetries and
        # cells geometries
        if self.lattice.symmetry_type in self.CASES_FOR_TRANSLATION[
            self.lattice.cells_type]:
            # Evaluate the new center of the lattice, if it has not been
            # translated yet, and apply the translation to the regions
            # and the lattice compound
            lattice_cmpd = self.__apply_lattice_elements_translation(
                lattice_cmpd, self.__evaluate_lattice_center())
            # Re-evaluate the lattice borders
            self.borders = build_compound_borders(lattice_cmpd)
        # Extract the lattice edges from the lattice compound to analyse
        self.lattice_edges = extract_sub_shapes(lattice_cmpd, ShapeType.EDGE)

    def __update_edge_face_association(
            self,
            edges_name_vs_faces: Dict[str, List[Any | Face]],
            subface: Face,
            edges: List[Any]) -> None:
        """
        Method that updates the given dictionary of edge names VS connected
        faces with the provided 'Face' object and the list of edges.
        For each edge, its ID is checked for its presence among the keys of
        the given dictionary: if so, the corresponding list is updated with
        the 'Face' object, otherwise a new entry is created. The entry has
        as key the edge name, as value a list with the GEOM edge object as
        first element, followed by the given 'Face' object.

        Parameters
        ----------
        edges_name_vs_faces : Dict[str, List[Any | Face]]
            A dictionary of edge names VS connected faces
        subface : Face
            A 'Face' object to associate the given edges to
        edges : List[Any]
            A list of GEOM edge objects each associated to the given face
        """
        for edge in edges:
            # Get the edge name
            edge_name = get_shape_name(edge)
            # Check if the edge ID is already stored, if so, append
            # the current face, otherwise add a new entry
            if edge_name in edges_name_vs_faces:
                edges_name_vs_faces[edge_name].append(subface)
            else:
                edges_name_vs_faces[edge_name] = [edge, subface]


def analyse_lattice(lattice: Lattice,
                    geom_type: GeometryType,
                    property_type: PropertyType) -> LatticeDataExtractor:
    """
    Function that performs the lattice analysis in order to extract the
    needed information about the regions, and the associated properties,
    the edges, the association between edges and faces connected to each
    of them, and the edges representing the lattice boundaries.

    Parameters
    ----------
    lattice : Lattice
        The instance of the 'Lattice' class storing the geometrical data
        about the lattice to analyse
    geom_type : GeometryType = GeometryType.TECHNOLOGICAL
        The type of geometry of the lattice cells to use in the analysis
    property_type : PropertyType = PropertyType.MATERIAL
        The type of property associated to the lattice regions to use in
        the analysis

    Returns
    -------
    LatticeDataExtractor
        Object collecting all the information about the geometry and the
        properties extracted from the lattice.
    """
    # Instantiate the class for extracting the geometric data from the lattice
    # according to the given type of geometry
    data_extractor = LatticeDataExtractor(lattice, geom_type)
    # Call its method for performing the analysis
    data_extractor.build_faces(property_type)
    edge_name_vs_faces = data_extractor.build_edges_and_faces_association()
    data_extractor.build_edges(edge_name_vs_faces)
    data_extractor.build_boundaries()
    data_extractor.print_log_analysis(edge_name_vs_faces)

    # Return the instance
    return data_extractor

def build_edge_id(edge: Any) -> str:
    """
    Function that builds an unique ID for a GEOM edge object in the geometry
    to process.
    This is performed by retrieving characteristic information abount the edge
    by means of the GEOM ``KindOfShape()`` function. This provides a list
    containing the type of shape and a series of parameters that describe the
    shape itself.
    According to the shape type we could have:

    - CIRCLE xc yc zc dx dy dz R
      (X-Y-Z center coordinates, X-Y-Z normal vector elements, circle radius)
    - ARC_CIRCLE xc yc zc dx dy dz R x1 y1 z1 x2 y2 z2
      (X-Y-Z center coordinates, X-Y-Z normal vector elements, arc radius,
      X-Y-Z coordinates of arc starting and ending points)
    - SEGMENT x1 y1 z1 x2 y2 z2
      (X-Y-Z coordinates of segment starting and ending points)

    If any shape type other than 'CIRCLE', 'ARC_CIRCLE' and 'SEGMENT' is
    provided, an exception is raised.

    Parameters
    ----------
    edge  : Any
            A GEOM edge object

    Returns
    -------
    A string representing the built unique ID for the edge.
    """
    # Get the information of the given GEOM shape object
    data = get_kind_of_shape(edge)
    # Check if the shape is one of the admitted edges
    if str(data[0]) not in ['CIRCLE', 'ARC_CIRCLE', 'SEGMENT']:
        raise ValueError(f"The shape, whose information is '{data}', is not "
                         "one of the admitted edges 'CIRCLE', 'ARC_CIRCLE', "
                         "'SEGMENT'")
    # Loop through all the other information about the edge and build
    # the ID by appending all the info
    return "EDGE_" + str(data[0]) + "_" + "_".join(f"{info:.6g}"
                                                   for info in data[1:])

def classify_lattice_edges(edges: List[Any]) -> Dict[str, Any]:
    """
    Function that classifies the given lattice edges (as GEOM objects)
    by building a dictionary with keys being a unique ID and values the
    corresponding GEOM edge object.
    Handled lattice cases are:
    - complete lattice;
    - boxed lattice with applied symmetry.

    Parameters
    ----------
    edges : List[Any]
            A list of the GEOM edge objects to be classified

    Returns
    -------
    A dictionary storing the edges IDs and the GEOM edges themselves.
    """
    # # Initialize the list of GEOM edges
    # edges = []
    # # Handle the case of a lattice with applied symmetry
    # if not lattice.symmetry_type == SymmetryType.FULL:
    #     edges = geompy.SubShapeAllSortedCentres(lattice.lattice_symm,
    #                                             EDGE_TYPE)
    # else:
    #     # Full lattice case
    #     edges = geompy.SubShapeAllSortedCentres(lattice.lattice_edges,
    #                                             EDGE_TYPE)
    ids_edges = {}
    for indx, edge in enumerate(edges):
        # Set the name of the GEOM edge object
        set_shape_name(edge, f"EDGE_{indx + 1}")
        # Build a unique ID for the edge
        edge_id = build_edge_id(edge)
        # Add an entry to the dictionary storing edges ID VS the
        # corresponding GEOM edge objects
        ids_edges[edge_id] = edge
    # Return the built dictionary
    return ids_edges
