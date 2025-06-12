"""
Module containing classes that deals with the extraction of the geometric
information from the lattice built in SALOME. The functionalities in this
module serve for preparing all the data for the output TDT file generation.
"""

import math

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple, Union, Self

from glow.generator.support import BoundaryType, GeometryType, \
    LatticeGeometryType, PropertyType, SymmetryType
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import ShapeType, \
    extract_sorted_sub_shapes, fuse_edges_in_wire, get_closed_free_boundary, \
    get_in_place, get_kind_of_shape, get_min_distance, get_point_coordinates, \
    get_shape_name, make_common, make_compound, make_face, make_fuse, \
    make_vertex, make_vertex_inside_face, make_vertex_on_curve, set_shape_name


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
        # Extract the face number from the corresponding GEOM name attribute
        self.no = int(get_shape_name(self.face).split('_')[1])
        # Build a point within the face
        self.inner_point = get_point_coordinates(
            make_vertex_inside_face(self.face))
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

    This object can be used to determine the position of one or two faces
    relative to an edge.
    The position is determined either by performing a scalar product of
    two vectors, or by observing the position of a point relative to a
    circle.

    Attributes
    ----------
    no      : int
              The edge number extracted from the GEOM edge object name
              attribute
    edge    : Any
              A GEOM object of type EDGE representing an edge
    data    : Any
              Characteristic data of the given GEOM edge object depending
              on its type
    kind    : str
              Type of edge; admitted types are: CIRCLE, ARC_CIRCLE, SEGMENT
    right   : Any
              The GEOM face object to the right of the current edge
    left    : Any
              The GEOM face object to the left of the current edge
    """
    def __init__(self, edge: Any, *faces: List[Face]):
        # Get the number of the edge directly from the name attribute of the
        # corresponding GEOM edge object
        try:
            self.no = int(get_shape_name(edge).split('_')[1])
        except:
            print(f"Error on edge: {get_shape_name(edge)}")
            raise

        self.edge  = edge
        # Store all the information of the given GEOM edge object
        self.data  = get_kind_of_shape(edge)
        # Convert the first piece of information, retrieved from the GEOM edge object,
        # into string
        self.kind  = str(self.data[0])

        # Initialize to 'None' both the right and the left GEOM face objects
        self.right = None
        self.left  = None
        # Loop through all the given GEOM face objects associated to the current edge.
        # This allows to define the faces on the right and those on the left wrt the
        # edge.
        for f in faces:
            # Add a face connected to the edge
            self.__add_face(f)

        # Check everything is allright
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
        return self.data[-6:-3]

    def __add_face(self, face: Face) -> None:
        """
        Method that allows to define whether the given GEOM face object,
        connected to the edge, is placed to the right or to left of the
        same oriented edge.

        This analysis is based on the creation of a point at a very small
        distance from the edge and in a known direction (left or right).
        GEOM's 'minDistance' function is then used to find out whether
        this point is inside or outside (positive result) one of the
        connected face.
        In case the point is inside the face, then this face is on the
        right of the edge, whereas, if outside, the face is considered on
        the left.

        Parameters
        ----------
        face    : Face
                  'Face' object whose position (right or left) relative
                  to the edge has to be determined
        epsilon : float
                  Margin small enough to place a point wrt the edge
        """
        # Handle the edge types differently
        point = self.__build_point_on_edge_normal()

        # The face identification is based on the distance between the point
        # and face to analyse. It differs between a 'SEGMENT'-type edge and
        # the 'CIRCLE' and 'ARC_CIRCLE' ones.
        if get_min_distance(point, face.face) > 0:
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
        Method
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
            # Extract the X-Y-Z coordinates of the circle center the arc b
            # elongs to
            (xc, yc, zc) = self.data[1:4] # xc yc zc dx dy dz R x1 y1 z1 x2 y2 z2
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
            (x1, y1, z1, x2, y2, _) = self.data[1:7] # x1 y1 z1 x2 y2 z2
            # Build a GEOM point on the segment, positioned at its middle
            pm = make_vertex_on_curve(self.edge, 0.5)
            # Get the X-Y-Z coordinates of the point
            coord = get_point_coordinates(pm)
            # Define the normal of the segment oriented to the left of the
            # segment
            n = ((y1 - y2), (x2 - x1))
            # Get the length of the normal vector
            l = math.sqrt(n[0]*n[0] + n[1]*n[1])
            # Get the coordinates of a point positioned at a distance epsilon
            # from the segment along its normal
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

    Attributes
    ----------
    border      : Any
                  A GEOM edge object representing a border of the lattice
    type        : BoundaryType
                  Providing the type of BCs
    type_geo    : LatticeGeometryType
                  Providing the lattice type of geometry
    ox          : float
                  X-coordinate of the border start point (origin)
    oy          : float
                  Y-coordinate of the border start point (origin)
    angle       : float
                  Angle (in degrees) of the border defined from its origin to
                  its second point, so that it is always positive
    edge_indxs  : List[int]
                  List of the indices of the edges the border is made of
    lattice_o   : Any
                  A vertex object representing the lattice origin
    tx          : float
                  The X-component of the border axis
    ty          : float
                  The Y-component of the border axis
    """
    def __init__(self, border: Any, bc_type: BoundaryType,
                 type_geo: LatticeGeometryType, lattice_o: Any) -> None:
        # Initialize instance attributes
        self.border     : Any = border
        # FIXME the BC type should not be set when instantiating the class,
        # but it depends on the 'LatticeGeometryType' value
        self.type       : BoundaryType = bc_type
        self.type_geo   : LatticeGeometryType = type_geo
        self.ox         : Union[float, None] = None
        self.oy         : Union[float, None] = None
        self.angle      : float = 0.0
        self.edge_indxs : List[int] = []
        self.lattice_o  : Any = lattice_o
        self.tx         : float = 0.0
        self.ty         : float = 0.0

    def get_bc_type_number(self) -> str:
        """
        Method that returns the index associated to the BC type in the
        corresponding instance attribute dictionary.

        Returns
        -------
        A string describing the type for the BC this class instance
        refers to.
        """
        return self.type.value

    def __build_border_characteristics(self, lx: float, ly: float) -> None:
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
        """
        # Only the 'AXIAL_SYMMETRY', 'TRANSLATION' and 'VOID type of BC are
        # handled
        # print("BC TYPE =", self.type)
        if self.type not in {BoundaryType.AXIAL_SYMMETRY,
                             BoundaryType.TRANSLATION}:
            # Raise an exception if a different BC type is provided
            raise ValueError(f"The specified {self.type} BC type has not "
                             "an implemented treatment.")
        # Retrieve all the information from the GEOM edge object representing
        # one of the lattice borders
        data = get_kind_of_shape(self.border)
        # Get the X-Y coordinates of the two points of the edge
        self.ox, self.oy, _, x2, y2 = data[1:6]
        # Calculate the X-Y distances between the start-end points
        dx = x2 - self.ox
        dy = y2 - self.oy
        # Calculate the angle of the border in degrees
        self.angle = math.degrees(math.atan2(y2-self.oy, x2-self.ox))

        # The border origin must be defined so that the angle between
        # the start-end points is positive. If not, the edge start point
        # is inverted.
        if self.angle < -EPSILON or math.isclose(self.angle, 180.0):
            self.angle = math.degrees(math.atan2(self.oy-y2, self.ox-x2))
            self.ox = x2
            self.oy = y2

        # Initialize the border axes to the edge start point
        self.tx = self.ox
        self.ty = self.oy

        # Assign the BC type depending on the lattice type of geometry. The
        # geometries identified by 'RECTANGLE_TRAN' and 'HEXAGON_TRAN' need
        # to re-evaluate the border axes.
        match self.type_geo:
            case (LatticeGeometryType.SYMMETRIES_TWO |
                  LatticeGeometryType.RECTANGLE_SYM |
                  LatticeGeometryType.RECTANGLE_EIGHT |
                  LatticeGeometryType.SA60 |
                  LatticeGeometryType.S30):
                # Identifying a border characterized by the 'AXIAL_SYMMETRY'
                # type of BC, which corresponds to the 'REFL' case in DRAGON5
                self.type = BoundaryType.AXIAL_SYMMETRY
            case LatticeGeometryType.RA60 | LatticeGeometryType.R120:
                # Identifying a border characterized by any of the 'ROTATION'
                # or 'TRANSLATION' types of BC, which correspond to the 'ROTA'
                # or 'TRAN' cases respectively in DRAGON5. The position of the
                # border wrt the lattice center guides the choice.
                if get_min_distance(self.lattice_o, self.border) > 1e-7:
                    self.type = BoundaryType.TRANSLATION
                else:
                    self.type = BoundaryType.ROTATION
            case LatticeGeometryType.RECTANGLE_TRAN:
                # The BC information for the case of a cartesian geometry
                # with TRAN BCs follows the axes definition below:
                #             M=3
                #        ************
                #        *          *
                #    M=2 *          * M=4
                #        ************
                #             M=1
                if math.isclose(math.sin(math.radians(self.angle)), 0.0):
                    # The sign of 'dx' discriminates between the M=1 (dx > 0)
                    # and M=3 (dx < 0)
                    self.tx = 0.0
                    self.ty = (dx/abs(dx)) * lx
                elif math.isclose(math.sin(math.radians(self.angle)), 1.0):
                    # The sign of 'dy' discriminates between the M=4 (dy > 0)
                    # and M=2 (dy < 0)
                    self.tx = -(dy/abs(dy)) * ly
                    self.ty = 0.0
                else:
                    raise RuntimeError(
                        f"The border has an angle of {self.angle}° which is "
                        "not admitted")
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
                    raise RuntimeError("The case of an hexagonal lattice "
                        "made by cells rotated by 0° is not handled.")
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
                    f"The {self.type_geo} lattice geometry type is not "
                    "currently handled.")

    def find_edges_on_border(self,
                             boundaries: Any,
                             edge_name_vs_no: Dict[str, int],
                             id_vs_edge: Dict[Any, str],
                             dimensions: Tuple[float, float]) -> None:
        """
        Method that finds and associates all the GEOM edges related to the
        lattice border this instance refers to.
        These edges are the ones that are connected to a single face only,
        as these are the edges on the borders of the lattice.

        Parameters
        ----------
        boundaries      : Any
                          A GEOM compound object made from the list of GEOM
                          edge objects representing the lattice boundaries
        edge_name_vs_no : Dict[str, int]
                          A dictionary of edge name VS its index
        id_vs_edge      : Dict[Any, str]
                          A dictionary of all the lattice edge IDs VS the
                          corresponding GEOM edge objects
        dimensions      : Tuple[float, float]
                          The X-Y characteristic dimensions of the cell/lattice
        """
        # Call the method for defining the border characteristics
        # (origin and angle | axes)
        self.__build_border_characteristics(*dimensions)
        # Check if the boundary edge origin and its orientation angle have
        # been defined correctly.
        if (self.ox is None or self.oy is None):
            raise AssertionError("No border line defined!")

        # Build a GEOM compound object by extracting from the given lattice
        # borders the shapes that coincide with the lattice border this class
        # instance refers to or those that could be a part of it
        subshapes = get_in_place(boundaries, self.border)

        # print("\nboundaries = ", geompy.KindOfShape(boundaries))
        # print("self.border = ", geompy.KindOfShape(self.border))
        # print("len(subshapes) = ", len(geompy.SubShapeAllSortedCentres(subshapes, EDGE_TYPE)))

        # Loop through all the border sub-shapes (as GEOM objects) of type edge
        for shape in extract_sorted_sub_shapes(subshapes, ShapeType.EDGE):
            # Extract the shape type from all the information about the current
            # one
            shape_type = get_kind_of_shape(shape)[0]
            # Check the retrieved shape is of type 'SEGMENT'
            # TODO add support also for ARC_CIRCLE
            if str(shape_type) == 'SEGMENT':
                # Build the edge ID
                edge_id = build_edge_id(shape)
                # Retrieve the GEOM edge object corresponding to the built ID
                edge_ref = id_vs_edge[edge_id]
                # Get the name attribute of the GEOM edge object
                name = get_shape_name(edge_ref)

                # print("@@@\nBorder", name, "ID", id, "no", edge_name_vs_no[name], "\n")
                # name = geompy.SubShapeName(shape, partition)
                # Append the edge number to the list
                # FIXME the name contains the number -> extract from that??
                self.edge_indxs.append(edge_name_vs_no[name])
            else:
                # Raise an exception if the found sub-shape type is not
                # 'SEGMENT'
                raise AssertionError("Only edges of type 'SEGMENT' can be "
                                     "contained in a lattice border! "
                                     f"(found {shape_type})")


class LatticeDataExtractor():
    # FIXME - to transform into a dataclass since its methods have to be called
    # soon after it is instantiated in order to get the needed info
    """
    Class that extracts the needed geometrical data from the given lattice
    to be used further on for generating the output TDT file.
    It determines the association of edges with connected faces and those
    edges that represents the lattice borders.

    Parameters
    ----------
    lattice   : Lattice
                The class storing the lattice information

    Attributes
    ----------
    lattice             : Lattice
                          The class storing the lattice information
    borders             : Any
                          The GEOM compound object grouping the edges
                          representing the lattice borders
    boundaries          : List[Boundary]
                          The list of 'Boundary' objects storing the lattice
                          borders information
    subfaces            : List[Any]
                          The list of 'Face' objects storing information about
                          each subface of the lattice given its cells
                          subdivision according to either the technological
                          geometry or its sectorization
    edge_name_vs_faces  : Dict[str, List[Any]]
                          A dictionary storing all the faces connected to the
                          edges, given their name
    edges               : List[Edge]
                          The list of 'Edge' objects storing information about
                          each edge of the lattice
    id_vs_edge          : Dict[str, Any]
                          A dictionary of the lattice edge ID VS the GEOM edge
                          object reference
    """
    def __init__(self, lattice: Lattice) -> None:
        # Initialize the instance attributes
        self.lattice: Lattice = lattice
        self.borders: List[Any] = []
        self.lattice_subfaces: List[Any] = []
        self.lattice_edges: List[Any] = []
        self.boundaries: List[Boundary] = []
        self.subfaces: List[Face] = []
        self.edge_name_vs_faces: Dict[str, List[Any]] = {}
        self.edges: List[Edge] = []
        # Extract the information about borders, faces and edges of the lattice
        # to be stored according to its symmetry
        self.__preprocess()
        # Associate each edge with an index and build a dictionary
        self.id_vs_edge: Dict[str, Any] = classify_lattice_edges(
            self.lattice_edges)

    def __preprocess(self):
        """
        Method that pre-initializes the information to be stored from the
        lattice according to its symmetry.

        Parameters
        ----------
        lattice : Lattice
                  The object storing the information to pre-process.
        """
        # Get the GEOM compound identifying either the full lattice of a part
        # of it, if a symmetry is applied
        lattice_cmpd = self.lattice.lattice_cmpd
        if not self.lattice.symmetry_type == SymmetryType.FULL:
            lattice_cmpd = self.lattice.lattice_symm
            # FIXME in case of cartesian-type lattice, it must be translated
            # so that the left-most corner coincides with the XYZ origin
        # Extract the lattice borders
        self.borders = self.__build_borders(lattice_cmpd)
        # Rebuild the lattice compound by performing a 'common' operation
        # with the lattice face identified by its borders
        lattice_face = make_face(self.borders)
        lattice_cmpd = make_common(lattice_cmpd, lattice_face)

        # Extract both the lattice subfaces and the unique edges from the
        # result of the common operation
        self.lattice_subfaces = extract_sorted_sub_shapes(
            lattice_cmpd, ShapeType.FACE)
        self.lattice_edges = extract_sorted_sub_shapes(
            lattice_cmpd, ShapeType.EDGE)
        # Sort both lists in terms of the distance of each element from the
        # lattice center
        self.lattice_subfaces.sort(
            key=lambda item: get_min_distance(item,
                                              self.lattice.lattice_center))
        self.lattice_edges.sort(
            key=lambda item: get_min_distance(item,
                                              self.lattice.lattice_center))

    def __build_borders(self, lattice_cmpd: Any) -> List[Any]:
        """
        FIXME To extract a function from this method
        Method that extracts the borders of the lattice from the GEOM
        compound object storing the faces and edges of the lattice.

        Parameters
        ----------
        lattice_cmpd  : Any
                        The lattice GEOM compound object

        Returns
        -------
        The list of GEOM edge objects that represent the lattice borders.
        """
        # Extract a list of closed boundaries from the lattice compound
        closed_boundaries = get_closed_free_boundary(lattice_cmpd)
        # Handle the case where more than a closed wire is extracted from
        # the lattice compound
        if len(closed_boundaries) > 1:
            # Build a face for each of the extracted closed wires after
            # fusing adjacent edges
            shapes =  []
            for wire in closed_boundaries:
                wire_mod = fuse_edges_in_wire(wire)
                shapes.append(make_face(wire_mod))
            # Fuse all the faces into a single shape
            shape = make_fuse(shapes)
            # Return the edges of the fused shape
            return extract_sorted_sub_shapes(shape, ShapeType.EDGE)

        # Suppress vertices internal to the edges of the wire
        borders_wire = fuse_edges_in_wire(closed_boundaries[0])
        # Extract the edge objects of the border wire
        return extract_sorted_sub_shapes(borders_wire, ShapeType.EDGE)

    def build_faces(self,
                    geom_type : GeometryType,
                    property_type : PropertyType) -> None:
        """
        Method that builds a list of 'Face' objects from the GEOM face objects,
        extracted from the lattice regions, and the property values.
        Each region in the lattice corresponds to a region in a cell,
        according to the given type of geometry (either technological or
        sectorized).
        Each region must be associated with a value for the given property
        type. If not the case, an exception is raised.

        Parameters
        ----------
        geom_type : GeometryType
            The type of geometry of the lattice cells identifying the regions
        property_type : PropertyType = PropertyType.MATERIAL
            The type of property associated to the lattice regions

        Raises
        ------
        RuntimeError:
            - If no properties are associated to a lattice region.
            - If no value for the given property type is associated to a
              lattice region.
        """
        # FIXME Lattice 'Region' objects are built only if the 'show()' method
        # is called. If not the case, the analysis could be performed with an
        # incorrect number of regions. Use the given 'geom_type' parameter
        # to rebuild the regions. See issue #16.
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

    def build_edges_and_faces_association(self) -> Dict[str, List[Face]]:
        """
        Method that constructs a list of 'Edge' objects from the GEOM edge
        extracted from the GEOM faces of the lattice.

        Given the constructed 'Edge' objects, a dictionary is built, its keys
        being the name of the edge; this contains a global index to identify
        the edge in the lattice.
        For each edge name as the key, the corresponding faces connected to it
        are stored as the dictionary values.

        Returns
        -------
        A dictionary of edge name VS the list of faces, as 'Face' objects,
        connected to the corresponding edge.
        """
        # Initialize the dictionary storing the edges names VS the list of
        # connected faces
        edges_name_vs_faces : Dict[str, List[Face]] = {}

        # Loop through all the lattice subfaces (Face objects)
        print("LEN SUBFACES", len(self.subfaces))
        for subface in self.subfaces:
            print(subface)
            # Loop through all the edges of the current subface
            for subface_edge, edge_id in subface.edge_vs_id.items():
                # Extract the corresponding GEOM edge object in a try-except
                # block. In case of an exception for keys not present, a
                # further analysis is performed
                try:
                    unique_edge = self.id_vs_edge[edge_id]
                except KeyError as exc:
                    # print("KEYERROR ->", id)
                    # If not found, it could be the case of a border edge
                    # between two adiacent cells. The dictionary stores the
                    # unique IDs only, which means only one of the two is
                    # present.
                    # Common edges belonging to different cells could not have
                    # the same orientation, i.e. starting and ending points
                    # are inverted.
                    # When looking for the edge by extracting it from the
                    # subfaces, the same border edge is retrieved two times
                    # but since the orientations are different, one cannot
                    # be found in the dictionary.

                    # Find the lattice edge, from the compound made of the
                    # unique edges, that corresponds to the one of the subface
                    # being analysed.
                    # This is needed because edges shared between adiacent
                    # cells are considered only once in the lattice edges
                    # compound with a specific orientation. When analysing the
                    # edges of the subface, the orientation is different, so
                    # it is needed to retrieve the corresponding one in the
                    # lattice which was used to build the dictionary of edges
                    # IDs VS GEOM edges.
                    edge = get_in_place(self.lattice_edges, subface_edge)
                    if not edge:
                        raise AssertionError("No corresponding edge in the "
                                            "lattice could be retrieved for "
                                            "the subface edge whose data is "
                                            f"{edge_id}") from exc
                    # Build the ID
                    edge_id = build_edge_id(edge)
                    # Extract the stored edge in the ID VS edge dictionary
                    unique_edge = self.id_vs_edge[edge_id]

                # Get the edge name containing a global index of all the unique
                # edges in the lattice
                edge_name = get_shape_name(unique_edge)
                # Check if the edge ID is already stored, if so, append the
                # current face, otherwise add a new entry
                if edge_name in edges_name_vs_faces:
                    # The corresponding entry in the dictionary is updated
                    # by appending the found 'Face' object
                    edges_name_vs_faces[edge_name].append(subface)
                else:
                    # Declare a new entry for the dictionary of edge names
                    # VS list of 'Face' objects. The entry value is
                    # initialized as a list containing the GEOM edge object
                    # and the found 'Face' object.
                    edges_name_vs_faces[edge_name] = [unique_edge, subface]
        # -----------------------------------------------------------
        # Build the edges in the internal data model ('Edge' objects)
        # -----------------------------------------------------------
        # Loop through all the lists of faces associated to each edge in the
        # lattice
        print("LEN edges_name_vs_faces:", len(edges_name_vs_faces.values()))
        for shapes in edges_name_vs_faces.values():
            # Instantiate an object of the 'Edge' class, given the list of
            # associated shapes. The '*' is for unpacking the list, since
            # the first element is the GEOM edge object, while the others
            # are the face ones.
            e = Edge(*shapes)
            # Append the just built 'Edge' object to the list of edges in
            # the lattice
            self.edges.append(e)

        # Return the dictionary of edge names VS the list of connected faces
        return edges_name_vs_faces

    def build_boundaries(self) -> None:
        """
        Method that constructs a list of boundary edges (connected to a single
        face). These edges are then sorted by boundary condition plane.
        """
        # No boundaries to extract if an 'ISOTROPIC' type of geometry, meaning
        # 'VOID' or 'ALBE 1.0' BCs in DRAGON5.
        if self.lattice.type_geo == LatticeGeometryType.ISOTROPIC:
            return
        # Initialize the list of 'Boundary' objects
        boundary_edges = []
        # Initialize a dictionary of edge name VS its global index
        edge_name_to_no = {}

        # Loop through all the stored 'Edge' objects
        print("LEN EDGES", len(self.edges))
        for edge in self.edges:
            # Check if the current edge has no faces on the left or on the
            # right, i.e. it represents a boundary edge
            if not edge.left or not edge.right:
                # Append the corresponding GEOM edge object to the list
                boundary_edges.append(edge.edge)
                # Get the name attribute from the GEOM edge object
                name = edge.edge.GetName()
                # Add the edge name-number entry to the dictionary
                edge_name_to_no[name] = edge.no

        # Build a compound from the list of boundary edges
        outline = make_compound(boundary_edges)

        # Loop through all the GEOM line objects representing the lattice edges
        print("LEN BORDERS:", len(self.borders))
        for border in self.borders:
            # Build an object of the 'Boundary' class
            boundary = Boundary(border=border,
                                bc_type=self.lattice.boundary_type,
                                type_geo=self.lattice.type_geo,
                                lattice_o=self.lattice.lattice_center)
            boundary.find_edges_on_border(
                outline, edge_name_to_no, self.id_vs_edge,
                (self.lattice.lx, self.lattice.ly))
            # Append the just built 'Boundary' object to the corresponding list
            # print("  --> ", boundary.tx, boundary.ty, boundary.angle)
            self.boundaries.append(boundary)

    def print_log_analysis(self, edge_name_vs_faces: Dict[str, List[Any]]):
        """
        Method that prints on the stdout the log of the data extraction from
        the lattice.

        Parameters
        ----------
        edge_name_vs_faces  : Dict[str, List[Any]]
                              A dictionary of the name of the lattice edges
                              VS the GEOM faces connected to it
        """
        # Displays the number of faces, the number of edges and the edges associated
        # with one face and two faces.
        n1 = n2 = n0 = 0
        for f in edge_name_vs_faces.values():
            if len(f) == 2:
                n1 += 1
            elif len(f) == 3:
                n2 += 1
            else:
                n0 += 1
        print("\t# subparts (faces) : ", len(self.subfaces))
        print("\t# edges found : ", len(self.edges))
        print("\t# edges w/two faces :", n2)
        print("\t# edges w/one faces :", n1)
        print("\t# edges w errors :", n0)


def analyse_lattice(lattice: Lattice,
                    geom_type : GeometryType,
                    property_type : PropertyType) -> LatticeDataExtractor:
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
    # Instantiate the class for extracting the data from the lattice and
    # analyse them
    data_extractor = LatticeDataExtractor(lattice)
    # Call its method for performing the analysis
    data_extractor.build_faces(geom_type, property_type)
    edge_name_vs_faces = data_extractor.build_edges_and_faces_association()
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
        # TODO Qui potrei costruire oggetto della classe 'Edge' passando indice
        # globale e oggetto GEOM, da cui estrarre ID

        # Set the name of the GEOM edge object
        set_shape_name(edge, f"EDGE_{indx + 1}")
        # Build a unique ID for the edge
        edge_id = build_edge_id(edge)
        # Add an entry to the dictionary storing edges ID VS the
        # corresponding GEOM edge objects
        ids_edges[edge_id] = edge
    # Return the built dictionary
    return ids_edges
