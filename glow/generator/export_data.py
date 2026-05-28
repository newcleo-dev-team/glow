"""
Module containing classes that support the export of the geometry layout by
providing a mean for storing the needed information about faces, edges and
boundaries.
"""
import math

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple, Self

from glow.geometry_layouts.cells import Region
from glow.support.types import EDGE_NAME_VS_TYPE, BoundaryType, EdgeType, \
    LayoutGeometryType, PropertyType
from glow.support.utility import check_shape_expected_types, \
    get_angle_between_points, get_id_from_name, get_id_from_shape, \
    is_vertex_on_edge
from glow.interface.geom_interface import ShapeType, \
    extract_sorted_sub_shapes, extract_sub_shapes, \
    get_in_place, get_kind_of_shape, get_point_coordinates, get_shape_name, \
    is_point_inside_shape, make_vertex, make_vertex_inside_face, \
    make_vertex_on_curve, set_shape_name


# Sufficiently small value used to determine face-edge connectivity by
# building a point on an edge and shifting to one side of the epsilon.
EPSILON = 1e-05


@dataclass(order=True)
class FaceData():
    """
    Class that provides a data representation for a face object of the layout,
    i.e. a ``Region`` object, which represents a calculation zone containing a
    list of properties.
    This dataclass can be ordered on the basis of the ``no`` attribute, which
    provides a global index for the faces in the layout.
    """
    region: Region
    """
    A ``Region`` object associating a GEOM face of the layout to its
    properties.
    """
    no: int
    """Global index of the face."""
    property_types: List[PropertyType] = field(default_factory=list)
    """List of property types that must be associated with the region."""
    inner_point: Tuple[float, float, float] = field(init=False)
    """XYZ coordinates of a point within the GEOM face object."""
    edge_vs_id: Dict[Any, str] = field(default_factory=dict, init=False)
    """
    Dictionary associating the GEOM edge objects of the current face object
    to the corresponding IDs, based on the edges geometrical characteristics.
    """
    sort_index: int = field(init=False, repr=False)
    """Indicating the attribute used to sort ``FaceData`` objects."""

    def __post_init__(self) -> None:
        """
        Method that is automatically run after the dataclass initialization
        for setting all the attributes that depends on others. In addition,
        the attribute that allows to order instances of this class on the
        basis of the ``no`` attribute is set as well.

        Raises
        ------
        RuntimeError
            If no properties or none of the requested properties or no value
            for the property have been assigned to the region of the current
            instance.
        """
        # Build a point inside the face
        self.inner_point = get_point_coordinates(
            make_vertex_inside_face(self.region))
        # Check the region has the required properties
        if not self.region.properties:
            raise RuntimeError(
                "No properties have been assigned for the region "
                f"'{self.region.name}' of {self}."
            )
        for p_type in self.property_types:
            try:
                value = self.region.properties[p_type]
            except KeyError:
                raise RuntimeError(
                    f"No '{p_type.name}' property type has been defined for "
                    f"the region '{self.region.name}' of {self}."
                )
            if not value:
                raise RuntimeError(
                    f"No value for the '{p_type.name}' property type has "
                    f"been defined for the region '{self.region.name}' of "
                    f"{self}."
                )
        # Define the attribute for sorting instances of this class
        self.sort_index = self.no
        # Extract the edges and associate an ID to each
        edges = extract_sub_shapes(self.region, ShapeType.EDGE)
        for edge in edges:
            # Build the ID of the edge, based on its geometric characteristics
            self.edge_vs_id[edge] = build_edge_id(edge)

    def __str__(self) -> str:
        return f"Region {self.no}, name={get_shape_name(self.region)}, " + \
               f"properties={self.region.properties}"


class EdgeData():
    """
    Class that provides the data representation for an edge object of the
    layout.
    It provides an global index allowing to uniquely identify the edge, and
    two ``FaceData`` attributes, allowing to identify the face on the left
    and on the right of the edge.
    The identification of the left/right faces is performed by building
    a point on the edge's normal so that it is sligtly on the left of the
    edge.

    Parameters
    ----------
    edge : Any
        The GEOM object of type EDGE representing an edge.
    faces : Tuple[FaceData, ...]
        Providing the ``FaceData`` objects the edge belongs to.

    Attributes
    ----------
    no : int
        The edge number extracted from the GEOM edge object name.
    edge : Any
        A GEOM object of type EDGE representing an edge.
    data : Any
        Characteristic data of the given GEOM edge object.
    kind : str
        Indicating the type of edge with admitted values being ``CIRCLE``,
        ``ARC_CIRCLE``, ``SEGMENT``.
    right : FaceData | None
        A ``FaceData`` object providing the information for the face to the
        right of the edge, or ``None`` if the edge does not have any face on
        its right.
    left : FaceData | None
        A ``FaceData`` object providing the information for the face to the
        left of the edge, or ``None`` if the edge does not have any face on
        its left.
    """
    def __init__(self, edge: Any, *faces: FaceData) -> None:
        # Store all the information of the given GEOM edge object
        self.data: List[Any] = get_kind_of_shape(edge)
        no_faces = len(faces)
        try:
            # Check the type of the received edge is correct
            check_shape_expected_types(edge, [ShapeType.EDGE])
            # Get the number of the edge directly from the name attribute of
            # the corresponding GEOM edge object
            self.no: int = get_id_from_shape(edge)
            # Check the number of associated faces is one or two
            if no_faces not in [1, 2]:
                raise RuntimeError(
                    f"An invalid number of face objecs ({no_faces}) is "
                    "associated to the edge."
                )
        except RuntimeError as e:
            raise RuntimeError(
                f"Error with 'EdgeData' whose data is: {self.data}"
            ) from e
        # Store the edge object
        self.edge: Any  = edge
        # Get the type of edge
        self.kind: EdgeType = EDGE_NAME_VS_TYPE[str(self.data[0])][0]
        # Initialize to 'None' both the right and the left 'FaceData' objects
        self.right: FaceData | None = None
        self.left: FaceData | None = None
        # Loop through all the given GEOM face objects associated to the
        # current edge. This allows to define the faces on the right and
        # those on the left wrt the edge.
        for f in faces:
            # Add a face connected to the edge
            self._determine_face_position(f)
        # Check that both right and left faces have been assigned, if the
        # edge is shared by two faces
        if no_faces > 1 and (not self.right or not self.left):
            raise RuntimeError(
                f"The edge no. {self.no} has 2 faces (no. "
                f"{[face.no for face in faces]}), but they have not be "
                "correctly assigned to the left and right attributes "
                f"(left = {self.left}, right = {self.right})."
            )

    def _build_point_on_edge_normal(self, epsilon: float = EPSILON) -> Any:
        """
        Method that builds a vertex object positioned at an infinitesimal
        distance from the edge object this instance refers to.
        Depending on the edge type, this point is built according to the
        following rules:
        - ``CIRCLE``: the point is positioned slightly on the left wrt the
          X-coordinate of the point laying on the right-most position of
          the circle (identified by the centre X-coordinate + the radius).
        - ``ARC_CIRCLE``: given a point positioned at the middle of the arc,
          another point is built on the vector connecting the arc centre and
          the first point. This second point has both its X-Y coordinates
          slightly scaled down by a reduction factor so that its distance
          from the arc centre is less than the radius.
        - ``SEGMENT``: the normalized left-oriented vector normal to the edge
          is calculated. Its X-Y components are used to determine a point
          slightly to the left of a point positioned at the middle of the
          segment.

        Parameters
        ----------
        epsilon : float
            Indicating a value small enough to determine a point to the left
            of the edge.

        Returns
        -------
        Any
            A vertex object representing a point slightly to the left of the
            middle point of the edge.
        """
        if self.kind == EdgeType.CIRCLE:
            # Extract the X-Y-Z coordinates of the circle centre
            (xc, yc, zc) = self.data[1:4]
            # Extract the circle radius
            radius = self.data[7]
            # Build a point (as a GEOM object) positioned at an infinitesimal
            # distance from the starting point of the circle
            return make_vertex((xc+radius-epsilon, yc, zc))
        if self.kind == EdgeType.ARC_CIRCLE:
            # Extract the X-Y-Z coordinates of the circle centre the arc
            # belongs to
            (xc, yc, zc) = self.data[1:4]
            # Extract the radius of the arc
            radius = self.data[7]
            # Build a GEOM point on the arc, positioned at its middle
            pm = make_vertex_on_curve(self.edge, 0.5)
            # Get the X-Y-Z coordinates of the point
            coord = get_point_coordinates(pm)
            # Build a 2D vector from the arc centre to the arc middle point,
            # which is slightly scaled down by a reduction factor so that
            # its length is less than the radius
            vm = (
                (coord[0] - xc)*(1.0 - epsilon/radius),
                (coord[1] - yc)*(1.0 - epsilon/radius)
            )
            # Build a GEOM point along the direction of the just built vector
            return make_vertex((xc + vm[0], yc + vm[1], zc))
        if self.kind == EdgeType.SEGMENT:
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
            (xm, ym) = (
                (coord[0] + epsilon*n[0]/l), (coord[1] + epsilon*n[1]/l)
            )
            # Build a GEOM point positioned at the calculated coordinates
            return make_vertex((xm, ym, z1))
        # If here, raise an exception as the edge has not a valid type
        raise ValueError(f"The kind of shape {self.kind} is not valid!")

    def _determine_face_position(self, face: FaceData) -> None:
        """
        Method that allows to define whether the given GEOM face object,
        connected to the edge, is placed to the right or to left of the
        edge object the instance refers to.

        This analysis is based on the creation of a point at a very small
        distance from the edge middle point so that it belongs to the
        left-oriented normal vector for the edge.
        The identification of the face position (left or right) relative
        to the edge depends on the edge nature:
        - for ``SEGMENT``-type edges, the given face is considered on the
          left of the edge if the point belongs to the face; on the contrary,
          we have a right face;
        - for ``ARC_CIRCLE`` and ``CIRCLE``-type edges, the criteria is the
          opposite, i.e. the face is considered on the right of the edge if
          the point belongs to the face; on the contrary, we have a left face.

        Parameters
        ----------
        face : FaceData
            The ``FaceData`` object whose position (right or left) relative
            to the edge has to be determined.
        epsilon : float
            Margin small enough to place a point wrt the edge.
        """
        # Build the point on the left of the edge and use it to identify the
        # face position relative to the edge
        if not is_point_inside_shape(
            self._build_point_on_edge_normal(), face.region
        ):
            if self.kind == EdgeType.SEGMENT:
                self.right = face
            else:
                self.left = face
        else:
            if self.kind == EdgeType.SEGMENT:
                self.left = face
            else:
                self.right = face

    def _get_origin(self) -> Tuple[float, float, float]:
        """
        Method that retrieves the 'origin' point (X-Y-Z coordinates) of
        the edge (i.e. the edge starting point) according to its type.
        This information is retrieved from the instance attribute storing
        the characteristic data of the edge.

        Returns
        -------
        Tuple[float, float, float]
            A tuple with the X-Y-Z coordinates of the edge starting point.
        """
        # Handle the point retrieval differently if the edge is a circle
        if self.kind == EdgeType.CIRCLE:
            # Get the X-Y-Z coordinates of the circle centre
            x1, y1, z1  = self.data[1:4]
            # Get the circle radius
            r = self.data[-1]
            # Return the point (cx + r, cy, cz) with cx, cy, cz the centre
            # of the circle and r its radius.
            return x1 + r, y1, z1
        # Handle the other edge types by returning the X-Y-Z coordinates of
        # the first point of the edge
        return tuple(self.data[-6:-3])

    def __str__(self) -> str:
        # Get the name of the current edge and the string representation of
        # the associated 'FaceData' objects on its right and on its left
        name = get_shape_name(self.edge)
        right = str(self.right) if self.right else "None"
        left = str(self.left) if self.left else "None"
        return (
            f"Edge {name} \n\torigin = {self._get_origin()},"
            f"\n\tleft face = {left}, \n\tright face = {right}"
        )

    # ------------------------------------------
    # Methods for comparing two 'Edge' instances
    # ------------------------------------------
    def __lt__(self, other: Self) -> bool:
        return self.no < other.no

    def __le__(self, other: Self) -> bool:
        return self.no <= other.no

    def __eq__(self, other: Self) -> bool:
        return self.no == other.no

    def __ne__(self, other: Self) -> bool:
        return self.no != other.no

    def __gt__(self, other: Self) -> bool:
        return self.no > other.no

    def __ge__(self, other: Self) -> bool:
        return self.no >= other.no


class BoundaryData:
    """
    Class that provides the data structure for a GEOM edge object being a
    border of the layout.

    Parameters
    ----------
    border : Any
        A GEOM edge object representing a border of the layout.
    type_geo : LatticeGeometryType
        Providing the layout type of geometry.
    layout_o : Any
        A vertex object representing the layout centre point.
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the layout the border refers to.

    Attributes
    ----------
    border : Any
        A GEOM edge object representing a border of the layout.
    type : BoundaryType
        Providing the type of BCs.
    angle : float
        Angle (in degrees) of the border defined from its origin to its
        second point, so that it is always positive.
    edge_indxs : List[int]
        List of the indices of the edges the border is made of.
    tx : float
        The X-component of the border axis.
    ty : float
        The Y-component of the border axis.
    """
    def __init__(
            self,
            border: Any,
            type_geo: LayoutGeometryType,
            layout_o: Any,
            dimensions: Tuple[float, float]
        ) -> None:
        # Check the received border is an edge
        try:
            check_shape_expected_types(border, [ShapeType.EDGE])
        except RuntimeError as e:
            raise RuntimeError(
                "Error while initializing the 'BoundaryData' instance."
            ) from e
        # Initialize instance attributes
        self.type : BoundaryType
        self.border : Any = border
        self.angle : float = 0.0
        self.edge_indxs : List[int] = []
        self.tx : float = 0.0
        self.ty : float = 0.0
        # Set the border characteristics in terms of type and border axis
        self._build_border_characteristics(*dimensions, layout_o, type_geo)

    def find_edges_on_border(
            self, boundaries: Any, id_vs_edge: Dict[str, Any]
        ) -> None:
        """
        Method that finds and associates all the GEOM edges related to the
        layout border this instance refers to.
        These edges are the ones that are connected to a single face only,
        as these are the edges on the borders of the layout.

        Parameters
        ----------
        boundaries : Any
            A GEOM compound object made from the list of GEOM edge objects
            belonging to the layout boundaries.
        id_vs_edge : Dict[Any, str]
            A dictionary of all the layout edge IDs VS the corresponding
            GEOM edge objects.

        Raises
        ------
        RuntimeError
            In case the boundary edge is not of type ``SEGMENT``.
        """
        # Loop through all the sub-edges that are part of the layout border
        # this class instance refers
        for shape in extract_sorted_sub_shapes(
            get_in_place(boundaries, self.border), ShapeType.EDGE):
            # Extract the shape type from all the information about the current
            # one
            shape_type = get_kind_of_shape(shape)[0]
            # Check the retrieved shape is of type 'SEGMENT'
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
                    f"layout border! (found {shape_type})"
                )

    def get_bc_type_number(self) -> int:
        """
        Method that returns the index associated to the BC type in the
        corresponding instance attribute dictionary.

        Returns
        -------
        int
            An integer describing the type for the BC this class instance
            refers to.
        """
        return self.type.value

    def _build_border_characteristics(
            self,
            lx: float,
            ly: float,
            layout_o: Any,
            type_geo: LayoutGeometryType
        ) -> None:
        """
        Method that defines the characteristics of a layout border that
        represents a boundary for the layout itself.
        These characteristics are defined in terms of the associated BC type
        (value of the ``BoundaryType`` enumeration), the border axes and its
        angle; these depend on the specific layout type of geometry (as
        value of the ``LatticeGeometryType`` enumeration).

        Parameters
        ----------
        lx : float
            The X-characteristic dimension of the layout.
        ly : float
            The Y-characteristic dimension of the layout.
        layout_o : Any
            Vertex object representing the layout centre point.
        """
        # Get the X-Y coordinates of the two end points of the edge
        x1, y1, _, x2, y2 = get_kind_of_shape(self.border)[1:6]
        # Calculate the X-Y distances between the start-end points
        dx = x2 - x1
        dy = y2 - y1
        # Calculate the angle of the border in degrees
        self.angle = get_angle_between_points(
            (x1, y1, 0.0), (x2, y2, 0.0), True
        )

        # The border origin must be defined so that the angle between
        # the start-end points is positive. If not, the edge start point
        # is inverted.
        if (
            self.angle < -EPSILON
            or math.isclose(self.angle, 180.0, abs_tol=EPSILON)
        ):
            self.angle = get_angle_between_points(
                (x2, y2, 0.0), (x1, y1, 0.0), True
            )
            x1 = x2
            y1 = y2
        # If the angle is close to zero, it can happen that it is written
        # with a minus sign which is unnecessary; hence, it is set to 0.0
        # explicitly
        if math.isclose(abs(self.angle), 0.0, abs_tol=EPSILON):
            self.angle = 0.0

        # Initialize the border axes to the edge start point
        self.tx = x1
        self.ty = y1

        # Assign the BC type depending on the layout type of geometry. The
        # geometries identified by 'RECTANGLE_TRAN' and 'HEXAGON_TRAN' need
        # to re-evaluate the border axes.
        match type_geo:
            case (LayoutGeometryType.SYMMETRIES_TWO |
                  LayoutGeometryType.RECTANGLE_SYM |
                  LayoutGeometryType.RECTANGLE_EIGHT |
                  LayoutGeometryType.SA60 |
                  LayoutGeometryType.S30):
                # Identifying a border characterized by the 'AXIAL_SYMMETRY'
                # type of BC, which corresponds to the 'REFL' case in DRAGON5
                self.type = BoundaryType.AXIAL_SYMMETRY
            case (LayoutGeometryType.RA60 |
                  LayoutGeometryType.R120 |
                  LayoutGeometryType.ROTATION):
                # Identifying a border characterized by any of the 'ROTATION'
                # or 'TRANSLATION' types of BC, which correspond to the 'ROTA'
                # or 'TRAN' cases respectively in DRAGON5. The position of the
                # border wrt the layout centre guides the choice.
                if is_vertex_on_edge(layout_o, self.border):
                    self.type = BoundaryType.ROTATION
                else:
                    self.type = BoundaryType.TRANSLATION
            case LayoutGeometryType.RECTANGLE_TRAN:
                # The BC information for the case of a cartesian geometry
                # with TRAN BCs follows the axes definition below:
                #              M=3 (0,-ly)
                #             ************
                #  M=2 (lx,0) *          * M=4 (-lx,0)
                #             ************
                #              M=1 (0,ly)
                if math.isclose(
                    math.sin(math.radians(self.angle)), 0.0, abs_tol=1e-6
                ):
                    # The sign of 'dx' discriminates between the M=1 (dx > 0)
                    # and M=3 (dx < 0)
                    self.tx = 0.0
                    self.ty = (dx/abs(dx)) * ly
                elif math.isclose(
                    math.sin(math.radians(self.angle)),
                    1.0,
                    abs_tol=EPSILON
                ):
                    # The sign of 'dy' discriminates between the M=4 (dy > 0)
                    # and M=2 (dy < 0)
                    self.tx = -(dy/abs(dy)) * lx
                    self.ty = 0.0
                else:
                    raise RuntimeError(
                        f"The border has an angle of {self.angle}° which is "
                        "not one of the admitted values (0°, 90°) for a "
                        "cartesian geometry with TRAN as BC."
                    )
                # Assign the BC type
                self.type = BoundaryType.TRANSLATION
            case LayoutGeometryType.HEXAGON_TRAN:
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
                if math.isclose(
                    math.sin(math.radians(self.angle)),
                    1.0,
                    abs_tol=EPSILON
                ):
                    raise RuntimeError(
                        "The border refers to a Y-oriented hexagon which "
                        "is not admitted for tracking."
                    )
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
                    f"The {type_geo} layout geometry type is not "
                    "currently handled."
                )


# -------------------------------------------------------------------------- #
#                                FUNCTIONS                                   #
# -------------------------------------------------------------------------- #

def build_edge_id(edge: Any) -> str:
    """
    Function that builds a unique ID for a GEOM edge object in the geometry
    layout to process.
    This is performed by retrieving characteristic information about the
    edge in terms of a list containing the type of shape and a series of
    parameters that describe the shape itself.
    According to the shape type we could have:

    - `CIRCLE xc yc zc dx dy dz R`
      (X-Y-Z centre coordinates, X-Y-Z normal vector elements, circle radius).
    - `ARC_CIRCLE xc yc zc dx dy dz R x1 y1 z1 x2 y2 z2`
      (X-Y-Z centre coordinates, X-Y-Z normal vector elements, arc radius,
      X-Y-Z coordinates of arc starting and ending points).
    - `SEGMENT x1 y1 z1 x2 y2 z2`
      (X-Y-Z coordinates of segment starting and ending points).

    If any shape type other than ``CIRCLE``, ``ARC_CIRCLE`` and ``SEGMENT``
    is provided, an exception is raised.

    Parameters
    ----------
    edge : Any
        A GEOM edge object to build an ID for.

    Returns
    -------
    str
        A string representing the unique ID for the edge.

    Raises
    ------
    RuntimeError
        If the type of the provided edge is not among the allowed ``CIRCLE``,
        ``ARC_CIRCLE`` and ``SEGMENT``.
    """
    # Get the information of the given GEOM shape object
    data = get_kind_of_shape(edge)
    # Check if the shape is one of the admitted edges
    if str(data[0]) not in ['CIRCLE', 'ARC_CIRCLE', 'SEGMENT']:
        raise RuntimeError(
            f"The shape, whose information is '{data}', is not one of the "
            "admitted edges types 'CIRCLE', 'ARC_CIRCLE', 'SEGMENT'")
    # Loop through all the other information about the edge and build
    # the ID by appending all the info
    return "EDGE_" + str(data[0]) + "_" + "_".join(f"{info:.6g}"
                                                   for info in data[1:])

def classify_layout_edges(edges: List[Any]) -> Dict[str, Any]:
    """
    Function that classifies the given layout edges (as GEOM objects)
    by building a dictionary with keys being a unique ID and values the
    corresponding GEOM edge object.
    The IDs are built from the geometric characteristics of the edges;
    for each edge its ``name`` internal attribute is set using a global
    index over all the layout's edges.

    Parameters
    ----------
    edges : List[Any]
        A list of the GEOM edge objects to be classified.

    Returns
    -------
    Dict[str, Any]
        A dictionary storing the edges IDs and the corresponding GEOM edges.

    Raises
    ------
    RuntimeError
        If any of the elements in the input list is not an edge of the
        allowed ``CIRCLE``, ``ARC_CIRCLE`` and ``SEGMENT`` types.
    """
    # Initialize the dictionary of GEOM edges
    ids_edges = {}
    try:
        for indx, edge in enumerate(edges):
            # Set the name of the GEOM edge object
            set_shape_name(edge, f"EDGE_{indx + 1}")
            # Add an entry to the dictionary storing edges ID VS the
            # corresponding GEOM edge objects
            ids_edges[build_edge_id(edge)] = edge
    except RuntimeError as e:
        raise RuntimeError(
            "Error while classifying the layout's edges.") from e
    # Return the built dictionary
    return ids_edges
