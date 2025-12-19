"""
Module containing utility functions to support the construction and analysis
of the geometry layouts.
"""
import math
import random

from types import CellType
from typing import Any, List, Tuple

from glow.interface.geom_interface import ShapeType, add_to_study, \
    extract_sorted_sub_shapes, extract_sub_shapes, fuse_edges_in_wire, \
    get_angle_between_shapes, get_basic_properties, get_closed_free_boundary, get_kind_of_shape, \
    get_min_distance, get_point_coordinates, get_selected_object, \
    get_shape_name, get_shape_type, is_gui_available, make_arc_edge, make_cdg, make_compound, \
    make_cut, make_edge, make_face, make_fuse, make_partition, \
    make_translation, make_vector_from_points, make_vertex
from glow.support.types import CELL_VS_SYMM_VS_TYP_GEO, LatticeGeometryType, \
    SymmetryType


# List of the RGB color codes taken by varying each RGB value with steps of 10
# in the range (0-256)
RGB_COLORS = [(r, g, b) for r in range(0, 256, 10)
                        for g in range(0, 256, 10)
                        for b in range(0, 256, 10)]


def are_same_shapes(shape1: Any, shape2: Any, shapes_type: ShapeType) -> bool:
    """
    Function that determines whether two shapes are the same based on their
    geometric properties (perimeter, area and volume), a cut operation, and
    the shape type.

    Parameters
    ----------
    shape1 : Any
        The first shape to compare.
    shape2 : Any
        The second shape to compare.
    shapes_type : ShapeType
        The type of sub-shapes to check for after the cut operation.

    Returns
    -------
    bool
        ``True`` if the two shapes are considered the same (i.e., they have
        the same geometric properties and the cut result does not contain
        any sub-shapes of the specified type), ``False`` otherwise.

    Raises
    ------
    RuntimeError
        If the two shapes are not of the same type.
    """
    # Check whether the two shapes have the same type
    if get_shape_type(shape1) != get_shape_type(shape2):
        raise RuntimeError("Shapes not compatible")
    # Check whether the two shapes has same perimeter, area and volume
    for bp1, bp2 in zip(get_basic_properties(shape1),
                        get_basic_properties(shape2)):
        if not math.isclose(bp1, bp2, abs_tol=1e-5):
            return False
    # Perform a cut operation
    cut = make_cut(shape1, shape2)
    # If they are the same shape, the cut result must be a compound that
    # have no shapes of the indicated type
    if get_shape_type(cut) != ShapeType.COMPOUND:
        return False
    return len(extract_sub_shapes(cut, shapes_type)) == 0


def build_arcs_for_rounded_corners(
        rounded_corners: List[Tuple[int, float]],
        center: Tuple[float, float, float],
        height: float,
        width: float) -> List[Any]:
    """
    Function that builds the arcs (as GEOM edge objects) that represent the
    rounded corners of a rectangle with given dimensions.
    The information about which corners are rounded and the corresponding
    radius is provided by the `rounded_corners` input list.

    Parameters
    ----------
    rounded_corners : List[Tuple[int, float]]
        List of tuples containing the rectangle corner ID and the
        corresponding radius.
    center : Tuple[float, float, float]
        The XYZ coordinates of the center of the rectangle.
    height : float
        The height of the rectangle.
    width : float
        The width of the rectangle.

    Returns
    -------
    List[Any]
        A list of GEOM edge objects representing the arcs of circle for
        the rounded corners on the borders of the rectangle.
    """
    arcs = []
    max_radius = min(width/2, height/2)
    for rc in rounded_corners:
        # Check the radius of the corner does not exceed the minimum
        # among the characteristic dimensions
        if rc[1] > max_radius:
            raise RuntimeError(
                f"The corner no. {rc[0]} has a radius of {rc[1]} that "
                f"exceeds the maximum allowed value of {max_radius}.")
        # Build the XY coordinates of the arc's start/end points and its
        # center
        match rc[0]:
            case 0:
                # Calculate the XY coordinates of the corner
                x0 = center[0] - width/2
                y0 = center[1] - height/2
                # Build the GEOM objects to create an arc
                xy1 = x0, y0 + rc[1]
                xy2 = x0 + rc[1], y0
                cxy = (x0 + rc[1], y0 + rc[1])
            case 1:
                # Calculate the XY coordinates of the corner
                x0 = center[0] + width/2
                y0 = center[1] - height/2
                # Build the GEOM objects to create an arc
                xy1 = x0 - rc[1], y0
                xy2 = x0, y0 + rc[1]
                cxy = (x0 - rc[1], y0 + rc[1])
            case 2:
                # Calculate the XY coordinates of the corner
                x0 = center[0] + width/2
                y0 = center[1] + height/2
                # Build the GEOM objects to create an arc
                xy1 = x0, y0 - rc[1]
                xy2 = x0 - rc[1], y0
                cxy = (x0 - rc[1], y0 - rc[1])
            case 3:
                # Calculate the XY coordinates of the corner
                x0 = center[0] - width/2
                y0 = center[1] + height/2
                # Build the GEOM objects to create an arc
                xy1 = x0 + rc[1], y0
                xy2 = x0, y0 - rc[1]
                cxy = (x0 + rc[1], y0 - rc[1])
        v1 = make_vertex((*xy1, 0.0))
        v2 = make_vertex((*xy2, 0.0))
        cv = make_vertex((*cxy, 0.0))
        # Add the arc edge to the list
        arcs.append(make_arc_edge(cv, v1, v2))
    # Return the built arcs
    return arcs


def build_compound_borders(cmpd: Any) -> List[Any]:
    """
    Function that extracts the external borders of the given compound object
    and returns them as a list of edge objects.

    Parameters
    ----------
    cmpd : Any
        The compound object whose borders to extract.

    Returns
    -------
    List[Any]
        A list of edge objects representing the given compound external
        boundary.

    Raises
    ------
    RuntimeError
        If no closed boundaries could be extracted from the given compound
        object.
    """
    # Extract a list of closed boundaries from the given compound object
    closed_boundaries = get_closed_free_boundary(
        make_partition([cmpd], [], ShapeType.COMPOUND))
    # Handle the case where more than a closed wire is extracted from
    # the compound
    if len(closed_boundaries) > 1:
        # Build a face for each of the extracted closed wires after
        # fusing adjacent edges
        shapes =  []
        for wire in closed_boundaries:
            face = make_face(fuse_edges_in_wire(wire))
            if get_shape_type(face) == ShapeType.COMPOUND:
                shapes.extend(extract_sub_shapes(face, ShapeType.FACE))
            else:
                shapes.append(face)
        # Fuse all the faces into a single shape
        shape = make_fuse(shapes)
        # Return the edges of the fused shape
        return extract_sub_shapes(shape, ShapeType.EDGE)

    # Suppress vertices internal to the edges of the wire
    borders_wire = fuse_edges_in_wire(closed_boundaries[0])
    # Initialize a list of list each containing the edges lying on the same
    # border
    groups_of_collinear_edges: List[List[Any]] = []
    # Loop through all the sorted edges
    for edge in extract_sorted_sub_shapes(borders_wire, ShapeType.EDGE):
        if str(get_kind_of_shape(edge)[0]) != 'SEGMENT':
            groups_of_collinear_edges.append([edge])
            continue
        # Check if the current edge belongs to any group of collinear edges
        for group in groups_of_collinear_edges:
            if is_collinear(edge, group):
                group.append(edge)
                break
        else:
            groups_of_collinear_edges.append([edge])
    # Build the borders edges as single edge objects
    border_edges = []
    for edge in groups_of_collinear_edges:
        if len(edge) > 1 and str(get_kind_of_shape(edge[0])[0]) == 'SEGMENT':
            # Get start-end vertices of the edges on the border
            v1 = extract_sorted_sub_shapes(edge[0], ShapeType.VERTEX)[0]
            v2 = extract_sorted_sub_shapes(edge[-1], ShapeType.VERTEX)[1]
            border_edges.append(make_edge(v1, v2))
        else:
            border_edges.append(edge[0])

    return border_edges


def build_contiguous_edges(vertices: List[Any]) -> List[Any]:
    """
    Function that builds a list of contiguous edge objects each sharing
    one of the vertices to form a closed path.

    Parameters
    ----------
    vertices : List[Any]
        List of vertex objects representing the start-end points of the
        contiguous edges.

    Returns
    -------
    List[Any]
        The list of contiguous edge objects.

    Raises
    ------
    RuntimeError
        If less than 3 vertices are provided.
    RuntimeError
        If two consecutive vertices in the list are the same, as no edge
        could be built from them.
    """
    # Check if the number of vertices is enough to build a closed path
    if len(vertices) < 3:
        raise RuntimeError("A closed path requires at least three points.")
    edges = []
    # Check the given vertices are contiguous
    for i in range(len(vertices)):
        if (math.isclose(get_min_distance(vertices[i],
                                          vertices[(i+1) % len(vertices)]),
                         0.0)):
            raise RuntimeError("It is not possible to build an edge from two "
                               "consecutive equal vertices.")
        # Build and return a list of the contiguous edges
        edges.append(make_edge(vertices[i], vertices[(i+1) % len(vertices)]))
    return edges


def build_z_axis_from_vertex(vertex: Any) -> Any:
    """
    Function that builds an axis originating from the given GEOM vertex and
    being perpendicular to the XY plane.

    Parameters
    ----------
    vertex : Any
        The GEOM vertex being the origin point of the Z-axis.

    Returns
    -------
    Any
        The GEOM edge being the Z-axis originating from the given vertex.
    """
    # Get the vertex coordinates
    v_x, v_y, _ = get_point_coordinates(vertex)
    # Return the vector
    return make_vector_from_points(vertex, make_vertex((v_x, v_y, 1.0)))


def check_shape_expected_types(shape: Any,
                               expected_types: List[ShapeType]) -> None:
    """
    Function that checks if the type of the given shape matches any of
    the expected types.

    Parameters
    ----------
    shape : Any
        The shape object whose type is to be checked.
    expected_types : List[ShapeType]
        A list of expected shape types.

    Raises
    ------
    RuntimeError
        If the type of the shape is not in the list of expected types.
    """
    type = get_shape_type(shape)
    if type not in expected_types:
        raise RuntimeError(
            f"The GEOM object has a non-valid '{type}' type: any of the "
            f"'{expected_types}' objects are expected.")


def check_type_geo_consistency(
        type_geo: LatticeGeometryType,
        cell_type: CellType,
        symmetry_type: SymmetryType
    ) -> None:
    """
    Function that checks if the given type of geometry is valid for the
    indicated type of cell and the type of symmetry.

    Parameters
    ----------
    type_geo : LatticeGeometryType
        The type of geometry of the lattice.
    cell_type : CellType
        The type of cell.
    symmetry_type : SymmetryType
        The type of symmetry.

    Raises
    ------
    RuntimeError
        If the given lattice type of geometry does not match with the
        indicated cell and symmetry types.
    """
    try:
        # Get the list of types of geometry available for the lattice
        types_geo = CELL_VS_SYMM_VS_TYP_GEO[cell_type][symmetry_type]
        if type_geo not in types_geo:
            raise KeyError
    except KeyError:
        raise RuntimeError(
            f"The given type of geometry '{type_geo}' is not compatible "
            f"with the indicated type of cell (i.e. '{cell_type}') and the "
            f"applied symmetry type '{symmetry_type}'. "
            f"Expected values are {types_geo}.")


def compute_point_by_reference(
        point: Any,
        ref_point: Any,
        new_ref_coords: Tuple[float, float, float]
    ) -> Tuple[float, float, float]:
    """
    Function that calculates the new coordinates of the given vertex object
    so that it keeps the same relative vector (i.e., same distance and
    direction) from a reference point which is moved in another position.
    The coordinates of the reference point before and after its translation
    are provided.

    Parameters
    ----------
    point : Any
        The vertex object whose new coordinates must be evaluated.
    ref_point : Any
        The reference vertex object the point is relative to.
    new_ref_coords : Tuple[float, float, float]
        The XYZ coordinates of the reference point after translation.

    Returns
    -------
    Tuple[float, float, float]
        The new coordinates of the given point so that it preserves its
        relative position with respect to its translated reference point.
    """
    # Get the point position wrt the previous point position
    relative_vctr = tuple(coord2 - coord1 for coord1, coord2 in zip(
        get_point_coordinates(ref_point),
        get_point_coordinates(point)))
    # Return the updated shape position relative to the new point position
    return tuple(coord1 + delta for coord1, delta in zip(
        new_ref_coords,
        relative_vctr))


def generate_unique_random_colors(
        no_colors: int) -> List[Tuple[int, int, int]]:
    """
    Function for generating a specified number of random unique RGB colors.

    Parameters
    ----------
    no_colors : int
        The number of RGB colors to generate.

    Returns
    -------
    List[Tuple[int, int, int]]
        A list of tuples, each providing the 3 integer values identifying an
        RGB color.

    Raises
    ------
    RuntimeError
        If requesting more colors than the ones available.
    """
    # Raise an exception if requesting more than the available number of RGB
    # codes
    no_available_colors = len(RGB_COLORS)
    if no_colors > no_available_colors:
        raise RuntimeError(
            f"The requested number ({no_colors}) of RGB color codes exceeds "
            f"the number of available ones ({no_available_colors}).")
    # Declare a seed for the random number generation, so to produce the
    # same set of colors
    random.seed(50)
    # Return 'no_colors'-number of unique colors
    return random.sample(RGB_COLORS, no_colors)


def get_angle_between_points(
        point1: Tuple[float, float, float],
        point2: Tuple[float, float, float]) -> float:
    """
    Function that, given two points, calculates the angle between the line
    connecting the two points and the X-axis.

    Parameters
    ----------
    point1 : Tuple[float, float, float]
        First point.
    point2: Tuple[float, float, float]
        Second point.

    Returns
    -------
    float
        The angle, in radians, between the line connecting the two points
        and the X-axis.
    """
    return math.atan2(point1[1] - point2[1], point1[0] - point2[0])


def get_id_from_name(name: str) -> int:
    """
    Function that extracts the index of the shape whose name is provided.
    The shape's name must be defined as ``<name>_<id>``.

    Parameters
    ----------
    name : str
        The name of the shape in the SALOME viewer.

    Raises
    ------
    RuntimeError
        If no integer index can be extracted from the given name.

    Returns
    -------
    int
        An integer being the global index associated to the shape whose name
        is given as input.
    """
    try:
        return int(name.split('_')[1])
    except:
        raise RuntimeError("No index could be retrieved for the given "
                           f"shape's name '{name}'.")


def get_id_from_shape(shape: Any) -> int:
    """
    Function that extracts the index of the given shape from its name.
    The shape's name must have been previously assigned as ``<name>_<id>``.

    Parameters
    ----------
    shape : Any
        The generic GEOM shape object to get its ID, if any.

    Raises
    ------
    RuntimeError
        If no name has been assigned to the shape, or no integer index
        can be extracted from the shape's name.

    Returns
    -------
    int
        An integer being the global index associated to the given shape.
    """
    # Get the number of the edge directly from the name attribute of the
    # corresponding GEOM edge object
    name = get_shape_name(shape)
    if not name:
        raise RuntimeError("No name has been assigned to the shape.")
    return get_id_from_name(name)


def is_collinear(edge: Any, collinear_edges: List[Any]) -> bool:
    """
    Function that determines whether the given edge is collinear with a group
    of collinear edges.
    Collinearity is determined by checking if the angle between the edge and
    any edge in the group is close to ``0`` or ``pi``, and if the minimum
    distance between the edge and the infinite axis defined by the group is
    close to zero.

    Parameters
    ----------
    edge : Any
        The edge to check for collinearity.
    group_of_collinear_edges : List[Any]
        A list of edges that are collinear.

    Returns
    -------
    bool
        True if the edge is collinear with the group, False otherwise.
    """
    for collinear_edge in collinear_edges:
        angle = get_angle_between_shapes(edge, collinear_edge)
        if (
            math.isclose(angle, 0.0, abs_tol=1e-5) or
            math.isclose(angle, math.pi, abs_tol=1e-5)
        ):
            distance = get_min_distance(
                edge, make_infinite_axis(collinear_edges[0]))
            if math.isclose(distance, 0.0, abs_tol=1e-5):
                return True
    return False


def make_infinite_axis(edge, length=1e6) -> Any:
    """
    Function that creates an infinite-like axis along the direction of the
    given edge.

    It computes the direction vector of the provided edge and extends it
    in both directions by a specified length, resulting in a much longer edge
    that simulates an infinite axis.

    Parameters
    ----------
    edge : Any
        The edge object from which to derive the axis direction.
    length : float, optional
        The distance to extend the axis in both directions from the edge's
        endpoints. Default is 1e6.

    Returns
    -------
    Any
        An edge object resulting from extending the given edge along its
        direction vector.

    Notes
    -----
    The function assumes that the edge lies in the XY plane (Z=0).
    """
    # Get start and end points of the edge
    v1, v2 = extract_sorted_sub_shapes(edge, ShapeType.VERTEX)
    p1 = get_point_coordinates(v1)
    p2 = get_point_coordinates(v2)
    # Compute direction vector
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dz = 0.0
    # Normalize direction
    norm = (dx**2 + dy**2 + dz**2)**0.5
    dx /= norm
    dy /= norm
    dz /= norm
    # Create extended points
    p_start = (p1[0] - dx * length, p1[1] - dy * length, p1[2] - dz * length)
    p_end = (p2[0] + dx * length, p2[1] + dy * length, p2[2] + dz * length)
    # Build and return the extended edge
    return make_edge(make_vertex(p_start), make_vertex(p_end))


def retrieve_selected_object(error_msg: str) -> Any:
    """
    Function that retrieves the geometrical object currently selected in
    the SALOME study. If more than one or none is selected, an exception
    with the given message is raised.

    Parameters
    ----------
    error_msg : str
        The error message to display when the incorrect number of shapes
        is selected.

    Returns
    -------
    Any
        The GEOM object currently selected in the SALOME study.
    """
    if not is_gui_available():
        raise RuntimeError(
            "This function can only be called from the SALOME GUI")
    shape = get_selected_object()
    if not shape:
        raise RuntimeError(error_msg)
    return shape


def sort_shapes_from_vertex(
        shapes: List[Any], vertex: Any, reverse: bool = False) -> List[Any]:
    """
    Function that sorts the given shapes based on their minimum distance from
    a given vertex.

    Parameters
    ----------
    shapes : List[Any]
        A list of shape objects to sort.
    vertex : Any
        The reference vertex used to compute distances.
    reverse : bool = False
        The sorting order. It defaults to ``False``

    Returns
    -------
    List[Any]
        A new list of shapes sorted (by default in ascending order) according
        to their minimum distance from the given vertex.
    """
    return sorted(
        shapes,
        key=lambda shape: get_min_distance(vertex, shape),
        reverse=reverse
    )


def sort_vertices_radially(vertices) -> List[Any]:
    """
    Function that sorts a list of vertices in radial order around their
    centroid.
    The function computes the centroid of the given vertices, then sorts them
    based on the angle each vertex makes with respect to the centroid.

    Parameters
    ----------
    vertices : list
        A list of vertex objects to be sorted.

    Returns
    -------
    List[Any]
        A list of vertex objects sorted in radial order around the centroid.
    """
    # Compute the centroid X-Y coordinates
    coords = [get_point_coordinates(v) for v in vertices]
    cx = sum(p[0] for p in coords) / len(coords)
    cy = sum(p[1] for p in coords) / len(coords)
    # Sort vertices according to the angle wrt to the centroid
    return sorted(
        vertices,
        key=lambda v: get_angle_between_points(get_point_coordinates(v),
                                               (cx, cy, 0.0))
    )


def translate_wrt_reference(
        shape: Any,
        old_ref_point: Any,
        new_ref_coords: Tuple[float, float, float]) -> Any:
    """
    Function that translates a geometric shape so that it keeps its relative
    distance from a reference point, whose previous position is provided as
    a vertex object, while its new position by means of its XYZ coordinates.

    Parameters
    ----------
    shape : Any
        The geometric shape to be translated.
    old_ref_point : Any
        A vertex object identifying the previous reference point of the shape.
    new_ref_coords : Tuple[float, float, float]
        The coordinates of the reference point for the shape after its
        translation.

    Returns
    -------
    Any
        The given shape translated so that it keeps its relative distance from
        the new reference point.
    """
    # Build a vertex identifying the shape CDG
    cdg = make_cdg(shape)
    # Get the shape's CDG coordinates so that the relative distance from
    # the new reference point is kept the same
    shape_to_center = compute_point_by_reference(
        cdg, old_ref_point, new_ref_coords)
    return make_translation(
        shape, make_vector_from_points(cdg, make_vertex(shape_to_center)))
