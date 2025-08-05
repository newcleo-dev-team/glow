"""
Module containing utility functions to support the construction and analysis
of the geometry layouts.
"""
import math
import random
from typing import Any, List, Tuple

from glow.interface.geom_interface import ShapeType, \
    extract_sorted_sub_shapes, extract_sub_shapes, fuse_edges_in_wire, \
    get_basic_properties, get_closed_free_boundary, get_min_distance, \
    get_point_coordinates, get_selected_object, get_shape_name, \
    get_shape_type, is_gui_available, make_cdg, make_cut, make_edge, \
    make_face, make_fuse, make_translation, make_vector_from_points, \
    make_vertex


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
    closed_boundaries = get_closed_free_boundary(cmpd)
    # Handle the case where more than a closed wire is extracted from
    # the compound
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
