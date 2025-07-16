"""
Module containing utility functions to support the construction and analysis
of the geometry layout.
"""
import math
from typing import Any, List, Tuple

from glow.interface.geom_interface import ShapeType, \
  extract_sorted_sub_shapes, extract_sub_shapes, fuse_edges_in_wire, \
  get_basic_properties, get_closed_free_boundary, get_inertia_matrix, \
  get_point_coordinates, get_selected_object, get_shape_type, make_cdg, make_cut, make_face, \
  make_fuse, make_translation, make_vector_from_points, make_vertex


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
    bool:
        `True` if the two shapes are considered the same (i.e., they have
        the same geometric properties and the cut result does not contain
        any sub-shapes of the specified type), `False` otherwise.

    Raises
    ------
    RuntimeError: If the two shapes are not of the same type.
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
    Function that extracts the external boundary of the given compound object
    and returns them as a list of edge objects.

    Parameters
    ----------
    cmpd : Any
        The compound object whose borders to extract

    Returns
    -------
    A list of edge objects representing the given compound external boundary.
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


def update_relative_pos(point: Any,
                        rel_point_pre: Any,
                        new_pos: Tuple[float, float, float]) -> Tuple:
    """
    Function that calculates the position of the given point which is
    relative to another one, which has to be moved into the given XYZ
    position.

    Parameters
    ----------
    point         : Any
                    The point whose new coordinates must be evaluated
    rel_point_pre : Any
                    The point the first parameter is relative to
    new_pos       : Tuple[float, float, float]
                    The XYZ coordinates of the moved point (2nd parameter)

    Returns
    -------
    A tuple with the XYZ coordinates of the given point so that it still
    keeps its relative position with the moved point.
    """
    # Get the point position wrt the previous point position
    relative_pos = tuple(coord1 - coord2 for coord1, coord2 in zip(
        get_point_coordinates(rel_point_pre),
        get_point_coordinates(point)))
    # Return the updated shape position relative to the new point position
    return tuple(coord1 - coord2 for coord1, coord2 in zip(
        new_pos,
        relative_pos))


def compare_compounds_rotation(
        cmpd1: Any, cmpd2: Any, valid_angles: List[float]) -> float:
    """
    Function for comparing the rotation angle of two given compounds. If they
    have a different angle, their difference is calculated and checked against
    the allowed values for the lattice.

    Parameters
    ----------
    cmpd1 : Any
        One of the geometrical compounds to compare
    cmpd2 : Any
        One of the geometrical compounds to compare
    valid_angles: List[float]
        List of admitted angles for the compound

    Returns
    -------
    The difference (in degrees) of the rotation angles of the two compounds,
    if they have different rotations; otherwise, the common rotation angle.
    """
    # Get the principal axis angle for both compounds
    angle1 = get_principal_axis_angle(cmpd1)
    angle2 = get_principal_axis_angle(cmpd2)
    # Check if the rotation angle is the same
    if angle1 != angle2:
        print("The rotation angle has changed")
        # Calculate the rotation angle difference
        d_angle = angle2 - angle1
        # Check if the rotation angle is one of the allowed ones
        if not abs(d_angle) in valid_angles:
            raise ValueError("The lattice compound has been rotated by "
                             f"an invalid angle of {d_angle}. The allowed "
                             f"ones are {valid_angles}°.")
        # Return the difference of the compounds rotation angles
        return d_angle
    # Return the common rotation angle of the compounds
    return angle1


def get_principal_axis_angle(shape: Any) -> float:
    """
    Function for getting the principal axis angle of the given geometrical
    shape.
    The rotation angle of the shape in the XY plane is given from the atan
    of the XY components of this axis.

    Parameters
    ----------
    shape : Any
            The geometrical shape in the XY plane whose rotation angle has
            to be determined

    Returns
    -------
    The rotation angle (in degrees) of the shape in the XY plane.
    """
    # Get the shape inertia matrix
    inertia_matrix = get_inertia_matrix(shape)
    # Extract the principal X-Y axes coordinates (2D shape)
    px = inertia_matrix[-3]
    py = inertia_matrix[-2]
    # Return the rotation angle in degrees
    return math.degrees(math.atan2(py, px))


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
        The geometric shape to be translated
    old_ref_point : Any
        A vertex object identifying the previous reference point of the shape
    new_ref_coords : Tuple[float, float, float]
        The coordinates of the reference point for the shape after its
        translation

    Returns
    -------
    The given shape translated so that it keeps its relative distance from
    the new reference point.
    """
    # Build a vertex identifying the shape CDG
    cdg = make_cdg(shape)
    # Get the shape's CDG coordinates so that the relative distance from
    # the new reference point is kept the same
    shape_to_center = update_relative_pos(cdg, old_ref_point, new_ref_coords)
    return make_translation(
        shape, make_vector_from_points(cdg, make_vertex(shape_to_center)))

def get_id_from_name(name: str) -> int:
    """
    Function that extracts the index of the shape whose name is provided.
    The shape's name must be defined as `<name>_<id>`.

    Parameters
    ----------
    name : str
        The name of the shape in the SALOME viewer

    Raises
    ------
    RuntimeError
        If no integer index can be extracted from the given name

    Returns
    -------
    An integer being the global index associated to the shape whose name
    is given as input.
    """
    try:
        return int(name.split('_')[1])
    except:
        raise RuntimeError("No index could be retrieved for the given "
                           f"shape's name '{name}'.")

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
    """
    shape = get_selected_object()
    if not shape:
        raise RuntimeError(error_msg)
    return shape
