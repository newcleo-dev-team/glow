"""
Module containing the class and functions for handling the construction of the
shape of the symmetry that corresponds to the type of layout (i.e. generic,
Cartesian and hexagonal).
"""
import math

from copy import deepcopy
from dataclasses import dataclass, field
from functools import wraps
from typing import Callable, Tuple

from glow.geometry_layouts.geometries import Rectangle, Surface, \
    build_parallelogram, build_regular_triangle, build_right_triangle, \
    build_right_triangle_from_catheti
from glow.interface.geom_interface import get_bounding_box
from glow.support.types import SymmetryType


@dataclass
class SymmetryDomain:
    """
    Dataclass that represents the geometric domain required to build the
    characteristic shape that is associated to a given symmetry type.

    This class stores the 2D shape identifying the full layout, its centre,
    and automatically computes the X-Y bounding extents of the layout shape.
    """
    origin: Tuple[float, float, float]
    """The XYZ coordinates of the full layout centre."""
    full_layout_shape: Surface
    """The `Surface` instance representing the 2D shape of the full layout."""
    x_bounds: Tuple[float, float] = field(init=False)
    """
    The min-max values of the X-bounding extent of the full layout shape.
    """
    y_bounds: Tuple[float, float] = field(init=False)
    """
    The min-max values of the Y-bounding extent of the full layout shape.
    """

    def __post_init__(self) -> None:
        """
        Method that is automatically run after the dataclass initialization
        for setting the X-Y bounds of the geometry context from the full
        layout shape.
        """
        x_min, x_max, y_min, y_max = get_bounding_box(self.full_layout_shape)
        self.x_bounds = (x_min, x_max)
        self.y_bounds = (y_min, y_max)


# -------------------------------------------------------------------------- #
#                                FUNCTIONS                                   #
# -------------------------------------------------------------------------- #

def build_cartesian_symmetry_shape(
        symmetry: SymmetryType,
        context: SymmetryDomain
    ) -> Surface:
    """
    Function that builds the geometric shape that corresponds to the given
    symmetry type and domain for a Cartesian-type layout.
    In addition to the symmetry types available for all the geometry layouts,
    this function supports those symmetries specific for geometry layouts
    based on a rectangular characteristic shape.

    Supported symmetries are:

    - FULL: returns the characteristic shape of the layout.
    - HALF: builds a rectangle representing half of the layout.
    - QUARTER: builds a rectangle representing one quarter of the layout.
    - EIGHTH: builds a right triangle representing one eighth of the layout.
    - DIAG: builds a right triangle representing the diagonal symmetry of
      the layout.

    Parameters
    ----------
    symmetry : SymmetryType
        The symmetry type for which the characteristic shape is built.
    domain : SymmetryDomain
        Instance providing the domain of the full layout in terms of
        XY-bounding extents, full shape and its centre.

    Returns
    -------
    Surface
        A geometric shape representing the requested symmetry.

    Raises
    ------
    RuntimeError
        If the indicated symmetry type is not among the supported ones for a
        Cartesian-type layout.
    """
    # Get the XY dimensions of the bounding box for the geometry layout
    x_min, x_max = context.x_bounds
    y_min, y_max = context.y_bounds
    # Get the coordinates of the layout's centre
    o_x, o_y, o_z = context.origin
    # Match the shape with the type of symmetry
    match symmetry:
        case SymmetryType.DIAG:
            return build_right_triangle_from_catheti(
                x_max - x_min, y_max - y_min, (x_min, y_min, o_z)
            )
        case SymmetryType.EIGHTH:
            return build_right_triangle_from_catheti(
                x_max - o_x, y_max - o_y, (o_x, o_y, o_z)
            )
        case _:
            return build_common_symmetry_shape(symmetry, context)


def build_common_symmetry_shape(
        symmetry: SymmetryType,
        domain: SymmetryDomain
    ) -> Surface:
    """
    Function that builds the geometric shape that corresponds to the given
    symmetry type and domain for a generic layout.
    This function supports only the symmetry types that are common to all
    layouts, indipendently from their characteristic shape.

    Supported symmetries are:

    - FULL: returns the characteristic shape of the layout.
    - HALF: builds a rectangle representing half of the layout.
    - QUARTER: builds a rectangle representing one quarter of the layout.

    Parameters
    ----------
    symmetry : SymmetryType
        The symmetry type for which the characteristic shape is built.
    domain : SymmetryDomain
        Instance providing the domain of the full layout in terms of
        XY-bounding extents, full shape and its centre.

    Returns
    -------
    Surface
        A geometric shape representing the requested symmetry.

    Raises
    ------
    RuntimeError
        If the indicated symmetry type is not supported.
    """
    # Get the XY dimensions of the bounding box for the geometry layout
    x_min, x_max = domain.x_bounds
    y_min, y_max = domain.y_bounds
    # Get the coordinates of the layout's centre
    o_x, o_y, o_z = domain.origin
    # Match the shape with the type of symmetry
    match symmetry:
        case SymmetryType.FULL:
            return deepcopy(domain.full_layout_shape)
        case SymmetryType.HALF:
            return Rectangle(
                (o_x + (x_max - o_x)/2, o_y, o_z),
                y_max - y_min,
                x_max - o_x
            )
        case SymmetryType.QUARTER:
            return Rectangle(
                (o_x + (x_max - o_x)/2, o_y + (y_max - o_y)/2, o_z),
                (y_max - y_min)/2,
                (x_max - x_min)/2
            )
        case _:
            raise RuntimeError(f"Symmetry '{symmetry}' not supported.")


def build_hex_symmetry_shape(
        symmetry: SymmetryType,
        domain: SymmetryDomain
    ) -> Surface:
    """
    Function that builds the geometric shape that corresponds to the given
    symmetry type and domain for a hexagonal-type layout.
    In addition to the symmetry types available for all the geometry layouts,
    this function supports those symmetries specific for geometry layouts
    based on a hexagonal characteristic shape.

    Supported symmetries are:

    - FULL: returns the characteristic shape of the layout.
    - HALF: builds a rectangle representing half of the layout.
    - QUARTER: builds a rectangle representing one quarter of the layout.
    - THIRD: builds a parallelogram representing one third of the layout.
    - SIXTH: builds a regular triangle representing one sixth of the layout.
    - TWELFTH: builds a right triangle representing one twelfth of the
      layout.

    Parameters
    ----------
    symmetry : SymmetryType
        The symmetry type for which the characteristic shape is built.
    domain : SymmetryDomain
        Instance providing the domain of the full layout in terms of
        XY-bounding extents, full shape and its centre.

    Returns
    -------
    Surface
        A geometric shape representing the requested symmetry.

    Raises
    ------
    RuntimeError
        If the indicated symmetry type is not among the supported ones for a
        hexagonal layout.
    """
    # Get the XY dimensions of the bounding box for the geometry layout
    x_min, x_max = domain.x_bounds
    y_min, y_max = domain.y_bounds
    # Get the coordinates of the layout's centre
    o_x, o_y, o_z = domain.origin
    # Get the characteristic dimension of the hexagon
    side_len = (x_max - x_min) / 2
    apothem = (y_max - y_min) / 2
    # Match the shape with the type of symmetry
    match symmetry:
        case SymmetryType.THIRD:
            return build_parallelogram(
                side_len,
                side_len,
                60,
                (-(o_x + side_len/2), -(o_y + apothem), o_z)
            )
        case SymmetryType.SIXTH:
            return build_regular_triangle(
                side_len,
                (-(o_x + side_len/2), -(o_y + apothem), o_z)
            )
        case SymmetryType.TWELFTH:
            return build_right_triangle(
                side_len, side_len*math.cos(math.pi/6), domain.origin
            )
        case _:
            return build_common_symmetry_shape(symmetry, domain)


def use_symmetry_logic(symmetry_logic: Callable):
    """
    Decorator that replaces the method for building the characteristic shape
    of the symmetry for the ``Fillable`` class and its subclasses.
    Depending on the nature of the layout (i.e. generic, Cartesian, hexagonal
    ), the method is decorated with this function which calls the specific
    function providing the desired behaviour.

    Parameters
    ----------
    symmetry_logic : Callable
        Function implementing the logic for building the shape of the symmetry
        according to the layout type. Admitted functions must receive the
        type of symmetry (as ``SymmetryType`` element) and the instance of
        ``SymmetryDomain`` providing the characteristics of the layout. They
        must return the shape of the symmetry as a ``Surface`` object.

    Returns
    -------
    Callable
        A decorator that overrides the method for building the shape of the
        symmetry in the ``Fillable`` class and in its subclasses.
    """
    def decorator(method):
        # Function wrapping the symmetry shape building method; it ignores
        # the 'self' parameter
        @wraps(method)
        def wrapper(
                _, symmetry: SymmetryType, context: SymmetryDomain
            ):
            return symmetry_logic(symmetry, context)
        return wrapper
    return decorator
