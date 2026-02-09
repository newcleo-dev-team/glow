"""
Module containing enumeration classes for expressing different types (e.g.
BCs, geometry, etc.) used throughout the source code.
"""
from enum import Enum
from typing import Dict, List, Tuple


class BoundaryType(Enum):
    """
    Enumeration for defining the types of the boundary conditions applied to
    the borders of the geometry layout.
    """
    VOID:int = 0
    """Indicating a vacuum + albedo BC."""
    REFL:int = 1
    """Indicating a specular reflection BC."""
    TRANSLATION: int = 2
    """Indicating a translation BC."""
    ROTATION: int = 3
    """Indicating a rotation BC."""
    AXIAL_SYMMETRY: int = 4
    """Indicating an axial symmetry BC."""
    CENTRAL_SYMMETRY: int = 5
    """Indicating a central symmetry BC."""


class EdgeType(Enum):
    """
    Enumeration for identifying each type of edge object with an index.
    """
    SEGMENT: int = 1
    """Identifying a segment-type edge."""
    CIRCLE: int = 2
    """ Identifying a circle-type edge."""
    ARC_CIRCLE: int = 3
    """Identifying an arc of circle-type edge."""


class GeometryType(Enum):
    """
    Enumeration for identifying the type of geometry layout.
    Available values are the technological geometry, which describes the
    layout in terms of the regions identified by the different materials, and
    the refined geometry.
    """
    TECHNOLOGICAL: int = 0
    """Identifying the technological geometry of the layout."""
    SECTORIZED: int = 1
    """Identifying the sectorized geometry of the layout."""


class LayoutGeometryType(Enum):
    """
    Enumeration for defining the geometry types of the layout as requested by
    the ``SALT:`` module of DRAGON5.
    Values higher than 2, i.e. ``ROTATION``, refer to geometry layouts for
    which a cycling tracking needs to be applied (`TSPC` type).
    The values from ``0`` to ``2`` included, refer to geometry layouts for
    which a uniform tracking needs to be applied (`TISO` type).
    """
    ISOTROPIC: int = 0
    """Generic geometry with vacuum or isotropic reflection."""
    SYMMETRIES_TWO: int = 1
    """Generic geometry with symmetries of two axis of angle pi/n, n>0."""
    ROTATION: int = 2
    """Generic geometry with rotation of angle 2*pi/n, n>1."""
    RECTANGLE_TRAN: int = 5
    """Rectangular geometry with translation on all sides."""
    RECTANGLE_SYM: int = 6
    """Rectangular geometry with symmetry on all sides."""
    RECTANGLE_EIGHT: int = 7
    """1/8 rectangular assembly with symmetries on all sides."""
    SA60: int = 8
    """Isosceles triangle geometry with symmetries on all sides."""
    HEXAGON_TRAN: int = 9
    """Hexagonal geometry with translations on all sides."""
    RA60: int = 10
    """Isosceles triangle geometry with rotation and translation."""
    R120: int = 11
    """Lozenge geometry with rotation and translation."""
    S30: int = 12
    """Triangle geometry identifying a symmetry of 1/12 of the layout."""


class SymmetryType(Enum):
    """
    Enumeration for defining the symmetry types of the layout.
    """
    FULL: int = 0
    """Identifying a complete layout."""
    HALF: int = 2
    """Identifying an half of the layout."""
    THIRD: int = 3
    """Identifying a third of the layout."""
    QUARTER: int = 4
    """Identifying a quarter of the layout."""
    SIXTH: int = 6
    """Identifying a sixth of the layout."""
    EIGHTH: int = 8
    """Identifying an eighth of the layout."""
    TWELFTH: int = 12
    """Identifying an twelfth of the layout."""
    DIAG: int = 13
    """Identifying an half of the layout along its diagonal."""


class LayoutType(Enum):
    """
    Enumeration for defining the geometric types of the layout.
    """
    RECT: int = 0
    """Identifying a cartesian (i.e. rectangular) layout."""
    HEX: int = 1
    """Identifying a hexagonal layout."""
    GENERIC: int = 2
    """Identifying a generic layout."""


class PropertyType(Enum):
    """
    Enumeration for defining the property types that can be associated
    to each region of the technological geometry layout.
    """
    MATERIAL: int = 0
    """Identifying the material property type."""


# Dictionary associating for each type of layout, the valid combinations of
# types of symmetry and layout types of geometry
LAYOUT_VS_SYMM_VS_TYP_GEO : Dict[
    LayoutType, Dict[SymmetryType, List[LayoutGeometryType]]] = {
    LayoutType.HEX: {
        SymmetryType.FULL: [
            LayoutGeometryType.ISOTROPIC, LayoutGeometryType.HEXAGON_TRAN],
        SymmetryType.THIRD: [
            LayoutGeometryType.ROTATION, LayoutGeometryType.R120],
        SymmetryType.SIXTH: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.ROTATION,
            LayoutGeometryType.SA60,
            LayoutGeometryType.RA60],
        SymmetryType.TWELFTH: [
            LayoutGeometryType.SYMMETRIES_TWO, LayoutGeometryType.S30],
    },
    LayoutType.RECT: {
        SymmetryType.FULL: [
            LayoutGeometryType.ISOTROPIC,
            LayoutGeometryType.RECTANGLE_TRAN,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.HALF: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.DIAG: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.QUARTER: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.EIGHTH: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_EIGHT],
    },
    LayoutType.GENERIC: {
        SymmetryType.FULL: [
            LayoutGeometryType.ISOTROPIC,
            LayoutGeometryType.HEXAGON_TRAN,
            LayoutGeometryType.RECTANGLE_TRAN,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.HALF: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.DIAG: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.QUARTER: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_SYM],
        SymmetryType.EIGHTH: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.RECTANGLE_EIGHT],
        SymmetryType.THIRD: [
            LayoutGeometryType.ROTATION, LayoutGeometryType.R120],
        SymmetryType.SIXTH: [
            LayoutGeometryType.SYMMETRIES_TWO,
            LayoutGeometryType.ROTATION,
            LayoutGeometryType.SA60,
            LayoutGeometryType.RA60],
        SymmetryType.TWELFTH: [
            LayoutGeometryType.SYMMETRIES_TWO, LayoutGeometryType.S30],
    },
}


# Dictionary of edges' type name VS a tuple containing the corresponding
# attribute of the 'EdgeType' enumeration and a descriptive string
EDGE_NAME_VS_TYPE : Dict[str, Tuple[EdgeType, str]] = {
    "SEGMENT"     : (EdgeType.SEGMENT, "line segment"),
    "CIRCLE"      : (EdgeType.CIRCLE, "circle"),
    "ARC_CIRCLE"  : (EdgeType.ARC_CIRCLE, "circular arc")
}
