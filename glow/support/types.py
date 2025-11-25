"""
Module containing enumeration classes for expressing different types (e.g.
BCs, geometry, etc.) used throughout the code.
"""
from enum import Enum
from typing import Dict, List, Tuple


class BoundaryType(Enum):
    """
    Enumeration for defining the types of the lattice's boundary condition
    applied to the borders of its geometry layout.
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
    Enumeration for defining either the lattice or the single cell
    geometry type.
    """
    TECHNOLOGICAL: int = 0
    """Identifying the cell/lattice technological geometry."""
    SECTORIZED: int = 1
    """Identifying the cell/lattice sectorized geometry."""


class LatticeGeometryType(Enum):
    """
    Enumeration for defining the lattice's geometry types. Values higher than
    2, i.e. ``ROTATION``, refer to geometry layouts for which a cycling
    tracking needs to be applied (`TSPC` type).
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
    """Triangle geometry identifying a symmetry of 1/12 of an assembly."""


class SymmetryType(Enum):
    """
    Enumeration for defining the lattice's symmetry types.
    """
    FULL: int = 0
    """Identifying a complete lattice."""
    HALF: int = 2
    """Identifying an half of the lattice."""
    THIRD: int = 3
    """Identifying a third of the lattice."""
    QUARTER: int = 4
    """Identifying a quarter of the lattice."""
    SIXTH: int = 6
    """Identifying a sixth of the lattice."""
    EIGHTH: int = 8
    """Identifying an eighth of the lattice."""
    TWELFTH: int = 12
    """Identifying an twelfth of the lattice."""


class CellType(Enum):
    """
    Enumeration for defining the geometric types of cells.
    """
    RECT: int = 0
    """Identifying a cartesian (i.e. rectangular) cell."""
    HEX: int = 1
    """Identifying a hexagonal cell."""


class PropertyType(Enum):
    """
    Enumeration for defining the property types that can be associated
    to each cell/lattice region.
    """
    MATERIAL: int = 0
    """Identifying the material property type."""
    MACRO: int = 1
    """Identifying the macro region a region belongs to."""


# Dictionary associating for each type of cells, the valid combinations of
# types of symmetry and lattice types of geometry
CELL_VS_SYMM_VS_TYP_GEO : Dict[
    CellType, Dict[SymmetryType, List[LatticeGeometryType]]] = {
    CellType.HEX : {
        SymmetryType.FULL : [
            LatticeGeometryType.ISOTROPIC, LatticeGeometryType.HEXAGON_TRAN],
        SymmetryType.THIRD : [
            LatticeGeometryType.ROTATION, LatticeGeometryType.R120],
        SymmetryType.SIXTH : [
            LatticeGeometryType.SYMMETRIES_TWO,
            LatticeGeometryType.ROTATION,
            LatticeGeometryType.SA60,
            LatticeGeometryType.RA60],
        SymmetryType.TWELFTH : [
            LatticeGeometryType.SYMMETRIES_TWO, LatticeGeometryType.S30],
    },
    CellType.RECT : {
        SymmetryType.FULL : [
            LatticeGeometryType.ISOTROPIC,
            LatticeGeometryType.RECTANGLE_TRAN,
            LatticeGeometryType.RECTANGLE_SYM],
        SymmetryType.HALF : [
            LatticeGeometryType.SYMMETRIES_TWO,
            LatticeGeometryType.RECTANGLE_SYM],
        SymmetryType.QUARTER : [
            LatticeGeometryType.SYMMETRIES_TWO,
            LatticeGeometryType.RECTANGLE_SYM],
        SymmetryType.EIGHTH : [
            LatticeGeometryType.SYMMETRIES_TWO,
            LatticeGeometryType.RECTANGLE_EIGHT],
    }
}


# Dictionary of edges' type name VS a tuple containing the corresponding
# attribute of the 'EdgeType' enumeration and a descriptive string
EDGE_NAME_VS_TYPE : Dict[str, Tuple[EdgeType, str]] = {
    "SEGMENT"     : (EdgeType.SEGMENT, "line segment"),
    "CIRCLE"      : (EdgeType.CIRCLE, "circle"),
    "ARC_CIRCLE"  : (EdgeType.ARC_CIRCLE, "circular arc")
}
