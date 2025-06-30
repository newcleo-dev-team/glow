"""
Module containing enumeration classes for expressing different types (e.g.
BCs, geometry, etc.) used throughout the code.
Utility functions are herein declared as well.
"""
import random

from enum import Enum
from typing import Dict, List, Tuple


class BoundaryType(Enum):
    """
    Enumeration for defining the lattice boundary condition types.

    Attributes
    ----------
    VOID              : int = 0
                        Indicating a vacuum + albedo BC
    REFL              : int = 1
                        Indicating a specular reflexion BC
    TRANSLATION       : int = 2
                        Indicating a translation BC
    ROTATION          : int = 3
                        Indicating a rotation BC
    AXIAL_SYMMETRY    : int = 4
                        Indicating an axial symmetry BC
    CENTRAL_SYMMETRY  : int = 5
                        Indicating a central symmetry BC
    """
    VOID              : int = 0
    REFL              : int = 1
    TRANSLATION       : int = 2
    ROTATION          : int = 3
    AXIAL_SYMMETRY    : int = 4
    CENTRAL_SYMMETRY  : int = 5


class GeometryType(Enum):
    """
    Enumeration for defining either the lattice or the single cell
    geometry type.

    Attributes
    ----------
    TECHNOLOGICAL : int = 0
        Identifying the cell/lattice technological geometry
    SECTORIZED    : int = 1
        Identifying the cell/lattice sectorized geometry
    """
    TECHNOLOGICAL : int = 0
    SECTORIZED    : int = 1


class LatticeGeometryType(Enum):
    """
    Enumeration for defining the lattice geometry types.

    Attributes
    ----------
    ISOTROPIC       : int = 0
                      Generic geometry with vacuum or isotropic reflexion
    SYMMETRIES_TWO  : int = 1
                      Generic geometry with symmetries of two axis of angle
                      pi/n, n>0
    ROTATION        : int = 2
                      Generic geometry with rotation of angle 2*pi/n, n>1
    RECTANGLE_TRAN  : int = 5
                      Rectangular geometry with translation on all sides
    RECTANGLE_SYM   : int = 6
                      Rectangular geometry with symmetry on all sides
    RECTANGLE_EIGHT : int = 7
                      1/8 assembly with symmetries on all sides
    SA60            : int = 8
                      Isosceles triangle geometry with symmetries on all
                      sides (SA60)
    HEXAGON_TRAN    : int = 9
                      Hexagonal geometry with translations on all sides
    RA60            : int = 10
                      Isosceles triangle geometry with RA60 rotation and
                      translation
    R120            : int = 11
                      Lozenge geometry with R120 rotation and translation
    S30             : int = 12
                      Triangle geometry identifying a symmetry of one twelfth
                      of an assembly
    """
    ISOTROPIC       : int = 0
    SYMMETRIES_TWO  : int = 1
    ROTATION        : int = 2
    RECTANGLE_TRAN  : int = 5
    RECTANGLE_SYM   : int = 6
    RECTANGLE_EIGHT : int = 7
    SA60            : int = 8
    HEXAGON_TRAN    : int = 9
    RA60            : int = 10
    R120            : int = 11
    S30             : int = 12


class SymmetryType(Enum):
    """
    Enumeration for defining the lattice symmetry types.

    Attributes
    ----------
    FULL    : int = 0
              Identifying a complete lattice
    HALF    : int = 2
              Identifying an half of the lattice
    THIRD   : int = 3
              Identifying a third of the lattice
    QUARTER : int = 4
              Identifying a quarter of the lattice
    SIXTH   : int = 6
              Identifying a sixth of the lattice
    EIGHTH  : int = 8
              Identifying an eighth of the lattice
    TWELFTH : int = 12
              Identifying an twelfth of the lattice
    """
    FULL    : int = 0
    HALF    : int = 2
    THIRD   : int = 3
    QUARTER : int = 4
    SIXTH   : int = 6
    EIGHTH  : int = 8
    TWELFTH : int = 12


class CellType(Enum):
    """
    Enumeration for defining the cell types.

    Attributes
    ----------
    RECT    : int = 0
              Identifying a cartesian (i.e. rectangular) cell
    HEX     : int = 1
              Identifying a hexagonal cell
    """
    RECT    : int = 0
    HEX     : int = 1


class PropertyType(Enum):
    """
    Enumeration for defining the property types that can be associated
    to a cell region.

    Attributes
    ----------
    MATERIAL  : int = 0
                Identifying the property for the cell region material
    """
    MATERIAL  : int = 0
    # FIXME add other types


# Dictionary associating the compatible BC type to the lattice type of
# geometry
TYPEGEO_VS_BC : Dict[LatticeGeometryType, List[BoundaryType]] = {
    LatticeGeometryType.ISOTROPIC        : [BoundaryType.VOID],
    LatticeGeometryType.SYMMETRIES_TWO   : [BoundaryType.AXIAL_SYMMETRY],
    LatticeGeometryType.RECTANGLE_TRAN   : [BoundaryType.TRANSLATION],
    LatticeGeometryType.RECTANGLE_SYM    : [BoundaryType.AXIAL_SYMMETRY],
    LatticeGeometryType.RECTANGLE_EIGHT  : [BoundaryType.AXIAL_SYMMETRY],
    LatticeGeometryType.SA60             : [BoundaryType.AXIAL_SYMMETRY],
    LatticeGeometryType.HEXAGON_TRAN     : [BoundaryType.TRANSLATION],
    LatticeGeometryType.RA60             : [BoundaryType.TRANSLATION,
                                            BoundaryType.ROTATION],
    LatticeGeometryType.R120             : [BoundaryType.TRANSLATION,
                                            BoundaryType.ROTATION],
    LatticeGeometryType.S30              : [BoundaryType.AXIAL_SYMMETRY]
}


# Dictionary associating for each type of cells, the valid combinations of
# types of symmetry and lattice types of geometry
CELL_VS_SYMM_VS_TYP_GEO : Dict[
    CellType, Dict[SymmetryType, List[LatticeGeometryType]]] = {
    CellType.HEX : {
        SymmetryType.FULL : [
            LatticeGeometryType.ISOTROPIC, LatticeGeometryType.HEXAGON_TRAN],
        SymmetryType.THIRD : [
            LatticeGeometryType.SYMMETRIES_TWO, LatticeGeometryType.R120],
        SymmetryType.SIXTH : [
            LatticeGeometryType.ISOTROPIC,
            LatticeGeometryType.SA60,
            LatticeGeometryType.RA60],
        SymmetryType.TWELFTH : [
            LatticeGeometryType.ISOTROPIC, LatticeGeometryType.S30],
    },
    CellType.RECT : {
        SymmetryType.FULL : [
            LatticeGeometryType.ISOTROPIC,
            LatticeGeometryType.RECTANGLE_TRAN],
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


RGB_COLORS = [(r, g, b) for r in range(0, 256, 10)
                        for g in range(0, 256, 10)
                        for b in range(0, 256, 10)]

def generate_unique_random_colors(
        no_colors: int) -> List[Tuple[int, int, int]]:
    """
    Function for generating a specified number of random unique RGB colors.

    Parameters
    ----------
    no_colors : int
                The number of RGB colors to generate

    Returns
    -------
    A list of tuples, each providing the 3 integer values identifying an
    RGB color.
    """
    # Declare a seed for the random number generation, so to produce the
    # same set of colors
    random.seed(50)
    # Return 'no_colors'-number of unique colors
    return random.sample(RGB_COLORS, no_colors)
