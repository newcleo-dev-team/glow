"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.lattice` module have a valid implementation.
"""
import math
import unittest

from typing import List

from glow.generator.support import *
from glow.geometry_layouts.cells import Cell, HexCell, RectCell
from glow.geometry_layouts.lattices import Lattice
from glow.geometry_layouts.utility import are_same_shapes
from glow.interface.geom_interface import *


class TestLattice(unittest.TestCase):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Lattice` class.

    This test suite provides common setup and a set of tests to ensure that
    implementations of `Lattice` correctly handle the operations performed
    on a lattice, whose type is either cartesian or hexagonal.

    Attributes
    ----------
    name : str
        The name of the lattice when displayed in the SALOME viewer.
    hex_cells : List[HexCell]
        The list of `HexCell` objects the hexagonal lattice to test is made
        of.
    rect_cells : List[RectCell]
        The list of `RectCell` objects the cartesian lattice to test is made
        of.
    o : Any
        A vertex object positioned at the origin of the XYZ space.
    lattice : Lattice
        The `Lattice` object to test, made either by cartesian or hexagonal
        cells.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `Lattice` class.
        It initializes the attributes common to all the tests.
        """
        self.name: str = "Lattice"
        self.lattice: Lattice | None = Lattice(name=self.name)
        self.o = make_vertex((0.0, 0.0, 0.0))
        hex_cell = HexCell(name="Hexagonal Cell")
        rect_cell = RectCell(name="Cartesian Cell")
        self.hex_cells: List[HexCell] = self.__set_up_hex_cells(hex_cell)
        self.rect_cells: List[RectCell] = self.__set_up_rect_cells(rect_cell)

    def test_init(self) -> None:
        """
        Method that tests the correct initialization of the `Lattice` class.
        """
        # No cells case
        self.__assess_lattice_init([], [], 0)
        # No cells case with center different from XYZ origin
        center = (1.0, 1.0, 0.0)
        self.lattice = Lattice(center=center)
        self.__assess_lattice_init([], [], 0, center)
        # Hexagonal cells case
        self.lattice = Lattice(self.hex_cells)
        self.__assess_lattice_init(self.hex_cells, [], 1)
        # Cartesian case (odd number of cells)
        self.rect_cells = self.__set_up_rect_cells(self.rect_cells[0])
        self.lattice = Lattice(self.rect_cells)
        self.__assess_lattice_init(self.rect_cells, [], 1)
        # Cartesian case (even number of cells)
        self.rect_cells = self.__set_up_rect_cells(self.rect_cells[0], True)
        self.lattice = Lattice(self.rect_cells)
        self.__assess_lattice_init(self.rect_cells, [], 1)

    def __assess_lattice_init(
            self,
            cells: List[Cell],
            box_layers: List[float],
            no_rings: int,
            center: Tuple[float] | None = None) -> None:
        """
        Method that assesses the correct initialization of a `Lattice`
        instance depending on whether any cells, any box layers, or center
        have been provided.

        Parameters
        ----------
        cells : List[Cell]
            List of cells in the lattice.
        box_layers : List[float]
            List of the thicknesses of the box enclosing the lattice.
        no_rings : int
            The number of rings of cells in the lattice.
        center : Tuple[float] | None = None
            The center of the lattice, if any.
        """
        if center is None:
            center = (0.0, 0.0, 0.0)
        self.assertTrue(
            are_same_shapes(
                self.lattice.lattice_center,
                make_vertex(center),
                ShapeType.VERTEX)
        )
        if cells:
            lx = cells[0].figure.lx
            ly = cells[0].figure.ly
            self.assertTrue(self.lattice.layers and all(
                layer for layer in self.lattice.layers))
            self.assertEqual(self.lattice.cells_type, cells[0].cell_type)
            self.assertEqual(
                self.lattice.cells_rot, math.degrees(cells[0].rotation))
            self.assertAlmostEqual(
                self.lattice.distance,
                max(math.dist(
                    get_point_coordinates(self.lattice.lattice_center),
                    get_point_coordinates(c.figure.o)) for c in cells)
            )
        else:
            lx = 0.0
            ly = 0.0
            self.assertFalse(self.lattice.layers and all(
                layer for layer in self.lattice.layers))
            self.assertIsNone(self.lattice.cells_type)
            self.assertIsNone(self.lattice.cells_rot)
            self.assertEqual(self.lattice.distance, 0.0)

        self.assertEqual(len(self.lattice.lattice_cells), len(cells))
        self.assertEqual(len(self.lattice.layers[0]), len(cells))
        self.assertEqual(self.lattice.name, self.name)
        self.assertIsNone(self.lattice.lattice_entry_id)
        self.assertEqual(self.lattice.rings_no, no_rings)
        self.assertEqual(self.lattice.type_geo, LatticeGeometryType.ISOTROPIC)
        self.assertEqual(self.lattice.symmetry_type, SymmetryType.FULL)
        self.assertEqual(self.lattice.lx, lx)
        self.assertEqual(self.lattice.ly, ly)
        self.assertEqual(len(self.lattice.box_layers), len(box_layers))
        self.assertIsNone(self.lattice.lattice_symm)
        self.assertEqual(len(self.lattice.regions), 0)
        self.assertEqual(self.lattice.displayed_geom,
                         GeometryType.TECHNOLOGICAL)
        self.assertFalse(self.lattice.is_update_needed)

        if box_layers:
            self.assertIsNotNone(self.lattice.lattice_box)
        else:
            self.assertIsNone(self.lattice.lattice_box)


    def __set_up_hex_cells(self, hex_cell: HexCell) -> List[HexCell]:
        """
        Method that builds a list of hexagonal cells with a central cell
        surrounded by six other cells.

        Parameters
        ----------
        hex_cell : HexCell
            The `HexCell` object representing an hexagonal cell.

        Returns
        -------
        List[HexCell]
            A list of hexagonal cells with a central one surrounded by six
            cells.
        """
        dx = 2*hex_cell.apothem*math.sin(math.radians(60))
        dy = 2*hex_cell.apothem*math.cos(math.radians(60))
        hex_cell.rotate(90.0)
        return [
            hex_cell,
            hex_cell.translate((dx, dy, 0)),
            hex_cell.translate((-dx, dy, 0)),
            hex_cell.translate((-dx, -dy, 0)),
            hex_cell.translate((dx, -dy, 0)),
            hex_cell.translate((0, 2*dy, 0)),
            hex_cell.translate((0, -2*dy, 0))
        ]

    def __set_up_rect_cells(self,
                            rect_cell: RectCell,
                            is_even: bool = False) -> List[RectCell]:
        """
        Method that builds a list of cartesian cells. Depending on the
        boolean flag `is_even`, the resulting list is made by a central
        cell surrounded by eight other cells, if `False`, or by cells
        without any central one replicating a pattern with an even number
        of cells.

        Parameters
        ----------
        rect_cell : RectCell
            The `RectCell` object representing a cartesian cell.
        is_even : bool
            Boolean flag indicating the kind of pattern of cells (either with
            an odd or even number of cells).

        Returns
        -------
        List[RectCell]
            A list of cartesian cells with a specific pattern.
        """
        dx = rect_cell.width
        dy = rect_cell.height
        if is_even:
            dx /= 2
            dy /= 2
            return [
                rect_cell.translate((dx, dy, 0)),
                rect_cell.translate((-dx, dy, 0)),
                rect_cell.translate((-dx, -dy, 0)),
                rect_cell.translate((dx, -dy, 0)),
                rect_cell.translate((2*dx, 0, 0)),
                rect_cell.translate((2*dx, dy, 0)),
                rect_cell.translate((2*dx, 2*dy, 0)),
                rect_cell.translate((dx, 2*dy, 0)),
                rect_cell.translate((-dx, 2*dy, 0)),
                rect_cell.translate((-2*dx, 2*dy, 0)),
                rect_cell.translate((-2*dx, dy, 0)),
                rect_cell.translate((-2*dx, 0, 0)),
                rect_cell.translate((-2*dx, -dy, 0)),
                rect_cell.translate((-2*dx, -2*dy, 0)),
                rect_cell.translate((-dx, -2*dy, 0)),
                rect_cell.translate((dx, -2*dy, 0)),
                rect_cell.translate((2*dx, -2*dy, 0)),
                rect_cell.translate((2*dx, -dy, 0)),
            ]
        return [
            rect_cell,
            rect_cell.translate((dx, 0, 0)),
            rect_cell.translate((dx, dy, 0)),
            rect_cell.translate((0, dy, 0)),
            rect_cell.translate((-dx, dy, 0)),
            rect_cell.translate((-dx, 0, 0)),
            rect_cell.translate((-dx, -dy, 0)),
            rect_cell.translate((0, -dy, 0)),
            rect_cell.translate((dx, -dy, 0))
        ]


if __name__ == "__main__":
    unittest.main()