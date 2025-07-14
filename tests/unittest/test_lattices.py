"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.lattice` module have a valid implementation.
"""
import math
import unittest

from typing import List

from glow.generator.support import *
from glow.geometry_layouts.cells import Cell, HexCell, RectCell
from glow.geometry_layouts.geometries import GenericSurface, Hexagon, Rectangle, Surface, build_hexagon
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
    box_layers : List[float]
        List storing the thickness for each layer the lattice's box is made of
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
        self.box_layers: List[float] = [0.075, 0.075]

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

    def test_build_lattice_box(self) -> None:
        """
        Method that tests the correct implementation of the method
        `build_lattice_box` of the `Lattice` class.
        """
        # Instantiate the lattice with the hexagonal cells
        self.lattice = Lattice(self.hex_cells)
        # Build the comparison figure
        box_face = self.__build_hex_box()
        # Build the lattice box
        self.lattice.build_lattice_box(self.box_layers)
        # Verify the correctness of the box surface
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_box.face,
                            box_face,
                            ShapeType.COMPOUND)
        )
        # Instantiate the lattice with the cartesian cells
        self.lattice = Lattice(self.rect_cells)
        # Build the comparison figure
        box_face = self.__build_rect_box()
        # Build the lattice box
        self.lattice.build_lattice_box(self.box_layers)
        # Verify the correctness of the box surface
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_box.face,
                            box_face,
                            ShapeType.COMPOUND)
        )
        # Instantiate the lattice with the cartesian cells (even cells number)
        self.lattice = Lattice(
            self.__set_up_rect_cells(self.rect_cells[0], True))
        # Build the comparison figure
        box_face = self.__build_rect_box()
        # Build the lattice box
        self.lattice.build_lattice_box(self.box_layers)
        # Verify the correctness of the box surface
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_box.face,
                            box_face,
                            ShapeType.COMPOUND)
        )

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

    def __build_hex_box(self) -> Any:
        """
        Method that builds a face object made by several encapsulated
        hexagon surfaces in number equal to the box layers.

        Returns
        -------
        Any
            A face object made by assembling several hexagons together.
        """
        box_surfaces: List[Hexagon] = []
        box_apothem = self.hex_cells[0].edge_length * (
            2 + math.sin(math.pi/6))
        if self.box_layers[0] > 0:
            box_layers = [0.0] + self.box_layers
        for thick in box_layers:
            box_apothem += thick
            box_surfaces.append(
                build_hexagon(
                    box_apothem,
                    get_point_coordinates(self.lattice.lattice_center)))
        # Perform a 'partition' operation to assemble the lattice box
        return make_partition(
            [hex.face for hex in box_surfaces], [], ShapeType.FACE)

    def __build_rect_box(self) -> Any:
        """
        Method that builds a face object made by several encapsulated
        rectangle surfaces in number equal to the box layers.

        Returns
        -------
        Any
            A face object made by assembling several rectangles together.
        """
        # Declare the starting dimensions of the box
        x_min, x_max, y_min, y_max = get_bounding_box(
            make_compound([cell.face for cell in self.lattice.lattice_cells]))
        height = y_max - y_min
        width = x_max - x_min
        box_surfaces: List[Rectangle] = []
        if self.box_layers[0] > 0:
            box_layers = [0.0] + self.box_layers
        # Build a rectangle for each layer
        for thick in box_layers:
            height += 2*thick
            width += 2*thick
            box_surfaces.append(
                Rectangle(get_point_coordinates(self.lattice.lattice_center),
                          height,
                          width))
        # Perform a 'partition' operation to assemble the lattice box
        return make_partition(
            [rect.face for rect in box_surfaces], [], ShapeType.FACE)

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
        dx = hex_cell.apothem
        dy = 3/2*hex_cell.edge_length
        hex_cell.rotate(90.0)
        return [
            hex_cell,
            hex_cell.translate((dx, dy, 0)),
            hex_cell.translate((-dx, dy, 0)),
            hex_cell.translate((-dx, -dy, 0)),
            hex_cell.translate((dx, -dy, 0)),
            hex_cell.translate((2*dx, 0, 0)),
            hex_cell.translate((-2*dx, 0, 0))
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