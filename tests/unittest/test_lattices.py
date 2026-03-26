"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.lattices` module have a valid implementation.
"""
import unittest

from math import isclose, pi, sin
from typing import List, Tuple

from glow.geometry_layouts.cells import CartesianCell, Cell, HexCell
from glow.geometry_layouts.geometries import Circle, Hexagon, Rectangle
from glow.geometry_layouts.lattices import CartesianLattice, HexLattice, \
    Lattice, compute_subdivision_points_on_borders, ensure_not_zero, \
    get_cell_at_centres
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_entities import wrap_shape
from glow.interface.geom_interface import ShapeType, get_angle_between_shapes, \
    get_bounding_box, make_compound, make_face, make_translation, \
    make_vector_from_points, make_vertex
from glow.support.types import PropertyType
from glow.support.utility import are_same_shapes, build_compound_borders, \
    build_z_axis_from_vertex
from tests.unittest.support_funcs import set_up_hex_cells, set_up_rect_cells
from tests.unittest.test_fillable_layouts import TestFillable


class TestCartesianLattice(TestFillable):
    """
    Test class for assessing the `CartesianLattice` class implementation.

    This test class contains unit tests for verifying the correct behaviour
    of the `CartesianLattice` class, including initialization and ring
    addition.
    Tests dealing with operations common to all `Fillable` subclasses are
    declared in the `TestFillable` class this class inherits from.

    Attributes
    ----------
    fillable : CartesianLattice
        The Cartesian lattice under test, whose class inherits from
        `Fillable`.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `CartesianLattice`
        class.
        """
        centres = [
            (0.5, 0.5, 0.0),
            (-0.5, 0.5, 0.0),
            (-0.5, -0.5, 0.0),
            (0.5, -0.5, 0.0)
        ]
        cells = [
            CartesianCell(
                centre, base_props={PropertyType.MATERIAL: "MAT"}
            ) for centre in centres
        ]
        self.fillable: CartesianLattice = CartesianLattice(
            cells, centre=(0.0, 0.0, 0.0), name="TestCartesianLattice"
        )

    def test_init(self) -> None:
        """
        Method that tests the `CartesianLattice` initialization with
        default parameters.
        """
        lattice = CartesianLattice()
        self.assertEqual(lattice.name, f"CartesianLattice_{id(lattice)}")
        self.assertEqual(lattice.layers, [])
        self.assertTrue(lattice.state.is_update_needed)
        self.assertTrue(
            are_same_shapes(
                lattice.o, make_vertex((0.0, 0.0, 0.0)), ShapeType.VERTEX
            )
        )

    def test_init_non_default(self) -> None:
        """
        Method that tests the `CartesianLattice` initialization with
        specified centre, cells and name.
        """
        centre = (5.0, 5.0, 0.0)
        cells = [
            CartesianCell((5.5, 5.5, 0.0)),
            CartesianCell((4.5, 5.5, 0.0)),
            CartesianCell((4.5, 4.5, 0.0)),
            CartesianCell((5.5, 4.5, 0.0)),
        ]
        lattice = CartesianLattice(cells, centre, "Lattice")
        self.assertIsNotNone(lattice.o)
        self.assertTrue(
            are_same_shapes(lattice.o, make_vertex(centre), ShapeType.VERTEX)
        )
        self.assertEqual(lattice.name, f"Lattice_{id(lattice)}")
        self.assertEqual(len(lattice.layers), 1)
        self.assertEqual(len(lattice.layers[0]), 4)
        for c1, c2 in zip(lattice.layers[0], cells):
            self.assertTrue(
                are_same_shapes(c1, c2, ShapeType.COMPOUND)
            )
        self.assertEqual(len(lattice.regions), 4)
        self.assertIsNotNone(lattice.shape)
        self.assertTrue(
            are_same_shapes(
                lattice.shape, Rectangle(centre, 2, 2), ShapeType.FACE
            )
        )

    def test_add_ring_of_cells_valid_ring(self) -> None:
        """
        Method that tests the `add_ring_of_cells` method with a valid ring
        index without specifying the layer index.
        It verifies that the ring of cells is correctly added to the lattice
        in a new layer.
        """
        # Declare a cell and add a ring of it in a new layer
        cell = CartesianCell((0.5, 0.5, 0.0))
        initial_layers = len(self.fillable.layers)
        self.fillable.add_ring_of_cells(cell, ring_index=1)
        # Verify the number of layers has increased by one and that the layer
        # contains four cells
        self.assertEqual(len(self.fillable.layers), initial_layers + 1)
        self.assertEqual(len(self.fillable.layers[-1]), 4)

        # Verify the case in which a new ring of cells is added
        self.fillable.add_ring_of_cells(cell, ring_index=2)
        # Verify the state has updated
        self.assertTrue(self.fillable.state.is_update_needed)
        self.assertEqual(len(self.fillable.layers), initial_layers + 2)
        self.assertEqual(len(self.fillable.layers[-1]), 12)
        # Verify that the shape has updated
        self.assertIsNotNone(self.fillable.shape)
        face = make_face(
            build_compound_borders(make_compound(self.fillable.get_regions()))
        )
        self.assertTrue(
            are_same_shapes(self.fillable.shape, face, ShapeType.FACE)
        )
        # Verify the characteristic dimensions has updated
        x_min, x_max, y_min, y_max = get_bounding_box(face)
        self.assertAlmostEqual(self.fillable.dimensions[0], x_max - x_min)
        self.assertAlmostEqual(self.fillable.dimensions[1], y_max - y_min)

    def test_add_ring_of_cells_invalid_indices(self) -> None:
        """
        Method that tests the `add_ring_of_cells` method when invalid values
        for the ring index or the layer one are specified.
        It verifies that a `ValueError` exception is raised.
        """
        # Verify the 'ValueError' exception is raised with zero ring index
        cell = CartesianCell((0.5, 0.5, 0.0))
        with self.assertRaises(ValueError):
            self.fillable.add_ring_of_cells(cell, ring_index=0)

        # Verify the 'ValueError' exception is raised with a layer index that
        # does not correspond to any present layer
        with self.assertRaises(ValueError):
            self.fillable.add_ring_of_cells(cell, 1, layer_index=99)

    def test_add_ring_of_cells_valid_layer(self) -> None:
        """
        Method that tests the `add_ring_of_cells` method with a valid layer
        index.
        It verifies that the ring of cells is correctly added to the lattice
        in the indicated layer.
        """
        # Declare a cell and add a ring of it in a specific layer
        cell = CartesianCell((0.5, 0.5, 0.0))
        initial_layers = len(self.fillable.layers)
        no_cells_in_layer_0 = len(self.fillable.layers[0])
        self.fillable.add_ring_of_cells(cell, ring_index=1, layer_index=0)
        # Verify the number of layers have not changed and that they contain
        # the added cells + the one present before
        self.assertEqual(len(self.fillable.layers), initial_layers)
        self.assertEqual(len(self.fillable.layers[0]), no_cells_in_layer_0+4)

    def test_add_rings_of_cells_valid_rings(self) -> None:
        """
        Method that tests the `add_rings_of_cells` method with valid
        parameters.
        It verifies that multiple rings of cells are correctly added to the
        lattice in a new layer.
        """
        # Declare a cell and add a ring of it
        cell = CartesianCell((0.5, 0.5, 0.0))
        initial_layers = len(self.fillable.layers)
        self.fillable.add_rings_of_cells(cell, no_rings=2, ring_index=1)
        # Verify the number of layers has increased by one
        self.assertEqual(len(self.fillable.layers), initial_layers + 1)
        # Verify the number of cells in the new layer is 16
        self.assertEqual(len(self.fillable.layers[-1]), 16)
        # Verify the state has updated
        self.assertTrue(self.fillable.state.is_update_needed)
        # Verify the shape is not None
        self.assertIsNotNone(self.fillable.shape)
        face = make_face(
            build_compound_borders(make_compound(self.fillable.get_regions()))
        )
        self.assertTrue(
            are_same_shapes(self.fillable.shape, face, ShapeType.FACE)
        )
        # Verify the characteristic dimensions has updated
        x_min, x_max, y_min, y_max = get_bounding_box(face)
        self.assertAlmostEqual(self.fillable.dimensions[0], x_max - x_min)
        self.assertAlmostEqual(self.fillable.dimensions[1], y_max - y_min)

    def test_add_rings_of_cells_invalid_indices(self) -> None:
        """
        Method that tests the `add_rings_of_cells` method when invalid values
        for the number of rings, the starting ring index, or the layer index
        are specified. It verifies that a `ValueError` exception is raised.
        """
        # Declare a cell and add a ring of it
        cell = CartesianCell((0.5, 0.5, 0.0))
        # Verify the 'ValueError' exception is raised with zero number of
        # rings
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(cell, no_rings=0)
        # Verify the 'ValueError' exception is raised with negative rings
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(cell, no_rings=-1)
        # Verify the 'ValueError' exception is raised with zero ring index
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(cell, no_rings=1, ring_index=0)
        # Verify the 'ValueError' exception is raised with invalid layer index
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(
                cell, no_rings=1, ring_index=1, layer_index=99
            )


class TestHexLattice(TestFillable):
    """
    Test class for assessing the HexLattice class implementation.

    This test class contains unit tests for verifying the correct behaviour
    of the `HexLattice` class, including initialization and ring addition.
    Tests dealing with operations common to all `Fillable` subclasses are
    declared in the `TestFillable` class this class inherits from.

    Attributes
    ----------
    fillable : HexLattice
        The hexagonal lattice under test, whose class inherits from
        `Fillable`.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `HexLattice` class.
        """
        # Declare a hexagonal cell rotated by 90° as the hexagonal lattice
        # expects cells are provided Y-oriented
        cell = HexCell(base_props={PropertyType.MATERIAL: "MAT"})
        cell.rotate(90)
        self.fillable: HexLattice = HexLattice(
            [cell], (0.0, 0.0, 0.0), "TestHexLattice"
        )

    def test_init(self) -> None:
        """
        Method that tests the `HexLattice` initialization with default
        parameters.
        """
        lattice = HexLattice()
        self.assertEqual(lattice.name, f"HexLattice_{id(lattice)}")
        self.assertEqual(lattice.layers, [])
        self.assertTrue(lattice.state.is_update_needed)
        self.assertTrue(
            are_same_shapes(
                lattice.o, make_vertex((0.0, 0.0, 0.0)), ShapeType.VERTEX
            )
        )

    def test_init_non_default(self) -> None:
        """
        Method that tests the `HexLattice` initialization with
        specified centre, cells and name.
        """
        centre = (5.0, 5.0, 0.0)
        cells = set_up_hex_cells(HexCell(centre))
        no_cells = len(cells)
        lattice = HexLattice(cells, centre, "Lattice")
        self.assertIsNotNone(lattice.o)
        self.assertTrue(
            are_same_shapes(lattice.o, make_vertex(centre), ShapeType.VERTEX)
        )
        self.assertEqual(lattice.name, f"Lattice_{id(lattice)}")
        self.assertEqual(len(lattice.layers), 1)
        self.assertEqual(len(lattice.layers[0]), no_cells)
        for c1, c2 in zip(lattice.layers[0], cells):
            self.assertTrue(
                are_same_shapes(c1, c2, ShapeType.COMPOUND)
            )
        self.assertEqual(len(lattice.regions), no_cells)
        self.assertIsNotNone(lattice.shape)
        face = make_face(
            build_compound_borders(make_compound(lattice.regions))
        )
        self.assertTrue(
            are_same_shapes(lattice.shape, face, ShapeType.FACE)
        )

    def test_add_ring_of_cells_valid_ring(self) -> None:
        """
        Method that tests the `add_ring_of_cells` method with a valid ring
        index without specifying the layer index.
        It verifies that the ring of cells is correctly added to the lattice
        in a new layer.
        """
        # Declare a cell and add a ring of it in a new layer
        cell = HexCell()
        cell.rotate(90)
        initial_layers = len(self.fillable.layers)
        self.fillable.add_ring_of_cells(cell, ring_index=1)
        # Verify the number of layers has increased by one and that the layer
        # contains four cells
        self.assertEqual(len(self.fillable.layers), initial_layers + 1)
        self.assertEqual(len(self.fillable.layers[-1]), 6)

        # Verify the case in which a new ring of cells is added
        self.fillable.add_ring_of_cells(cell, ring_index=2)
        # Verify the state has updated
        self.assertTrue(self.fillable.state.is_update_needed)
        self.assertEqual(len(self.fillable.layers), initial_layers + 2)
        self.assertEqual(len(self.fillable.layers[-1]), 12)
        # Verify that the shape has updated
        self.assertIsNotNone(self.fillable.shape)
        _, _, y_min, y_max = get_bounding_box(
            make_compound(self.fillable.get_regions())
        )
        apothem = (y_max - y_min) / 2
        shape = Hexagon(edge_length=apothem / sin(pi/3))
        self.assertTrue(
            are_same_shapes(self.fillable.shape, shape, ShapeType.FACE)
        )
        # Verify the characteristic dimensions has updated
        x_min, x_max, y_min, y_max = get_bounding_box(shape)
        self.assertTrue(
            isclose(
                self.fillable.dimensions[0],
                (x_max - x_min) / 2,
                abs_tol=1e-6
            )
        )
        self.assertTrue(
            isclose(
                self.fillable.dimensions[1],
                (y_max - y_min) / 2,
                abs_tol=1e-6
            )
        )

    def test_add_ring_of_cells_invalid_indices(self) -> None:
        """
        Method that tests the `add_ring_of_cells` method when invalid values
        for the ring index or the layer one are specified.
        It verifies that a `ValueError` exception is raised.
        """
        # Verify the 'ValueError' exception is raised with zero ring index
        cell = CartesianCell((0.5, 0.5, 0.0))
        with self.assertRaises(ValueError):
            self.fillable.add_ring_of_cells(cell, ring_index=0)

        # Verify the 'ValueError' exception is raised with a layer index that
        # does not correspond to any present layer
        with self.assertRaises(ValueError):
            self.fillable.add_ring_of_cells(cell, 1, layer_index=99)

    def test_add_ring_of_cells_valid_layer(self) -> None:
        """
        Method that tests the `add_ring_of_cells` method with a valid layer
        index.
        It verifies that the ring of cells is correctly added to the lattice
        in the indicated layer.
        """
        # Declare a cell and add a ring of it in a specific layer
        cell = CartesianCell((0.5, 0.5, 0.0))
        initial_layers = len(self.fillable.layers)
        no_cells_in_layer_0 = len(self.fillable.layers[0])
        self.fillable.add_ring_of_cells(cell, ring_index=1, layer_index=0)
        # Verify the number of layers have not changed and that they contain
        # the added cells + the one present before
        self.assertEqual(len(self.fillable.layers), initial_layers)
        self.assertEqual(len(self.fillable.layers[0]), no_cells_in_layer_0+6)

    def test_add_rings_of_cells_valid_rings(self) -> None:
        """
        Method that tests the `add_rings_of_cells` method with valid
        parameters.
        It verifies that multiple rings of cells are correctly added to the
        lattice in a new layer.
        """
        # Declare a cell and add a ring of it
        cell = CartesianCell((0.5, 0.5, 0.0))
        initial_layers = len(self.fillable.layers)
        self.fillable.add_rings_of_cells(cell, no_rings=2, ring_index=1)
        # Verify the number of layers has increased by one
        self.assertEqual(len(self.fillable.layers), initial_layers + 1)
        # Verify the number of cells in the new layer is 16
        self.assertEqual(len(self.fillable.layers[-1]), 18)
        # Verify the state has updated
        self.assertTrue(self.fillable.state.is_update_needed)
        # Verify the shape is not None
        self.assertIsNotNone(self.fillable.shape)
        _, _, y_min, y_max = get_bounding_box(
            make_compound(self.fillable.get_regions())
        )
        apothem = (y_max - y_min) / 2
        shape = Hexagon(edge_length=apothem / sin(pi/3))
        self.assertTrue(
            are_same_shapes(self.fillable.shape, shape, ShapeType.FACE)
        )
        # Verify the characteristic dimensions has updated
        x_min, x_max, y_min, y_max = get_bounding_box(shape)
        self.assertTrue(
            isclose(
                self.fillable.dimensions[0],
                (x_max - x_min) / 2,
                abs_tol=1e-6
            )
        )
        self.assertTrue(
            isclose(
                self.fillable.dimensions[1],
                (y_max - y_min) / 2,
                abs_tol=1e-6
            )
        )

    def test_add_rings_of_cells_invalid_indices(self) -> None:
        """
        Method that tests the `add_rings_of_cells` method when invalid values
        for the number of rings, the starting ring index, or the layer index
        are specified. It verifies that a `ValueError` exception is raised.
        """
        # Declare a cell and add a ring of it
        cell = CartesianCell((0.5, 0.5, 0.0))
        # Verify the 'ValueError' exception is raised with zero number of
        # rings
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(cell, no_rings=0)
        # Verify the 'ValueError' exception is raised with negative rings
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(cell, no_rings=-1)
        # Verify the 'ValueError' exception is raised with zero ring index
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(cell, no_rings=1, ring_index=0)
        # Verify the 'ValueError' exception is raised with invalid layer index
        with self.assertRaises(ValueError):
            self.fillable.add_rings_of_cells(
                cell, no_rings=1, ring_index=1, layer_index=99
            )


class TestLattice(TestFillable):
    """
    Test class for assessing the `Lattice` class implementation.

    This test class contains unit tests for verifying the correct behaviour
    of the `Lattice` class, including initialization and fillables addition.
    Tests dealing with operations common to all `Fillable` subclasses are
    declared in the `TestFillable` class this class inherits from.

    Attributes
    ----------
    fillable : Lattice
        The lattice under test, whose class inherits from `Fillable`.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `Lattice` class.
        """
        cell = CartesianCell(base_props={PropertyType.MATERIAL: "MAT"})
        self.fillable: Lattice = Lattice(
            cells=set_up_rect_cells(cell, True),
            centre=(0.0, 0.0, 0.0),
            name="TestLattice"
        )

    def test_init(self) -> None:
        """
        Method that tests the `Lattice` initialization with default
        parameters.
        """
        lattice = Lattice()
        self.assertEqual(lattice.name, f"Lattice_{id(lattice)}")
        self.assertEqual(lattice.layers, [])
        self.assertTrue(lattice.state.is_update_needed)
        self.assertTrue(
            are_same_shapes(
                lattice.o, make_vertex((0.0, 0.0, 0.0)), ShapeType.VERTEX
            )
        )

    def test_init_non_default(self) -> None:
        """
        Method that tests the `Lattice` initialization with specified centre,
        cells and name.
        """
        centre = (5.0, 5.0, 0.0)
        cells = [
            CartesianCell((5.5, 5.5, 0.0)),
            CartesianCell((4.5, 5.5, 0.0)),
            CartesianCell((4.5, 4.5, 0.0)),
            CartesianCell((5.5, 4.5, 0.0)),
        ]
        lattice = Lattice(cells, centre, "Lattice")
        self.assertIsNotNone(lattice.o)
        self.assertTrue(
            are_same_shapes(lattice.o, make_vertex(centre), ShapeType.VERTEX)
        )
        self.assertEqual(lattice.name, f"Lattice_{id(lattice)}")
        self.assertEqual(len(lattice.layers), 1)
        self.assertEqual(len(lattice.layers[0]), 4)
        for c1, c2 in zip(lattice.layers[0], cells):
            self.assertTrue(
                are_same_shapes(c1, c2, ShapeType.COMPOUND)
            )
        self.assertEqual(len(lattice.regions), 4)
        self.assertIsNotNone(lattice.shape)
        self.assertTrue(
            are_same_shapes(
                lattice.shape, Rectangle(centre, 2, 2), ShapeType.FACE
            )
        )

    def test_add(self) -> None:
        """
        Method that tests the `add` method by checking if the lattice's shape
        and the characteristic dimensions are correctly updated.
        """
        # Store the lattice's characteristic shape and dimensions
        pre_shape = self.fillable.shape
        pre_dimensions = self.fillable.dimensions
        # Add a circular region within the lattice's borders
        self.fillable.add(Region(Circle()))
        # Verify the shape and the characteristic dimensions have not changed
        self.assertTrue(
            are_same_shapes(self.fillable.shape, pre_shape, ShapeType.FACE)
        )
        self.assertEqual(self.fillable.dimensions, pre_dimensions)

        # Add a cell outside the lattice's borders
        self.fillable.add(CartesianCell((2.0, 2.0, 0.0)))
        # Verify the shape and the characteristic dimensions have changed
        self.assertFalse(
            are_same_shapes(self.fillable.shape, pre_shape, ShapeType.FACE)
        )
        self.assertNotEqual(self.fillable.dimensions, pre_dimensions)
        self.assertEqual(
            self.fillable.dimensions, tuple([i + 0.5 for i in pre_dimensions])
        )


class TestComputeSubdivisionPointsOnBorders(unittest.TestCase):
    """
    Test class for the `compute_subdivision_points_on_borders` function of
    the `glow.geometry_layouts.lattices` module.

    This test class verifies the correct generation of evenly spaced vertices
    along the borders of a geometric surface.
    """
    def test_returned_list_of_points(self) -> None:
        """
        Method that tests that the function returns a list of tuples with the
        correct number of points.
        """
        rect = Rectangle((0.0, 0.0, 0.0), 10.0, 10.0)
        points = compute_subdivision_points_on_borders(4, rect)
        self.assertIsInstance(points, List)
        self.assertTrue(len(points), 16)
        for point in points:
            self.assertIsInstance(point, Tuple)
            self.assertEqual(len(point), 3)

    def test_returned_correct_points_coordinates(self) -> None:
        """
        Method that tests that the function returns a list of points having
        the expected coordinates.
        """
        rect = Rectangle((0.0, 0.0, 0.0), 10.0, 10.0)
        no_vrtcs = 2
        points = compute_subdivision_points_on_borders(no_vrtcs, rect)
        expected_coords = [
            (-5.0, -5.0),
            (0.0, -5.0),
            (5.0, -5.0),
            (5.0, 0.0),
            (5.0, 5.0),
            (0.0, 5.0),
            (-5.0, 5.0),
            (-5.0, 0.0)
        ]
        for p in points:
            found = False
            for ep in expected_coords:
                if all(
                    [isclose(p[i], ep[i], abs_tol=1e-6) for i in range(2)]
                ):
                    found = True
                    break
            self.assertTrue(found)


class TestEnsureNotZero(unittest.TestCase):
    """
    Test class for the `ensure_not_zero` function of the
    `glow.geometry_layouts.lattices` module.

    This test class verifies the correct validation of non-zero values
    and appropriate exception raising.
    """
    def test_raises_error_when_zero(self) -> None:
        """
        Method that tests that `ValueError` is raised when the nput value
        is zero.
        """
        with self.assertRaises(ValueError):
            ensure_not_zero(0)

    def test_no_error(self) -> None:
        """
        Method that tests that no error is raised for values different from
        zero.
        """
        try:
            ensure_not_zero(5)
            ensure_not_zero(-5)
        except ValueError:
            self.fail(
                "Function 'ensure_not_zero' raised 'ValueError' unexpectedly "
                "for value different from zero."
            )

    def test_custom_message(self) -> None:
        """
        Method that tests that a custom error message is properly used when
        the value is zero.
        """
        custom_msg = "Custom error message"
        with self.assertRaisesRegex(ValueError, custom_msg):
            ensure_not_zero(0, custom_msg)


class TestGetCellAtCentres(unittest.TestCase):
    """
    Test class for the `get_cell_at_centres` function of the
    `glow.geometry_layouts.lattices` module.

    This test class verifies the correct placement and cloning of cells
    at specified centre positions.
    """
    def test_returns_list_of_cells(self) -> None:
        """
        Method that tests that the function returns a list of `Cell` objects.
        """
        # Create a cell with rectangular shape
        rect = Rectangle((0.0, 0.0, 0.0), 5.0, 5.0)
        cell = Cell(shape=rect)
        centres = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
        z_axis = wrap_shape(
            build_z_axis_from_vertex(make_vertex((0.0, 0.0, 0.0)))
        )
        # Get the list of cells
        result = get_cell_at_centres(cell, centres, 0.0, z_axis)

        # Verify the correctness of the returned output list
        self.assertIsInstance(result, List)
        self.assertEqual(len(result), len(centres))
        for cell, pos in zip(result, centres):
            self.assertIsInstance(cell, Cell)
            self.assertTrue(
                are_same_shapes(cell.o, make_vertex(pos), ShapeType.VERTEX)
            )

    def test_returns_rotated_cells(self) -> None:
        """
        Method that tests that the returned cells are rotated around the
        given axis.
        """
        # Create a cell with rectangular shape
        rect = Rectangle((0.0, 0.0, 0.0), 5.0, 5.0)
        cell = Cell(shape=rect)
        centres = [(0.0, 0.0, 0.0), (5.0, 5.0, 0.0), (10.0, 10.0, 0.0)]
        # Build the rotation axis in a given position
        axis_centre = (1, 1, 0)
        axis_centre_vrtx = make_vertex(axis_centre)
        z_axis = wrap_shape(build_z_axis_from_vertex(axis_centre_vrtx))
        # Build reference vector before the rotation (from cell centre to
        # axis centre)
        ref_vect = make_vector_from_points(cell.o, axis_centre_vrtx)
        rot_angle = 30.0

        # Get the list of cells
        result = get_cell_at_centres(cell, centres, rot_angle, z_axis)

        # Verify the correct rotation happened
        for c in result:
            # Translate the original axis in the cell's centre, if needed
            ref_vect_tr = ref_vect
            if not are_same_shapes(cell.o, c.o, ShapeType.VERTEX):
                ref_vect_tr = make_translation(
                    ref_vect, make_vector_from_points(cell.o, c.o)
                )
            # Build reference vector for each cell after the rotation
            ref_vect2 = make_vector_from_points(c.o, axis_centre_vrtx)
            # Check the correct rotation happened
            self.assertEqual(c.rot_angle, rot_angle)
            self.assertTrue(
                isclose(
                    get_angle_between_shapes(ref_vect_tr, ref_vect2),
                    rot_angle,
                    abs_tol=1e-6
                )
            )
