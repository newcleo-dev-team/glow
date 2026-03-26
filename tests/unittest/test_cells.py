"""
Module containing unittest class to assess that the classes of the
`glow.geometry_layouts.cells` module have a valid implementation.
"""
import math
from typing import List

from glow.geometry_layouts.cells import CartesianCell, Cell, HexCell
from glow.geometry_layouts.geometries import Circle, Hexagon, Rectangle
from glow.geometry_layouts.lattices import Lattice
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_kind_of_shape, get_shape_type, make_compound, make_vertex
from glow.support.types import GeometryType, PropertyType
from glow.support.utility import are_same_shapes
from tests.unittest.test_fillable_layouts import TestFillable


class TestCell(TestFillable):
    """
    Test case for verifying the construction operations of the `Cell` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Cell` class can be correctly instantiated and that it properly
    handles sectorisation operations.
    Tests dealing with operations common to all `Fillable` subclasses are
    declared in the `TestFillable` class this class inherits from.

    Attributes
    ----------
    fillable : Cell
        The `Cell` object under test, whose class inherits from `Fillable`.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `Cell` class.
        """
        super().setUp()
        self.fillable: Cell = Cell(
            shape=Circle(center=(0.0, 0.0, 0.0), radius=4.0),
            base_props={PropertyType.MATERIAL: "MAT_CELL"},
            name="TestCell"
        )

    def test_init(self) -> None:
        """
        Method that tests the `Cell` initialisation with default parameters.
        """
        # Instantiate a 'Cell' object
        shape = Circle()
        cell = Cell(shape)
        # Verify the correct initialisation
        self.assertIsNotNone(cell.name)
        self.assertEqual(cell.name, f"Cell_{id(cell)}")
        self.assertEqual(len(cell.layers), 1)
        self.assertEqual(len(cell.regions), 1)
        self.assertEqual(len(cell.layers[0]), 1)
        self.assertIsInstance(cell.layers[0][0], Region)
        self.assertIsNotNone(cell.shape)
        self.assertTrue(
            are_same_shapes(cell.shape, shape, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(cell.layers[0][0], shape, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(cell.regions[0], shape, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(
                cell.geom_obj, make_compound([shape]), ShapeType.COMPOUND
            )
        )
        self.assertIsNone(cell.layers[0][0].properties)
        self.assertEqual(cell.dimensions, shape.dimensions)

    def test_init_with_properties(self) -> None:
        """
        Method that tests the `Cell` initialisation with defined properties.
        """
        # Instantiate a 'Cell' object with properties
        props = {PropertyType.MATERIAL: "MAT"}
        cell = Cell(Circle(), props)
        # Verify the correct assignment of properties to the single region
        # available
        self.assertEqual(len(cell.regions), 1)
        self.assertEqual(cell.regions[0].properties, props)

    def test_sectorize_inputs(self) -> None:
        """
        Method that tests the `sectorize` method with valid sectors and
        starting angles and with unexpected `kwargs` which raises a
        `TypeError` exception.
        """
        # Verify the sectorisation with valid arguments produces the expected
        # edges
        sectors_no = [2]
        angles = [0.0]
        self.fillable.sectorize(sectors_no, angles)
        self.assertIn(GeometryType.SECTORIZED, self.fillable.geometry_maps)
        self.assertEqual(
            len(
                extract_sub_shapes(
                    self.fillable.geometry_maps[GeometryType.SECTORIZED],
                    ShapeType.EDGE
                )
            ),
            sectors_no[0]
        )

        # Verify the sectorisation raises 'TypeError' with invalid kwargs
        with self.assertRaises(TypeError):
            self.fillable.sectorize(sectors_no, angles, unexpected_param=True)

    def test_sectorize_with_nested_fillable(self) -> None:
        """
        Method that tests the sectorization raises a `RuntimeError` exception
        when the cell contains nested `Fillable` objects.
        """
        # Add a 'Lattice' to the cell
        self.fillable.add(
            Lattice(cells=[self.fillable])
        )
        # Verify the sectorisation raises an exception for a cell with nested
        # layouts
        with self.assertRaises(RuntimeError):
            self.fillable.sectorize([2], [0.0])


class TestCartesianCell(TestFillable):
    """
    Test case for verifying the construction operations of the `CartesianCell`
    class.

    This test suite provides common setup and a set of tests to ensure that
    the `CartesianCell` class can be correctly instantiated and that it
    properly handles sectorisation operations.
    Tests dealing with operations common to all `Fillable` subclasses are
    declared in the `TestFillable` class this class inherits from.

    Attributes
    ----------
    fillable : CartesianCell
        The `CartesianCell` object under test, whose class inherits from
        `Fillable`.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `CartesianCell`
        class.
        """
        super().setUp()
        self.fillable: CartesianCell = CartesianCell(
            center=(0.0, 0.0, 0.0),
            width_height=(2.0, 2.0),
            base_props={PropertyType.MATERIAL: "MAT_CELL"},
            name="TestCartesianCell"
        )

    def test_init(self) -> None:
        """
        Method that tests the `CartesianCell` initialisation with default
        parameters.
        """
        rect = Rectangle()
        # Instantiate a 'CartesianCell' object
        cell = CartesianCell()
        # Verify the correct initialisation
        self.assertIsNotNone(cell.name)
        self.assertEqual(cell.name, f"Cartesian_Cell_{id(cell)}")
        self.assertEqual(len(cell.layers), 1)
        self.assertEqual(len(cell.regions), 1)
        self.assertEqual(len(cell.layers[0]), 1)
        self.assertIsInstance(cell.layers[0][0], Region)
        self.assertIsNotNone(cell.shape)
        self.assertIsInstance(self.fillable.shape, Rectangle)
        self.assertTrue(
            are_same_shapes(
                cell.o, make_vertex((0.0, 0.0, 0.0)), ShapeType.VERTEX
            )
        )
        self.assertTrue(
            are_same_shapes(cell.shape, rect, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(cell.layers[0][0], rect, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(cell.regions[0], rect, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(
                cell.geom_obj, make_compound([rect]), ShapeType.COMPOUND
            )
        )
        self.assertIsNone(cell.layers[0][0].properties)
        self.assertEqual(cell.dimensions, rect.dimensions)

    def test_init_with_properties(self) -> None:
        """
        Method that tests the `CartesianCell` initialisation with defined
        properties.
        """
        # Instantiate a 'CartesianCell' object with properties
        props = {PropertyType.MATERIAL: "MAT"}
        cell = CartesianCell(base_props=props)
        # Verify the correct assignment of properties to the single region
        # available
        self.assertEqual(len(cell.regions), 1)
        self.assertEqual(cell.regions[0].properties, props)

    def test_init_with_rounded_corners(self) -> None:
        """
        Method that tests the `CartesianCell` initialisation with rounded
        corners.
        """
        # Instantiate a 'CartesianCell' with rounded corners
        curv_radius = 0.2
        rounded_corners = [(0, curv_radius), (1, curv_radius)]
        cell = CartesianCell(
            center=(0.0, 0.0, 0.0),
            width_height=(2.0, 2.0),
            rounded_corners=rounded_corners,
            name="RoundedCell"
        )

        # Assess the correct creation of the arc borders
        border_arcs = [
            e for e in cell.shape.borders
            if str(get_kind_of_shape(e)[0]) == "ARC_CIRCLE"
        ]
        self.assertEqual(len(border_arcs), len(rounded_corners))
        lx, ly = cell.dimensions
        for arc in border_arcs:
            # Get the geometric information about each arc
            data = get_kind_of_shape(arc)
            self.assertTrue(get_shape_type(arc) == ShapeType.EDGE)
            self.assertTrue(
                math.isclose(data[7], curv_radius, abs_tol=1e-6)
            )
        self.assertTrue(
            are_same_shapes(
                cell.shape,
                Rectangle(
                    height=ly, width=lx, rounded_corners=rounded_corners
                ),
                ShapeType.FACE
            )
        )

    def test_sectorize_inputs(self) -> None:
        """
        Method that tests the `sectorize` method with valid sectors and
        starting angles and with unexpected `kwargs` which does not have
        any effect.
        """
        # Verify the sectorisation with valid arguments produces the expected
        # number of edges
        sectors_no = [2]
        angles = [0.0]
        self.fillable.sectorize(sectors_no, angles)
        self.assertIn(GeometryType.SECTORIZED, self.fillable.geometry_maps)
        self.assertEqual(
            len(
                extract_sub_shapes(
                    self.fillable.geometry_maps[GeometryType.SECTORIZED],
                    ShapeType.EDGE
                )
            ),
            sectors_no[0]
        )

        # Verify the sectorisation with invalid kwargs does not have any
        # effect
        self.fillable.sectorize(sectors_no, angles, unexpected_param=True)
        self.assertEqual(
            len(
                extract_sub_shapes(
                    self.fillable.geometry_maps[GeometryType.SECTORIZED],
                    ShapeType.EDGE
                )
            ),
            sectors_no[0]
        )

    def test_sectorize_valid_windmill(self) -> None:
        """
        Method that tests the `sectorize` method when the `windmill` flag
        is specified.
        """
        # Verify the sectorisation with the right combination of sectors
        # numbers and starting angles (8 and 16 sectors combinations)
        sectors_no = [8]
        for val in [22.5 + j * 45 for j in range(sectors_no[0])]:
            angles = [val]
            self.__assess_valid_windmill_sect(sectors_no, angles)
        sectors_no = [16]
        for val in [0.0 + j * 45 for j in range(sectors_no[0])]:
            angles = [val]
            self.__assess_valid_windmill_sect(sectors_no, angles)

        # Apply the sectorisation with windmill flag set to False
        sectors_no = [8]
        angles = [22.5]
        self.fillable.sectorize(sectors_no, angles, windmill=False)
        self.assertIn(GeometryType.SECTORIZED, self.fillable.geometry_maps)
        self.assertEqual(
            len(
                extract_sub_shapes(
                    self.fillable.geometry_maps[GeometryType.SECTORIZED],
                    ShapeType.EDGE
                )
            ),
            sectors_no[0]
        )

    def test_sectorize_windmill_invalid_config(self) -> None:
        """
        Method that tests the `sectorize` method with the `windmill` parameter
        is `True`, but it is specified a combination of number of sectors and
        starting angles that is unsupported for the windmill sectorisation.
        """
        sectors_no = [4]
        angles = [22.5]
        self.fillable.sectorize(sectors_no, angles, windmill=True)
        # Verify no windmill edges have been created
        self.assertIn(GeometryType.SECTORIZED, self.fillable.geometry_maps)
        self.assertEqual(
            len(
                extract_sub_shapes(
                    self.fillable.geometry_maps[GeometryType.SECTORIZED],
                    ShapeType.EDGE
                )
            ),
            sectors_no[0]
        )

    def __assess_valid_windmill_sect(
            self, sectors_no: List[int], angles: List[float]
        ) -> None:
        """
        Method that verifies the correct application of the sectorisation
        with the creation of the four additional sectorisation edges
        representative of the windmill option.
        """
        self.fillable.sectorize(sectors_no, angles, windmill=True)
        # Verify that the resulting edges include the four edges of the
        # windmill sectorisation
        self.assertIn(GeometryType.SECTORIZED, self.fillable.geometry_maps)
        self.assertEqual(
            len(
                extract_sub_shapes(
                    self.fillable.geometry_maps[GeometryType.SECTORIZED],
                    ShapeType.EDGE
                )
            ),
            sectors_no[0] + 4
        )


class TestHexCell(TestFillable):
    """
    Test case for verifying the construction operations of the `HexCell`
    class.

    This test suite provides common setup and a set of tests to ensure that
    the `HexCell` class can be correctly instantiated and that it properly
    handles sectorisation operations.
    Tests dealing with operations common to all `Fillable` subclasses are
    declared in the `TestFillable` class this class inherits from.

    Attributes
    ----------
    fillable : HexCell
        The `HexCell` object under test, whose class inherits from `Fillable`.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the `HexCell` class.
        """
        super().setUp()
        self.fillable: HexCell = HexCell(
            center=(0.0, 0.0, 0.0),
            side=1.0,
            base_props={PropertyType.MATERIAL: "MAT_CELL"},
            name="TestHexCell"
        )

    def test_init(self) -> None:
        """
        Method that tests the `HexCell` initialisation with default
        parameters.
        """
        hex = Hexagon()
        # Instantiate a 'HexCell' object
        cell = HexCell()
        # Verify the correct initialisation
        self.assertIsNotNone(cell.name)
        self.assertEqual(cell.name, f"Hexagonal_Cell_{id(cell)}")
        self.assertEqual(len(cell.layers), 1)
        self.assertEqual(len(cell.regions), 1)
        self.assertEqual(len(cell.layers[0]), 1)
        self.assertIsInstance(cell.layers[0][0], Region)
        self.assertIsNotNone(cell.shape)
        self.assertIsInstance(self.fillable.shape, Hexagon)
        self.assertTrue(
            are_same_shapes(
                cell.o, make_vertex((0.0, 0.0, 0.0)), ShapeType.VERTEX
            )
        )
        self.assertTrue(
            are_same_shapes(cell.shape, hex, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(cell.layers[0][0], hex, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(cell.regions[0], hex, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(
                cell.geom_obj, make_compound([hex]), ShapeType.COMPOUND
            )
        )
        self.assertIsNone(cell.layers[0][0].properties)
        self.assertEqual(cell.dimensions, hex.dimensions)

    def test_init_with_properties(self) -> None:
        """
        Method that tests the `HexCell` initialisation with defined
        properties.
        """
        # Instantiate a 'HexCell' object with properties
        props = {PropertyType.MATERIAL: "MAT"}
        cell = HexCell(base_props=props)
        # Verify the correct assignment of properties to the single region
        # available
        self.assertEqual(len(cell.regions), 1)
        self.assertEqual(cell.regions[0].properties, props)

    def test_sectorize_inputs(self) -> None:
        """
        Method that tests the `sectorize` method with valid sectors and
        starting angles and with unexpected `kwargs` which does not have
        any effect.
        """
        # Verify the sectorisation with valid arguments produces the expected
        # number of edges
        sectors_no = [6]
        angles = [0.0]
        self.fillable.sectorize(sectors_no, angles)
        self.assertIn(GeometryType.SECTORIZED, self.fillable.geometry_maps)
        self.assertEqual(
            len(
                extract_sub_shapes(
                    self.fillable.geometry_maps[GeometryType.SECTORIZED],
                    ShapeType.EDGE
                )
            ),
            sectors_no[0]
        )
