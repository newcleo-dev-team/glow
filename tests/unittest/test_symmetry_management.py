"""
Module containing unittest class to assess that the symmetry management
functions of the module `glow.geometry_layouts.symmetry_management.py`
have a valid implementation.
"""
import math
import unittest

from typing import Tuple

from glow.geometry_layouts.geometries import Rectangle, Hexagon, Surface
from glow.geometry_layouts.symmetry_management import SymmetryDomain, \
    build_cartesian_symmetry_shape, build_hex_symmetry_shape
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_basic_properties, get_point_coordinates, make_vertex
from glow.support.types import SymmetryType
from glow.support.utility import are_same_shapes


class TestSymmetryDomain(unittest.TestCase):
    """
    Test case for verifying the `SymmetryDomain` dataclass implementation.

    Attributes
    ----------
    rect_shape : Rectangle
        Rectangular shape used in tests.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for `SymmetryDomain` tests.
        """
        self.rect_shape = Rectangle((0.0, 0.0, 0.0), 4.0, 4.0)

    def test_init_hex_domain(self) -> None:
        """
        Method that tests the initialisation with a hexagonal shape.
        """
        origin = (1.0, 1.0, 0.0)
        hex_shape = Hexagon(origin, edge_length=2.0)
        domain = SymmetryDomain(
            origin=origin, full_layout_shape=hex_shape
        )

        # Verify the correct initialisation
        self.assertEqual(domain.origin, origin)
        self.assertEqual(domain.full_layout_shape, hex_shape)
        self.assertIsNotNone(domain.x_bounds)
        self.assertIsNotNone(domain.y_bounds)
        self.assertTrue(
            all(
                math.isclose(a, b, abs_tol=1e-6)
                for a, b in zip(
                    domain.x_bounds,
                    (
                        origin[0] - hex_shape.dimensions[0],
                        origin[0] + hex_shape.dimensions[0]
                    )
                )
            )
        )
        self.assertTrue(
            all(
                math.isclose(a, b, abs_tol=1e-6)
                for a, b in zip(
                    domain.y_bounds,
                    (
                        origin[1] - hex_shape.dimensions[1],
                        origin[1] + hex_shape.dimensions[1]
                    )
                )
            )
        )

    def test_init_rect_domain(self) -> None:
        """
        Method that tests the initialisation with a rectangular shape.
        """
        origin = (1.0, 1.0, 0.0)
        rect_shape = Rectangle(origin, 4.0, 4.0)
        domain = SymmetryDomain(
            origin=origin, full_layout_shape=rect_shape
        )

        # Verify the correct initialisation
        self.assertEqual(domain.origin, origin)
        self.assertEqual(domain.full_layout_shape, rect_shape)
        self.assertIsNotNone(domain.x_bounds)
        self.assertIsNotNone(domain.y_bounds)
        self.assertTrue(
            all(
                math.isclose(a, b, abs_tol=1e-6)
                for a, b in zip(
                    domain.x_bounds,
                    (
                        origin[0] - rect_shape.dimensions[0]/2,
                        origin[0] + rect_shape.dimensions[0]/2
                    )
                )
            )
        )
        self.assertTrue(
            all(
                math.isclose(a, b, abs_tol=1e-6)
                for a, b in zip(
                    domain.y_bounds,
                    (
                        origin[1] - rect_shape.dimensions[1]/2,
                        origin[1] + rect_shape.dimensions[1]/2
                    )
                )
            )
        )


class TestBuildCartesianSymmetryShape(unittest.TestCase):
    """
    Test case for verifying the construction of symmetry shapes for
    Cartesian-type layouts.

    Attributes
    ----------
    centre : Tuple[float, float, float]
        The XYZ coordinates of the centre of the full layout shape.
    context : SymmetryDomain
        The dataclass storing the data used for building the symmetry shape.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for Cartesian symmetry shape
        building functions.
        """
        self.centre: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        self.context: SymmetryDomain = SymmetryDomain(
            origin=self.centre,
            full_layout_shape=Rectangle(self.centre, 4.0, 4.0)
        )

    def test_build_cartesian_symmetry_full(self) -> None:
        """
        Method that tests building a FULL symmetry shape for a Cartesian
        layout. It verifies that the result is a copy of the full layout
        shape.
        """
        result = build_cartesian_symmetry_shape(
            SymmetryType.FULL, self.context
        )
        # Verify the returned shape is the full layout
        self.assertTrue(
            are_same_shapes(
                result, self.context.full_layout_shape, ShapeType.FACE
            )
        )

    def test_build_cartesian_symmetry_half(self) -> None:
        """
        Method that tests building a HALF symmetry shape for a Cartesian
        layout. It verifies that the result is half of the source rectangle
        by checking the shape's centre and the area.
        """
        result = build_cartesian_symmetry_shape(
            SymmetryType.HALF, self.context
        )
        # Verify the returned shape is the half layout
        self.__assess_cartesian_symmetry(
            result,
            (
                self.centre[0],
                self.centre[1]
                - self.context.full_layout_shape.dimensions[1] / 2,
                0.0
            ),
            2
        )
        expected_shape_o = (
            self.centre[0] + self.context.full_layout_shape.dimensions[0] / 4,
            self.centre[1],
            0.0
        )
        self.assertTrue(
            are_same_shapes(
                result.o, make_vertex(expected_shape_o), ShapeType.VERTEX
            )
        )

    def test_build_cartesian_symmetry_quarter(self) -> None:
        """
        Method that tests building a QUARTER symmetry shape for a Cartesian
        layout. It verifies that the result is a quarter of the source
        rectangle by checking the shape's centre and the area.
        """
        result = build_cartesian_symmetry_shape(
            SymmetryType.QUARTER, self.context
        )
        # Verify the returned shape is a quarter of the layout
        self.__assess_cartesian_symmetry(result, self.centre, 4)
        expected_shape_o = (
            self.centre[0] + self.context.full_layout_shape.dimensions[0] / 4,
            self.centre[1] + self.context.full_layout_shape.dimensions[1] / 4,
            0.0
        )
        self.assertTrue(
            are_same_shapes(
                result.o, make_vertex(expected_shape_o), ShapeType.VERTEX
            )
        )

    def test_build_cartesian_symmetry_eighth(self) -> None:
        """
        Method that tests building an EIGHTH symmetry shape for a Cartesian
        layout. It verifies that the result is an eighth of the source
        rectangle by checking the shape's left corner and the area.
        """
        result = build_cartesian_symmetry_shape(
            SymmetryType.EIGHTH, self.context
        )
        # Verify the returned shape is an eighth of the layout
        self.__assess_cartesian_symmetry(result, self.centre, 8)

    def test_build_cartesian_symmetry_diag(self) -> None:
        """
        Method that tests building a DIAG symmetry shape for a Cartesian
        layout. It verifies that the result is half of the source rectangle
        by checking the shape's left corner and the area.
        """
        result = build_cartesian_symmetry_shape(
            SymmetryType.DIAG, self.context
        )
        # Verify the returned shape is half of the layout along the diagonal
        self.__assess_cartesian_symmetry(
            result,
            (
                self.centre[0]
                - self.context.full_layout_shape.dimensions[0] / 2,
                self.centre[1]
                - self.context.full_layout_shape.dimensions[1] / 2,
                0.0
            ),
            2
        )

    def test_build_cartesian_symmetry_unsupported_type(self) -> None:
        """
        Method that tests that an unsupported symmetry type for Cartesian
        layouts raises a `RuntimeError`.
        """
        with self.assertRaises(RuntimeError):
            build_cartesian_symmetry_shape(SymmetryType.THIRD, self.context)

    def __assess_cartesian_symmetry(
            self,
            result: Surface,
            expected_lower_left: Tuple[float, float, float],
            area_factor: float
        ) -> None:
        """
        Method that assesses the correctness of a computed Cartesian symmetry
        shape by verifying that:

        - the position of the lower-left corner match the expected given
          coordinates;
        - the area of the intersection between the full rectangle and the
          shape of the symmetry matches (within a tolerance) with the area of
          the full rectangle scaled by the given factor.

        Parameters
        ----------
        result : Surface
            The geometric surface representing the shape of the analysed
            symmetry. It used to extract a portion of the rectangular full
            layout.
        expected_lower_left : Tuple[float, float, float]
            The expected coordinates of the lower-left vertex of the portion
            of the full rectangle.
        area_factor : float
            A factor for scaling the area of the full rectangle.
        """
        # Get the XYZ coordinates of the vertices of the symmetry shape
        coords = [
            get_point_coordinates(p)
            for p in extract_sub_shapes(result, ShapeType.VERTEX)
        ]
        # Retrieve the lower-left corner
        lower_left = min(coords, key=lambda c: (c[0], c[1], c[2]))
        # Verify the correct positioning of the corner and the value of the
        # resulting area
        self.assertTrue(
            are_same_shapes(
                make_vertex(lower_left),
                make_vertex(expected_lower_left),
                ShapeType.VERTEX
            )
        )
        self.assertAlmostEqual(
            get_basic_properties(result)[1],
            get_basic_properties(self.context.full_layout_shape)[1]
            / area_factor
        )


class TestBuildHexSymmetryShape(unittest.TestCase):
    """
    Test case for verifying the construction of symmetry shapes for
    hexagonal-type layouts.

    Attributes
    ----------
    centre : Tuple[float, float, float]
        The XYZ coordinates of the centre of the full layout shape.
    context : SymmetryDomain
        The dataclass storing the data used for building the symmetry shape.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for hexagonal symmetry shape
        building functions.
        """
        self.centre: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        self.context = SymmetryDomain(
            origin=self.centre,
            full_layout_shape=Hexagon(self.centre, 4)
        )

    def test_build_hex_symmetry_full(self) -> None:
        """
        Method that tests building a FULL symmetry shape for a hexagonal
        layout. It verifies that the result is a copy of the full layout
        shape.
        """
        result = build_hex_symmetry_shape(SymmetryType.FULL, self.context)
        # Verify the returned shape is the full layout
        self.assertTrue(
            are_same_shapes(
                result, self.context.full_layout_shape, ShapeType.FACE
            )
        )

    def test_build_hex_symmetry_half(self) -> None:
        """
        Method that tests building a HALF symmetry shape for a hexagonal
        layout. It verifies that the result is half of the source hexagon
        by checking the shape's lower-left corner and the area.
        """
        result = build_hex_symmetry_shape(SymmetryType.HALF, self.context)
        # Verify the returned shape is the half layout
        self.__assess_hex_symmetry(
            result,
            (
                self.centre[0],
                self.centre[1]
                - round(self.context.full_layout_shape.dimensions[1], 6),
                0.0
            ),
            2
        )

    def test_build_hex_symmetry_quarter(self) -> None:
        """
        Method that tests building a QUARTER symmetry shape for a hexagonal
        layout. It verifies that the result is a quarter of the source
        hexagon by checking the shape's lower-left corner and the area.
        """
        result = build_hex_symmetry_shape(SymmetryType.QUARTER, self.context)
        # Verify the returned shape is a quarter of the layout
        self.__assess_hex_symmetry(result, self.centre, 4)

    def test_build_hex_symmetry_third(self) -> None:
        """
        Method that tests building a THIRD symmetry shape for a hexagonal
        layout. It verifies that the result is a third of the source
        hexagon by checking the shape's lower-left corner and the area.
        """
        result = build_hex_symmetry_shape(SymmetryType.THIRD, self.context)
        # Verify the returned shape is a third of the layout
        self.__assess_hex_symmetry(
            result,
            (
                self.centre[0]
                - self.context.full_layout_shape.dimensions[0] / 2,
                self.centre[1]
                - round(self.context.full_layout_shape.dimensions[1], 6),
                0.0
            ),
            3
        )

    def test_build_hex_symmetry_sixth(self) -> None:
        """
        Method that tests building a SIXTH symmetry shape for a hexagonal
        layout. It verifies that the result is a sixth of the source
        hexagon by checking the shape's lower-left corner and the area.
        """
        result = build_hex_symmetry_shape(SymmetryType.SIXTH, self.context)
        # Verify the returned shape is a sixth of the layout
        self.__assess_hex_symmetry(
            result,
            (
                self.centre[0]
                - self.context.full_layout_shape.dimensions[0] / 2,
                self.centre[1]
                - round(self.context.full_layout_shape.dimensions[1], 6),
                0.0
            ),
            6
        )

    def test_build_hex_symmetry_twelfth(self) -> None:
        """
        Method that tests building a TWELFTH symmetry shape for a hexagonal
        layout. It verifies that the result is a twelfth of the source
        hexagon by checking the shape's lower-left corner and the area.
        """
        result = build_hex_symmetry_shape(SymmetryType.TWELFTH, self.context)
        # Verify the returned shape is a twelfth of the layout
        self.__assess_hex_symmetry(result, self.centre, 12)

    def test_build_hex_symmetry_unsupported_type(self) -> None:
        """
        Method that tests that an unsupported symmetry type raises a
        `RuntimeError`.
        """
        # Try to build with an unsupported symmetry type
        with self.assertRaises(RuntimeError):
            build_hex_symmetry_shape(SymmetryType.DIAG, self.context)

    def __assess_hex_symmetry(
            self,
            result: Surface,
            expected_lower_left: Tuple[float, float, float],
            area_factor: float
        ) -> None:
        """
        Method that assesses the correctness of a computed hexagonal symmetry
        shape by verifying that:

        - the position of the lower-left corner match the expected given
          coordinates;
        - the area of the intersection between the full hexagon and the shape
          of the symmetry matches (within a tolerance) with the area of the
          full hexagon scaled by the given factor.

        Parameters
        ----------
        result : Surface
            The geometric surface representing the shape of the analysed
            symmetry. It used to extract a portion of the hexagonal full
            layout.
        expected_lower_left : Tuple[float, float, float]
            The expected coordinates of the lower-left vertex of the portion
            of the full hexagon.
        area_factor : float
            A factor for scaling the area of the full hexagon.
        """
        # Get the XYZ coordinates of the vertices of the symmetry shape
        coords = [
            get_point_coordinates(p)
            for p in extract_sub_shapes(result, ShapeType.VERTEX)
        ]
        # Retrieve the lower-left corner
        lower_left = min(coords, key=lambda c: (c[0], c[1], c[2]))
        # Verify the correct positioning of the corner and the value of the
        # resulting area
        self.assertTrue(
            are_same_shapes(
                make_vertex(lower_left),
                make_vertex(expected_lower_left),
                ShapeType.VERTEX
            )
        )
        self.assertAlmostEqual(
            get_basic_properties(self.context.full_layout_shape * result)[1],
            get_basic_properties(self.context.full_layout_shape)[1]
            / area_factor,
            5
        )
