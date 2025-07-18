"""
Module containing unittest classes to assess that the functions of the
`glow.geometry_layouts.utility` module have a valid implementation.
"""
import math
import unittest

from glow.geometry_layouts.geometries import Rectangle
from glow.support.utility import *
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_min_distance, make_cdg, make_circle, make_compound, make_face, \
    make_partition, make_vertex


class TestUtilityFunctions(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the functions
    declared in the `utility.py` module.

    Attributes
    ----------
    center : Tuple[float, float, float]
        A tuple representing the XYZ coordinates of the origin.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the functions declared
        in the `utility.py` module.
        It initializes the attributes common to all the tests.
        """
        self.center = (0.0, 0.0, 0.0)

    def test_are_same_shapes(self) -> None:
        """
        Method that tests the implementation of the function `are_same_shapes`
        declared in the `utility.py` module.
        """
        # Declare the shapes to compare
        shape1 = make_face(
            [make_circle(make_vertex(self.center), None, 1.0)]
        )
        shape2 = make_face(
            [make_circle(make_vertex(self.center), None, 2.0)]
        )
        shape3 = make_face(
            [make_circle(make_vertex(self.center), None, 1.0)]
        )
        # Verify an exception is raised when comparing shapes having different
        # types
        with self.assertRaises(RuntimeError):
            are_same_shapes(
                shape1,
                make_compound(extract_sub_shapes(shape2, ShapeType.EDGE)),
                ShapeType.FACE)
        # Verify the returned value matches the shapes to compare
        self.assertFalse(are_same_shapes(shape1, shape2, ShapeType.FACE))
        self.assertTrue(are_same_shapes(shape1, shape3, ShapeType.FACE))
        self.assertTrue(
            are_same_shapes(
                make_compound(extract_sub_shapes(shape1, ShapeType.EDGE)),
                make_compound(extract_sub_shapes(shape3, ShapeType.EDGE)),
                ShapeType.EDGE)
        )

    def test_build_compound_borders(self) -> None:
        """
        Method that tests the implementation of the function
        `build_compound_borders` declared in the `utility.py` module.
        """
        # Declare the reference shape
        outer_shape = Rectangle(height=2, width=2)
        shape = make_partition(
            [outer_shape.face],
            [
                Rectangle(height=1, width=1).face,
                make_circle(make_vertex(self.center), None, 0.4),
                make_circle(make_vertex(self.center), None, 0.2)
            ],
            ShapeType.FACE
        )
        # Get the borders of the compound
        borders = build_compound_borders(shape)
        # Verify the correct borders extraction
        for border in borders:
            found = False
            for b in outer_shape.borders:
                if are_same_shapes(border, b, ShapeType.EDGE):
                    found = True
                    break
            self.assertTrue(found)

    def test_check_shape_expected_types(self) -> None:
        """
        Method that tests the implementation of the function
        `check_shape_expected_types` declared in the `utility.py` module.
        """
        # Build shapes of different types
        vrtx_shape = make_vertex(self.center)
        edge_shape = make_circle(make_vertex(self.center), None, 0.4)
        face_shape = make_face([edge_shape])
        cmpd_shape = make_partition(
            [Rectangle(height=2, width=2).face],
            [
                Rectangle(height=1, width=1).face,
                edge_shape,
                make_circle(make_vertex(self.center), None, 0.2)
            ],
            ShapeType.FACE
        )
        # Verify the types of the shapes
        try:
            check_shape_expected_types(face_shape, [ShapeType.FACE])
            check_shape_expected_types(cmpd_shape, [ShapeType.COMPOUND])
            check_shape_expected_types(edge_shape, [ShapeType.EDGE])
            check_shape_expected_types(vrtx_shape, [ShapeType.VERTEX])
        except:
            self.fail("Test failed as an exception was raised.")
        # Verify the exception is raised if the shape does not have any of
        # the indicated types
        with self.assertRaises(RuntimeError):
            check_shape_expected_types(face_shape, [ShapeType.EDGE])

    def test_compute_point_by_reference(self) -> None:
        """
        Method that tests the implementation of the function
        `compute_point_by_reference` declared in the `utility.py` module.
        """
        # Declare the points
        ref_point = make_vertex(self.center)
        point_b = make_vertex((1.0, 1.0, 0.0))
        ref_coords2 = (0.5, 0.5, 0.0)
        # Calculate the new coordinates of the point
        coords_b2 = compute_point_by_reference(
            point_b, ref_point, ref_coords2)
        # Verify the new point has kept the same distance from the moved
        # reference point
        self.assertTrue(
            math.isclose(
                math.sqrt(sum((xyz2 - xyz1)*(xyz2 - xyz1)
                              for xyz1, xyz2 in zip(ref_coords2, coords_b2))),
                get_min_distance(ref_point, point_b)
            )
        )

    def test_generate_unique_random_colors(self) -> None:
        """
        Method that tests the implementation of the function
        `generate_unique_random_colors` declared in the `utility.py`
        module.
        """
        # Generate the indicated number of colors
        no_colors = 10
        colors = generate_unique_random_colors(no_colors)
        # Verify the correct number of colors have been generated and that
        # they are different
        self.assertEqual(len(colors), no_colors)
        self.assertEqual(len(colors), len(set(colors)))
        # Verify the colors have been extracted from the list of available
        # colors
        for color in colors:
            self.assertIn(color, RGB_COLORS)
        # Verify an exception is raised when requesting more than the
        # available number of colors
        with self.assertRaises(RuntimeError):
            generate_unique_random_colors(len(RGB_COLORS) + 1)

    def test_get_id_from_name(self) -> None:
        """
        Method that tests the implementation of the function
        `get_id_from_name` declared in the `utility.py` module.
        """
        # Verify the exception is raised if providing a name in an incorrect
        # format
        with self.assertRaises(RuntimeError):
            get_id_from_name("Cell")
        with self.assertRaises(RuntimeError):
            get_id_from_name("Cell1")
        # Verify the ID index can be retrieved
        self.assertEqual(get_id_from_name("Cell_01"), 1)
        self.assertEqual(get_id_from_name("Cell_001"), 1)

    def test_retrieve_selected_object(self) -> None:
        """
        Method that tests the implementation of the function
        `retrieve_selected_object` declared in the `utility.py` module.
        """
        # Verify the exception is raised as the function is called from
        # outside the SALOME GUI
        with self.assertRaises(RuntimeError):
            retrieve_selected_object("")

    def test_translate_wrt_reference(self) -> None:
        """
        Method that tests the implementation of the function
        `translate_wrt_reference` declared in the `utility.py` module.
        """
        # Declare the shape to test; it can be seen as a portion of a 1x1
        # rectangle centered in the XYZ origin
        shape = Rectangle(center=(0.75, 0.75, 0.0),
                          height=0.25,
                          width=0.25).face
        ref_pnt = make_vertex(self.center)
        new_ref_coords = (1.0, 1.0, 0.0)
        translated_shape = translate_wrt_reference(
            shape, ref_pnt, new_ref_coords)
        # Verify the shape has been traslated so that its relative distance
        # from the translated reference point is kept
        self.assertTrue(
            math.isclose(
                get_min_distance(
                    make_vertex(new_ref_coords), make_cdg(translated_shape)),
                get_min_distance(ref_pnt, make_cdg(shape))
            )
        )
