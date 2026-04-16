"""
Module containing unittest classes to assess that the functions of the
`glow.support.utility` module have a valid implementation.
"""
import math
import unittest

from typing import Sequence

from glow.geometry_layouts.geometries import Circle, Rectangle, Surface
from glow.interface.geom_entities import GeomWrapper, wrap_shape
from glow.support.utility import *
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_min_distance, make_cdg, make_circle, make_compound, make_face, \
    make_partition, make_vertex, set_shape_name


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
                ShapeType.FACE
            )
        # Verify the returned value matches the shapes to compare
        self.assertFalse(are_same_shapes(shape1, shape2, ShapeType.FACE))
        self.assertTrue(are_same_shapes(shape1, shape3, ShapeType.FACE))
        self.assertTrue(
            are_same_shapes(
                make_compound(extract_sub_shapes(shape1, ShapeType.EDGE)),
                make_compound(extract_sub_shapes(shape3, ShapeType.EDGE)),
                ShapeType.EDGE
            )
        )

    def test_build_arcs_for_rounded_corners(self) -> None:
        """
        Method that tests the implementation of the function
        `build_arcs_for_rounded_corners` declared in the `utility.py` module.
        """
        # Verify an exception is raised when any of the radii of the corners
        # is greater than the maximum allowed radius
        with self.assertRaises(RuntimeError):
            build_arcs_for_rounded_corners(
                [(0, 1), (1, 2.5)], (0.0, 0.0, 0.0), 1, 1
            )
        # Verify the function builds the correct arc edges
        rounded_corners = [(0, 0.5), (1, 0.75)]
        width = 3
        height = 3
        corner_arcs = [
            make_arc_center(
                make_vertex((-1.0, -1.0, 0.0)),
                make_vertex((-width/2, -1.0, 0.0)),
                make_vertex((-1.0, -height/2, 0.0))
            ),
            make_arc_center(
                make_vertex((0.75, -0.75, 0.0)),
                make_vertex((0.75, -height/2, 0.0)),
                make_vertex((width/2, -0.75, 0.0))
            )
        ]
        arcs = build_arcs_for_rounded_corners(
            rounded_corners, (0.0, 0.0, 0.0), height, width
        )
        self.assertEqual(len(arcs), len(corner_arcs))
        for a1, a2 in zip(arcs, corner_arcs):
            self.assertTrue(are_same_shapes(a1, a2, ShapeType.EDGE))


    def test_build_compound_borders(self) -> None:
        """
        Method that tests the implementation of the function
        `build_compound_borders` declared in the `utility.py` module.
        """
        # Declare the reference shape
        outer_shape = Rectangle(height=2, width=2)
        shape = make_partition(
            [outer_shape],
            [
                Rectangle(height=1, width=1),
                make_circle(make_vertex(self.center), None, 0.4),
                make_circle(make_vertex(self.center), None, 0.2)
            ],
            ShapeType.FACE
        )
        # Get the borders of the compound
        self.__assess_borders(outer_shape, wrap_shape(shape))

        # Build a compound whose borders are cut by circles in the corners
        shape = make_partition(
            [shape],
            [
                Circle((-1, 1, 0), 0.5),
                Circle((-1, -1, 0), 0.5),
                Circle((1, 1, 0), 0.5),
                Circle((1, -1, 0), 0.5)
            ],
            ShapeType.FACE
        )
        # Verify the borders are the same as the previous case
        self.__assess_borders(outer_shape, wrap_shape(shape))

    def test_build_contiguous_edges(self) -> None:
        """
        Method that tests the implementation of the function
        `build_contiguous_edges` declared in the `utility.py` module.
        """
        # Build a list of vertices
        vertices = [
            make_vertex(coords) for coords in [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 0.5, 0.0)
            ]
        ]
        # Build the contiguous edges
        edges = build_contiguous_edges(vertices)
        # Check the first vertex of the first edge coincides with the second
        # vertex of the last edge
        self.assertTrue(
            are_same_shapes(
                extract_sub_shapes(edges[0], ShapeType.VERTEX)[0],
                extract_sub_shapes(edges[-1], ShapeType.VERTEX)[1],
                ShapeType.VERTEX
            )
        )
        # Verify an exception is raised if providing less than 2 vertices or
        # if two consecutive vertices coincide
        with self.assertRaises(RuntimeError):
            _ = build_contiguous_edges(vertices[:2])
        with self.assertRaises(RuntimeError):
            # Add the same vertex two times consecutively
            vertices.insert(1, vertices[0])
            _ = build_contiguous_edges(vertices[:2])

    def test_build_z_axis_from_vertex(self) -> None:
        """
        Method that tests the implementation of the function
        `build_z_axis_from_vertex` declared in the `utility.py` module.
        """
        # Declare the vertices of the axis
        o = (1, 1, 0.0)
        axis_vertices = [make_vertex(o), make_vertex((*o[:2], 1.0))]
        # Build the Z-axis from the origin point
        z_axis = build_z_axis_from_vertex(make_vertex(o))
        # Verify the axis vertices are the expected ones
        z_axis_vertices = extract_sub_shapes(z_axis, ShapeType.VERTEX)
        self.__assess_shapes_in_list_coincide(
            axis_vertices, z_axis_vertices, ShapeType.VERTEX
        )

    def test_build_subdvision_vertices_on_edge(self) -> None:
        """
        Method that tests the implementation of the function
        `build_subdvision_vertices_on_edge` declared in the `utility.py`
        module.
        """
        # Declare the subdivision points for an edge parallel to the X-axis
        p1 = make_vertex((0.0, 0.0, 0.0))
        p2 = make_vertex((2.0, 0.0, 0.0))
        subdiv_pnts = [
            p1,
            make_vertex((0.5, 0.0, 0.0)),
            make_vertex((1.0, 0.0, 0.0)),
            make_vertex((1.5, 0.0, 0.0)),
        ]
        edge = make_edge(p1, p2)
        # Build the subdivision points on the edge
        vertices = build_subdvision_vertices_on_edge(4, edge)
        # Verify the two list of vertices coincides
        self.__assess_shapes_in_list_coincide(
            vertices, subdiv_pnts, ShapeType.VERTEX
        )

        # Build the subdivision points on the edge from a starting position
        vertices = build_subdvision_vertices_on_edge(4, edge, 0.25)
        subdiv_pnts = [
            make_vertex((0.5, 0.0, 0.0)),
            make_vertex((0.875, 0.0, 0.0)),
            make_vertex((1.25, 0.0, 0.0)),
            make_vertex((1.625, 0.0, 0.0)),
        ]
        # Verify the two list of vertices coincides
        self.__assess_shapes_in_list_coincide(
            vertices, subdiv_pnts, ShapeType.VERTEX
        )
        # Verify an exception is raised if providing a value outside the
        # admitted range
        with self.assertRaises(ValueError):
            _ = build_subdvision_vertices_on_edge(4, edge, 1.0)
            _ = build_subdvision_vertices_on_edge(4, edge, 1.25)


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
            [Rectangle(height=2, width=2)],
            [
                Rectangle(height=1, width=1),
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

    def test_check_type_geo_consistency(self) -> None:
        """
        Method that tests the implementation of the function
        `check_type_geo_consistency` declared in the `utility.py` module.
        """
        # Verify the exception is raised when a wrong combination of values
        # is provided
        with self.assertRaises(RuntimeError):
            check_type_geo_consistency(
                LayoutGeometryType.HEXAGON_TRAN,
                LayoutType.RECT,
                SymmetryType.DIAG
            )
        # Verify the correct combination of values does not raise an exception
        try:
            check_type_geo_consistency(
                LayoutGeometryType.ROTATION,
                LayoutType.HEX,
                SymmetryType.THIRD
            )
        except:
            self.fail("Test failed as an exception was raised.")

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
                math.sqrt(
                    sum((xyz2 - xyz1)*(xyz2 - xyz1)
                        for xyz1, xyz2 in zip(ref_coords2, coords_b2))
                ),
                get_min_distance(ref_point, point_b)
            )
        )

    def test_flatten_list(self) -> None:
        """
        Method that tests the implementation of the function
        `flatten_list` declared in the `utility.py` module.
        """
        # Verify the function against a flat list
        self.assertEqual(list(flatten_list([1, 2, 3])), [1, 2, 3])
        # Verify the function against a nested list
        self.assertEqual(
            list(flatten_list([1, [2, [3, [4]]]])), [1, 2, 3, 4]
        )
        # Verify the function against a nested list with None values
        self.assertEqual(
            list(flatten_list([1, None, [2, None, 3], None])),
            [1, 2, 3]
        )
        # Verify the function against a nested list with mixed type values
        self.assertEqual(
            list(flatten_list([1, 'a', [True, [3.14, ['x']]]])),
            [1, 'a', True, 3.14, 'x']
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

    def test_get_angle_between_points(self) -> None:
        """
        Method that tests the implementation of the function
        `get_angle_between_points` declared in the `utility.py` module.
        """
        # Verify the function outputs the correct angles
        self.assertAlmostEqual(
            get_angle_between_points((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
            0.0
        )
        self.assertAlmostEqual(
            get_angle_between_points((0.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
            math.pi/4
        )
        self.assertAlmostEqual(
            get_angle_between_points((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0)),
            math.pi
        )
        self.assertAlmostEqual(
            get_angle_between_points((0.0, 0.0, 0.0), (1.0, 1.0, 0.0), True),
            45.0
        )

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

    def test_get_id_from_shape(self) -> None:
        """
        Method that tests the implementation of the function
        `get_id_from_shape` declared in the `utility.py` module.
        """
        # Verify the exception is raised if providing a shape whose name is
        # in an incorrect format
        shape = make_circle(make_vertex(self.center), None, 1.0)
        set_shape_name(shape, "Cell")
        with self.assertRaises(RuntimeError):
            get_id_from_shape(shape)
        set_shape_name(shape, "Cell1")
        with self.assertRaises(RuntimeError):
            get_id_from_shape(shape)
        # Verify the ID index can be retrieved
        set_shape_name(shape, "Cell_01")
        self.assertEqual(get_id_from_shape(shape), 1)
        set_shape_name(shape, "Cell_001")
        self.assertEqual(get_id_from_shape(shape), 1)

    def test_get_vertex_polar_position(self) -> None:
        """
        Method that tests the implementation of the function
        `get_vertex_polar_position` declared in the `utility.py` module.
        """
        # Verify the returned angle and distance is (-pi, 0.0) for a point
        # that coincides with the origin
        self.assertEqual(
            get_vertex_polar_position(
                make_vertex((0.0, 0.0, 0.0)), make_vertex((0.0, 0.0, 0.0))
            ),
            (-math.pi, 0.0)
        )
        # Verify angle and distance for CCW ordering
        self.__assess_values_in_sequences(
            get_vertex_polar_position(
                make_vertex((0.0, 1.0, 0.0)),
                make_vertex((0.0, 0.0, 0.0)),
                ref_vec=(1, 0),
                is_cw=False
            ),
            (math.pi/2, 1.0)
        )
        # Verify angle and distance for CW ordering
        self.__assess_values_in_sequences(
            get_vertex_polar_position(
                make_vertex((0.0, 1.0, 0.0)),
                make_vertex((0.0, 0.0, 0.0)),
                ref_vec=(1, 0),
                is_cw=True
            ),
            (3*math.pi/2, 1.0)
        )
        # Verify the exception is raised when indicating (0, 0) the reference
        # vector
        with self.assertRaises(ValueError):
            _ = get_vertex_polar_position(
                make_vertex((0.0, 1.0, 0.0)),
                make_vertex((0.0, 0.0, 0.0)),
                ref_vec=(0, 0),
                is_cw=True
            )

    def test_get_vertices_on_edges(self) -> None:
        """
        Method that tests the implementation of the function
        `get_vertices_on_edges` declared in the `utility.py` module.
        """
        ref_vertex = make_vertex(self.center)
        # Define the reference edges
        ref_edges = Rectangle().borders
        # Define vertices on the reference edges
        vertices_on_ref_edges = [
            make_vertex((0.5, 0.207107, 0.0)),
            make_vertex((0.207107, 0.5, 0.0)),
            make_vertex((-0.207107, 0.5, 0.0)),
            make_vertex((-0.5, 0.207107, 0.0)),
            make_vertex((-0.5, -0.207107, 0.0)),
            make_vertex((-0.207107, -0.5, 0.0)),
            make_vertex((0.207107, -0.5, 0.0)),
            make_vertex((0.5, -0.207107, 0.0)),
        ]
        result = get_vertices_on_edges(
            build_contiguous_edges(vertices_on_ref_edges),
            ref_edges,
            ref_vertex
        )
        # Verify the returned vertices coincides with those on the reference
        # edges
        self.__assess_shapes_in_list_coincide(
            result, vertices_on_ref_edges, ShapeType.VERTEX
        )
        # Verify the returned list is empy if the edges to check do not have
        # vertices on the reference ones
        rect = Rectangle()
        rect.rotate(45)
        result = get_vertices_on_edges(
            rect.borders, ref_edges, ref_vertex
        )
        self.__assess_shapes_in_list_coincide(result, [], ShapeType.VERTEX)


    def test_retrieve_selected_object(self) -> None:
        """
        Method that tests the implementation of the function
        `retrieve_selected_object` declared in the `utility.py` module.
        """
        # Verify the exception is raised as the function is called from
        # outside the SALOME GUI
        with self.assertRaises(RuntimeError):
            _ = retrieve_selected_object("")

    def test_sort_shapes_from_vertex(self) -> None:
        """
        Method that tests the implementation of the function
        `sort_shapes_from_vertex` declared in the `utility.py` module.
        """
        o = make_vertex(self.center)
        # Declare shapes at different distances from the XYZ origin
        shapes = [
            Circle((1.5, 0.0, 0.0)),
            Circle((1.25, 1.5, 0.0)),
            Circle((2.0, 2.0, 0.0)),
            Rectangle((0.0, 0.0, 0.0))
        ]
        ascending_ordering = [
            shapes[3],
            shapes[0],
            shapes[1],
            shapes[2]
        ]
        descending_ordering = ascending_ordering[::-1]
        # Verify the returned list is in ascending order
        self.__assess_shapes_in_list_coincide(
            sort_shapes_from_vertex(shapes, o),
            ascending_ordering,
            ShapeType.FACE
        )
        # Verify the returned list is in descending order
        self.__assess_shapes_in_list_coincide(
            sort_shapes_from_vertex(shapes, o, True),
            descending_ordering,
            ShapeType.FACE
        )

    def test_translate_wrt_reference(self) -> None:
        """
        Method that tests the implementation of the function
        `translate_wrt_reference` declared in the `utility.py` module.
        """
        # Declare the shape to test; it can be seen as a portion of a 1x1
        # rectangle centered in the XYZ origin
        shape = Rectangle(
            center=(0.75, 0.75, 0.0), height=0.25, width=0.25
        )
        ref_pnt = make_vertex(self.center)
        new_ref_coords = (1.0, 1.0, 0.0)
        translated_shape = translate_wrt_reference(
            shape, ref_pnt, new_ref_coords
        )
        # Verify the shape has been traslated so that its relative distance
        # from the translated reference point is kept
        self.assertTrue(
            math.isclose(
                get_min_distance(
                    make_vertex(new_ref_coords), make_cdg(translated_shape)
                ),
                get_min_distance(ref_pnt, make_cdg(shape))
            )
        )

    def __assess_borders(self, outer_shape: Surface, shape: GeomWrapper):
        """
        Method that validates that all borders extracted from the given shape
        are present in the surface which represents the characteristic shape.
        This method computes borders for `shape` using
        ``build_compound_borders``.
        For each extracted border, it checks whether an equivalent border
        exists among the attribute `outer_shape.borders`, using the function
        ``are_same_shapes``. The test fails if any border is not found.

        Parameters
        ----------
        outer_shape : Any
            The ``Surface`` object expected to contain all borders.
        shape : GeomWrapper
            The shape from which borders will be extracted.
        """
        borders = build_compound_borders(shape)
        # Verify the correct borders extraction
        for border in borders:
            found = False
            for b in outer_shape.borders:
                if are_same_shapes(border, b, ShapeType.EDGE):
                    found = True
                    break
            self.assertTrue(found)

    def __assess_shapes_in_list_coincide(
            self,
            shapes_1: List[Any],
            shapes_2: List[Any],
            shapes_type: ShapeType
        ) -> None:
        """
        Method that validates that corresponding shapes in the two given lists
        of GEOM objects, whose type is provided as third parameter, coincide.

        Parameters
        ----------
        shapes_1 : List[Any]
            First list of GEOM shapes.
        shapes_2 : list of Any
            Second list of GEOM shapes.
        shapes_type : ShapeType
            The specific shape type to use when determining shape equivalence.
        """
        # Verify the two lists have the same size
        self.assertEqual(len(shapes_1), len(shapes_2))
        # Verify corresponding shapes of the two lists are the same
        for s_1, s_2 in zip(shapes_1, shapes_2):
            self.assertTrue(are_same_shapes(s_1, s_2, shapes_type))

    def __assess_values_in_sequences(
            self, seq_1: Sequence, seq_2: Sequence
        ) -> None:
        """
        Method that validats that the corresponding values in the two given
        sequences are the same.

        Parameters
        ----------
        seq_1 : Sequence
            First sequence of values.
        seq_2 : Sequence
            Second sequence of values.
        """
        # Verify the two sequences have the same size
        self.assertEqual(len(seq_1), len(seq_2))
        # Verify corresponding values of the two sequences are the same
        for s_1, s_2 in zip(seq_1, seq_2):
            self.assertTrue(math.isclose(s_1, s_2, abs_tol=1e-6))
