"""
Module containing unittest classes to assess that the classes and functions
of the `glow.interface.geom_interface` module have a valid implementation.
"""
import unittest

from glow.geometry_layouts.geometries import Circle, Rectangle
from glow.interface.geom_entities import GeomWrapper, Compound, Edge, Face, \
    Vertex, wrap_shape
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_point_coordinates, get_shape_name, make_common, make_compound, \
    make_cut, make_fuse, make_partition, make_vertex
from glow.support.utility import are_same_shapes


class TestGeomWrapperInit(unittest.TestCase):
    """
    Test case for verifying the correct initialisation of the `GeomWrapper`
    class.
    """
    def test_init_with_none(self):
        """
        Test the initialisation with `None`.
        """
        wrapper = GeomWrapper(None, [ShapeType.VERTEX])
        self.assertIsNone(wrapper.geom_obj)

    def test_init_with_valid_vertex(self):
        """
        Test the initialisation with a valid vertex object.
        """
        vertex_geom = make_vertex((0, 0, 0))
        wrapper = GeomWrapper(vertex_geom, [ShapeType.VERTEX])
        self.assertEqual(wrapper.geom_obj, vertex_geom)

    def test_init_default_name(self):
        """
        Test the initialisation sets default empty name.
        """
        vertex_geom = make_vertex((0, 0, 0))
        wrapper = GeomWrapper(vertex_geom, [ShapeType.VERTEX])
        self.assertEqual(wrapper.name, "")

    def test_init_with_invalid_type_raises_error(self):
        """
        Test the initialisation with an invalid type raises `RuntimeError`.
        """
        vertex_geom = make_vertex((0, 0, 0))
        with self.assertRaises(RuntimeError):
            GeomWrapper(vertex_geom, [ShapeType.EDGE])


class TestGeomWrapperProperties(unittest.TestCase):
    """
    Test case for verifying the correct behaviour of the getters and setters
    for the properties of the `GeomWrapper` class.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment.
        """
        self.vertex_geom = make_vertex((0, 0, 0))
        self.wrapper = GeomWrapper(self.vertex_geom, [ShapeType.VERTEX])

    def test_geom_obj_getter(self):
        """
        Test `geom_obj` property getter returns the stored GEOM object.
        """
        self.assertEqual(self.wrapper.geom_obj, self.vertex_geom)

    def test_geom_obj_setter_valid(self):
        """
        Test setting `geom_obj` to a valid GEOM object of the expected type.
        """
        new_vertex = make_vertex((1, 2, 3))
        self.wrapper.geom_obj = new_vertex
        self.assertEqual(self.wrapper.geom_obj, new_vertex)

    def test_geom_obj_setter_none(self):
        """
        Test setting `geom_obj` to `None`.
        """
        self.wrapper.geom_obj = None
        self.assertIsNone(self.wrapper.geom_obj)

    def test_geom_obj_setter_invalid_type(self):
        """
        Test setting `geom_obj` to an invalid type raises `RuntimeError`.
        """
        circle = Circle()
        with self.assertRaises(RuntimeError):
            self.wrapper.geom_obj = circle.geom_obj

    def test_name_getter_default(self) -> None:
        """
        Test `name` property getter returns empty string by default.
        """
        self.assertEqual(self.wrapper.name, "")

    def test_name_setter(self) -> None:
        """
        Test setting `name` property.
        """
        self.wrapper.name = "Geometry"
        self.assertEqual(self.wrapper.name, "Geometry")

    def test_name_setter_updates_geom_object(self):
        """
        Test setting `name` property also updates the name of the wrapped
        GEOM object.
        """
        self.wrapper.name = "MyGeometry"
        self.assertEqual(self.wrapper.name, "MyGeometry")
        # Verify the name is actually set on the GEOM object
        self.assertEqual(get_shape_name(self.wrapper.geom_obj), "MyGeometry")

    def test_name_setter_with_none_geom_obj(self):
        """
        Test setting `name` property when `geom_obj` is `None`.
        """
        self.wrapper.geom_obj = None
        self.wrapper.name = "NoGeom"
        self.assertEqual(self.wrapper.name, "NoGeom")

    def test_name_setter_empty_string(self):
        """
        Test setting `name` property to empty string.
        """
        self.wrapper.name = "TestName"
        self.wrapper.name = ""
        self.assertEqual(self.wrapper.name, "")

    def test_name_persistence(self):
        """
        Test `name` property persists after updating the `geom_obj` instance
        attribute.
        """
        self.wrapper.name = "Name"
        new_vertex = make_vertex((1, 1, 1))
        self.wrapper.geom_obj = new_vertex
        self.assertEqual(self.wrapper.name, "Name")


class TestGeomWrapperOperators(unittest.TestCase):
    """
    Test case for verifying the correct behaviour of the arithmetic operators
    of the `GeomWrapper` class.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment.
        """
        self.shape_1 = Circle(radius=2)
        self.shape_2 = Circle((1.0, 1.0, 0.0))
        self.wrapper_1 = Face(self.shape_1.geom_obj)
        self.wrapper_2 = Face(self.shape_2.geom_obj)

    def test_add_single_operand(self):
        """
        Test the `__add__` method (`+` operator) with a single operand.
        """
        result = self.wrapper_1 + self.wrapper_2
        cmpr_shape = extract_sub_shapes(
            make_fuse([self.shape_1.geom_obj, self.shape_2.geom_obj]),
            ShapeType.FACE
        )[0]
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(cmpr_shape, result, ShapeType.FACE),
            ShapeType.FACE
        )

    def test_add_sequence_operands(self):
        """
        Test the `__add__` method (`+` operator) with a sequence of operands.
        """
        result = self.wrapper_1 + [self.wrapper_2]
        cmpr_shape = extract_sub_shapes(
            make_fuse([self.shape_1.geom_obj, self.shape_2.geom_obj]),
            ShapeType.FACE
        )[0]
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(cmpr_shape, result, ShapeType.FACE)
        )

    def test_add_returns_new_instance(self):
        """
        Test the `__add__` method (`+` operator) returns a new instance
        without modifying operands.
        """
        result = self.wrapper_1 + self.wrapper_2
        self.assertIsNot(result, self.wrapper_1)
        self.assertIsNot(result, self.wrapper_2)
        self.assertIsNotNone(self.wrapper_1.geom_obj)
        self.assertIsNotNone(self.wrapper_2.geom_obj)

    def test_sub_single_operand(self):
        """
        Test the `__sub__` method (`-` operator) with a single operand.
        """
        result = self.wrapper_1 - self.wrapper_2
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(
                make_cut(self.shape_1, self.shape_2),
                result,
                ShapeType.FACE
            )
        )

    def test_sub_returns_new_instance(self):
        """
        Test the `__sub__` method (`-` operator) returns a new instance
        without modifying operands.
        """
        result = self.wrapper_1 - self.wrapper_2
        self.assertIsNot(result, self.wrapper_1)
        self.assertIsNot(result, self.wrapper_2)
        self.assertIsNotNone(self.wrapper_1.geom_obj)
        self.assertIsNotNone(self.wrapper_2.geom_obj)

    def test_mul_single_operand(self):
        """
        Test the `__mul__` method (`*` operator) with a single operand.
        """
        result = self.wrapper_1 * self.wrapper_2
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(
                make_common(self.shape_1, self.shape_2),
                result,
                ShapeType.FACE
            )
        )

    def test_mul_returns_new_instance(self):
        """
        Test the `__mul__` method (`*` operator) returns a new instance
        without modifying operands.
        """
        result = self.wrapper_1 * self.wrapper_2
        self.assertIsNot(result, self.wrapper_1)
        self.assertIsNot(result, self.wrapper_2)
        self.assertIsNotNone(self.wrapper_1.geom_obj)
        self.assertIsNotNone(self.wrapper_2.geom_obj)

    def test_truediv_single_operand(self):
        """
        Test the `__truediv__` method (`/` operator) with a single operand.
        """
        result = self.wrapper_1 / self.wrapper_2
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(
                make_partition(
                    [self.shape_1], [self.shape_2], ShapeType.FACE
                ),
                result,
                ShapeType.FACE
            )
        )

    def test_truediv_sequence_operands(self):
        """
        Test the `__truediv__` method (`/` operator) with a sequence of
        operands.
        """
        circle_medium = Circle(radius=1.5)
        wrapper_medium = wrap_shape(circle_medium.geom_obj)
        result = self.wrapper_1 / [self.wrapper_2, wrapper_medium]
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(
                make_partition(
                    [self.shape_1],
                    [self.shape_2, wrapper_medium],
                    ShapeType.FACE
                ),
                result,
                ShapeType.FACE
            )
        )

    def test_truediv_returns_new_instance(self):
        """
        Test the `__truediv__` method (`/` operator) returns a new instance
        without modifying operands.
        """
        result = self.wrapper_1 / self.wrapper_2
        self.assertIsNot(result, self.wrapper_1)
        self.assertIsNot(result, self.wrapper_2)
        self.assertIsNotNone(self.wrapper_1.geom_obj)
        self.assertIsNotNone(self.wrapper_2.geom_obj)

    def test_floordiv_single_operand(self):
        """
        Test the `__floordiv__` method (`//` operator) with a single operand.
        """
        result = self.wrapper_1 // self.wrapper_2
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(
                make_partition(
                    [self.shape_1, self.shape_2], [], ShapeType.FACE
                ),
                result,
                ShapeType.FACE
            )
        )

    def test_floordiv_sequence_operands(self):
        """
        Test the `__floordiv__` method (`//` operator) with a sequence of
        operands.
        """
        circle_medium = Circle(radius=1.5)
        wrapper_medium = wrap_shape(circle_medium.geom_obj)
        result = self.wrapper_1 // [self.wrapper_2, wrapper_medium]
        self.assertIsInstance(result, GeomWrapper)
        self.assertIsNotNone(result.geom_obj)
        self.assertTrue(
            are_same_shapes(
                make_partition(
                    [self.shape_1, self.shape_2, wrapper_medium],
                    [],
                    ShapeType.FACE
                ),
                result,
                ShapeType.FACE
            )
        )

    def test_floordiv_returns_new_instance(self):
        """
        Test the `__floordiv__` method (`//` operator) returns a new instance
        without modifying operands.
        """
        result = self.wrapper_1 // self.wrapper_2
        self.assertIsNot(result, self.wrapper_1)
        self.assertIsNot(result, self.wrapper_2)
        self.assertIsNotNone(self.wrapper_1.geom_obj)
        self.assertIsNotNone(self.wrapper_2.geom_obj)


class TestGeomWrapperGetattr(unittest.TestCase):
    """
    Test case for verifying the correct behaviour of the `__getattr__`
    delegation.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment.
        """
        self.vertex_geom = make_vertex((1, 2, 3))
        self.wrapper = GeomWrapper(self.vertex_geom, [ShapeType.VERTEX])

    def test_getattr_delegates_to_geom_obj(self):
        """
        Test the `__getattr__` delegates attribute access to the wrapped GEOM
        object.
        """
        # Get the vertex coordinates by accessing to wrapped GEOM Object
        # through the wrapper class instance
        coords = get_point_coordinates(self.wrapper)
        self.assertEqual(coords, (1, 2, 3))

    def test_getattr_raises_without_geom_obj(self):
        """
        Test the `__getattr__` raises `AttributeError` when trying to access
        attributes before the `_geom_obj` attribute is set.
        """
        wrapper = GeomWrapper.__new__(GeomWrapper)
        with self.assertRaises(AttributeError):
            _ = wrapper.some_nonexistent_attribute

        with self.assertRaises(AttributeError):
            print(self.wrapper.some_nonexistent_attribute)


class TestGeomWrapperRepr(unittest.TestCase):
    """
    Test case for verifying the correct behaviour of the `__repr__` method.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment.
        """
        self.vertex_geom = make_vertex((0, 0, 0))
        self.wrapper = GeomWrapper(self.vertex_geom, [ShapeType.VERTEX])

    def test_repr_contains_class_name(self):
        """
        Test the `__repr__` method contains the class name.
        """
        repr_str = repr(self.wrapper)
        self.assertIn("GeomWrapper", repr_str)

    def test_repr_contains_shape_name(self):
        """
        Test the `__repr__` method contains the shape name of the GEOM object.
        """
        self.wrapper.name = "TestShape"
        self.assertIn("TestShape", repr(self.wrapper))

    def test_repr_with_unnamed_shape(self):
        """
        Test the `__repr__` method handles unnamed shapes gracefully.
        """
        repr_str = repr(self.wrapper)
        self.assertIn("GeomWrapper", repr_str)
        # Should contain either the name or "Unnamed" fallback
        self.assertTrue(
            "Unnamed" in repr_str or "GeomWrapper" in repr_str
        )

    def test_repr_format(self):
        """
        Test the `__repr__` method returns properly formatted string.
        """
        self.wrapper.name = "GeometryShape"
        repr_str = repr(self.wrapper)
        self.assertIn("<", repr_str)
        self.assertIn(">", repr_str)
        self.assertIn(":", repr_str)


class TestWrapShape(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the function
    `wrap_shape` in the `geom_entities.py` module.
    """
    def test_wrap_vertex(self):
        """
        Test wrapping a vertex shape.
        """
        vertex_geom = make_vertex((0, 0, 0))
        result = wrap_shape(vertex_geom)
        self.assertIsInstance(result, Vertex)
        self.assertEqual(result.geom_obj, vertex_geom)

    def test_wrap_edge(self):
        """
        Test wrapping an edge shape.
        """
        circle = Circle()
        edge_geom = circle.geom_obj
        edges = extract_sub_shapes(edge_geom, ShapeType.EDGE)
        if edges:
            result = wrap_shape(edges[0])
            self.assertIsInstance(result, Edge)
            self.assertEqual(result.geom_obj, edges[0])

    def test_wrap_face(self):
        """
        Test wrapping a face shape.
        """
        rectangle = Rectangle()
        face_geom = rectangle.geom_obj
        faces = extract_sub_shapes(face_geom, ShapeType.FACE)
        if faces:
            result = wrap_shape(faces[0])
            self.assertIsInstance(result, Face)
            self.assertEqual(result.geom_obj, faces[0])

    def test_wrap_compound(self):
        """
        Test wrapping a compound shape.
        """
        compound_geom = make_compound([Circle(), Circle(radius=2)])
        result = wrap_shape(compound_geom)
        self.assertIsInstance(result, Compound)
        self.assertEqual(result.geom_obj, compound_geom)

    def test_wrap_invalid_object(self):
        """
        Test wrapping invalid object raises `TypeError`.
        """
        invalid_obj = "not a geom object"
        with self.assertRaises(TypeError):
            wrap_shape(invalid_obj)
