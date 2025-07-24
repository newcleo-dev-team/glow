"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.geometries` module have a valid implementation.
"""
import math
import unittest

from abc import ABC
from copy import deepcopy
from typing import Callable

from glow.support.utility import are_same_shapes
from glow.interface.geom_interface import *
from glow.geometry_layouts.geometries import GenericSurface, Hexagon, \
    Rectangle, Surface, Circle, build_hexagon
from support_funcs import build_hex_geom_elements, make_ref_vectors


class TestSurface(ABC, unittest.TestCase):
    """
    Abstract base test case for verifying the geometric operations and
    visualization capabilities of the `Surface` subclasses.

    This test suite provides common setup and a set of tests to ensure that
    implementations of `Surface` correctly handle rotation, translation, and
    visualization of their geometric elements (i.e. faces, borders, and
    vertices).
    Tests methods declared herein are skipped as no `Surface` subclass
    is instantiated.

    Attributes
    ----------
    normal_vect : Tuple[float, float, float]
        The normal vector to the surface, defaulting to (0, 0, 1).
    o : Any
        The origin vertex of the surface.
    face : Any | None
        The face associated with the surface, defined in subclasses instances.
    name : str
        The name of the surface.
    surf : Surface | None
        The surface instance under test, defined in subclasses instances as
        an object of the `Surface` subclasses.
    rotation : float
        The initial rotation angle of the surface, in degrees.
    rotation_angle : float
        The angle by which to rotate the surface during tests.
    z_axis : Any
        The axis vector used for rotation tests.
    center_after_transl : Tuple[float, float, float]
        The target center position of the surface after translation.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the subclasses of the
        `Surface` class.
        It initializes the common geometric characteristics for a SALOME
        surface while setting to `None` the ones specific to the geometric
        surface being tested.
        """
        # Setup the common geometric characteristics used for testing the
        # subclasses of 'Surface'
        center: Tuple[float, float, float] = (0, 0, 0)
        self.normal_vect: Tuple[float, float, float] = (0, 0, 1)
        self.o = make_vertex(center)
        self.face: Any | None = None
        self.name: str = "Surface"
        self.surf: Surface | None = None
        self.rotation: float = 0.0
        self.rotation_angle: float = 90.0
        self.z_axis: Any = make_vector_from_points(
            self.o, make_vertex((center[0], center[1], 1)))
        self.center_after_transl: Tuple[float, float, float] = (1, 1, 0)

    def test_rotate(self) -> None:
        """
        Method that tests the rotation functionality of the `Surface` subclass
        object.

        This test verifies that the `rotate` method of the `Surface` instance
        is applied correctly by comparing the rotation angle with the one
        between a reference vector before and after the rotation.

        Notes
        -----
        This test is skipped if run by the `TestSurface` class as it requires
        the proper instantiation of any of the `Surface` subclass the `surf`
        attribute refers to.
        """
        # Check the surface istantiation
        self.__check_surface_setup()
        # Build reference vectors before the rotation
        face_ref_vect, out_circle_ref_vect = make_ref_vectors(self.surf)
        # Rotate the surface
        self.surf.rotate(self.rotation_angle)
        # Build reference vectors after the rotation
        face_ref_vect2, out_circle_ref_vect2 = make_ref_vectors(self.surf)
        # Check the correct rotation happened
        self.__assert_rotation(
            face_ref_vect,
            out_circle_ref_vect,
            face_ref_vect2,
            out_circle_ref_vect2)

    def test_rotate_from_axis(self) -> None:
        """
        Method that tests the rotation functionality around a specified axis
        of the `Surface` subclass object.

        This test verifies that the `rotate_from_axis` method of the `Surface`
        instance is applied correctly by comparing the rotation angle with the
        one between a reference vector before and after the rotation.

        Notes
        -----
        This test is skipped if run by the `TestSurface` class as it requires
        the proper instantiation of any of the `Surface` subclass the `surf`
        attribute refers to.
        """
        # Check the surface istantiation
        self.__check_surface_setup()
        # Build reference vectors before the rotation
        face_ref_vect, out_circle_ref_vect = make_ref_vectors(self.surf)
        # Rotate the surface, using the declared axis
        self.surf.rotate_from_axis(self.rotation_angle, self.z_axis)
        # Build reference vectors after the rotation
        face_ref_vect2, out_circle_ref_vect2 = make_ref_vectors(self.surf)
        # Check the correct rotation happened
        self.__assert_rotation(
            face_ref_vect,
            out_circle_ref_vect,
            face_ref_vect2,
            out_circle_ref_vect2)

    def test_translation(self) -> None:
        """
        Method that tests the translation of a `Surface` subclass object and
        verifies the correctness of the resulting position of its geometric
        elements.

        Notes
        -----
        This test is skipped if run by the `TestSurface` class as it requires
        the proper instantiation of any of the `Surface` subclass the `surf`
        attribute refers to.
        """
        # Check the surface istantiation
        self.__check_surface_setup()
        # Store the elements for assessing the correct translation
        new_pos_vrtx = make_vertex(self.center_after_transl)
        cdg_pre = make_cdg(self.surf.face)
        distance = get_min_distance(self.surf.o, new_pos_vrtx)
        surf_org = deepcopy(self.surf)

        # Translate the surface
        self.surf.translate(self.center_after_transl)
        # Check the correct translation happened
        self.assertTrue(
            are_same_shapes(self.surf.o, new_pos_vrtx, ShapeType.VERTEX))
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre, make_cdg(self.surf.face)),
                distance))
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre, make_cdg(self.surf.out_circle)),
                distance))
        for v1, v2 in zip(surf_org.vertices, self.surf.vertices):
            self.assertTrue(
                math.isclose(get_min_distance(v1, v2), distance))
        for b1, b2 in zip(surf_org.borders, self.surf.borders):
            self.assertTrue(
                math.isclose(
                    get_min_distance(
                        make_vertex_on_curve(b1, 0.0),
                        make_vertex_on_curve(b2, 0.0)),
                    distance))

    def test_show_face(self) -> None:
        """
        Method that tests the `show_face` method of the `Surface` subclass
        object.

        It ensures that the face is correctly displayed in SALOME by:
        - verifying that SALOME assigns an entry ID to the face object;
        - checking that the geometric object the entry ID corresponds to
          matches the expected shape.

        Notes
        -----
        This test is skipped if run by the `TestSurface` class as it requires
        the proper instantiation of any of the `Surface` subclass the `surf`
        attribute refers to.
        """
        # Check the surface istantiation
        self.__check_surface_setup()
        self.surf.show_face()
        self.assertTrue(self.surf.face_entry_id is not None)
        self.assertTrue(are_same_shapes(
            self.surf.face,
            get_object_from_id(self.surf.face_entry_id),
            ShapeType.FACE))

    def test_show_borders(self) -> None:
        """
        Method that tests the `show_borders` method of the `Surface` subclass
        object.

        It ensures that the surface's borders are correctly displayed in
        SALOME by:
        - verifying that trying to show the borders without displaying the
          surface's face first raises an exception;
        - verifying that SALOME assigns an entry ID to each edge object being
          a border of the surface.

        Notes
        -----
        This test is skipped if run by the `TestSurface` class as it requires
        the proper instantiation of any of the `Surface` subclass the `surf`
        attribute refers to.
        """
        # Check the surface istantiation
        self.__check_surface_setup()
        # Test the 'Surface' borders display
        self.__test_show_surf_elements(self.surf.show_borders)

    def test_show_edges_and_borders(self) -> None:
        """
        Method that tests the `show_edges_and_vertices` method of the
        `Surface` subclass object.

        It ensures that the surface's borders and vertices are correctly
        displayed in SALOME by:
        - verifying that trying to show the borders and vertices without
          displaying the surface's face first raises an exception;
        - verifying that SALOME assigns an entry ID to each edge and vertex
          object being shown.

        Notes
        -----
        This test is skipped if run by the `TestSurface` class as it requires
        the proper instantiation of any of the `Surface` subclass the `surf`
        attribute refers to.
        """
        # Check the surface istantiation
        self.__check_surface_setup()
        # Test the 'Surface' edges and vertices display
        self.__test_show_surf_elements(self.surf.show_edges_and_vertices)
        for v in self.surf.vertices:
            _ = get_id_from_object(v)

    def __assert_rotation(self,
                          face_ref_vect: Any,
                          out_circle_ref_vect: Any,
                          face_ref_vect2: Any,
                          out_circle_ref_vect2: Any) -> None:
        """
        Method that asserts that the surface rotation and the angle between
        the given reference vectors match the expected rotation angle.

        Parameters
        ----------
        face_ref_vect : Any
            The reference vector on the face before rotation.
        out_circle_ref_vect : Any
            The reference vector on the outer circle before rotation.
        face_ref_vect2 : Any
            The reference vector on the face after rotation.
        out_circle_ref_vect2 : Any
            The reference vector on the outer circle after rotation.
        """
        self.assertEqual(self.surf.rotation,
                         math.radians(self.rotation_angle))
        self.assertTrue(
            math.isclose(
                get_angle_between_shapes(face_ref_vect, face_ref_vect2),
                self.rotation_angle))
        self.assertTrue(
            math.isclose(
                get_angle_between_shapes(out_circle_ref_vect,
                                        out_circle_ref_vect2),
                self.rotation_angle))

    def __check_surface_setup(self) -> None:
        """
        Method that checks whether the test class has the `surf` attribute
        correctly initialized. If the 'surf' attribute is `None`, the test
        this method is called in is skipped by providing an explanatory
        message.

        This method is called within each method of the class `TestSurface`
        that operates on the subclass of the `Surface` class the attribute
        `surf` is set to.
        """
        if self.surf is None:
            self.skipTest("TestSurface: attribute 'surf' not initialized. "
                          "This test is only valid in subclasses.")

    def __test_show_surf_elements(self, show_func: Callable) -> None:
        """
        Method that tests the behaviour when displaying the geometric elements
        of a `Surface` instance by means of the given `Callable` object.

        This method verifies that an exception is raised when attempting to
        show the geometric elements of the surface without having shown its
        corresponding face first. This is performed by running the provided
        `Callable` object, i.e. the `Surface` method to test.
        Then, the correct procedure is applied and it is checked that an entry
        ID is correctly associated to each of the shown borders of the
        surface.

        Parameters
        ----------
        show_func : Callable
            A function to be called for displaying the surface elements.
        """
        # Verify the exception is raised when incorrectly showing the
        # geometric elements of the 'Surface' object
        with self.assertRaises(RuntimeError):
            if self.surf.face_entry_id:
                self.surf.face_entry_id = None
            show_func()
        # Correctly show the surface's elements indicated by the callable
        self.surf.show_face()
        show_func()
        # Check if the borders' IDs in the SALOME study can be retrieved
        for b in self.surf.borders:
            _ = get_id_from_object(b)


class TestCircle(TestSurface):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Circle` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Circle` class can be correctly instantiated and that it properly
    handles rotation, translation, and visualization of their geometric
    elements (i.e. faces, borders, and vertices).
    Tests dealing with the above-mentioned operations are declared in the
    `TestSurface` class this class inherits from. They are run here, as this
    class declares the `surf` attribute.

    In addition, there are tests to check whether the operations for updating
    the `Circle` object's face are correctly handled.

    Attributes
    ----------
    In addition to the attributes declared in the `TestSurface` superclass,
    there are the following ones:

    radius : float
          The radius of the circle.
    name : str
        The name of the circle's face when added in the SALOME study.
    surf : Circle
        The `Surface` subclass that describes the geometric characteristics
        of a circle.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for a `Circle` class.
        It initializes the common geometric characteristics for a SALOME
        surface and the ones specific for describing a circle.
        The `surf` attribute is assigned to an instance of the `Circle` class
        so that the test methods can be run and addressed to the correct
        geometric object.
        """
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.radius = 5.0
        self.name = "Circle"
        self.surf: Circle = Circle(
            center=get_point_coordinates(self.o),
            normal_vect=self.normal_vect,
            radius=self.radius,
            name=self.name
        )

    def test_circle_init(self) -> None:
        """
        Method that tests the initialization of the `Circle` object by
        verifying it is correctly instantiated with the provided center,
        normal vector, radius, and name.
        It also checks that the geometric properties and associated shapes
        (vertex, border, face) are correctly set and match the expected
        reference objects.
        """
        # Declare a circle object with the same geometric characteristics for
        # comparison purposes
        circle = make_circle(
            center=self.o, axis=None, radius=self.radius)
        # Check the correct instantiation
        self.assertEqual(self.surf.radius, self.radius)
        self.assertEqual(self.surf.name, self.name)
        self.assertTrue(
            math.isclose(get_min_distance(self.surf.o, self.o), 0.0))
        self.assertEqual(len(self.surf.vertices), 1)
        self.assertEqual(len(self.surf.borders), 1)
        self.assertTrue(
            math.isclose(get_min_distance(self.surf.vertices[0], self.o), 0.0)
        )
        self.assertTrue(
            are_same_shapes(self.surf.borders[0], circle, ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(self.surf.out_circle, circle, ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(
                self.surf.face, make_face([circle]), ShapeType.FACE)
        )
        self.assertEqual(self.surf.lx, self.radius)
        self.assertEqual(self.surf.ly, self.radius)
        self.assertEqual(self.surf.rotation, self.rotation)
        self.assertEqual(self.surf.face_entry_id, None)

    def test_circle_build_borders(self) -> None:
        """
        Method that tests the implementation of the `_build_borders` method
        for a `Circle` class.
        It is checked whether the returned object is:
        - a list of objects;
        - it contains only one element;
        - the object is an EDGE having the geometric charateristics used to
          initialize the corresponding `Circle` instance.
        """
        # Build the list of borders of the surface
        borders = self.surf._build_borders()
        # Assess the correct creation of the borders
        self.assertTrue(isinstance(borders, List))
        self.assertEqual(len(borders), 1)
        border = borders[0]
        self.assertTrue(get_shape_type(border) == ShapeType.EDGE)
        self.assertTrue(
            str(get_kind_of_shape(border)[0]) == 'CIRCLE')
        self.assertTrue(
            are_same_shapes(make_cdg(border), self.o, ShapeType.VERTEX))
        self.assertTrue(
            math.isclose(
                get_min_distance(self.o, border),
                self.radius)
        )

    def test_circle_update_from_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `Circle` class.
        A face object representing a circle is built and used to update the
        one of a `Circle` instance.
        Afterwards, the geometric characteristics of the `Circle` object are
        checked to assess they matches the ones of the new face.
        """
        # Instantiate the 'Circle' object without initializing its attributes
        c = Circle.__new__(Circle)
        # Setup the face to update the 'Circle' object with
        face = make_face([make_circle(self.o, None, self.radius)])
        # Update the face and the attributes of the 'Circle' object
        c.update_from_face(face)

        # Check the instance has been updated correctly
        self.assertTrue(math.isclose(c.radius, self.radius))
        self.assertTrue(
            math.isclose(get_min_distance(c.o, make_cdg(face)), 0.0)
        )
        self.assertEqual(len(c.vertices), 1)
        self.assertEqual(len(c.borders), 1)
        self.assertTrue(
            are_same_shapes(c.borders[0],
                            extract_sub_shapes(face, ShapeType.EDGE)[0],
                            ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(c.face, face, ShapeType.FACE)
        )
        self.assertTrue(math.isclose(c.lx, self.radius))
        self.assertTrue(math.isclose(c.ly, self.radius))

    def test_circle_update_from_non_valid_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `Circle` class when providing an invalid face object.
        It verifies that an exception is correctly raised when updating the
        `Circle` instance with an object resulting from the partition of two
        shapes. This object has type `ShapeType.COMPOUND` which is not a valid
        option for the `update_from_face` method.
        """
        # Instantiate the 'Circle' object without initializing its attributes
        c = Circle.__new__(Circle)
        # Setup the non-valid object to update the 'Circle' object with
        face = make_partition(
            [make_face(make_circle(self.o, None, self.radius))],
            [make_face(make_circle(self.o, None, self.radius-1))],
            ShapeType.FACE)
        # Verify the exception is raised when updating the 'Circle' object
        with self.assertRaises(RuntimeError):
            c.update_from_face(face)


class TestRectangle(TestSurface):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Rectangle` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Rectangle` class can be correctly instantiated and that it properly
    handles rotation, translation, and visualization of their geometric
    elements (i.e. faces, borders, and vertices).
    Tests dealing with the above-mentioned operations are declared in the
    `TestSurface` class this class inherits from. They are run here, as this
    class declares the `surf` attribute.

    In addition, there are tests to check whether the operations for updating
    the `Rectangle` object's face are correctly handled.

    Attributes
    ----------
    In addition to the attributes declared in the `TestSurface` superclass,
    there are the following ones:

    height : float
        The height of the rectangle.
    width : float
        The width of the rectangle.
    rounded_corners : List[Tuple[int, float]]
        Indicating the corner index and the curvature radius
    name : str
        The name of the rectangle's face when added in the SALOME study.
    surf : Rectangle
        The `Surface` subclass that describes the geometric characteristics
        of a rectangle.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for a `Rectangle` class.
        It initializes the common geometric characteristics for a SALOME
        surface and the ones specific for describing a rectangle.
        The `surf` attribute is assigned to an instance of the `Rectangle`
        class so that the test methods can be run and addressed to the
        correct geometric object.
        """
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.height = 1.0
        self.width = 1.0
        self.curv_radius = 0.2
        self.rounded_corners = [(i, self.curv_radius) for i in range(4)]
        self.name = "Rectangle"
        self.surf: Rectangle = Rectangle(
            center=get_point_coordinates(self.o),
            height=self.height,
            width=self.width,
            name=self.name
        )

    def test_rect_init(self) -> None:
        """
        Method that tests the initialization of the `Rectangle` object by
        verifying it is correctly instantiated with the provided center,
        height, width, and name.
        It also checks that the geometric properties and associated shapes
        (vertex, border, face) are correctly set and match the expected
        reference objects.
        """
        # Build a rectangular shape with the same geometric characteristics
        # for comparison purposes
        vertices, rectangle_edges = self.__build_rect_geom_elements(
            self.width, self.height
        )
        ref_circle = make_circle(
            self.o,
            None,
            math.sqrt(self.width*self.width/4 + self.height*self.height/4))
        # Check the correct instantiation
        self.__assess_instantiation(vertices, rectangle_edges, self.o)
        self.assertTrue(
            are_same_shapes(self.surf.out_circle, ref_circle, ShapeType.EDGE)
        )
        self.assertEqual(self.surf.rotation, self.rotation)
        self.assertEqual(self.surf.face_entry_id, None)
        self.assertEqual(self.surf.name, self.name)

    def test_build_borders(self) -> None:
        """
        Method that tests the implementation of the `_build_borders` method
        for a `Rectangle` class.
        It is checked whether the returned object is:
        - a list of objects;
        - it contains four elements;
        - each element is an EDGE of type `SEGMENT`;
        - each element has the geometric charateristics used to initialize
          the corresponding `Rectangle` instance.
        """
        # Build the list of borders of the surface
        borders = self.surf._build_borders()
        # Assess the correct creation of the borders
        self.assertTrue(isinstance(borders, List))
        self.assertEqual(len(borders), 4)
        for border in borders:
            self.assertTrue(get_shape_type(border) == ShapeType.EDGE)
            self.assertTrue(
                str(get_kind_of_shape(border)[0]) == 'SEGMENT')
            self.assertTrue(
                math.isclose(
                    get_min_distance(make_cdg(border), self.o),
                    self.surf.lx/2) or
                math.isclose(
                    get_min_distance(make_cdg(border), self.o),
                    self.surf.ly/2)
            )
            border_len = get_basic_properties(border)[0]
            self.assertTrue(
                math.isclose(border_len, self.surf.lx) or
                math.isclose(border_len, self.surf.ly)
            )

    def test_build_borders_with_rounded_corners(self) -> None:
        """
        Method that tests the implementation of the private method
        `__build_borders_with_rounded_corners` for a `Rectangle` class.
        It is checked whether the returned object is:
        - a list of objects;
        - it contains four elements;
        - each element is an EDGE of type `ARC_CIRCLE`;
        - each element has its vertices and radius set accordingly.
        """
        # Build the arcs representing the rounded corners
        arcs = self.surf._Rectangle__build_borders_with_rounded_corners(
            self.rounded_corners,
            get_point_coordinates(self.o),
            self.surf.ly,
            self.surf.lx
        )
        # Assess the correct creation of the arc borders
        self.assertTrue(isinstance(arcs, List))
        self.assertEqual(len(arcs), 4)
        for arc in arcs:
            # Get the geometric information about each arc
            data = get_kind_of_shape(arc)
            self.assertTrue(get_shape_type(arc) == ShapeType.EDGE)
            self.assertTrue(str(data[0]) == 'ARC_CIRCLE')
            for v in extract_sub_shapes(arc, ShapeType.VERTEX):
                self.assertTrue(
                    math.isclose(
                        get_min_distance(v, self.o),
                        math.sqrt(
                            self.surf.lx/2*self.surf.lx/2 + (
                                self.surf.ly/2 - self.curv_radius)*(
                                    self.surf.ly/2 - self.curv_radius))) or
                    math.isclose(
                        get_min_distance(v, self.o),
                        math.sqrt(
                            (self.surf.lx/2 - self.curv_radius)*(
                                self.surf.lx/2 - self.curv_radius) +
                            self.surf.ly/2*self.surf.ly/2))
                )
            self.assertTrue(
                math.isclose(data[7], self.curv_radius)
            )

    def test_update_from_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `Rectangle` class.
        A face object representing a rectangle is built and used to update the
        one of a `Rectangle` instance.
        Afterwards, the geometric characteristics of the `Rectangle` object
        are checked to assess they matches the ones of the new face.
        """
        # Dimensions of the new rectangle
        self.width = 2.0
        self.height = 1.0
        # Instantiate the 'Rectangle' object without initializing its
        # attributes
        self.surf = Rectangle.__new__(Rectangle)
        # Setup the rectangular shape to update the 'Rectangle' object with
        vertices, edges = self.__build_rect_geom_elements(self.width,
                                                          self.height)
        face = make_face(edges)
        # Update the face and the attributes of the 'Rectangle' object
        self.surf.update_from_face(face)

        # Check the instance has been updated correctly
        self.__assess_instantiation(vertices, edges, make_cdg(face))

    def test_update_from_non_valid_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `Rectangle` class when providing an invalid face object.
        It verifies that an exception is correctly raised when updating the
        `Rectangle` instance with an object resulting from the partition of
        two shapes (having an invalid type `ShapeType.COMPOUND`) or a
        non-rectangular shape.
        """
        # Setup the non-valid object to update the 'Rectangle' object with
        face = make_partition(
            [make_face(make_circle(self.o, None, self.width))],
            [make_face(make_circle(self.o, None, self.width-0.1))],
            ShapeType.FACE)
        # Verify the exception is raised when updating the 'Rectangle' object
        with self.assertRaises(RuntimeError):
            self.surf.update_from_face(face)
        with self.assertRaises(RuntimeError):
            self.surf.update_from_face(
                make_face(make_circle(self.o, None, self.width)))

    def __assess_instantiation(
            self, vertices: List[Any], edges: List[Any], center: Any) -> None:
        """
        Method that verifies whether the `Rectangle` class has been correctly
        instantiated by checking that its geometric characteristics match the
        ones used to initialize the instance.

        Parameters
        ----------
        vertices : List[Any]
            The list of vertex objects of the reference rectangle.
        edges : List[Any]
            The list of edge objects of the reference rectangle.
        center : Any
            The vertex object being the center of the reference rectangle.
        """
        self.assertTrue(
            math.isclose(self.surf.lx, self.width, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(self.surf.ly, self.height, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(get_min_distance(self.surf.o, center), 0.0))
        self.assertEqual(len(self.surf.vertices), 4)
        self.assertEqual(len(self.surf.borders), 4)
        for v_rect, v_ref in zip(self.surf.vertices, vertices):
            self.assertTrue(
                math.isclose(get_min_distance(v_rect, v_ref), 0.0)
        )
        for b_rect, b_ref in zip(self.surf.borders, edges):
            self.assertTrue(
                are_same_shapes(b_rect, b_ref, ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(
                self.surf.face, make_face(edges), ShapeType.FACE)
        )

    def __build_rect_geom_elements(
            self, width: float, height: float) -> Tuple[List[Any], List[Any]]:
        """
        Method that, given the width and height, builds the vertex and edge
        objects that represent a rectangle.

        Parameters
        ----------
        width : float
            The width of the rectangle.
        height : float
            The height of the rectangle.

        Returns
        -------
        Tuple[List[Any], List[Any]]
            A tuple with the list of vertices and edges of the rectangle.
        """
        vertices = [
            make_vertex((-width/2, -height/2, 0)),
            make_vertex((width/2, -height/2, 0)),
            make_vertex((width/2, height/2, 0)),
            make_vertex((-width/2, height/2, 0)),
        ]
        edges = [
            make_edge(vertices[i], vertices[(i+1) % 4]) for i in range(4)
        ]
        return vertices, edges


class TestHexagon(TestSurface):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Hexagon` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Hexagon` class can be correctly instantiated and that it properly
    handles rotation, translation, and visualization of their geometric
    elements (i.e. faces, borders, and vertices).
    Tests dealing with the above-mentioned operations are declared in the
    `TestSurface` class this class inherits from. They are run here, as this
    class declares the `surf` attribute.

    In addition, there are tests to check whether the operations for updating
    the `Hexagon` object's face are correctly handled.

    Attributes
    ----------
    In addition to the attributes declared in the `TestSurface` superclass,
    there are the following ones:

    edge_length : float
        The length of the hexagon's edge.
    apothem : float
        The length of the hexagon's apothem.
    name : str
        The name of the hexagon's face when added in the SALOME study.
    surf : Hexagon
        The `Surface` subclass that describes the geometric characteristics
        of a hexagon.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for a `Hexagon` class.
        It initializes the common geometric characteristics for a SALOME
        surface and the ones specific for describing a hexagon.
        The `surf` attribute is assigned to an instance of the `Hexagon`
        class so that the test methods can be run and addressed to the
        correct geometric object.
        """
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.edge_length = 1.0
        self.apothem = self.edge_length * math.sin(math.pi/3)
        self.name = "Hexagon"
        self.surf: Hexagon = Hexagon(
            center=get_point_coordinates(self.o),
            edge_length=self.edge_length,
            name=self.name
        )

    def test_hex_init(self) -> None:
        """
        Method that tests the initialization of the `Hexagon` object by
        verifying it is correctly instantiated with the provided center,
        edge length, and name.
        It also checks that the geometric properties and associated shapes
        (vertex, border, face) are correctly set and match the expected
        reference objects.
        """
        # Build a hexagonal shape with the same geometric characteristics
        # for comparison purposes
        vertices, edges = build_hex_geom_elements(self.o, self.edge_length)
        ref_circle = make_circle(self.o, None, self.edge_length)
        # Check the correct instantiation
        self.__assess_instantiation(vertices, edges, self.o)
        self.assertTrue(
            are_same_shapes(self.surf.out_circle, ref_circle, ShapeType.EDGE)
        )
        self.assertEqual(self.surf.rotation, self.rotation)
        self.assertEqual(self.surf.face_entry_id, None)
        self.assertEqual(self.surf.name, self.name)

    def test_build_borders(self) -> None:
        """
        Method that tests the implementation of the `_build_borders` method
        for a `Hexagon` class.
        It is checked whether the returned object is:
        - a list of objects;
        - it contains six elements;
        - each element is an EDGE of type `SEGMENT`;
        - each element has the geometric charateristics used to initialize
          the corresponding `Hexagon` instance.
        """
        # Build the list of borders of the surface
        borders = self.surf._build_borders()
        # Assess the correct creation of the borders
        self.assertTrue(isinstance(borders, List))
        self.assertEqual(len(borders), 6)
        for border in borders:
            self.assertTrue(get_shape_type(border) == ShapeType.EDGE)
            self.assertTrue(
                str(get_kind_of_shape(border)[0]) == 'SEGMENT')
            self.assertTrue(
                math.isclose(
                    get_min_distance(make_cdg(border), self.o),
                    self.surf.ly)
            )
            border_len = get_basic_properties(border)[0]
            self.assertTrue(
                math.isclose(border_len, self.surf.lx)
            )

    def test_update_from_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `Hexagon` class.
        A face object representing a hexagon is built and used to update the
        one of a `Hexagon` instance.
        Afterwards, the geometric characteristics of the `Hexagon` object
        are checked to assess they matches the ones of the new face.
        """
        # Dimensions of the new hexagon
        self.edge_length = 2.0
        self.apothem = self.edge_length * math.sin(math.pi/3)
        # Instantiate the 'Hexagon' object without initializing its
        # attributes
        self.surf = Hexagon.__new__(Hexagon)
        # Setup the hexagonal shape to update the 'Hexagon' object with
        vertices, edges = build_hex_geom_elements(self.o, self.edge_length)
        face = make_face(edges)

        # Update the face and the attributes of the 'Hexagon' object
        self.surf.update_from_face(face)

        # Check the instance has been updated correctly
        self.__assess_instantiation(vertices, edges, make_cdg(face))

    def test_update_from_non_valid_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `Hexagon` class when providing an invalid face object.
        It verifies that an exception is correctly raised when updating the
        `Hexagon` instance with an object resulting from the partition of
        two shapes (having an invalid type `ShapeType.COMPOUND`) or a
        non-hexagonal shape.
        """
        # Setup the non-valid object to update the 'Hexagon' object with
        face = make_partition(
            [make_face(make_circle(self.o, None, self.edge_length))],
            [make_face(make_circle(self.o, None, self.edge_length-0.1))],
            ShapeType.FACE)
        # Verify the exception is raised when updating the 'Hexagon' object
        with self.assertRaises(RuntimeError):
            self.surf.update_from_face(face)
        with self.assertRaises(RuntimeError):
            self.surf.update_from_face(
                make_face(make_circle(self.o, None, self.edge_length)))

    def __assess_instantiation(
            self, vertices: List[Any], edges: List[Any], center: Any) -> None:
        """
        Method that verifies whether the `Hexagon` class has been correctly
        instantiated by checking that its geometric characteristics match the
        ones used to initialize the instance.

        Parameters
        ----------
        vertices : List[Any]
            The list of vertex objects of the reference hexagon.
        edges : List[Any]
            The list of edge objects of the reference hexagon.
        center : Any
            The vertex object being the center of the reference hexagon.
        """
        self.assertTrue(
            math.isclose(self.surf.lx, self.edge_length, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(self.surf.ly, self.apothem, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(get_min_distance(self.surf.o, center), 0.0))
        self.assertEqual(len(self.surf.vertices), 6)
        self.assertEqual(len(self.surf.borders), 6)
        for v_hex, v_ref in zip(self.surf.vertices, vertices):
            self.assertTrue(
                math.isclose(get_min_distance(v_hex, v_ref), 0.0)
        )
        for b_hex, b_ref in zip(self.surf.borders, edges):
            self.assertTrue(
                are_same_shapes(b_hex, b_ref, ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(
                self.surf.face, make_face(edges), ShapeType.FACE)
        )


class TestGenericSurface(TestSurface):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `GenericSurface` class.

    This test suite provides common setup and a set of tests to ensure that
    the `GenericSurface` class can be correctly instantiated and that it
    properly handles rotation, translation, and visualization of their
    geometric elements (i.e. faces, borders, and vertices).
    Tests dealing with the above-mentioned operations are declared in the
    `TestSurface` class this class inherits from. They are run here, as this
    class declares the `surf` attribute.

    In addition, there are tests to check whether the operations for updating
    the `GenericSurface` object's face are correctly handled.

    Attributes
    ----------
    In addition to the attributes declared in the `TestSurface` superclass,
    there are the following ones:

    face : Any
        The face object to build the `GenericSurface` instance from.
    name : str
        The name of the generic surface when added in the SALOME study.
    surf : GenericSurface
        The `Surface` subclass that describes the geometric characteristics
        of a generic surface.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for a `GenericSurface` class.
        It initializes the common geometric characteristics for a SALOME
        surface and the ones specific for a generic surface.
        The `surf` attribute is assigned to an instance of the
        `GenericSurface` class so that the test methods can be run and
        addressed to the correct geometric object.
        """
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.face = self.__build_generic_face()
        self.name = "GenericSurface"
        self.surf: GenericSurface = GenericSurface(
            face=self.face,
            name=self.name
        )

    def test_init(self) -> None:
        """
        Method that tests the initialization of the `GenericSurface` object
        by verifying it is correctly instantiated with the provided face and
        name.
        It also checks that the geometric properties and associated shapes
        (vertex, border, face) are correctly set and match the expected
        reference objects.
        """
        # Extract the vertices and edges from the generic surface
        vertices = extract_sub_shapes(self.face, ShapeType.VERTEX)
        edges = extract_sub_shapes(self.face, ShapeType.EDGE)
        # Check the correct instantiation
        self.__assess_instantiation(vertices, edges, self.o)
        self.assertEqual(self.surf.rotation, self.rotation)
        self.assertEqual(self.surf.face_entry_id, None)
        self.assertEqual(self.surf.name, self.name)

    def test_build_borders(self) -> None:
        """
        Method that tests the implementation of the `_build_borders` method
        for a `GenericSurface` class.
        It is checked whether the returned object is:
        - a list of objects;
        - it contains six elements;
        - each element is an EDGE of type `SEGMENT`, `ARC_CIRCLE` or `CIRCLE`;
        - each element has the geometric charateristics used to initialize
          the corresponding `GenericSurface` instance.
        """
        # Build the list of borders of the surface
        borders = self.surf._build_borders()
        # Assess the correct creation of the borders
        self.assertTrue(isinstance(borders, List))
        for border in borders:
            self.assertTrue(get_shape_type(border) == ShapeType.EDGE)
            self.assertTrue(
                str(get_kind_of_shape(border)[0]) in [
                    'SEGMENT', 'ARC_CIRCLE', 'CIRCLE'])

    def test_update_from_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `GenericSurface` class.
        A generic face object is built and used to update the one of a
        `GenericSurface` instance.
        Afterwards, the geometric characteristics of the `GenericSurface`
        object are checked to assess they matches the ones of the new face.
        """
        # Instantiate the 'GenericSurface' object without initializing its
        # attributes
        self.surf = GenericSurface.__new__(GenericSurface)
        # Setup the shape to update the 'GenericSurface' object with
        self.face = make_face(make_face([make_circle(self.o, None, 5.0)]))

        # Update the face and the attributes of the 'GenericSurface' object
        self.surf.update_from_face(self.face)

        # Check the instance has been updated correctly
        self.__assess_instantiation(
            extract_sub_shapes(self.face, ShapeType.VERTEX),
            extract_sub_shapes(self.face, ShapeType.EDGE),
            make_cdg(self.face))

    def test_update_from_non_valid_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `GenericSurface` class when providing an invalid face object.
        It verifies that an exception is correctly raised when updating the
        `GenericSurface` instance with an object that has an invalid type,
        i.e. neither `ShapeType.COMPOUND`, nor `ShapeType.FACE`.
        """
        # Verify the exception is raised when updating the 'GenericSurface'
        # object
        with self.assertRaises(RuntimeError):
            self.surf.update_from_face(make_circle(self.o, None, 1.0))

    def __assess_instantiation(
            self, vertices: List[Any], edges: List[Any], center: Any) -> None:
        """
        Method that verifies whether the `GenericSurface` class has been
        correctly instantiated by checking that its geometric characteristics
        match the ones used to initialize the instance.

        Parameters
        ----------
        vertices : List[Any]
            The list of vertex objects of the reference surface.
        edges : List[Any]
            The list of edge objects of the reference surface.
        center : Any
            The vertex object being the center of the reference surface.
        """
        # Extract the extension of the surface's bounding box
        b_box = get_bounding_box(self.face)
        lx = (b_box[1] - b_box[0]) / 2
        ly = (b_box[3] - b_box[2]) / 2
        # Check the attributes of 'GenericSurface'
        self.assertTrue(
            math.isclose(self.surf.lx, lx, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(self.surf.ly, ly, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(
                get_min_distance(self.surf.o, center), 0.0, abs_tol=1e-6)
        )
        self.assertEqual(len(self.surf.vertices), len(vertices))
        self.assertEqual(len(self.surf.borders), len(edges))
        self.assertTrue(
            are_same_shapes(self.surf.out_circle,
                            make_circle(self.o, None, max(lx, ly)),
                            ShapeType.EDGE)
        )
        for v_surf, v_ref in zip(self.surf.vertices, vertices):
            self.assertTrue(
                math.isclose(get_min_distance(v_surf, v_ref), 0.0)
        )
        for b_surf, b_ref in zip(self.surf.borders, edges):
            self.assertTrue(
                are_same_shapes(b_surf, b_ref, ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(
                self.surf.face, self.face, ShapeType.FACE)
        )

    def __build_generic_face(self) -> Any:
        """
        Method that builds a generic shape made by two hexagons one into the
        other.

        Return
        ------
        Any
            A GEOM object made by partitioning a hexagon face with another
            hexagon with smaller dimensions.
        """
        _, edges1 = build_hex_geom_elements(self.o, 1.0)
        _, edges2 = build_hex_geom_elements(self.o, 2.0)
        return make_partition(
            [make_face(edges2)], [make_face(edges1)], ShapeType.FACE)


class TestHexagonBuilder(unittest.TestCase):
    """
    Test case for verifying the correct functioning of the `build_hexagon`
    function that allows to build a `Hexagon` instance from the apothem of
    the hexagonal shape.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the `build_hexagon`
        function.
        """
        self.apothem = 1.0
        self.center = (1.0, 1.0, 0.0)
        self.edge_length = self.apothem / math.sin(math.pi/3)
        self.o = make_vertex(self.center)

    def test_build_hexagon(self) -> None:
        """
        Method that verifies the `build_hexagon` function by checking if the
        attributes of the built `Hexagon` instance are correctly set.
        """
        # Build the 'Hexagon' instance
        hex = build_hexagon(self.apothem, self.center)
        # Check the correct instantiation
        self.assertTrue(
            math.isclose(hex.lx, self.edge_length, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(hex.ly, self.apothem, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(get_min_distance(hex.o, self.o), 0.0, abs_tol=1e-6)
        )
        self.assertEqual(len(hex.vertices), 6)
        self.assertEqual(len(hex.borders), 6)
        self.assertTrue(
            all(
                math.isclose(
                    get_min_distance(v, self.o),
                    self.edge_length,
                    abs_tol=1e-6) for v in hex.vertices
            )
        )
        self.assertTrue(
            all(
                math.isclose(
                    round(get_basic_properties(b)[0], 6),
                    self.edge_length,
                    abs_tol=1e-6) for b in hex.borders
            )
        )


if __name__ == "__main__":
    unittest.main()
