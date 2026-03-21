"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.geometries` module have a valid implementation.
"""
import math
import unittest

from copy import deepcopy

from glow.interface.geom_entities import Edge, Face, wrap_shape
from glow.support.utility import are_same_shapes
from glow.interface.geom_interface import *
from glow.geometry_layouts.geometries import Surface, Hexagon, Rectangle, \
    Surface, Circle, build_hexagon_from_apothem, build_parallelogram, \
    build_regular_triangle, build_right_triangle, \
    build_right_triangle_from_catheti
from tests.unittest.support_funcs import build_hex_geom_elements, \
    make_ref_vectors


class TestSurface(unittest.TestCase):
    """
    Base test case for verifying the geometric operations and visualization
    capabilities of the `Surface` class and of its subclasses.

    This test suite provides common setup and a set of tests to ensure that
    implementations of `Surface` correctly handle rotation, scaling,
    translation, and visualization of their geometric elements.

    Attributes
    ----------
    o : Any
        The origin vertex of the surface.
    name : str
        The name of the surface.
    surface : Surface
        The surface instance under test. It is build by cutting two hexagons
        with different sizes, placed in the same position.
    rotation : float
        The initial rotation angle of the surface, in degrees.
    rotation_angle : float
        The angle by which to rotate the surface during tests.
    z_axis : Any
        The axis vector used for rotation tests.
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
        # self.normal_vect: Tuple[float, float, float] = (0, 0, 1)
        self.o = make_vertex(center)
        # self.face: Any | None = None
        self.name: str = "Surface"
        self.surface: Surface = Surface(
            Hexagon(edge_length=2) - Hexagon()
        )
        self.rotation: float = 0.0
        self.rotation_angle: float = 90.0
        self.z_axis: Any = make_vector_from_points(
            self.o, make_vertex((center[0], center[1], 1)))

    def test_init_no_face(self) -> None:
        """
        Method that tests the initialisation of the `Surface` class in case
        no face is provided.
        """
        # Skip the test if run from a subclass of 'Surface'
        self.__skip_if_subclass()
        # Test the initialisation of the 'Surface' class without any face and
        # centre
        surface = Surface(None)
        self.assertTrue(
            are_same_shapes(surface.o, self.o, ShapeType.VERTEX)
        )
        self.assertEqual(surface.borders, [])
        self.assertEqual(surface.name, "Surface")
        self.assertIsNone(surface.geom_obj)
        self.assertIsNone(surface.entry_id)
        self.assertEqual(surface.dimensions, (0.0, 0.0))
        self.assertEqual(surface.rot_angle, 0.0)

    def test_init_wrong_face(self) -> None:
        """
        Method that tests the initialisation of the `Surface` class in case
        a GEOM object other than a face is provided.
        """
        # Skip the test if run from a subclass of 'Surface'
        self.__skip_if_subclass()
        # Test the exception is raised when initialising the 'Surface' class
        # with a compound
        with self.assertRaises(RuntimeError):
            Surface(Hexagon(edge_length=2) // Hexagon())

    def test_init_face(self) -> None:
        """
        Method that tests the initialisation of the `Surface` class in case
        a valid face is provided.
        """
        # Skip the test if run from a subclass of 'Surface'
        self.__skip_if_subclass()
        # Test the initialisation of the 'Surface' class with a valid face and
        # no centre
        face_side = 2
        face = Hexagon(edge_length=face_side) - Hexagon()
        surface = Surface(face)
        self.assertTrue(
            are_same_shapes(surface.o, self.o, ShapeType.VERTEX)
        )
        for b, e in zip(
            surface.borders,
            [wrap_shape(e) for e in extract_sub_shapes(face, ShapeType.EDGE)]
        ):
            self.assertTrue(are_same_shapes(b, e, ShapeType.EDGE))

        self.assertEqual(surface.name, "Surface")
        self.assertIsNotNone(surface.geom_obj)
        self.assertTrue(
            are_same_shapes(surface.geom_obj, face.geom_obj, ShapeType.FACE)
        )
        self.assertIsNone(surface.entry_id)
        self.assertAlmostEqual(surface.dimensions[0], 2*face_side)
        self.assertTrue(
            math.isclose(
                surface.dimensions[1],
                2 * face_side * math.sin(math.pi/3),
                abs_tol=1e-6
            )
        )
        self.assertEqual(surface.rot_angle, 0.0)
        # Test the initialisation of the 'Surface' class with a valid face and
        # centre
        centre = (1.0, 1.0, 0.0)
        trasl_vec = make_vector(centre)
        trasl_face = make_translation(face, trasl_vec)
        surface = Surface(wrap_shape(trasl_face), centre)
        self.assertTrue(
            are_same_shapes(
                surface.geom_obj, trasl_face, ShapeType.FACE
            )
        )
        self.assertTrue(
            are_same_shapes(surface.o, make_vertex(centre), ShapeType.VERTEX)
        )

    def test_rotate(self) -> None:
        """
        Method that tests the rotation functionality of the `Surface` class.

        This test verifies that the `rotate` method of the `Surface` instance
        is applied correctly by comparing the rotation angle with the one
        between a reference vector before and after the rotation.
        """
        # Build reference vector before the rotation
        face_ref_vect = make_ref_vectors(self.surface)
        # Rotate the surface
        self.surface.rotate(self.rotation_angle)
        # Build reference vector after the rotation
        face_ref_vect2 = make_ref_vectors(self.surface)
        # Check the correct rotation happened
        self.__assert_rotation(
            face_ref_vect,
            face_ref_vect2,
        )

        # Check no rotation happened if the rotation angle is zero
        self.surface.rotate(0.0)
        self.assertAlmostEqual(self.surface.rot_angle, self.rotation_angle)

    def test_rotate_from_axis(self) -> None:
        """
        Method that tests the rotation functionality around a specified axis
        of the `Surface` object.

        This test verifies that the `rotate` method of the `Surface` instance
        is applied correctly by comparing the rotation angle with the
        one between a reference vector before and after the rotation.
        """
        # Build reference vector before the rotation
        face_ref_vect = make_ref_vectors(self.surface)
        # Rotate the surface, using the declared axis
        self.surface.rotate(self.rotation_angle, self.z_axis)
        # Build reference vector after the rotation
        face_ref_vect2 = make_ref_vectors(self.surface)
        # Check the correct rotation happened
        self.__assert_rotation(
            face_ref_vect,
            face_ref_vect2,
        )

    def test_scale(self) -> None:
        """
        Method that tests the scaling of a `Surface` object by checking the
        area of the surface before and after the operation.
        """
        scale_factor = 2
        area_before = get_basic_properties(self.surface)[1]
        # Scale the surface wrt to its centre
        self.surface.scale(scale_factor)
        # Check the area has increased by the scaling factor^2
        self.assertTrue(
            math.isclose(
                area_before*scale_factor*scale_factor,
                get_basic_properties(self.surface)[1],
                abs_tol=1e-6
            )
        )

        # Verify the exception is raised if an invalid scaling factor is
        # provided
        with self.assertRaises(ValueError):
            self.surface.scale(0.0)
        with self.assertRaises(ValueError):
            self.surface.scale(-2.0)

    def test_show(self) -> None:
        """
        Method that tests the `show` method of the `Surface` class.

        It ensures that the face is correctly displayed in SALOME by:

        - verifying that SALOME assigns an entry ID to the face object;
        - checking that the geometric object the entry ID corresponds to
          matches the expected shape.
        """
        self.surface.show()
        self.assertTrue(self.surface.entry_id is not None)
        self.assertTrue(
            are_same_shapes(
                self.surface,
                get_object_from_id(self.surface.entry_id),
                ShapeType.FACE
            )
        )

    def test_translation(self) -> None:
        """
        Method that tests the translation of a `Surface` object and verifies
        the correctness of the resulting position of its geometric elements.
        """
        # Store the elements for assessing the correct translation
        center_after_transl = (1, 1, 0)
        new_pos_vrtx = make_vertex(center_after_transl)
        cdg_pre = make_cdg(self.surface)
        distance = get_min_distance(self.surface.o, new_pos_vrtx)
        surf_org = deepcopy(self.surface)

        # Translate the surface
        self.surface.translate(center_after_transl)
        # Check the correct translation happened
        self.assertTrue(
            are_same_shapes(self.surface.o, new_pos_vrtx, ShapeType.VERTEX)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre, make_cdg(self.surface)),
                distance,
                abs_tol=1e-6
            )
        )
        for b1, b2 in zip(surf_org.borders, self.surface.borders):
            self.assertTrue(
                math.isclose(
                    get_min_distance(
                        make_vertex_on_curve(b1, 0.0),
                        make_vertex_on_curve(b2, 0.0)
                    ),
                    distance,
                    abs_tol=1e-6
                )
            )

        # Check no translation happened if the new centre coincides with the
        # current one
        self.surface.translate(center_after_transl)
        self.assertTrue(
            are_same_shapes(self.surface.o, new_pos_vrtx, ShapeType.VERTEX)
        )

    def test_update(self) -> None:
        """
        Method that tests the update of a `Surface` object with another face.
        """
        # Skip the test if run from a subclass of 'Surface'
        self.__skip_if_subclass()
        # Verify the exception is raised when providing an invalid shape
        with self.assertRaises(RuntimeError):
            self.surface.update(
                Edge(make_edge(self.o, make_vertex((1.0, 0.0, 0.0))))
            )
        # Update the surface with another one
        new_layout = Circle()
        self.surface.update(new_layout)
        # Verify that the geometric elements has changed correctly
        self.assertTrue(
            are_same_shapes(self.surface, new_layout, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(self.surface.o, new_layout.o, ShapeType.VERTEX)
        )
        for b1, b2 in zip(self.surface.borders, new_layout.borders):
            self.assertTrue(
                are_same_shapes(b1, b2, ShapeType.EDGE)
            )
        for d1, d2 in zip(self.surface.dimensions, new_layout.dimensions):
            self.assertAlmostEqual(d1/2, d2)

    def __assert_rotation(
            self,
            face_ref_vect: Any,
            face_ref_vect2: Any,
        ) -> None:
        """
        Method that asserts that the surface rotation and the angle between
        the given reference vectors match the expected rotation angle.

        Parameters
        ----------
        face_ref_vect : Any
            The reference vector on the face before rotation.
        face_ref_vect2 : Any
            The reference vector on the face after rotation.
        """
        self.assertEqual(self.surface.rot_angle, self.rotation_angle)
        self.assertTrue(
            math.isclose(
                get_angle_between_shapes(face_ref_vect, face_ref_vect2),
                self.rotation_angle
            )
        )

    def __skip_if_subclass(self) -> None:
        """
        Method that checks whether the current test class is `TestSurface`
        and skips the test that runs this method if the class is a subclass.
        """
        if self.__class__ is not TestSurface:
            self.skipTest(
                f"{self.__class__.__name__}. This test is only valid in "
                "superclass 'TestSurface'."
            )


class TestCircle(TestSurface):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Circle` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Circle` class can be correctly instantiated and that it properly
    handles rotation, scaling, translation, and visualization of their
    geometric elements, as well as their update.
    Tests dealing with the above-mentioned operations are declared in the
    `TestSurface` class this class inherits from.

    Attributes
    ----------
    In addition to the attributes declared in the `TestSurface` superclass,
    there are the following ones:

    radius : float
          The radius of the circle.
    name : str
        The name of the circle's face when added in the SALOME study.
    surface : Circle
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
        self.surface: Circle = Circle(
            center=get_point_coordinates(self.o),
            radius=self.radius,
            name=self.name
        )

    def test_circle_init(self) -> None:
        """
        Method that tests the initialization of the `Circle` object by
        verifying it is correctly instantiated with the provided center,
        normal vector, radius, and name.
        It also checks that the geometric properties and associated shapes
        are correctly set and match the expected reference objects.
        """
        # Declare a circle object with the same geometric characteristics for
        # comparison purposes
        circle = make_circle(center=self.o, axis=None, radius=self.radius)
        # Check the correct instantiation
        self.assertEqual(self.surface.radius, self.radius)
        self.assertEqual(self.surface.name, self.name)
        self.assertTrue(
            are_same_shapes(self.surface.o, self.o, ShapeType.VERTEX)
        )
        self.assertEqual(len(self.surface.borders), 1)
        self.assertTrue(
            are_same_shapes(self.surface.borders[0], circle, ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(self.surface, make_face([circle]), ShapeType.FACE)
        )
        self.assertEqual(self.surface.dimensions[0], self.radius)
        self.assertEqual(self.surface.dimensions[1], self.radius)
        self.assertEqual(self.surface.rot_angle, self.rotation)
        self.assertEqual(self.surface.entry_id, None)

    def test_update(self) -> None:
        """
        Method that tests the implementation of the `update` method for a
        `Circle` class.
        A face object representing a circle is built and used to update the
        one of a `Circle` instance.
        Afterwards, the geometric characteristics of the `Circle` object are
        checked to assess they matches the ones of the new face.
        """
        # Setup the face to update the 'Circle' object with
        face = make_face([make_circle(self.o, None, self.radius)])
        # Update the face and the attributes of the 'Circle' object
        self.surface.update(Face(face))

        # Check the instance has been updated correctly
        self.assertTrue(
            math.isclose(self.surface.radius, self.radius, abs_tol=1e-6)
        )
        self.assertTrue(
            are_same_shapes(self.surface.o, make_cdg(face), ShapeType.VERTEX)
        )
        self.assertEqual(len(self.surface.borders), 1)
        self.assertTrue(
            are_same_shapes(
                self.surface.borders[0],
                extract_sub_shapes(face, ShapeType.EDGE)[0],
                ShapeType.EDGE
            )
        )
        self.assertTrue(are_same_shapes(self.surface, face, ShapeType.FACE))
        self.assertTrue(math.isclose(self.surface.dimensions[0], self.radius))
        self.assertTrue(math.isclose(self.surface.dimensions[1], self.radius))

    def test_update_from_non_valid_face(self) -> None:
        """
        Method that tests the implementation of the `update` method for a
        `Circle` class when providing an invalid face object.
        It verifies that an exception is correctly raised when updating the
        `Circle` instance with:

        - an object resulting from the partition of two shapes. This object
          has type `ShapeType.COMPOUND` which is not a valid option for the
          `update` method.
        - a face object resulting from cutting a circle with another one.
          Despite the correct `ShapeType.FACE` type, the shape does not
          represent a circle, hence it is not admitted.
        """
        # Instantiate the 'Circle' object without initializing its attributes
        c = Circle.__new__(Circle)
        # Setup the non-valid object to update the 'Circle' object with
        face = make_partition(
            [make_face(make_circle(self.o, None, self.radius))],
            [make_face(make_circle(self.o, None, self.radius-1))],
            ShapeType.FACE
        )
        # Verify the exception is raised when updating the 'Circle' object
        with self.assertRaises(RuntimeError):
            c.update(wrap_shape(face))

        # Build the face as the result of cutting a circle with another one,
        # thus resulting in a face with two borders
        face = make_cut(
            make_face(make_circle(self.o, None, self.radius)),
            make_face(make_circle(self.o, None, self.radius-1))
        )
        # Verify the exception is raised when updating the 'Circle' object
        with self.assertRaises(RuntimeError):
            c.update(wrap_shape(face))

        # Build the face as the result of cutting a circle with one having
        # a different centre, thus resulting in a face which is not a circular
        # face
        face = make_cut(
            make_face(make_circle(self.o, None, self.radius)),
            make_face(
                make_circle(make_vertex((0.5, 0.5, 0)), None, self.radius-1)
            )
        )
        # Verify the exception is raised when updating the 'Circle' object
        with self.assertRaises(RuntimeError):
            c.update(wrap_shape(face))


class TestRectangle(TestSurface):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Rectangle` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Rectangle` class can be correctly instantiated and that it properly
    handles rotation, scaling, translation, and visualization of their
    geometric elements, as well as their update.
    Tests dealing with the above-mentioned operations are declared in the
    `TestSurface` class this class inherits from.

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
    surface : Rectangle
        The `Surface` subclass that describes the geometric characteristics
        of a rectangle.
    curv_radius : float
        The curvature radius of each rounded corner.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for a `Rectangle` class.
        It initializes the common geometric characteristics for a SALOME
        surface and the ones specific for describing a rectangle.
        The `surface` attribute is assigned to an instance of the `Rectangle`
        class so that the test methods can be run and addressed to the
        correct geometric object.
        """
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.height: float = 1.0
        self.width: float = 1.0
        self.curv_radius: float = 0.2
        self.rounded_corners: List[Tuple[int, float]] = [
            (i, self.curv_radius) for i in range(4)
        ]
        self.name = "Rectangle"
        self.surface: Rectangle = Rectangle(
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
        are correctly set and match the expected reference objects.
        """
        # Build a rectangular shape with the same geometric characteristics
        # for comparison purposes
        rectangle_edges = self.__build_rect_geom_elements(
            self.width, self.height
        )
        # Check the correct instantiation
        self.__assess_instantiation(rectangle_edges, self.o)
        self.assertEqual(self.surface.rot_angle, self.rotation)
        self.assertEqual(self.surface.entry_id, None)
        self.assertEqual(self.surface.name, self.name)

    def test_init_with_rounded_corners(self) -> None:
        """
        Method that tests the initialization of the `Rectangle` object when
        rounded corners are specified.
        It is checked whether the borders of the `Rectangle` object contain
        edges of type `ARC_CIRCLE` with same curvature radius and centre.
        """
        # Instantiate the 'Rectangle' class with rounded corners
        self.surface = Rectangle(
            center=get_point_coordinates(self.o),
            height=self.height,
            width=self.width,
            rounded_corners=self.rounded_corners,
            name=self.name
        )
        # Assess the correct creation of the arc borders
        border_arcs = [
            e for e in self.surface.borders
            if str(get_kind_of_shape(e)[0]) == "ARC_CIRCLE"
        ]
        self.assertEqual(len(border_arcs), len(self.rounded_corners))
        lx, ly = self.surface.dimensions
        for arc in border_arcs:
            # Get the geometric information about each arc
            data = get_kind_of_shape(arc)
            self.assertTrue(get_shape_type(arc) == ShapeType.EDGE)
            for v in extract_sub_shapes(arc, ShapeType.VERTEX):
                self.assertTrue(
                    math.isclose(
                        get_min_distance(v, self.o),
                        math.sqrt(lx**2/4 + (ly/2 - self.curv_radius)**2)
                    ) or
                    math.isclose(
                        get_min_distance(v, self.o),
                        math.sqrt((lx/2 - self.curv_radius)**2 + ly**2/4)
                    )
                )
            self.assertTrue(
                math.isclose(data[7], self.curv_radius, abs_tol=1e-6)
            )

    def test_update(self) -> None:
        """
        Method that tests the implementation of the `update` method for a
        `Rectangle` class.
        A face object representing a rectangle is built and used to update the
        one of a `Rectangle` instance.
        Afterwards, the geometric characteristics of the `Rectangle` object
        are checked to assess they matches the ones of the new face.
        """
        # Dimensions of the new rectangle
        self.width = 2.0
        self.height = 1.0
        # Setup the rectangular shape to update the 'Rectangle' object with
        edges = self.__build_rect_geom_elements(self.width, self.height)
        face = make_face(edges)
        # Update the face and the attributes of the 'Rectangle' object
        self.surface.update(Face(face))

        # Check the instance has been updated correctly
        self.__assess_instantiation(edges, make_cdg(face))

    def test_update_from_non_valid_face(self) -> None:
        """
        Method that tests the implementation of the `update_from_face` method
        for a `Rectangle` class when providing an invalid face object.

        It verifies that an exception is correctly raised when updating the
        `Rectangle` instance with:

        - an object resulting from the partition of two shapes. This object
          has type `ShapeType.COMPOUND` which is not a valid option for the
          `update` method.
        - a face object not having a rectangular shape.
        """
        # Setup the non-valid object to update the 'Rectangle' object with
        face = make_partition(
            [make_face(make_circle(self.o, None, self.width))],
            [make_face(make_circle(self.o, None, self.width-0.1))],
            ShapeType.FACE
        )
        # Verify the exception is raised when updating the 'Rectangle' object
        with self.assertRaises(RuntimeError):
            self.surface.update(wrap_shape(face))
        with self.assertRaises(RuntimeError):
            self.surface.update(
                Face(make_face(make_circle(self.o, None, self.width)))
            )

    def __assess_instantiation(self, edges: List[Any], centre: Any) -> None:
        """
        Method that verifies whether the `Rectangle` class has been correctly
        instantiated by checking that its geometric characteristics match the
        ones used to initialize the instance.

        Parameters
        ----------
        edges : List[Any]
            The list of edge objects of the reference rectangle.
        center : Any
            The vertex object being the center of the reference rectangle.
        """
        self.assertTrue(
            math.isclose(self.surface.dimensions[0], self.width, abs_tol=1e-6)
        )
        self.assertTrue(
            math.isclose(
                self.surface.dimensions[1], self.height, abs_tol=1e-6
            )
        )
        self.assertTrue(
            are_same_shapes(self.surface.o, centre, ShapeType.VERTEX)
        )
        self.assertEqual(len(self.surface.borders), 4)
        for b_rect, b_ref in zip(self.surface.borders, edges):
            self.assertTrue(are_same_shapes(b_rect, b_ref, ShapeType.EDGE))
        self.assertTrue(
            are_same_shapes(self.surface, make_face(edges), ShapeType.FACE)
        )

    def __build_rect_geom_elements(
            self, width: float, height: float
        ) -> List[Any]:
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
        List[Any]
            The list of edges of the rectangle.
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
        return edges


class TestHexagon(TestSurface):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Hexagon` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Hexagon` class can be correctly instantiated and that it properly
    handles rotation, scaling, translation, and visualization of their
    geometric elements, as well as their update.
    Tests dealing with the above-mentioned operations are declared in the
    `TestSurface` class this class inherits from.

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
    surface : Hexagon
        The `Surface` subclass that describes the geometric characteristics
        of a hexagon.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for a `Hexagon` class.
        It initializes the common geometric characteristics for a SALOME
        surface and the ones specific for describing a hexagon.
        The `surface` attribute is assigned to an instance of the `Hexagon`
        class so that the test methods can be run and addressed to the
        correct geometric object.
        """
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.edge_length: float = 1.0
        self.apothem: float = self.edge_length * math.sin(math.pi/3)
        self.name = "Hexagon"
        self.surface: Hexagon = Hexagon(
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
        are correctly set and match the expected reference objects.
        """
        # Build a hexagonal shape with the same geometric characteristics
        # for comparison purposes
        edges = build_hex_geom_elements(self.o, self.edge_length)
        # Check the correct instantiation
        self.__assess_instantiation(edges, self.o)
        self.assertEqual(self.surface.rot_angle, self.rotation)
        self.assertEqual(self.surface.entry_id, None)
        self.assertEqual(self.surface.name, self.name)

    def test_update(self) -> None:
        """
        Method that tests the implementation of the `update` method
        for a `Hexagon` class.
        A face object representing a hexagon is built and used to update the
        one of a `Hexagon` instance.
        Afterwards, the geometric characteristics of the `Hexagon` object
        are checked to assess they matches the ones of the new face.
        """
        # Dimensions of the new hexagon
        self.edge_length = 2.0
        self.apothem = self.edge_length * math.sin(math.pi/3)
        # Setup the hexagonal shape to update the 'Hexagon' object with
        edges = build_hex_geom_elements(self.o, self.edge_length)
        face = make_face(edges)

        # Update the face and the attributes of the 'Hexagon' object
        self.surface.update(Face(face))

        # Check the instance has been updated correctly
        self.__assess_instantiation(edges, make_cdg(face))

    def test_update_from_non_valid_face(self) -> None:
        """
        Method that tests the implementation of the `update` method
        for a `Hexagon` class when providing an invalid face object.
        It verifies that an exception is correctly raised when updating the
        `Hexagon` instance with:

        - an object resulting from the partition of two shapes. This object
          has type `ShapeType.COMPOUND` which is not a valid option for the
          `update` method.
        - a face object not having a hexagonal shape.
        """
        # Setup the non-valid object to update the 'Hexagon' object with
        face = make_partition(
            [make_face(make_circle(self.o, None, self.edge_length))],
            [make_face(make_circle(self.o, None, self.edge_length-0.1))],
            ShapeType.FACE
        )
        # Verify the exception is raised when updating the 'Hexagon' object
        with self.assertRaises(RuntimeError):
            self.surface.update(wrap_shape(face))
        with self.assertRaises(RuntimeError):
            self.surface.update(
                Face(make_face(make_circle(self.o, None, self.edge_length)))
            )

    def __assess_instantiation(
            self, edges: List[Any], centre: Any) -> None:
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
        centre : Any
            The vertex object being the centre of the reference hexagon.
        """
        self.assertTrue(
            math.isclose(
                self.surface.dimensions[0], self.edge_length, abs_tol=1e-6
            )
        )
        self.assertTrue(
            math.isclose(
                self.surface.dimensions[1], self.apothem, abs_tol=1e-6
            )
        )
        self.assertTrue(
            are_same_shapes(self.surface.o, centre, ShapeType.VERTEX)
        )
        self.assertEqual(len(self.surface.borders), 6)
        for b_hex, b_ref in zip(self.surface.borders, edges):
            self.assertTrue(are_same_shapes(b_hex, b_ref, ShapeType.EDGE))
        self.assertTrue(
            are_same_shapes(self.surface, make_face(edges), ShapeType.FACE)
        )


class TestHexagonBuilder(unittest.TestCase):
    """
    Test case for verifying the correct functioning of the function
    `build_hexagon_from_apothem` that allows to build a `Hexagon` instance
    from the apothem of the hexagonal shape.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the function
        `build_hexagon_from_apothem`.
        """
        self.apothem = 1.0
        self.center = (1.0, 1.0, 0.0)
        self.edge_length = self.apothem / math.sin(math.pi/3)
        self.o = make_vertex(self.center)

    def test_build_hexagon_from_apothem(self) -> None:
        """
        Method that verifies the `build_hexagon` function by checking if the
        attributes of the built `Hexagon` instance are correctly set.
        """
        # Build the 'Hexagon' instance
        hex = build_hexagon_from_apothem(self.apothem, self.center)
        # Check the correct instantiation
        self.assertTrue(
            math.isclose(hex.dimensions[0], self.edge_length, abs_tol=1e-6))
        self.assertTrue(
            math.isclose(hex.dimensions[1], self.apothem, abs_tol=1e-6))
        self.assertTrue(
            are_same_shapes(hex.o, self.o, ShapeType.VERTEX)
        )
        self.assertEqual(len(hex.borders), 6)
        self.assertTrue(
            all(
                math.isclose(
                    round(get_basic_properties(b)[0], 6),
                    self.edge_length,
                    abs_tol=1e-6
                ) for b in hex.borders
            )
        )


class TestBuildSurfaceFunctions(unittest.TestCase):
    """
    Test case for verifying the correct functioning of the functions
    of the module `geometries.py` that allow to build a `Surface` instance
    representing specific shapes.
    """
    def test_build_parallelogram(self) -> None:
        """
        Method that verifies the `build_parallelogram` function by checking
        if the attributes of the built `Surface` instance are correctly set
        and represent a parallelogram.
        """
        # Input geometry data
        side_x = 10.0
        side_y = 5.0
        angle_deg = 30.0
        left_corner = (0.0, 0.0, 0.0)
        # Calculate the expected area
        angle_rad = math.radians(angle_deg)
        expected_area = side_x * side_y * math.sin(angle_rad)

        # Build the surface for the parallelogram
        surface = build_parallelogram(
            side_x=side_x,
            side_y=side_y,
            left_corner_angle=angle_deg,
            left_corner=left_corner
        )

        # Validate number of vertices
        vertices = extract_sub_shapes(surface, ShapeType.VERTEX)
        self.assertEqual(len(vertices), 4)

        # Calculate expected coordinates of vertices
        prj_x = side_y * math.cos(angle_rad)
        prj_y = side_y * math.sin(angle_rad)
        expected_coords = [
            (0.0, 0.0, 0.0),
            (side_x, 0.0, 0.0),
            (side_x + prj_x, prj_y, 0.0),
            (prj_x, prj_y, 0.0),
        ]
        # Validate coordinates of vertices
        for c_expected in expected_coords:
            self.assertTrue(
                any(
                    math.isclose(c[0], c_expected[0], abs_tol=1e-6) and
                    math.isclose(c[1], c_expected[1], abs_tol=1e-6) and
                    math.isclose(c[2], c_expected[2], abs_tol=1e-6)
                    for c in [get_point_coordinates(v) for v in vertices]
                )
            )

        # Validate the parallelogram area
        self.assertTrue(
            math.isclose(
                get_basic_properties(surface)[1],
                expected_area,
                abs_tol=1e-6
            )
        )

    def test_build_right_triangle(self) -> None:
        """
        Method that verifies the `build_right_triangle` function by checking
        if the attributes of the built `Surface` instance are correctly set
        and represent a right triangle.
        """
        # Input geometry data
        hyp = 10.0
        left_cat = 6.0
        left_corner = (0.0, 0.0, 0.0)
        # Compute expected geometric quantities, i.e. the right cathetus, its
        # projection on the hypothenuse, the height relative to the hypotenuse
        right_cat = math.sqrt(hyp*hyp - left_cat*left_cat)
        proj = (right_cat*right_cat) / hyp
        height = math.sqrt(right_cat*right_cat - proj*proj)
        # Calculate the expected vertices and area of the right triangle
        expected_vertices = [
            left_corner, (hyp, 0.0, 0.0), (hyp - proj, height, 0.0)
        ]
        expected_area = 0.5 * left_cat * right_cat

        # Build the right triangle
        surface = build_right_triangle(
            hypotenuse=hyp,
            cathetus=left_cat,
            left_corner=left_corner
        )
        # Verify the correctness of the generated surface
        self.__assess_triangle(expected_vertices, expected_area, surface)

    def test_build_right_triangle_from_catheti(self) -> None:
        """
        Method that verifies the `build_right_triangle_from_catheti` function
        by checking if the attributes of the built `Surface` instance are
        correctly set and represent a right triangle.
        """
        # Input geometry data
        horiz_cat = 10.0
        vert_cat = 6.0
        left_corner = (0.0, 0.0, 0.0)
        # Calculate the expected vertices and area of the right triangle
        expected_vertices = [
            left_corner, (horiz_cat, 0.0, 0.0), (horiz_cat, vert_cat, 0.0)
        ]
        expected_area = 0.5 * horiz_cat * vert_cat

        # Build the right triangle
        surface = build_right_triangle_from_catheti(
            cathetus_x=horiz_cat,
            cathetus_y=vert_cat,
            left_corner=left_corner
        )
        # Verify the correctness of the generated surface
        self.__assess_triangle(expected_vertices, expected_area, surface)

    def test_build_regular_triangle(self) -> None:
        """
        Method that verifies the `build_regular_triangle` function by
        checking if the attributes of the built `Surface` instance are
        correctly set and represent a regular triangle.
        """
        # Input geometry data
        side = 10.0
        left_corner = (0.0, 0.0, 0.0)
        height = side*math.sin(math.pi/3)
        # Calculate the expected vertices and area of the regular triangle
        expected_vertices = [
            left_corner, (side, 0.0, 0.0), (side/2, height, 0.0)
        ]
        expected_area = 0.5 * side * height

        # Build the right triangle
        surface = build_regular_triangle(
            side_length=side,
            left_corner=left_corner
        )
        # Verify the correctness of the generated surface
        self.__assess_triangle(expected_vertices, expected_area, surface)

    def __assess_triangle(
            self,
            expected_vertices: List[Tuple[float]],
            expected_area: float,
            surface: Surface
        ) -> None:
        """
        Method that verifies if the given `Surface` object represent a
        triangular shape by checking the number of vertices and their
        coordinates against the expected points, as well as the area
        against the expected value.

        Parameters
        ----------
        expected_vertices : List[Tuple[float]]
            The list of the expected XYZ coordinates of the triangle.
        expected_area : float
            The value of the expected are of the triangle.
        surface : Surface
            The `Surface` object representing the triangle to check.
        """
        # Verify number of vertices and their coordinates
        vertices = extract_sub_shapes(surface, ShapeType.VERTEX)
        self.assertEqual(len(vertices), 3)
        for ev in expected_vertices:
            self.assertTrue(
                any(
                    math.isclose(c[0], ev[0], abs_tol=1e-6) and
                    math.isclose(c[1], ev[1], abs_tol=1e-6) and
                    math.isclose(c[2], ev[2], abs_tol=1e-6)
                    for c in [get_point_coordinates(v) for v in vertices]
                )
            )

        # Verify the surface's area coincides with the expected one
        computed_area = get_basic_properties(surface)[1]
        self.assertTrue(
            math.isclose(computed_area, expected_area, abs_tol=1e-6)
        )


if __name__ == "__main__":
    unittest.main()
