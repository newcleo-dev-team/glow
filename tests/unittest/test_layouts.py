"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.layouts` module have a valid implementation.
"""
import unittest
import math
from typing import Tuple, Any, Dict

from glow.interface.geom_interface import extract_sorted_sub_shapes
from glow.geometry_layouts.geometries import Circle, Rectangle, \
    build_right_triangle_from_catheti
from glow.geometry_layouts.layouts import Layout, Region, LayoutState, \
    associate_colors_to_regions, build_compound_regions, \
    get_unique_values_for_property, is_layout_contained, DEFAULT_REGION_COLOR
from glow.support.types import GeometryType, SymmetryType, PropertyType
from glow.interface.geom_entities import Face, Vertex, Edge
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    extract_sorted_sub_shapes, get_angle_between_shapes, get_min_distance, \
    make_common, make_compound, make_cut, make_fuse, make_partition, \
    make_vector, make_vertex, make_translation, make_vector_from_points, \
    get_basic_properties, make_cdg, get_object_from_id, make_vertex_on_curve
from glow.support.utility import are_same_shapes, build_compound_borders


class ConcreteLayout(Layout):
    """
    Class that mocks a fully concrete implementation of the `Layout`
    abstract class. It is intended for testing purposes only, hence
    no real implementation of the abstract method is provided.
    """
    def rotate(self, angle: float, axis: Edge | None = None) -> None:
        pass

    def scale(self, factor: float, origin: Vertex | None = None) -> None:
        pass

    def show(self, *args: Any) -> None:
        pass

    def translate(self, new_cntr: Tuple[float, float, float]) -> None:
        pass

    def update(self, layout: Face) -> None:
        pass

# -------------------------------------------------------------------------- #
#                               Test Suites
# -------------------------------------------------------------------------- #

class TestLayout(unittest.TestCase):
    """
    Class for testing the abstract `Layout` class behaviour.
    """
    def test_abstract_instantiation(self) -> None:
        """
        Method that tests that the `Layout` class cannot be instantiated
        directly.
        """
        with self.assertRaises(TypeError):
            Layout()

    def test_concrete_instantiation(self) -> None:
        """
        Method that tests the default initialisation of a class that mocks
        a concrete class which derives from the `Layout` abstract class.
        """
        layout = ConcreteLayout()

        # Verify default properties
        self.assertEqual(layout.dimensions, (0.0, 0.0))
        self.assertIsNone(layout.entry_id)
        self.assertEqual(layout.rot_angle, 0.0)

        # Verify that 'o' is successfully initialized as a 'Vertex' at the
        # XYZ origin
        self.assertTrue(
            are_same_shapes(
                layout.o, make_vertex((0.0, 0.0, 0.0)),ShapeType.VERTEX
            )
        )


class TestLayoutState(unittest.TestCase):
    """
    Class for testing the `LayoutState` dataclass behaviour.
    """
    def test_init(self) -> None:
        """
        Method that tests the initialisation of a `LayoutState` instance
        with or without default values.
        """
        state = LayoutState()
        # Verify the initialisation with default values
        self.assertEqual(state.displayed_geom, GeometryType.TECHNOLOGICAL)
        self.assertFalse(state.is_update_needed)
        self.assertEqual(state.symmetry_type, SymmetryType.FULL)

        state = LayoutState(
            displayed_geom=GeometryType.SECTORIZED,
            is_update_needed=False,
            symmetry_type=SymmetryType.DIAG
        )
        # Verify the correct initialisation of its attributes
        self.assertEqual(state.displayed_geom, GeometryType.SECTORIZED)
        self.assertFalse(state.is_update_needed)
        self.assertEqual(state.symmetry_type, SymmetryType.DIAG)


class TestRegion(unittest.TestCase):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Region` class.

    This test suite provides common setup and a set of tests to ensure that
    the `Region` class can be correctly instantiated and that it properly
    handles rotation, scaling, translation, and visualization of their
    geometric elements, as well as their update.

    Attributes
    ----------
    face : Any
        The GEOM face object representing the surface the `Region` refers to.
    name : str
        The name of the `Region` instance.
    o : Vertex
        The `Vertex` object representing the CDG of the surface the `Region`
        refers to.
    properties : Dict[PropertyType, str]
        The properties associated with the `Region` instance.
    region : Region
        The `Region` instance under test.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for a `Region` class.
        """
        # Build the GEOM face the `Region` refers to
        self.face: Any = Rectangle().geom_obj
        self.name: str = "Region"
        self.o : Vertex = Vertex(make_cdg(self.face))
        self.properties: Dict[PropertyType, str] = {
            PropertyType.MATERIAL: "MAT"
        }
        self.region: Region = Region(self.face, self.name, self.properties)

    def test_init(self) -> None:
        """
        Method that tests the initialisation of a `Region` instance with or
        without default values.
        """
        region = Region(self.face)
        # Verify the initialisation with default values
        self.assertEqual(region.color, DEFAULT_REGION_COLOR)
        self.assertEqual(region.region_id, id(region))
        self.assertEqual(region.name, f"Region_{id(region)}")
        self.assertIsNone(region.properties)
        self.assertIsNone(region.entry_id)
        self.assertTrue(
            are_same_shapes(region.o, self.o, ShapeType.VERTEX)
        )

        # Verify the correct initialisation of its attributes
        self.assertTrue(
            are_same_shapes(self.region.geom_obj, self.face, ShapeType.FACE)
        )
        self.assertEqual(self.region.name, self.name)
        self.assertEqual(self.region.properties, self.properties)

        # Verify an exception is raised when instantiating a 'Region' without
        # a GEOM face or by providing a GEOM object not being a face
        with self.assertRaises(RuntimeError):
            Region(None)
        with self.assertRaises(RuntimeError):
            Region(make_compound(self.face))

    def test_clone(self) -> None:
        """
        Method that tests the implementation of the `clone` method for a
        `Region` class. It verifies the cloned region shares the same ID,
        properties and colour of the source region.
        """
        # Modify source region attributes to check modifications are kept in
        # the cloned region
        self.region.color = (255, 0, 0)
        self.region.name = "Source region"
        self.region.entry_id = "0:1:1"

        # Clone the source region
        cloned_region = self.region.clone()
        # Verify the two regions are distinct objects sharing the same
        # attributes
        self.assertIsNot(cloned_region, self.region)
        self.assertEqual(cloned_region.color, self.region.color)
        self.assertEqual(cloned_region.name, self.region.name)
        self.assertEqual(cloned_region.region_id, self.region.region_id)
        self.assertIsNone(cloned_region.entry_id)
        self.assertIsNot(cloned_region.properties, self.region.properties)
        for p1, p2 in zip(cloned_region.properties, self.region.properties):
            self.assertEqual(p1, p2)

    def test_reset_region_color(self) -> None:
        """
        Method that tests the implementation of the `reset_region_color`
        method for a `Region` class. It verifies that the region's colour
        is set back to the default value.
        """
        # Modify the region 'color' attribute
        self.region.color = (255, 0, 0)
        # Verify the region's colour is the default one after resetting it
        self.region.reset_region_color()
        self.assertEqual(self.region.color, DEFAULT_REGION_COLOR)

    def test_rotate(self) -> None:
        """
        Method that tests the implementation of the `rotate` method for a
        `Region` class. It verifies that the region's shape is correctly
        rotated wrt to the X-axis.
        """
        # Build reference vector before the rotation
        borders = build_compound_borders(self.region)
        face_ref_vect = make_vector_from_points(
            self.region.o, make_vertex_on_curve(borders[0], 0.0)
        )
        # Rotate the surface
        rot_angle = 90.0
        self.region.rotate(rot_angle)
        # Build reference vector after the rotation
        borders = build_compound_borders(self.region)
        face_ref_vect2 = make_vector_from_points(
            self.region.o, make_vertex_on_curve(borders[0], 0.0)
        )
        # Check the correct rotation happened
        self.assertEqual(self.region.rot_angle, rot_angle)
        self.assertTrue(
            math.isclose(
                get_angle_between_shapes(face_ref_vect, face_ref_vect2),
                rot_angle,
                abs_tol=1e-6
            )
        )

        # Check no rotation happened if the rotation angle is zero
        self.region.rotate(0.0)
        self.assertAlmostEqual(self.region.rot_angle, rot_angle)

    def test_rotate_from_axis(self) -> None:
        """
        Method that tests the implementation of the `rotate` method for a
        `Region` class when an axis of rotation is specified.
        """
        # Build the rotation axis in a given position
        axis_centre = (1, 1, 0)
        axis_centre_vrtx = make_vertex(axis_centre)
        axis = make_vector_from_points(
            axis_centre_vrtx, make_vertex((*axis_centre[:2], 1))
        )
        # Build reference vector before the rotation (from region centre to
        # axis centre)
        ref_vect = make_vector_from_points(self.region.o, axis_centre_vrtx)

        # Rotate the region around the specified axis
        rot_angle = 90.0
        self.region.rotate(rot_angle, axis)
        # Build reference vector after the rotation
        ref_vect2 = make_vector_from_points(self.region.o, axis_centre_vrtx)
        # Check the correct rotation happened
        self.assertEqual(self.region.rot_angle, rot_angle)
        self.assertTrue(
            math.isclose(
                get_angle_between_shapes(ref_vect, ref_vect2),
                rot_angle,
                abs_tol=1e-6
            )
        )

    def test_scale(self) -> None:
        """
        Method that tests the scaling of a `Region` object by checking the
        area of the surface before and after the operation.
        """
        scale_factor = 2
        area_before = get_basic_properties(self.region)[1]
        # Scale the surface wrt to its centre
        self.region.scale(scale_factor)
        # Check the area has increased by the scaling factor^2
        self.assertTrue(
            math.isclose(
                area_before*scale_factor*scale_factor,
                get_basic_properties(self.region)[1],
                abs_tol=1e-6
            )
        )

        # Verify the exception is raised if an invalid scaling factor is
        # provided
        with self.assertRaises(ValueError):
            self.region.scale(0.0)
        with self.assertRaises(ValueError):
            self.region.scale(-2.0)

    def test_set_region_color(self) -> None:
        """
        Method that tests the implementation of the `set_region_color` method
        for a `Region` class. It verifies the `color` attribute is correctly
        set.
        """
        # Test valid color
        new_color = (255, 0, 0)
        self.region.set_region_color(new_color)
        self.assertEqual(self.region.color, new_color)

        # Test invalid colors
        invalid_inputs = [(300, 0, 0), (0, -1, 0), (0, 0), (0, 0, 0, 0)]
        for c in invalid_inputs:
            with self.subTest(color=c):
                with self.assertRaises(ValueError):
                    self.region.set_region_color(c)

    def test_show(self) -> None:
        """
        Method that tests the implementation of the `show` method for a
        `Region` class.
        It ensures that the face is correctly displayed in SALOME by:

        - verifying that SALOME assigns an entry ID to the face object;
        - checking that the geometric object the entry ID corresponds to
          matches the expected shape.
        - verifying the entry ID changes after showing the region multiple
          times.
        - verifying an exception is raised if calling the method with any
          arguments.
        """
        self.region.show()
        self.assertIsNotNone(self.region.entry_id)
        self.assertTrue(
            are_same_shapes(
                self.region,
                get_object_from_id(self.region.entry_id),
                ShapeType.FACE
            )
        )

        # Test show with existing entry
        old_id = self.region.entry_id
        self.region.show()
        self.assertNotEqual(self.region.entry_id, old_id)

        # Test show with args
        with self.assertRaises(ValueError):
            self.region.show(1)

    def test_translate(self) -> None:
        """
        Method that tests the implementation of the `translate` method for a
        `Region` class. It verifies the correctness of the resulting position
        of its geometric elements.
        """
        # Store the elements for assessing the correct translation
        center_after_transl = (1, 1, 0)
        new_pos_vrtx = make_vertex(center_after_transl)
        cdg_pre = make_cdg(self.region)
        distance = get_min_distance(self.region.o, new_pos_vrtx)

        # Translate the surface
        self.region.translate(center_after_transl)
        # Check the correct translation happened
        self.assertTrue(
            are_same_shapes(self.region.o, new_pos_vrtx, ShapeType.VERTEX)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre, make_cdg(self.region)),
                distance,
                abs_tol=1e-6
            )
        )

        # Check no translation happened if the new centre coincides with the
        # current one
        self.region.translate(center_after_transl)
        self.assertTrue(
            are_same_shapes(self.region.o, new_pos_vrtx, ShapeType.VERTEX)
        )

    def test_boolean_subtraction(self) -> None:
        """
        Method that tests the cut functionality with which the `-` operator
        has been overloaded.
        """
        # Cut the region with a rectangular face representing a hole
        hole = Region(Rectangle(height=0.5, width=0.5), name="Hole")
        # Resulting region
        result = self.region - hole
        # Verify the result of the cut operation
        cut = make_cut(self.region, hole)
        self.assertTrue(
            are_same_shapes(cut, result, ShapeType.FACE)
        )
        self.assertAlmostEqual(
            get_basic_properties(result)[1],
            get_basic_properties(self.region)[1]
            - get_basic_properties(hole)[1]
        )
        self.assertEqual(result.properties, self.properties)

    def test_boolean_addition(self) -> None:
        """
        Method that tests the fuse functionality with which the `+` operator
        has been overloaded.
        """
        # Verify fuse operation of regions with identical properties
        r1 = Region(
            Rectangle((1, 0, 0)), properties={PropertyType.MATERIAL: "MAT"}
        )
        fused = self.region + r1
        self.assertIsInstance(fused, Region)
        self.assertTrue(
            are_same_shapes(
                fused,
                extract_sub_shapes(
                    make_fuse([self.region, r1]), ShapeType.FACE
                )[0],
                ShapeType.FACE
            )
        )
        self.assertAlmostEqual(
            get_basic_properties(fused)[1],
            get_basic_properties(self.region)[1]
            + get_basic_properties(r1)[1]
        )
        self.assertEqual(fused.properties[PropertyType.MATERIAL], "MAT")

        # Verify fuse operation of multiple regions with identical properties
        r2 = Region(
            Rectangle((1, 1, 0)), properties={PropertyType.MATERIAL: "MAT"}
        )
        fused = self.region + [r1, r2]
        self.assertIsInstance(fused, Region)
        self.assertTrue(
            are_same_shapes(
                fused,
                extract_sub_shapes(
                    make_fuse([self.region, r1, r2]), ShapeType.FACE
                )[0],
                ShapeType.FACE
            )
        )
        self.assertAlmostEqual(
            get_basic_properties(fused)[1],
            get_basic_properties(self.region)[1]
            + get_basic_properties(r1)[1]
            + get_basic_properties(r2)[1]
        )
        self.assertEqual(fused.properties[PropertyType.MATERIAL], "MAT")

        # 3. Failure branch: Property mismatch
        r3 = Region(r1.geom_obj, properties={PropertyType.MATERIAL: "MAT2"})
        with self.assertRaises(RuntimeError):
            _ = r1 + r3

    def test_boolean_common(self) -> None:
        """
        Method that tests the common functionality with which the `*` operator
        has been overloaded.
        """
        # Get the common part between the region and a rectangular surface
        # with (0.5, 0.5, 0) as centre, and two properties
        r1 = Region(
            Rectangle((0.5, 0.5, 0)),
            properties={
                PropertyType.MACRO: "MAC001", PropertyType.MATERIAL: "MAT1"
            }
        )
        # Resulting region
        result = self.region * r1
        # Verify the result of the common operation
        self.assertTrue(
            are_same_shapes(
                make_common(self.region, r1),
                result,
                ShapeType.FACE
            )
        )
        self.assertAlmostEqual(
            get_basic_properties(result)[1], 0.5*0.5
        )
        self.assertEqual(result.properties, r1.properties)

class TestLayoutFunctions(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the functions in
    the `layouts.py` module.

    Attributes
    ----------
    r1 : Region
        `Region` object for a rectangular surface with 'MAT1' material.
    r2 : Region
        `Region` object for a circular surface with 'MAT2' material.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the functions of the
        `layouts.py` module.
        """
        # Declare two regions with different values for the 'MATERIAL'
        # property
        self.r1 = Region(
            Rectangle(), properties={PropertyType.MATERIAL: "MAT1"}
        )
        self.r2 = Region(
            Circle(radius=0.25), properties={PropertyType.MATERIAL: "MAT2"}
        )

    def test_associate_colors_to_regions(self) -> None:
        """
        Method that verifies the `associate_colors_to_regions` function of the
        `layouts.py` module.
        """
        ref_regions = [self.r1, self.r2]
        # Verify that the indicated regions have two different colours
        associate_colors_to_regions(PropertyType.MATERIAL, ref_regions)
        self.assertNotEqual(self.r1.color, DEFAULT_REGION_COLOR)
        self.assertNotEqual(self.r2.color, DEFAULT_REGION_COLOR)
        self.assertNotEqual(self.r1.color, self.r2.color)

        # Verify that regions colours are set to the default value if the
        # function is called without specifying any property type
        associate_colors_to_regions(None, ref_regions)
        self.assertEqual(self.r1.color, DEFAULT_REGION_COLOR)
        self.assertEqual(self.r2.color, DEFAULT_REGION_COLOR)

    def test_build_compound_regions(self) -> None:
        """
        Method that verifies the `build_compound_regions` function of the
        `layouts.py` module.
        """
        ref_regions = [self.r1, self.r2]
        # Declare a compound made by the intersection of the two regions with
        # a triangular shape resulting in a halved compound
        cmpd = make_partition(ref_regions, [], ShapeType.FACE)
        cmpd = make_common(
            cmpd,
            build_right_triangle_from_catheti(1.0, 1.0, (-0.5, -0.5, 0.0))
        )
        # Get the faces for comparison purposes
        cmpd_faces = extract_sorted_sub_shapes(cmpd, ShapeType.FACE)

        # Build the compound regions
        cmpd_regions = build_compound_regions(cmpd, ref_regions)
        # Verify that the regions corresponds to the compound faces and that
        # they share the same properties of the source regions
        for cmpd_r in cmpd_regions:
            found = False
            for cmpd_f in cmpd_faces:
                if are_same_shapes(cmpd_r, cmpd_f, ShapeType.FACE):
                    found = True
                    break
            self.assertTrue(found)

        for cmpd_r in cmpd_regions:
            found = False
            for r in ref_regions:
                if cmpd_r.properties == r.properties:
                    found = True
                    break
            self.assertTrue(found)

        # Verify an exception is raised when any of the compound faces does
        # not have a corresponding region
        cmpd = make_translation(cmpd, make_vector((10, 10, 0)))
        with self.assertRaises(RuntimeError):
            _ = build_compound_regions(cmpd, ref_regions)

        # Verify an exception is raised when the compound does not have any
        # face
        with self.assertRaises(RuntimeError):
            _ = build_compound_regions(make_compound([]), ref_regions)

    def test_get_unique_values_for_property(self) -> None:
        """
        Method that verifies the `get_unique_values_for_property` function
        of the `layouts.py` module.
        """
        # Declare a third region with one of the property values already used
        r3 = Region(Rectangle(), properties={PropertyType.MATERIAL: "MAT1"})
        ref_regions = [self.r1, self.r2, r3]
        ref_mat_values = ["MAT1", "MAT2"]

        # Verify that the returned list of properties contains unique names
        values = get_unique_values_for_property(
            PropertyType.MATERIAL, ref_regions
        )
        self.assertEqual(len(values), 2)
        for r in ref_regions:
            found = False
            for mat in ref_mat_values:
                if r.properties[PropertyType.MATERIAL] == mat:
                    found = True
                    break
            self.assertTrue(found)

        # Verify an exception is raised when getting the unique property
        # values for a missing property type
        with self.assertRaises(RuntimeError):
            _ = get_unique_values_for_property(
                PropertyType.MACRO, ref_regions
            )
        # Verify all the regions' colours have been set to red
        for r in ref_regions:
            self.assertEqual(r.color, (255, 0, 0))

    def test_is_layout_contained(self) -> None:
        """
        Method that verifies the `is_layout_contained` function of the
        `layouts.py` module.
        """
        # Declare a surface greater than the stored regions
        container = Rectangle(width=10, height=10)

        # Verify that the surface contains any of the stored regions
        for r in [self.r1, self.r2]:
            self.assertTrue(
                is_layout_contained(container, r)
            )

        # Verify that a region is not contained if its bounding box is not
        # within that of the container
        r = self.r1.clone()
        r.translate((11, 11, 0))
        self.assertFalse(is_layout_contained(container, r))


if __name__ == "__main__":
    unittest.main()