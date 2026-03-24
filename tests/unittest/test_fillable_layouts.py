"""
Module containing unittest class to assess that the `Fillable` class of the
`glow.geometry_layouts.fillable_layouts` module have a valid implementation.
"""
import math
import unittest

from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.fillable_layouts import Fillable
from glow.geometry_layouts.geometries import Circle
from glow.geometry_layouts.lattices import CartesianLattice
from glow.geometry_layouts.layouts import DEFAULT_REGION_COLOR, Region
from glow.interface.geom_entities import Compound, Vertex
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_angle_between_shapes, get_basic_properties, get_id_from_object, \
    get_min_distance, get_object_from_id, get_shape_type, make_cdg, \
    make_common, make_compound, make_cut, make_vector_from_points, \
    make_vertex, make_vertex_on_curve
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.support.utility import are_same_shapes, build_compound_borders
from glow.interface.geom_entities import wrap_shape
from glow.interface.geom_interface import make_circle
from tests.unittest.support_funcs import capture_output


class TestFillable(unittest.TestCase):
    """
    Base test case for verifying the construction and transformation
    operations, as well as the visualization capabilities of the `Fillable`
    subclasses.

    This test suite provides common setup and a set of tests common to all
    the tests that ensure the correct behaviour of the sublclasses of
    `Fillable`.

    Parameters
    ----------
    fillable : Fillable
        The instance of the `Fillable` subclasses to test. It is redeclared
        in the test cases for each of the `Fillable` subclasses.
    o : Vertex
        A `Vertex` object representing a point in the XYZ origin.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the `Fillable` class
        in the `fillable_layouts.py` module, and of its subclasses.
        """
        self.fillable: Fillable | None = None
        self.o : Vertex = Vertex(make_vertex((0.0, 0.0, 0.0)))

    def test_abstract_instantiation(self) -> None:
        """
        Method that tests that the `Fillable` class cannot be instantiated
        directly.
        """
        with self.assertRaises(TypeError):
            Fillable()
    def test_add_region_without_position(self) -> None:
        """
        Method that tests adding a `Region` to a `Fillable` without specifying
        a position. The region should be placed at the `Fillable`'s centre.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create a region
        layout = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})
        initial_layers_no = len(self.fillable.layers)
        # Add the region without position
        self.fillable.add(layout)

        # Verify the region was added to a new layer
        self.assertEqual(len(self.fillable.layers), initial_layers_no + 1)
        self.assertEqual(len(self.fillable.layers[-1]), 1)
        # Verify the region has been added to the fillable's centre
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[-1][0].o,
                self.fillable.o,
                ShapeType.VERTEX
            )
        )

        # Verify the update flag is set
        self.assertTrue(self.fillable.state.is_update_needed)

    def test_add_region_with_position(self) -> None:
        """
        Method that tests adding a Region to a Fillable with a specific
        position. The region should be translated to the given position.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region' to add
        layout = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})
        position = (1.0, 2.0, 0.0)
        initial_layers_no = len(self.fillable.layers)
        # Add the region with a specific position
        self.fillable.add(layout, position=position)

        # Verify the region was added to a new layer
        self.assertEqual(len(self.fillable.layers), initial_layers_no + 1)
        self.assertEqual(len(self.fillable.layers[-1]), 1)
        self.assertTrue(self.fillable.state.is_update_needed)
        # Verify the region has been added to the indicated position
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[-1][0].o,
                make_vertex(position),
                ShapeType.VERTEX
            )
        )

    def test_add_to_existing_layer(self) -> None:
        """
        Method that tests adding multiple regions to the same layer.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region's to add
        layout1 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})
        layout2 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT2"})
        initial_layers_no = len(self.fillable.layers)
        initial_layer_0_no = len(self.fillable.layers[0])

        # Add first region to an existing layer
        self.fillable.add(layout1, layer_index=0)
        self.assertEqual(len(self.fillable.layers), initial_layers_no)
        self.assertEqual(len(self.fillable.layers[0]), initial_layer_0_no + 1)

        # Add second region to the same layer
        self.fillable.add(layout2, layer_index=0)
        self.assertEqual(len(self.fillable.layers), initial_layers_no)
        self.assertEqual(len(self.fillable.layers[0]), initial_layer_0_no + 2)

    def test_add_invalid_layer_index(self) -> None:
        """
        Method that tests that adding a region with an invalid layer index
        raises a ValueError.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region' to add
        layout = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})

        # Try to add to a non-existent layer with invalid index
        with self.assertRaises(ValueError):
            self.fillable.add(layout, layer_index=5)

    def test_add_region_name_includes_parent(self) -> None:
        """
        Method that tests that the added region's name includes the parent
        `Fillable`'s name.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region' to add
        layout = Region(
            Circle(),
            name="test_region",
            properties={PropertyType.MATERIAL: "MAT"}
        )
        original_name = layout.name

        # Add the region
        self.fillable.add(layout)
        # Verify the region name was updated
        added_region = self.fillable.layers[0][0]
        self.assertIn(self.fillable.name, added_region.name)
        self.assertIn(original_name, added_region.name)
    def test_add_region_without_position(self) -> None:
        """
        Method that tests adding a `Region` to a `Fillable` without specifying
        a position. The region should be placed at the `Fillable`'s centre.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create a region
        layout = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})
        initial_layers_no = len(self.fillable.layers)
        # Add the region without position
        self.fillable.add(layout)

        # Verify the region was added to a new layer
        self.assertEqual(len(self.fillable.layers), initial_layers_no + 1)
        self.assertEqual(len(self.fillable.layers[-1]), 1)
        # Verify the region has been added to the fillable's centre
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[-1][0].o,
                self.fillable.o,
                ShapeType.VERTEX
            )
        )

        # Verify the update flag is set
        self.assertTrue(self.fillable.state.is_update_needed)

    def test_add_region_with_position(self) -> None:
        """
        Method that tests adding a Region to a Fillable with a specific
        position. The region should be translated to the given position.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region' to add
        layout = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})
        position = (1.0, 2.0, 0.0)
        initial_layers_no = len(self.fillable.layers)
        # Add the region with a specific position
        self.fillable.add(layout, position=position)

        # Verify the region was added to a new layer
        self.assertEqual(len(self.fillable.layers), initial_layers_no + 1)
        self.assertEqual(len(self.fillable.layers[-1]), 1)
        self.assertTrue(self.fillable.state.is_update_needed)
        # Verify the region has been added to the indicated position
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[-1][0].o,
                make_vertex(position),
                ShapeType.VERTEX
            )
        )

    def test_add_to_existing_layer(self) -> None:
        """
        Method that tests adding multiple regions to the same layer.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region's to add
        layout1 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})
        layout2 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT2"})
        initial_layers_no = len(self.fillable.layers)
        initial_layer_0_no = len(self.fillable.layers[0])

        # Add first region to an existing layer
        self.fillable.add(layout1, layer_index=0)
        self.assertEqual(len(self.fillable.layers), initial_layers_no)
        self.assertEqual(len(self.fillable.layers[0]), initial_layer_0_no + 1)

        # Add second region to the same layer
        self.fillable.add(layout2, layer_index=0)
        self.assertEqual(len(self.fillable.layers), initial_layers_no)
        self.assertEqual(len(self.fillable.layers[0]), initial_layer_0_no + 2)

    def test_add_invalid_layer_index(self) -> None:
        """
        Method that tests that adding a region with an invalid layer index
        raises a ValueError.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region' to add
        layout = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})

        # Try to add to a non-existent layer with invalid index
        with self.assertRaises(ValueError):
            self.fillable.add(layout, layer_index=5)

    def test_add_region_name_includes_parent(self) -> None:
        """
        Method that tests that the added region's name includes the parent
        `Fillable`'s name.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare the 'Region' to add
        layout = Region(
            Circle(),
            name="test_region",
            properties={PropertyType.MATERIAL: "MAT"}
        )
        original_name = layout.name

        # Add the region
        self.fillable.add(layout)
        # Verify the region name was updated
        added_region = self.fillable.layers[-1][0]
        self.assertIn(self.fillable.name, added_region.name)
        self.assertIn(original_name, added_region.name)

    def test_apply_symmetry(self) -> None:
        """
        Method that tests applying different symmetry types to a `Fillable`
        instance.
        The symmetry map should be populated and the state should be updated.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Test the FULL, HALF, QUARTER symmetry types
        self.__assess_symmetry(SymmetryType.FULL)
        self.__assess_symmetry(SymmetryType.HALF)
        self.__assess_symmetry(SymmetryType.QUARTER)

    def test_apply_symmetry_updates_geom_obj(self) -> None:
        """
        Method that tests that applying symmetry updates the GEOM object
        if it is None before the operation.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Ensure geom_obj is None initially and the state needs an update
        self.fillable.geom_obj = None
        self.fillable.state.is_update_needed = True
        # Apply symmetry
        self.fillable.apply_symmetry(SymmetryType.FULL)

        # Verify geom_obj was updated
        self.assertIsNotNone(self.fillable.geom_obj)
        self.assertEqual(self.fillable.state.symmetry_type, SymmetryType.FULL)

    def test_apply_symmetry_respects_rotation(self) -> None:
        """
        Method that tests that applying symmetry respects the rotation angle
        of the Fillable instance by rotating the symmetry shape accordingly.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Set a non-zero rotation angle
        self.fillable.rot_angle = 45.0
        # Apply symmetry
        self.fillable.apply_symmetry(SymmetryType.HALF)

        # Verify the symmetry was applied
        self.assertIn(SymmetryType.HALF, self.fillable.symmetry_map)
        self.assertEqual(self.fillable.state.symmetry_type, SymmetryType.HALF)
        # Verify the rotation angle is still set
        self.assertEqual(
            self.fillable.symmetry_map[SymmetryType.HALF].rot_angle, 45.0
        )

    def test_clone(self) -> None:
        """
        Method that tests the implementation of the `clone` method for a
        `Fillable` subclass.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Clone the source fillable
        cloned_fillable = self.fillable.clone()
        # Verify the two fillables are distinct objects
        self.assertIsNot(cloned_fillable, self.fillable)

    def test_get_geometry_map(self) -> None:
        """
        Method that tests the implementation of the `get_geometry_map` method
        for a `Fillable` subclass.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()

        # Verify the exception is raised when no geometry map is stored
        with self.assertRaises(RuntimeError):
            self.fillable.get_geometry_map(GeometryType.SECTORIZED)

        # Create a mock compound for the 'SECTORIZED' geometry map
        mock_geom = Compound(
            make_compound([
                make_circle(make_vertex((0.0, 0.0, 0.0)), None, 1.0),
                make_circle(make_vertex((0.0, 0.0, 0.0)), None, 2.0)
            ])
        )
        self.fillable.geometry_maps[GeometryType.SECTORIZED] = mock_geom
        # Get the 'SECTORIZED' geometry map
        result = self.fillable.get_geometry_map(GeometryType.SECTORIZED)
        # Verify the returned geometry map is the one stored
        self.assertIs(result, mock_geom)

    def test_get_geometry_map_nested_layers(self) -> None:
        """
        Method that tests the `get_geometry_map` method when no geometry map
        is stored, so it is built from those of the nested `Fillable` objects
        in layers.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create a fillable with a nested structure
        cell = CartesianCell()
        cell.sectorize([8], [0])
        lattice = CartesianLattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        self.fillable = CartesianCell(width_height=(4, 4))
        self.fillable.add(lattice)

        # Try to get a geometry map that still doesn't exist so to trigger
        # its construction from nested fillables
        result = self.fillable.get_geometry_map(GeometryType.SECTORIZED)
        # Verify the result is a Compound of 72 edges (8 edges per 9 cells)
        self.assertIsInstance(result, Compound)
        self.assertEqual(
            len(extract_sub_shapes(result, ShapeType.EDGE)), 72
        )

    def test_get_geometry_map_empty(self) -> None:
        """
        Method that tests that `get_geometry_map` raises RuntimeError when
        the resulting compound has no edges (empty mapping).
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Ensure no geometry map is stored for the requested type
        self.fillable.restore()
        # Verify that an exception is raised as no geometry map can be
        # retrieved
        with self.assertRaises(RuntimeError):
            self.fillable.get_geometry_map(GeometryType.SECTORIZED)

    def test_get_regions(self) -> None:
        """
        Method that tests the `get_regions` method for `Fillable` subclasses.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Ensure the fillable has exactly three regions (characteristic
        # surface + two circular regions)
        self.fillable.restore()
        layout1 = Region(
            Circle(radius=2), properties={PropertyType.MATERIAL: "MAT1"}
        )
        layout2 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT2"})
        self.fillable.add(layout1)
        self.fillable.add(layout2)

        # Get regions from the fillable
        regions = self.fillable.get_regions()

        # Verify that all regions are returned
        self.assertEqual(len(regions), 3)
        self.assertTrue(
            are_same_shapes(regions[0], self.fillable.shape, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(regions[1], layout1, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(regions[2], layout2, ShapeType.FACE)
        )

    def test_get_regions_nested_fillables(self) -> None:
        """
        Method that tests the `get_regions` method when nested `Fillable`
        objects are present in the layers.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create a fillable with a nested structure
        cell = CartesianCell()
        lattice = CartesianLattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        self.fillable = CartesianCell(width_height=(4, 4))
        self.fillable.add(lattice)

        # Get regions from the parent fillable
        regions = self.fillable.get_regions()

        # Verify that all nested regions are returned
        self.assertEqual(len(regions), 10)

    def test_get_regions_with_symmetry(self) -> None:
        """
        Method that tests the `get_regions_with_symmetry` method for
        `Fillable` subclasses. Different symmetry types are applied and
        checked that the number of regions coincides with the number of
        faces of the compound being the common part between the full
        geometry layout and the shape of the symmetry.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Apply the 'QUARTER' common symmetry type
        symm_type = SymmetryType.QUARTER
        self.fillable.apply_symmetry(symm_type)
        self.fillable.update_hierarchical_structure()

        # Verify that the number of regions coincides with the number of
        # faces
        common = (
            wrap_shape(self.fillable.geom_obj)
            * self.fillable.symmetry_map[symm_type]
        )
        no_faces = len(
            extract_sub_shapes(make_compound([common]), ShapeType.FACE)
        )
        regions = self.fillable.get_regions_with_symmetry(symm_type)
        self.assertEqual(len(regions), no_faces)

        # Verify that the returned number of regions coincides with those
        # of the full layout when indicating a 'FULL' symmetry type
        no_faces = len(
            extract_sub_shapes(make_compound(self.fillable), ShapeType.FACE)
        )
        regions = self.fillable.get_regions_with_symmetry(SymmetryType.FULL)
        self.assertEqual(len(regions), no_faces)

        # Verify the exception is raised when the indicated symmetry type is
        # not present
        with self.assertRaises(RuntimeError):
            self.fillable.get_regions_with_symmetry(SymmetryType.DIAG)

    def test_print_region_info(self) -> None:
        """
        Method that tests the `get_regions_with_symmetry` method for
        `Fillable` subclasses. It verifies that:

        - the exception is raised when no shape is provided and nothing
          is selected in the study;
        - the output about a region of the fillable is correctly printed to
          the stdout;
        - in case of a region of the fillable without properties, the output
          does not contain any property information.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Verify the exception is raised when no shape is provided and nothing
        # is selected in the study
        with self.assertRaises(RuntimeError):
            self.fillable.print_region_info()

        # Declare a reference region and add it to the geometry layout
        region = Region(
            Circle(),
            "Region_1",
            {PropertyType.MACRO: "MAC001", PropertyType.MATERIAL: "MAT1"}
        )
        self.fillable.add(region)
        # Update the hierarchical structure to populate regions
        self.fillable.update_hierarchical_structure()
        # Verify output was printed
        captured = capture_output(
            self.fillable.print_region_info, region.geom_obj
        )
        self.assertIn(
            f"Properties of '{self.fillable.name}_{region.name}':", captured
        )
        self.assertIn("   MATERIAL: MAT1\n", captured)
        self.assertIn("   MACRO: MAC001\n", captured)

        # Clear the properties of the region of the fillable (last inserted
        # object)
        self.fillable.layers[-1][-1].properties.clear()
        # Verify the output does not indicate any property
        captured = capture_output(
            self.fillable.print_region_info, region.geom_obj
        )
        self.assertIn(
            f"Properties of '{self.fillable.name}_{region.name}':", captured
        )
        self.assertIn("   No associated properties.", captured)

    def test_print_region_info_invalid_shape(self) -> None:
        """
        Method that tests the `print_region_info` method when an shape that
        does not correspond to any region in the layout is provided.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create an invalid shape that is not in the layout
        invalid_shape = make_circle(make_vertex((5.0, 5.0, 0.0)), None, 1.0)
        # Verify the exception is raised for the invalid shape
        with self.assertRaises(RuntimeError):
            self.fillable.print_region_info(invalid_shape)

    def test_print_region_info_nested_fillable(self) -> None:
        """
        Method that tests the `print_region_info` method with nested
        `Fillable` objects in the hierarchical structure.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create a nested structure with regions
        cell = CartesianCell()
        region = Region(
            Circle(), "NestedRegion", {PropertyType.MATERIAL: "MAT1"}
        )
        cell.add(region)
        # Add the nested cell to the main fillable
        self.fillable.add(cell)
        # Update the hierarchical structure
        self.fillable.update_hierarchical_structure()

        # Verify output was printed with correct tree path information
        captured = capture_output(
            self.fillable.print_region_info, region.geom_obj
        )
        last_fillable_layer = len(self.fillable.layers)-1
        pos_in_fillable_layer = len(
            self.fillable.layers[last_fillable_layer]) - 1
        last_cell_layer = len(cell.layers) - 1
        pos_in_cell_layer = len(cell.layers[last_cell_layer]) - 1
        self.assertIn(
            "Access Code: <root-instance-name>.layers["
            f"{last_fillable_layer}][{pos_in_fillable_layer}]"
            f".layers[{last_cell_layer}][{pos_in_cell_layer}]",
            captured
        )
        self.assertIn(
            f"Properties of '{cell.name}_{region.name}':", captured
        )
        self.assertIn("   MATERIAL: MAT1\n", captured)

    def test_restore(self) -> None:
        """
        Method that tests the `restore` method for `Fillable` subclasses.
        It verifies that the fillable's layout contains one layer with only
        the region from the characteristic shape, while all the mappings
        are cleared
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Add a cell with a sectorised layout and apply a symmetry so that
        # the corresponding mappings are not empty
        cell = CartesianCell()
        cell.sectorize([8], [0])
        self.fillable.add(cell)
        self.fillable.apply_symmetry(SymmetryType.QUARTER)
        self.assertTrue(
            get_shape_type(
                self.fillable.get_geometry_map(GeometryType.SECTORIZED)
            ) == ShapeType.COMPOUND
        )
        self.assertTrue(
            get_shape_type(
                self.fillable.symmetry_map[SymmetryType.QUARTER]
            ) == ShapeType.FACE
        )
        # Get the fillable's shape
        shape = self.fillable.shape

        # Restore the fillable's layout
        self.fillable.restore()
        # Verify that the fillable has been correctly restored
        self.assertEqual(self.fillable.geometry_maps, {})
        self.assertEqual(self.fillable.symmetry_map, {})
        self.assertEqual(len(self.fillable.regions), 1)
        self.assertEqual(len(self.fillable.layers), 1)
        self.assertEqual(len(self.fillable.layers[0]), 1)
        self.assertTrue(
            are_same_shapes(
                self.fillable.regions[0], shape, ShapeType.FACE
            )
        )
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[0][0], shape, ShapeType.FACE
            )
        )

    def test_rotate(self) -> None:
        """
        Method that tests the implementation of the `rotate` method for
        `Fillable` subclassess. It verifies that all the layouts in the
        fillable's layers are correctly rotated wrt to the X-axis.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Build reference vector before the rotation
        borders = build_compound_borders(self.fillable)
        face_ref_vect = make_vector_from_points(
            self.fillable.o, make_vertex_on_curve(borders[0], 0.0)
        )
        # Rotate the fillable
        rot_angle = 90.0
        self.fillable.rotate(rot_angle)
        # Build reference vector after the rotation
        borders = build_compound_borders(self.fillable)
        face_ref_vect2 = make_vector_from_points(
            self.fillable.o, make_vertex_on_curve(borders[0], 0.0)
        )
        # Check the correct rotation happened
        self.assertEqual(self.fillable.rot_angle, rot_angle)
        self.assertTrue(
            math.isclose(
                get_angle_between_shapes(face_ref_vect, face_ref_vect2),
                rot_angle,
                abs_tol=1e-6
            )
        )

        # Check no rotation happened if the rotation angle is zero
        self.fillable.rotate(0.0)
        self.assertAlmostEqual(self.fillable.rot_angle, rot_angle)


    def test_rotate_from_axis(self) -> None:
        """
        Method that tests the implementation of the `rotate` method for
        `Fillable` subclasses when an axis of rotation is specified.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Build the rotation axis in a given position
        axis_centre = (1, 1, 0)
        axis_centre_vrtx = make_vertex(axis_centre)
        axis = make_vector_from_points(
            axis_centre_vrtx, make_vertex((*axis_centre[:2], 1))
        )
        # Build reference vector before the rotation (from fillable centre to
        # axis centre)
        ref_vect = make_vector_from_points(self.fillable.o, axis_centre_vrtx)

        # Rotate the fillable around the specified axis
        rot_angle = 90.0
        self.fillable.rotate(rot_angle, axis)
        # Build reference vector after the rotation
        ref_vect2 = make_vector_from_points(self.fillable.o, axis_centre_vrtx)
        # Check the correct rotation happened
        self.assertEqual(self.fillable.rot_angle, rot_angle)
        self.assertTrue(
            math.isclose(
                get_angle_between_shapes(ref_vect, ref_vect2),
                rot_angle,
                abs_tol=1e-6
            )
        )

    def test_scale(self) -> None:
        """
        Method that tests the scaling of `Fillable` subclasses instances by
        checking the area of the characteristic surface before and after the
        operation.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        scale_factor = 2
        area_before = get_basic_properties(self.fillable)[1]
        # Scale the fillable wrt to its centre
        self.fillable.scale(scale_factor)
        # Check the area has increased by the scaling factor^2
        self.assertTrue(
            math.isclose(
                area_before*scale_factor*scale_factor,
                get_basic_properties(self.fillable)[1],
                abs_tol=1e-6
            )
        )

        # Verify the exception is raised if an invalid scaling factor is
        # provided
        with self.assertRaises(ValueError):
            self.fillable.scale(0.0)
        with self.assertRaises(ValueError):
            self.fillable.scale(-2.0)

    def test_set_region_properties_exceptions(self) -> None:
        """
        Method that tests the `set_region_properties` method for `Fillable`
        subclasses. It verifies that:

        - the exception is raised when no shape is provided and nothing
          is selected in the study;
        - the exception is raised when the given shape does not have a
          corresponding region in the fillable's hierarchical tree;
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Verify the exception is raised when no shape is provided and nothing
        # is selected in the study
        with self.assertRaises(RuntimeError):
            self.fillable.set_region_properties({})

        # Create an invalid shape that is not in the layout
        invalid_shape = Circle((5.0, 5.0, 0.0))
        # Verify the exception is raised for the invalid shape
        with self.assertRaises(RuntimeError):
            self.fillable.set_region_properties({}, invalid_shape)

    def test_set_region_properties(self) -> None:
        """
        Method that tests the `set_region_properties` method for `Fillable`
        subclasses. It verifies that:

        - the properties of the region of the fillable are correctly updated;
        - in case of a region of the fillable without properties, its
          properties are correctly set;
        - properties are not modified if an empty dictionary is given.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Declare a reference region and add it to the geometry layout
        region = Region(
            Circle(),
            "Region_1",
            {PropertyType.MACRO: "MAC001", PropertyType.MATERIAL: "MAT1"}
        )
        self.fillable.add(region)
        # Update the hierarchical structure to populate regions
        self.fillable.update_hierarchical_structure()
        # Verify that only the indicated property is updated
        updated_prop = "UPDATED_MAT"
        self.fillable.set_region_properties(
            {PropertyType.MATERIAL: updated_prop}, region.geom_obj
        )
        self.assertEqual(
            self.fillable.layers[-1][-1].properties[PropertyType.MATERIAL],
            updated_prop
        )

        # Clear the properties of the tested region of the fillable
        self.fillable.layers[-1][-1].properties.clear()
        # Verify that the properties of the fillable's region are set
        updated_mac = "MAC002"
        self.fillable.set_region_properties(
            {
                PropertyType.MATERIAL: updated_prop,
                PropertyType.MACRO: updated_mac
            },
            region.geom_obj
        )
        self.assertEqual(
            self.fillable.layers[-1][-1].properties[PropertyType.MATERIAL],
            updated_prop
        )
        self.assertEqual(
            self.fillable.layers[-1][-1].properties[PropertyType.MACRO],
            updated_mac
        )

        # Verify properties are not modified if an empty dictionary is given
        self.fillable.set_region_properties({}, region.geom_obj)
        self.assertNotEqual(self.fillable.layers[-1][-1].properties, {})
        self.assertEqual(
            self.fillable.layers[-1][-1].properties[PropertyType.MATERIAL],
            updated_prop
        )
        self.assertEqual(
            self.fillable.layers[-1][-1].properties[PropertyType.MACRO],
            updated_mac
        )

    def test_show_exceptions(self) -> None:
        """
        Method that tests the `show` method for `Fillable` subclasses. It
        verifies that the exception is raised:

        - in case of an invalid argument type;
        - if any region does not have a value for the indicated property;
        - if no compound of edges is associated to the indicated geometry
          type.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()

        # Verify the invalid argument case
        with self.assertRaises(RuntimeError):
            self.fillable.show(SymmetryType.QUARTER)

        # Verify the case for a missing property value for a region
        self.fillable.add(Region(Circle()))
        with self.assertRaises(RuntimeError):
            self.fillable.show(PropertyType.MATERIAL)

        # Verify the case for a missing compound of edges for the given
        # geometry type
        with self.assertRaises(RuntimeError):
            self.fillable.show(GeometryType.SECTORIZED)

    def test_show_no_settings(self) -> None:
        """
        Method that tests the `show` method for `Fillable` subclasses without
        any visualisation settings. It verifies that:

        - the layout is correctly updated before display;
        - the GEOM object is added to the study;
        - regions are added to the study.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Add a region to the fillable
        region = Region(
            Circle((0.5, 0.5, 0.0), radius=0.01),
            "TestRegion",
            {PropertyType.MATERIAL: "MAT1"}
        )
        self.fillable.add(region)
        # Store the fillable's GEOM object before updating the structure
        pre_geom_obj = self.fillable.geom_obj

        # Call show with default parameters
        self.fillable.show()

        # Verify the fillable's GEOM object has changed
        self.assertIsNotNone(self.fillable.geom_obj)
        self.assertFalse(
            are_same_shapes(
                self.fillable.geom_obj, pre_geom_obj, ShapeType.COMPOUND
            )
        )
        # Verify the entry_id has been set and the corresponding GEOM object
        # can be retrieved
        self.assertIsNotNone(self.fillable.entry_id)
        self.assertTrue(
            are_same_shapes(
                get_object_from_id(self.fillable.entry_id),
                self.fillable.geom_obj,
                ShapeType.COMPOUND
            )
        )
        # Verify regions have been added to the study as children of the
        # fillable
        for region in self.fillable.regions:
            self.assertIsNotNone(region.entry_id)
            self.assertIn(self.fillable.entry_id, region.entry_id)
            self.assertIsNotNone(get_object_from_id(region.entry_id))

        # Verify the state indicates no update is needed
        self.assertFalse(self.fillable.state.is_update_needed)
        # Verify the displayed geometry type is TECHNOLOGICAL
        self.assertEqual(
            self.fillable.state.displayed_geom, GeometryType.TECHNOLOGICAL
        )

    def test_show_with_property_type(self) -> None:
        """
        Method that tests the `show` method with a `PropertyType` argument
        for enabling the regions' colour mapping.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Store the names of the materials prior to the update
        unique_mats = set()
        for region in self.fillable.get_regions():
            unique_mats.add(region.properties[PropertyType.MATERIAL])

        # Add regions with different materials
        region1 = Region(
            Circle(), "Region1", {PropertyType.MATERIAL: "MAT_REG_1"}
        )
        region2 = Region(
            Circle(radius=2), "Region2", {PropertyType.MATERIAL: "MAT_REG_2"}
        )
        self.fillable.add(region2)
        self.fillable.add(region1)

        # Call show with 'PropertyType.MATERIAL'
        self.fillable.show(PropertyType.MATERIAL)

        # Verify the layout was updated
        self.assertFalse(self.fillable.state.is_update_needed)
        # Verify regions have been added to the study as child of the
        # fillable, and each with a colour different from the default one
        unique_colors = set()
        for region in self.fillable.regions:
            self.assertIsNotNone(region.entry_id)
            self.assertIsNotNone(get_object_from_id(region.entry_id))
            self.assertIn(self.fillable.entry_id, region.entry_id)
            self.assertNotEqual(region.color, DEFAULT_REGION_COLOR)
            unique_colors.add(region.color)
        # Verify the number of colors has increased by 2
        self.assertTrue(len(unique_colors) - len(unique_mats), 2)

    def test_show_with_geometry_type(self) -> None:
        """
        Method that tests the `show` method with a `GeometryType` argument
        for enabling the visualisation of the corresponding edges.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Build a cell, sectorise it and add it to the fillable
        cell = CartesianCell()
        cell.sectorize([8], [0])
        self.fillable.add(cell)

        # Call show with 'SECTORIZED' geometry type
        self.fillable.show(GeometryType.SECTORIZED)

        # Verify the compound of edges has been added to the study as child
        # the fillable
        sect_cmpd_id = get_id_from_object(
            self.fillable.geometry_maps[GeometryType.SECTORIZED]
        )
        self.assertIsNotNone(sect_cmpd_id)
        self.assertIn(self.fillable.entry_id, sect_cmpd_id)
        # Verify the displayed geometry type is set correctly
        self.assertEqual(
            self.fillable.state.displayed_geom, GeometryType.SECTORIZED
        )

    def test_show_with_symmetry(self) -> None:
        """
        Method that tests the `show` method after applying a symmetry.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Apply 'QUARTER' symmetry
        symmetry_type = SymmetryType.QUARTER
        self.fillable.apply_symmetry(symmetry_type)

        # Call show method
        self.fillable.show()

        # Verify regions were collected for the given symmetry
        common = make_common(
            self.fillable, self.fillable.symmetry_map[symmetry_type]
        )
        symm_faces = extract_sub_shapes(
            make_compound([common]), ShapeType.FACE
        )
        self.assertEqual(len(symm_faces), len(self.fillable.regions))
        for f in symm_faces:
            found = False
            for r in self.fillable.regions:
                if are_same_shapes(f, r, ShapeType.FACE):
                    found = True
                    break
            self.assertTrue(found)

        # Verify the state symmetry type is set
        self.assertEqual(
            self.fillable.state.symmetry_type, SymmetryType.QUARTER
        )

    def test_show_clears_previous_view(self) -> None:
        """
        Method that tests that the `show` method clears the previous view
        before displaying new content. This is tested by verifying that two
        successive calls of the `show` method produces two values for the
        `entry_id` attribute.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()

        # Call show multiple times
        self.fillable.show()
        first_entry_id = self.fillable.entry_id
        self.fillable.show()
        second_entry_id = self.fillable.entry_id

        # Both entry IDs should be set and different and no GEOM object is
        # associated to the first entry ID
        self.assertIsNotNone(first_entry_id)
        self.assertIsNotNone(second_entry_id)
        self.assertNotEqual(first_entry_id, second_entry_id)
        self.assertIsNone(get_object_from_id(first_entry_id))
        self.assertIsNotNone(get_object_from_id(second_entry_id))

    def test_translate(self) -> None:
        """
        Method that tests the `translate` method for `Fillable` subclasses.
        It verifies the correctness of the resulting position of its
        geometric elements.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Store the elements for assessing the correct translation
        center_after_transl = (1, 1, 0)
        new_pos_vrtx = make_vertex(center_after_transl)
        cdg_pre = make_cdg(self.fillable)
        distance = get_min_distance(self.fillable.o, new_pos_vrtx)
        layouts_distances_1 = [
            [get_min_distance(layout.o, self.fillable.o) for layout in layer]
            for layer in self.fillable.layers
        ]

        # Translate the surface
        self.fillable.translate(center_after_transl)
        # Check the correct translation happened
        self.assertTrue(
            are_same_shapes(self.fillable.o, new_pos_vrtx, ShapeType.VERTEX)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre, make_cdg(self.fillable)),
                distance,
                abs_tol=1e-6
            )
        )
        # Verify the layouts in the fillable's layers keep the same relative
        # distance wrt the fillable's centre
        layouts_distances_2 = [
            [get_min_distance(layout.o, self.fillable.o) for layout in layer]
            for layer in self.fillable.layers
        ]
        for layer_1, layer_2 in zip(layouts_distances_1, layouts_distances_2):
            for layout_dist_1, layout_dist_2 in zip(layer_1, layer_2):
                self.assertAlmostEqual(layout_dist_1, layout_dist_2)

        # Check no translation happened if the new centre coincides with the
        # current one
        self.fillable.translate(center_after_transl)
        self.assertTrue(
            are_same_shapes(self.fillable.o, new_pos_vrtx, ShapeType.VERTEX)
        )

    def test_update(self) -> None:
        """
        Method that tests the implementation of the `update` method for
        `Fillable` subclasses.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Update the fillable with another one
        new_layout = Circle()
        self.fillable.update(new_layout)
        # Verify that the geometric elements has changed correctly
        self.assertTrue(
            are_same_shapes(
                self.fillable, make_compound([new_layout]), ShapeType.COMPOUND
            )
        )
        self.assertTrue(
            are_same_shapes(self.fillable.o, new_layout.o, ShapeType.VERTEX)
        )

    def test_update_hierarchical_structure_no_update_needed(self) -> None:
        """
        Method that tests that `update_hierarchical_structure` method returns
        immediately without making changes when `is_update_needed` is False.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Set the state to not need updating
        self.fillable.state.is_update_needed = False
        initial_geom_obj = self.fillable.geom_obj

        # Call update_hierarchical_structure
        self.fillable.update_hierarchical_structure()

        # Verify nothing changed
        self.assertIs(self.fillable.geom_obj, initial_geom_obj)

    def test_update_hierarchical_structure_single_layer(self) -> None:
        """
        Method that tests the `update_hierarchical_structure` method for
        `Fillable` subclasses with a single layer containing only regions.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Restore the fillable
        self.fillable.restore()
        # Add regions to the fillable
        layout1 = Region(
            Circle((0.5, 0.5, 0.0), radius=0.1),
            properties={PropertyType.MATERIAL: "MAT1"}
        )
        layout2 = Region(
            Circle(radius=0.5), properties={PropertyType.MATERIAL: "MAT2"}
        )
        self.fillable.add(layout1, layer_index=0)
        self.fillable.add(layout2, layer_index=0)
        # Ensure update is needed
        self.fillable.state.is_update_needed = True

        # Update the hierarchical structure
        self.fillable.update_hierarchical_structure()

        # Verify the state is no longer marked as needing update
        self.assertFalse(self.fillable.state.is_update_needed)
        # Verify the GEOM object has been updated
        self.assertIsNotNone(self.fillable.geom_obj)
        # Verify the first layer has not been cut by the two regions
        self.assertFalse(
            are_same_shapes(
                self.fillable.layers[0][0],
                (self.fillable.shape - layout1) - layout2,
                ShapeType.FACE
            )
        )
        # Verify regions are collected correctly
        self.assertEqual(
            len(self.fillable.get_regions()),
            len(extract_sub_shapes(self.fillable.geom_obj, ShapeType.FACE))
        )

    def test_update_hierarchical_structure_multiple_layers(self) -> None:
        """
        Method that tests the `update_hierarchical_structure` method for
        `Fillable` subclasses with multiple layers that need to be collapsed.
        Layers only contain regions.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create regions with different sizes to ensure overlapping
        layout1 = Region(
            Circle(radius=0.5), properties={PropertyType.MATERIAL: "MAT1"}
        )
        layout2 = Region(
            Circle(radius=0.25), properties={PropertyType.MATERIAL: "MAT2"}
        )
        # Add to different layers
        self.fillable.add(layout1)
        self.fillable.add(layout2)
        # Ensure update is needed
        self.fillable.state.is_update_needed = True

        # Update the hierarchical structure
        pre_geom_obj = self.fillable.geom_obj
        self.fillable.update_hierarchical_structure()

        # Verify the state is no longer marked as needing update
        self.assertFalse(self.fillable.state.is_update_needed)
        # Verify the layers have been cut
        self.assertTrue(
            are_same_shapes(
                make_compound(self.fillable.layers[-3]),
                make_compound([make_cut(pre_geom_obj, layout1)]),
                ShapeType.COMPOUND
            )
        )
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[-2][0],
                layout1 - layout2,
                ShapeType.FACE
            )
        )

        # Verify the GEOM object has been updated
        self.assertIsNotNone(self.fillable.geom_obj)

    def test_update_hierarchical_structure_with_nested_fillable(self) -> None:
        """
        Method that tests the `update_hierarchical_structure` method for
        `Fillable` subclasses with multiple layers that need to be collapsed.
        Nested `Fillable` objects are contained in the layers.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create a nested fillable (CartesianCell)
        nested_fillable = CartesianCell(width_height=(0.4, 0.4))
        region = Region(
            Circle(radius=0.05), properties={PropertyType.MATERIAL: "MAT1"}
        )
        nested_fillable.add(region)
        # Add the nested fillable to this fillable
        self.fillable.add(nested_fillable)
        # Ensure update is needed
        self.fillable.state.is_update_needed = True

        # Update the hierarchical structure
        pre_geom_obj = self.fillable.geom_obj
        self.fillable.update_hierarchical_structure()

        # Verify the state is no longer marked as needing update
        self.assertFalse(self.fillable.state.is_update_needed)
        # Verify nested fillable was also updated
        self.assertFalse(nested_fillable.state.is_update_needed)
        # Verify the GEOM object has been updated
        self.assertIsNotNone(self.fillable.geom_obj)
        # Verify the layers have been cut
        self.assertTrue(
            are_same_shapes(
                make_compound(self.fillable.layers[-2]),
                make_compound([make_cut(pre_geom_obj, nested_fillable)]),
                ShapeType.COMPOUND
            )
        )
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[-1][0].layers[0][0],
                nested_fillable.layers[0][0] - region,
                ShapeType.FACE
            )
        )
        # Verify regions are collected from nested structure
        self.assertEqual(
            len(self.fillable.get_regions()),
            len(extract_sub_shapes(self.fillable.geom_obj, ShapeType.FACE))
        )

    def test_update_hierarchical_structure_collapse_layers(self) -> None:
        """
        Method that tests the `update_hierarchical_structure` method for
        `Fillable` subclasses when the `collapse_layers` parameter is set to
        `True`, which reduces the fillable's hierarchical tree to a single
        layer.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Add multiple regions to different layers
        layout1 = Region(
            Circle(radius=0.5), properties={PropertyType.MATERIAL: "MAT1"}
        )
        layout2 = Region(
            Circle(radius=0.4), properties={PropertyType.MATERIAL: "MAT2"}
        )
        layout3 = Region(
            Circle(radius=0.3), properties={PropertyType.MATERIAL: "MAT3"}
        )
        self.fillable.add(layout1)
        self.fillable.add(layout2)
        self.fillable.add(layout3)
        # Ensure update is needed
        self.fillable.state.is_update_needed = True

        # Update the hierarchical structure with collapse_layers=True
        pre_geom_obj = self.fillable.geom_obj
        self.fillable.update_hierarchical_structure(collapse_layers=True)

        # Verify the layers have been collapsed to a single layer
        self.assertEqual(len(self.fillable.layers), 1)
        # Verify the state is no longer marked as needing update
        self.assertFalse(self.fillable.state.is_update_needed)
        # Verify the layers have been cut
        self.assertTrue(
            are_same_shapes(
                make_compound(self.fillable.layers[0][0]),
                make_compound([make_cut(pre_geom_obj, layout1)]),
                ShapeType.COMPOUND
            )
        )
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[0][1],
                (layout1 - layout2) - layout3,
                ShapeType.FACE
            )
        )
        self.assertTrue(
            are_same_shapes(
                self.fillable.layers[0][2], layout2 - layout3, ShapeType.FACE
            )
        )

    def test_update_hierarchical_structure_geometry_maps(self) -> None:
        """
        Method that tests that the `update_hierarchical_structure` method with
        `collapse_layers=True` preserves geometry maps from nested `Fillable`
        objects.
        """
        # Skip the test if run from this class
        self.__skip_if_superclass()
        # Create a nested fillable with a sectorized geometry map
        region = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})
        nested_cell = CartesianCell(width_height=(0.4, 0.4))
        nested_cell.add(region)
        nested_cell.sectorize([8], [0])
        # Add the nested fillable to this fillable
        self.fillable.add(nested_cell)
        # Ensure update is needed
        self.fillable.state.is_update_needed = True

        # Update the hierarchical structure with collapse_layers=True
        self.fillable.update_hierarchical_structure(collapse_layers=True)

        # Verify geometry maps have been preserved
        self.assertIn(GeometryType.SECTORIZED, self.fillable.geometry_maps)
        self.assertTrue(
            are_same_shapes(
                self.fillable.geometry_maps[GeometryType.SECTORIZED],
                nested_cell.geometry_maps[GeometryType.SECTORIZED],
                ShapeType.COMPOUND
            )
        )
        # Verify the state is no longer marked as needing update
        self.assertFalse(self.fillable.state.is_update_needed)

    def __assess_symmetry(self, symmetry_type: SymmetryType) -> None:
        """
        Method that checks the proper application of the given symmetry type.

        Parameters
        ----------
        symmetry_type : SymmetryType
            The type of symmetry to apply and verify.
        """
        # Apply the symmetry
        self.fillable.apply_symmetry(symmetry_type)
        # Verify the symmetry was applied
        self.assertIn(symmetry_type, self.fillable.symmetry_map)
        self.assertEqual(self.fillable.state.symmetry_type, symmetry_type)
        # Verify the symmetry shape is not empty and represents a face object
        symm_shape = self.fillable.symmetry_map[symmetry_type]
        self.assertIsNotNone(symm_shape.geom_obj)
        self.assertTrue(get_shape_type(symm_shape) == ShapeType.FACE)

    def __skip_if_superclass(self) -> None:
        """
        Method that checks whether the current test class is `TestSurface`
        and skips the test that runs this method if the class is
        `TestSurface`.
        """
        if self.__class__ is TestFillable:
            self.skipTest(
                f"{self.__class__.__name__}. This test is only valid in "
                "subclasses of 'TestFillable'."
            )


if __name__ == "__main__":
    unittest.main()
