"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.cells` module have a valid implementation.
"""
from abc import ABC
from copy import deepcopy
import math
from typing import Any, Dict, Tuple, Union
import unittest

from glow.generator.support import CellType, GeometryType, PropertyType
from glow.geometry_layouts.cells import Cell, HexCell, RectCell, Region, \
    get_region_info
from glow.geometry_layouts.geometries import Hexagon, Rectangle, Surface
from glow.geometry_layouts.utility import are_same_shapes
from glow.interface.geom_interface import *
from support_funcs import build_cell_ref_vectors


class TestRegion(unittest.TestCase):
    """
    Test case for verifying the correct storage of information about a
    generic surface being a region in a cell/lattice.

    Attributes
    ----------
    face : Any
        The face object the `Region` instance refers to.
    cdg : Union[Any, None]
        The vertex object representing the CDG of the region.
    face_entry_id : Union[str, None]
        The entry ID of the region's face when displayed in the SALOME viewer.
    inner_point : Any
        A vertex object representing a point within the region's boundaries.
    name : str
        The name associated to the region when displayed in the SALOME viewer.
    properties : Dict[PropertyType, str]
        A dictionary associating a value to each of the `PropertyType`
        attributed to the region.
    color : Tuple[int, int, int]
        A tuple providing the RGB code to colour the region's face when
        displayed in the SALOME viewer.
    region : Region
        An instance of the `Region` dataclass to test.
    """
    def setUp(self):
        """
        Method that sets up the test environment for a `Region` dataclass by
        initializing the geometric characteristics to be assigned to the
        cell's/lattice's region.
        """
        # Setup the geometric characteristics used for testing a 'Region'
        # instance
        self.face: Any = make_face(
            [make_circle(make_vertex((0.0, 0.0, 0.0)), None, 1.0)]
        )
        self.cdg: Any = make_cdg(self.face)
        self.face_entry_id: Union[str, None] = None
        self.inner_point: Any = make_vertex_inside_face(self.face)
        self.name: str = "Region"
        self.properties: Dict[PropertyType, str] = {
            PropertyType.MATERIAL, "FUEL"
        }
        self.color: Tuple[int, int, int] = (167, 167, 167)
        # Instantiate the 'Region' object to test
        self.region = Region(
            face=self.face,
            face_entry_id=self.face_entry_id,
            inner_point=self.inner_point,
            name=self.name,
            properties=self.properties)

    def test_init(self) -> None:
        """
        Method that tests the initialization of a `Region` object by
        verifying it is correctly instantiated with the provided data.
        """
        # Check the correct instantiation
        self.assertTrue(
            are_same_shapes(self.region.face, self.face, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(self.region.cdg, self.cdg, ShapeType.VERTEX)
        )
        self.assertEqual(self.region.face_entry_id, self.face_entry_id)
        self.assertTrue(
            are_same_shapes(
                self.region.inner_point, self.inner_point, ShapeType.VERTEX)
        )
        self.assertEqual(self.region.name, self.name)
        self.assertEqual(self.region.properties, self.properties)
        self.assertEqual(self.region.color, self.color)

    def test_calculate_cdg(self) -> None:
        """
        Method that tests the implementation of the private method
        `__calculate_cdg` of the `Region` class.
        """
        # Check the returned object is the region's CDG
        cdg = self.region._Region__calculate_cdg()
        self.assertTrue(
            are_same_shapes(cdg, self.cdg, ShapeType.VERTEX)
        )
        self.assertTrue(
            are_same_shapes(cdg, make_cdg(self.region.face), ShapeType.VERTEX)
        )
        # Check the method returns None, if no face is defined
        region = Region.__new__(Region)
        self.assertEqual(region._Region__calculate_cdg(), None)

    def test_reset_color_to_default(self) -> None:
        """
        Method that tests the implementation of the `reset_color_to_default`
        method of the `Region` class.
        """
        # Set the colour associated to the region to a value different from
        # the default one
        self.region.color = (255, 0, 0)
        # Reset the colour and assert it returned to the default value
        self.region.reset_color_to_default()
        self.assertEqual(self.region.color, self.color)

    def test_set_property_color(self) -> None:
        """
        Method that tests the implementation of the `set_property_color`
        method of the `Region` class.
        """
        # Change the color associated to the region
        new_color = (255, 0, 0)
        self.region.set_property_color(new_color)
        # Assess the change happened successfully
        self.assertEqual(self.region.color, new_color)


class TestCell(ABC, unittest.TestCase):
    """
    Abstract base test case for verifying the geometric operations and
    visualization capabilities of the `Cell` subclasses.

    This test suite provides common setup and a set of tests to ensure that
    implementations of `Cell` correctly handle the operations performed on a
    cell, whose type is either cartesian or hexagonal.
    Tests methods declared herein are skipped as the `Cell` class is abstract.
    They are run in subclasses of `TestCell` which tests the subclasses of
    `Cell`.

    Attributes
    ----------
    name : str
        The name of the cell when displayed in the SALOME viewer.
    cell : Cell | None
        The `Cell` superclass to test. It is initialized to `None` to skip
        tests acting on its subclasses.
    o : Any
        A vertex object positioned at the origin of the XYZ space.
    surf : Surface
        The `Surface` object providing the characteristic shape for a cell.
    sect_opts : Tuple[List[int], List[float]]
        Providing the sectorization options as the number of sectors and the
        starting angle for each cell-centered region.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the subclasses of the
        `Cell` class.
        It initializes the common geometric characteristics for a cell
        while setting to `None` the ones that are specific to a cartesian
        or a hexagonal cell.
        """
        self.name: str = "Cell"
        self.surf: Surface | None = None
        self.cell: Cell | None = None
        self.o = make_vertex((0.0, 0.0, 0.0))
        self.sect_opts = ([], [])

    def test_initialize_geometry(self) -> None:
        """
        Method that tests the implementation of the private method
        `__initialize_geometry` of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Initialize the cell's geometry layout
        self.cell._Cell__initialize_geometry(self.surf.face)
        # Verify the correct geometry layout assignement
        self.__assess_cell_face_initialization()

    def test_initialize_region_dicts(self) -> None:
        """
        Method that tests the implementation of the private method
        `__initialize_region_dicts` of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Initialize the cell's dictionaries
        self.cell._Cell__initialize_region_dicts()
        # Verify the correct initialization of the dictionaries
        self.__assess_cell_dict_initialization()

    def test_extract_subfaces_from_face(self) -> None:
        """
        Method that tests the implementation of the private method
        `__extract_subfaces_from_face` of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Extract the subfaces from a generic face object
        face = make_partition(
            [make_face(make_circle(self.o, None, 2.0))],
            [make_face(make_circle(self.o, None, 1.0))],
            ShapeType.FACE)
        subfaces = self.cell._Cell__extract_subfaces_from_face(face)
        # Verify the correct extraction
        self.assertEqual(len(subfaces), 2)
        for i in range(len(subfaces)):
            self.assertTrue(get_shape_type(subfaces[i]) == ShapeType.FACE)
        for i in range(len(subfaces)-1):
            self.assertTrue(
                get_min_distance(subfaces[i], self.o) <
                get_min_distance(subfaces[i+1], self.o)
            )

    def test_extract_subfaces(self) -> None:
        """
        Method that tests the implementation of the method `extract_subfaces`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Extract the cell's subfaces
        subfaces = self.cell.extract_subfaces()
        # Verify the correct extraction
        self.assertTrue(len(subfaces) >= 1)
        self.assertTrue(
            are_same_shapes(self.cell.face,
                            make_partition(subfaces, [], ShapeType.FACE),
                            ShapeType.FACE)
        )

    def test_get_centered_circles(self) -> None:
        """
        Method that tests the implementation of the method
        `get_centered_circles` of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Configure the cell
        self.cell.add_circle(0.1)
        self.cell.add_circle(0.1, (0.1, 0.0, 0.0))
        self.cell.add_circle(0.1, (-0.1, 0.0, 0.0))
        self.cell.add_circle(0.1, (0.0, 0.1, 0.0))
        self.cell.add_circle(0.1, (0.0, -0.1, 0.0))
        self.cell.add_circle(0.2)
        self.cell.add_circle(0.3)

        # Retrieve the cell-centered circles only
        centered_circles = self.cell.get_centered_circles()
        # Verify only the cell-centered 'Circle' objects have been retrieved
        self.assertEqual(len(centered_circles), 3)
        for cc in centered_circles:
            self.assertTrue(cc.radius in [0.1, 0.2, 0.3])
            self.assertTrue(
                get_min_distance(cc.o, self.cell.figure.o) < 1e-5
            )

    def test_get_regions_info(self) -> None:
        """
        Method that tests the implementation of the method `get_regions_info`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Verify an exception is raised when calling the method without
        # selecting any region
        with self.assertRaises(RuntimeError):
            self.cell.get_regions_info()
        # Test the inner implementation by verifying an exception is raised
        # when providing a region not beloging to the cell
        with self.assertRaises(RuntimeError):
            self.cell.add_circle(0.1)
            self.cell._Cell__build_regions()
            get_region_info(self.cell.figure.face, self.cell.regions)
        region_face = make_face([make_circle(self.o, None, 0.1)])
        # Check the message returns the name of the region and the associated
        # properties, if any are assigned
        self.__assess_get_regions_info(region_face,
                                       "No associated properties.")
        # Test the result message when a property is set
        self.cell.set_region_property(
            PropertyType.MATERIAL, "MAT", region_face)
        self.cell._Cell__build_regions()
        # Check the correct assignment
        self.__assess_get_regions_info(
            region_face, f"{PropertyType.MATERIAL.name}: MAT")

    def test_initialize_cell(self) -> None:
        """
        Method that tests the implementation of the private method
        `__initialize_cell` of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Initialize the cell's dictionaries
        self.cell._Cell__initialize_cell()
        # Verify the correct initialization
        self.__assess_cell_face_initialization()
        self.__assess_cell_dict_initialization()

    @unittest.skip
    def test_initialize_specific_cell(self) -> None:
        """
        Method for testing the initialization of the instance attributes
        that are specific for the subclasses of the `Cell` class.
        The implementation is specific to the subclasses of `TestCell`.
        """

    @unittest.skip
    def test_check_radius_vs_cell_dim(self) -> None:
        """
        Method for testing if the radius of a circle exceeds the cell's
        dimensions.
        The implementation is specific to the subclasses of `TestCell`.
        """

    def test_add_circle(self) -> None:
        """
        Method that tests the implementation of the method `add_circle`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Add a circle to the cell's geometric layout
        self.__assess_add_circle(None, 0.1, len(self.cell.inner_circles))
        # Add a circle to the cell's geometric layout in a different position
        self.__assess_add_circle(
            (0.1, 0.0, 0.0), 0.1, len(self.cell.inner_circles))
        # Test exceptions are raised when adding a circle that exceeds the
        # cell's dimensions or one that has already been added
        with self.assertRaises(RuntimeError):
            self.cell.add_circle(1.0)
        with self.assertRaises(RuntimeError):
            self.cell.add_circle(0.1)

        # Test the circle addition when a sectorization is present
        self.cell.sectorize(self.sect_opts[0][0:2], self.sect_opts[1][0:2])
        self.__assess_add_circle(
            None, 0.2, len(self.cell.inner_circles))

    def test_remove_circle(self) -> None:
        """
        Method that tests the implementation of the method `remove_circle`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Circle's data
        radius = 0.1
        position = (0.1, 0.0, 0.0)
        # Test exceptions are raised when using the method incorrectly, i.e.
        # when the cell does not have any circle or when trying to remove a
        # non existent one
        with self.assertRaises(RuntimeError):
            self.cell.remove_circle(0.2)
        self.cell.add_circle(radius)
        self.cell.add_circle(radius, position)
        with self.assertRaises(RuntimeError):
            self.cell.remove_circle(0.2)

        # Remove an existent circle in the cell's center
        self.__assess_remove_circle(
            None, radius, len(self.cell.inner_circles))
        # Remove an existent circle not in the cell's center
        self.__assess_remove_circle(
            position, radius, len(self.cell.inner_circles))
        # Verify the cell's face returned to the one of its figure
        self.assertTrue(
            are_same_shapes(
                self.cell.face, self.cell.figure.face, ShapeType.FACE)
        )

        # Test the circle removal when a sectorization is present
        self.cell.add_circle(radius)
        self.cell.add_circle(radius, position)
        self.cell.sectorize(self.sect_opts[0][0:2], self.sect_opts[1][0:2])
        self.__assess_remove_circle(
            None, radius, len(self.cell.inner_circles))
        self.__assess_remove_circle(
            position, radius, len(self.cell.inner_circles))
        # Verify the cell's face returned to the one of its figure
        self.assertTrue(
            are_same_shapes(
                self.cell.face, self.cell.figure.face, ShapeType.FACE)
        )

    def test_restore(self) -> None:
        """
        Method that tests the implementation of the method `restore` of
        the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Add a circle
        self.cell.add_circle(0.2)
        # Restore the cell
        self.cell.restore()
        # Assess the correct restore happened
        self.__assess_cell_face_initialization()
        self.__assess_cell_dict_initialization()

    def test_rotate(self) -> None:
        """
        Method that tests the implementation of the method `rotate` of
        the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Store the initial rotation angle
        rotation_0 = self.cell.rotation
        rot_angle = 90.0
        # Build reference vectors for each cell's face object before and
        # after the rotation
        ref_vectors1 = build_cell_ref_vectors(self.cell)
        self.cell.rotate(rot_angle)
        ref_vectors2 = build_cell_ref_vectors(self.cell)
        # Check the correct rotation happened by assessing that the angle
        # between the reference vectors is the rotation angle
        for v1, v2 in zip(ref_vectors1, ref_vectors2):
            self.assertTrue(
                math.isclose(
                    get_angle_between_shapes(v1, v2),
                    rot_angle)
            )
        # Check the 'rotation' attribute has updated
        self.assertTrue(
            math.isclose(
                self.cell.rotation, rotation_0 + math.radians(rot_angle))
        )

    def test_rotate_from_axis(self) -> None:
        """
        Method that tests the implementation of the method `rotate_from_axis`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Store the initial rotation angle and the rotation axis
        rotation_0 = self.cell.rotation
        rot_angle = 90.0
        center = get_point_coordinates(self.cell.figure.o)
        rot_axis = make_vector_from_points(
            self.o, make_vertex((center[0], center[1], 1)))
        # Build reference vectors for each cell's face object before and
        # after the rotation
        ref_vectors1 = build_cell_ref_vectors(self.cell)
        self.cell.rotate_from_axis(rot_angle, rot_axis)
        ref_vectors2 = build_cell_ref_vectors(self.cell)
        # Check the correct rotation happened by assessing that the angle
        # between the reference vectors is the rotation angle
        for v1, v2 in zip(ref_vectors1, ref_vectors2):
            self.assertTrue(
                math.isclose(
                    get_angle_between_shapes(v1, v2),
                    rot_angle)
            )
        # Check the 'rotation' attribute has updated
        self.assertTrue(
            math.isclose(
                self.cell.rotation, rotation_0 + math.radians(rot_angle))
        )

    def test_set_properties(self) -> None:
        """
        Method that tests the implementation of the method `set_properties`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Add a circle
        self.cell.add_circle(0.2)
        # Check the regions of the cell's technological geometry do not have
        # any property associated
        for region in self.cell.tech_geom_props:
            self.assertEqual(self.cell.tech_geom_props[region], {})

        # Assign the MATERIAL property to the cell's regions
        materials = {PropertyType.MATERIAL: ['MAT1', 'MAT2']}
        self.cell.set_properties(materials)

        # Verify the correct assignment of the materials to the regions
        # according to their distance from the cell's center
        regions = sorted(
            extract_sub_shapes(self.cell.face, ShapeType.FACE),
            key=lambda subface: get_min_distance(self.cell.figure.o, subface)
        )
        for i, (r1, r2) in enumerate(zip(regions, self.cell.tech_geom_props)):
            if are_same_shapes(r1, r2, ShapeType.FACE):
                self.assertEqual(
                    self.cell.tech_geom_props[r2][PropertyType.MATERIAL],
                    materials[PropertyType.MATERIAL][i]
                )

    def test_set_region_property(self) -> None:
        """
        Method that tests the implementation of the method
        `set_region_property` of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        property_value = (PropertyType.MATERIAL, "MAT")
        # Verify an exception is raised when calling the method without
        # selecting any region or providing a non-existent one
        with self.assertRaises(RuntimeError):
            self.cell.set_region_property(*property_value)
        with self.assertRaises(RuntimeError):
            self.cell.add_circle(0.1)
            self.cell.set_region_property(*property_value,
                                          self.cell.figure.face)

        # Call the method providing the region to update
        region_face = make_face([make_circle(self.o, None, 0.1)])
        self.cell.set_region_property(*property_value, region_face)
        # Check the correct assignment
        for region in self.cell.tech_geom_props:
            if are_same_shapes(region, region_face, ShapeType.FACE):
                self.assertTrue(
                    property_value[0] in self.cell.tech_geom_props[region]
                )
                self.assertEqual(
                    self.cell.tech_geom_props[region][property_value[0]],
                    property_value[1]
                )
                break

    def test_show(self) -> None:
        """
        Method that tests the implementation of the method `show` of the
        `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Add circles
        self.cell.add_circle(0.1)
        self.cell.add_circle(0.2)
        # Apply the sectorization
        self.cell.sectorize(*self.sect_opts)
        # Check the cell's face has not been displayed yet
        self.assertIsNone(self.cell.face_entry_id)

        # Verify cell's regions are correctly shown
        self.__assess_show(None, GeometryType.TECHNOLOGICAL)
        self.__assess_show(None, GeometryType.SECTORIZED)
        # Verify an exception is raised when trying to show regions according
        # to a property type without assigning values to regions
        with self.assertRaises(RuntimeError):
            self.__assess_show(PropertyType.MATERIAL,
                               GeometryType.TECHNOLOGICAL)

        # Apply values for the PropertyType.MATERIAL to regions
        self.cell.set_properties(
            {PropertyType.MATERIAL: ['MAT1', 'MAT2', 'MAT1']})
        # Verify cell's regions are shown according to the property values
        self.__assess_show(PropertyType.MATERIAL, GeometryType.TECHNOLOGICAL)
        self.__assess_show(PropertyType.MATERIAL, GeometryType.SECTORIZED)

    def test_translate(self) -> None:
        """
        Method that tests the implementation of the method `translate`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Store the elements for assessing the correct translation
        new_cntr = (1.0, 1.0, 0.0)
        new_cntr_vrtx = make_vertex(new_cntr)
        cdg_pre = make_cdg(self.cell.face)
        distance = get_min_distance(self.cell.figure.o, new_cntr_vrtx)
        cell_org = deepcopy(self.cell)

        # Translate the cell
        self.cell = self.cell.translate(new_cntr)
        # Check the correct translation happened
        self.assertTrue(
            are_same_shapes(
                self.cell.figure.o, new_cntr_vrtx, ShapeType.VERTEX)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre, make_cdg(self.cell.face)),
                distance)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre,
                                 make_cdg(self.cell.figure.out_circle)),
                distance)
        )
        for ic1, ic2 in zip(cell_org.inner_circles, self.cell.inner_circles):
            self.assertTrue(
                math.isclose(
                    get_min_distance(ic1.o, ic2.o), distance)
            )
        for tf1, tf2 in zip(cell_org.tech_geom_props,
                            self.cell.tech_geom_props):
            self.assertTrue(
                math.isclose(
                    get_min_distance(make_cdg(tf1), make_cdg(tf2)),
                    distance)
            )
        if self.cell.sectorized_face:
            self.assertTrue(
                math.isclose(
                    get_min_distance(
                        make_cdg(cell_org.sectorized_face),
                        make_cdg(self.cell.sectorized_face)),
                    distance)
            )

    def test_update_geometry(self) -> None:
        """
        Method that tests the implementation of the method `update_geometry`
        of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Verify an exception is raised when calling the method without
        # selecting any region
        with self.assertRaises(RuntimeError):
            self.cell.update_geometry()

    def test_update_geometry_from_face(self) -> None:
        """
        Method that tests the implementation of the method
        `update_geometry_from_face` of the `Cell` class.
        """
        # Check the cell istantiation
        self.__check_cell_setup()
        # Configure the cell
        self.cell.add_circle(0.1)
        self.cell.add_circle(0.3)
        self.cell.set_properties(
            {PropertyType.MATERIAL: ['MAT1', 'MAT2', 'MAT1']})
        # Build the face to update the cell's sectorized one with
        face = make_partition(
            [self.cell.face],
            [make_circle(self.o, None, 0.15),
             make_circle(self.o, None, 0.2),
             make_circle(self.o, None, 0.25)] + [
                make_edge(self.o, v) for v in self.cell.figure.vertices],
            ShapeType.FACE)
        # The cell's sectorized face should be None at this point
        self.assertIsNone(self.cell.sectorized_face)
        self.cell.update_geometry_from_face(GeometryType.SECTORIZED, face)
        # Verify the cell's sectorized face has been updated
        self.assertIsNotNone(self.cell.sectorized_face)
        for e1 in extract_sub_shapes(self.cell.sectorized_face,
                                     ShapeType.EDGE):
            found = False
            for e2 in extract_sub_shapes(face, ShapeType.EDGE):
                if are_same_shapes(e1, e2, ShapeType.EDGE):
                    found = True
                    break
            self.assertTrue(found)
        # Restore the cell geometry layout
        self.cell.restore()
        self.cell.add_circle(0.1)
        self.cell.add_circle(0.3)
        self.cell.set_properties(
            {PropertyType.MATERIAL: ['MAT1', 'MAT2', 'MAT1']})
        self.cell.sectorize(*self.sect_opts)

        # Update the cell's technological geometry with the same face
        face = make_partition(
            [self.cell.face],
            [make_circle(self.o, None, 0.15),
             make_circle(self.o, None, 0.2),
             make_circle(self.o, None, 0.25)],
            ShapeType.FACE)
        self.cell.update_geometry_from_face(GeometryType.TECHNOLOGICAL, face)
        for e1 in extract_sub_shapes(self.cell.face, ShapeType.EDGE):
            found = False
            for e2 in extract_sub_shapes(face, ShapeType.EDGE):
                if are_same_shapes(e1, e2, ShapeType.EDGE):
                    found = True
                    break
            self.assertTrue(found)
        face = make_partition(
            [self.cell.sectorized_face] + [face], [], ShapeType.FACE)
        for e1 in extract_sub_shapes(self.cell.sectorized_face,
                                     ShapeType.EDGE):
            found = False
            for e2 in extract_sub_shapes(face, ShapeType.EDGE):
                if are_same_shapes(e1, e2, ShapeType.EDGE):
                    found = True
                    break
            self.assertTrue(found)

        # Check the dictionaries of properties and sectorization options have
        # updated
        for f in extract_sub_shapes(self.cell.face, ShapeType.FACE):
            found = False
            for region in self.cell.tech_geom_props:
                if are_same_shapes(f, region, ShapeType.FACE):
                    found = True
                    self.assertTrue(
                        self.cell.tech_geom_props[region][
                            PropertyType.MATERIAL])
                    break
            self.assertTrue(found)
        for f in self.cell._Cell__extract_cell_centered_faces():
            found = False
            for region in self.cell.tech_geom_sect_opts:
                if are_same_shapes(f, region, ShapeType.FACE):
                    found = True
                    self.assertTrue(
                        self.cell.tech_geom_sect_opts[region])
                    break
            self.assertTrue(found)

    def __assess_add_circle(
            self,
            pos: Tuple[float, float, float] | None,
            radius: float,
            n_circles: int) -> None:
        """
        Method that assesses the operation of adding a circle to the face
        of a cell.

        Parameters
        ----------
        pos : Tuple[float, float, float] | None
            The XYZ coordinates of the circle to add, if any.
        radius : float
            The radius of the circle to add.
        n_circles : int
            An index indicating the current number of added circles.
        """
        # Add the circle
        self.cell.add_circle(radius, pos)
        # Verify the correct addition of a circle to the cell
        self.assertEqual(len(self.cell.inner_circles), n_circles+1)
        for i in range(len(self.cell.inner_circles)-1):
            self.assertTrue(
                self.cell.inner_circles[i].radius <=
                    self.cell.inner_circles[i+1].radius
            )
        if pos is None:
            pos = (0.0, 0.0, 0.0)
        self.assertTrue(
            self.__check_circle_presence(pos, radius)
        )
        for sf in extract_sub_shapes(self.cell.face, ShapeType.FACE):
            self.assertTrue(
                self.__check_face_in_dict(sf, self.cell.tech_geom_props)
            )
        # Check the dictionary for the sectorization options has updated
        if self.cell.sectorized_face is not None:
            for sf in self.cell._Cell__extract_cell_centered_faces():
                self.assertTrue(
                    self.__check_face_in_dict(sf,
                                              self.cell.tech_geom_sect_opts)
                )

    def __assess_cell_dict_initialization(self) -> None:
        """
        Method that tests that the cell's dictionaries are correctly
        initialized.
        """
        cell_regions = self.cell.extract_subfaces()
        self.assertEqual(self.cell.tech_geom_props, {
            region: {} for region in cell_regions
        })
        self.assertEqual(self.cell.tech_geom_sect_opts, {
            region: (1, 0) for region in cell_regions
        })

    def __assess_cell_face_initialization(self) -> None:
        """
        Method that tests that the cell's face is correctly set to the one
        of its corresponding `Surface`'s one.
        """
        self.assertTrue(
            are_same_shapes(self.cell.face, self.surf.face, ShapeType.FACE)
        )
        self.assertTrue(len(self.cell.inner_circles) == 0)
        self.assertIsNone(self.cell.sectorized_face)

    def __assess_get_regions_info(
            self, region_face: Any, str_to_find: str) -> None:
        """
        Method that assesses whether the right message is shown when
        displaying information about a specified region of the cell.

        Parameters
        ----------
        region_face : Any
            The cell's region to retrieve information from.
        str_to_find: str
            The informative string to find for in the one being produced.
        """
        # Retrieve the informative message about a region
        mssg = get_region_info(region_face, self.cell.regions)
        for region in self.cell.regions:
            if are_same_shapes(region.face, region_face, ShapeType.FACE):
                self.assertTrue(region.name in mssg)
                self.assertTrue(str_to_find in mssg)
                return

    def __assess_remove_circle(
            self,
            pos: Tuple[float, float, float] | None,
            radius: float,
            n_circles: int) -> None:
        """
        Method that assesses the operation of removing a circle from the face
        of a cell.

        Parameters
        ----------
        pos : Tuple[float, float, float] | None
            The XYZ coordinates of the circle to remove, if any.
        radius : float
            The radius of the circle to remove.
        n_circles : int
            An index indicating the current number of present circles.
        """
        # Remove the circle
        self.cell.remove_circle(radius, pos)
        # Verify the correct removal of a circle from the cell
        self.assertEqual(len(self.cell.inner_circles), n_circles-1)
        for i in range(len(self.cell.inner_circles)-1):
            self.assertTrue(
                self.cell.inner_circles[i].radius <=
                    self.cell.inner_circles[i+1].radius
            )
        if pos is None:
            pos = (0.0, 0.0, 0.0)
        self.assertFalse(
            self.__check_circle_presence(pos, radius)
        )
        for i, sf in enumerate(
                extract_sub_shapes(
                    make_compound([self.cell.face]), ShapeType.FACE)):
            self.assertTrue(
                self.__check_face_in_dict(sf, self.cell.tech_geom_props)
            )
        self.assertEqual(i+1, len(self.cell.tech_geom_props))
        # Check the dictionary for the sectorization options has updated
        if self.cell.sectorized_face is not None:
            for i, sf in enumerate(
                    extract_sub_shapes(
                        make_compound(
                            self.cell._Cell__extract_cell_centered_faces()),
                        ShapeType.FACE)):
                self.assertTrue(
                    self.__check_face_in_dict(sf,
                                              self.cell.tech_geom_sect_opts)
                )
            self.assertEqual(i+1, len(self.cell.tech_geom_sect_opts))

    def __assess_show(
            self,
            property: PropertyType | None,
            geometry: GeometryType) -> None:
        """
        Method that assesses the correct display of cell's regions according
        to whether regions should be colored by property values and to the
        geometry type, either `TECHNOLOGICAL` or `SECTORIZED`.

        Parameters
        ----------
        property : PropertyType | None
            The `PropertyType` according to which regions should be colored;
            if `None`, a default color must be used.
        geometry : GeometryType
            The `GeometryType` indicating whether regions are built from the
            technological or the sectorized geometry.
        """
        # Display the cell's technological geometry in the SALOME viewer
        self.cell.show(property, geometry)
        # Verify an entry ID has been associated to the cell's face
        self.assertIsNotNone(self.cell.face_entry_id)
        # Verify the regions corresponding to the geometry have been built
        self.assertTrue(len(self.cell.regions) > 0)
        if geometry == GeometryType.TECHNOLOGICAL:
            self.assertEqual(len(self.cell.regions),
                             len(self.cell.tech_geom_props.keys()))
            for region in self.cell.regions:
                found = False
                for tr in self.cell.tech_geom_props:
                    if are_same_shapes(region.face, tr, ShapeType.FACE):
                        found = True
                        break
                self.assertTrue(found)
        elif geometry == GeometryType.SECTORIZED:
            for region in self.cell.regions:
                found = False
                for tr in extract_sub_shapes(self.cell.sectorized_face,
                                             ShapeType.FACE):
                    if are_same_shapes(region.face, tr, ShapeType.FACE):
                        found = True
                        break
                self.assertTrue(found)
        # Verify if regions have the same default color or the ones having
        # the same property value share the same color
        if property is None:
            default_color = (167, 167, 167)
            self.assertTrue(
                all(region.color == default_color
                    for region in self.cell.regions)
            )
        else:
            prop_val_vs_color = {}
            for region in self.cell.regions:
                value = region.properties[property]
                if value in prop_val_vs_color:
                    self.assertEqual(region.color, prop_val_vs_color[value])
                else:
                    prop_val_vs_color[value] = region.color
        # Verify the regions have been displayed by checking their entry ID
        # has been defined so that they are children of the cell's face
        self.assertTrue(
            all(
                region.face_entry_id and region.face_entry_id.startswith(
                    self.cell.face_entry_id + ':')
                for region in self.cell.regions)
        )
        # Verify the attribute indicating the displayed geometry
        self.assertEqual(self.cell.displayed_geom, geometry)

    def __check_cell_setup(self) -> None:
        """
        Method that checks whether the test class has the `cell` attribute
        correctly initialized. If the `cell` attribute is `None`, the test
        this method is called in is skipped by providing an explanatory
        message.

        This method is called within each method of the class `TestCell`
        that operates on the subclass of the `Cell` class the attribute
        `cell` is set to.
        """
        if self.cell is None:
            self.skipTest("TestCell: attribute 'cell' not initialized. "
                          "This test is only valid in subclasses.")

    def __check_circle_presence(
            self, center: Tuple[float, float, float], radius: float) -> bool:
        """
        Method that looks for a circle among the edges extracted from the
        `face` attribute of a `Cell` instance.
        To be identified as the circle to find, the edge can be either a
        circle or an arc of circle with the same center position and radius
        given as inputs.

        Parameters
        ----------
        center: Tuple[float, float, float]
            The XYZ coordinates of the center of the circle to look for.
        radius: float
            The radius of the circle to look for.

        Returns
        -------
        bool:
            `True` if an edge with the given characteristics is found,
            `False` otherwise.
        """
        for edge in extract_sub_shapes(self.cell.face, ShapeType.EDGE):
            data = get_kind_of_shape(edge)
            if str(data[0]) in ["CIRCLE", "ARC_CIRCLE"]:
                # Verify the center and radius coincide with the given ones
                if tuple(data[1:4]) == center and data[7] == radius:
                    return True
        else:
            return False

    def __check_face_in_dict(self, face: Any, cell_dict: Dict) -> bool:
        """
        Method that checks if the given face object is present among the
        ones stored in the keys of the given dictionary.

        Parameters
        ----------
        face: Any
            The face object to look for.
        cell_dict: Dict
            The dictionary where to look for the face.

        Returns
        -------
        bool
            `True` if a face is found among the keys of the dictionary,
            `False` otherwise.
        """
        for f in cell_dict:
            if are_same_shapes(face, f, ShapeType.FACE):
                return True
        else:
            return False

class TestRectCell(TestCell):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `RectCell` class.

    This test suite provides common setup and a set of tests to ensure that
    implementations of `RectCell` correctly handle the operations performed
    on a cartesian cell, based on a `Rectangle` surface.
    Tests dealing with the cell's operations are declared in the `TestCell`
    class this class inherits from. They are run here, as this class declares
    the `cell` attribute.

    In addition, there are tests specific for a cartesian-type cell.

    Attributes
    ----------
    In addition to the attributes declared in the `TestCell` superclass,
    there are the following ones:

    name : str
        The name of the cell when displayed in the SALOME viewer.
    surf : Rectangle
        The `Surface` subclass providing the characteristic shape for a
        cartesian-type cell.
    cell : RectCell
        The `Cell` subclass that describes a cartesian-type cell.
    sect_opts : Tuple[List[int], List[float]]
        Providing the sectorization options as the number of sectors and the
        starting angle for each cell-centered region.
    height : float
        The height of the cartesian cell.
    width : float
        The width of the cartesian cell.
    """
    def setUp(self):
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.name = "Cartesian Cell"
        self.height = 1.0
        self.width = 1.0
        self.surf: Rectangle = Rectangle(
            center=get_point_coordinates(self.o),
            height=self.height,
            width=self.width
        )
        self.cell: RectCell = RectCell(
            center=get_point_coordinates(self.o),
            height_x_width=(self.height, self.width),
            name=self.name
        )
        self.sect_opts: Tuple[List[int], List[float]] = (
            [1, 4, 8], [0, 0, 22.5])

    def test_initialize_specific_cell(self) -> None:
        """
        Method for testing the initialization of the instance attributes
        that are specific for the `RectCell` class.
        """
        self.cell._initialize_specific_cell()
        self.assertIsNotNone(self.cell.windmills)
        self.assertEqual(len(self.cell.windmills), 0)

    def test_check_radius_vs_cell_dim(self) -> None:
        """
        Method for testing if a circle can be added to the cell by verifying
        that the circle's radius does not exceed the cell's dimensions.
        """
        # Test the exception is raised when the radius is not valid
        with self.assertRaises(RuntimeError):
            self.cell._check_radius_vs_cell_dim(1.0)
        # Test the check passes with a valid radius
        self.cell._check_radius_vs_cell_dim(0.1)

    def test_init(self) -> None:
        """
        Method that tests the correct initialization of the `RectCell` class.
        """
        self.assertEqual(self.cell.height, self.height)
        self.assertEqual(self.cell.width, self.width)
        self.assertEqual(self.cell.figure.lx, self.width)
        self.assertEqual(self.cell.figure.ly, self.height)
        self.assertEqual(self.cell.cell_type, CellType.RECT)

    def test_sectorize(self) -> None:
        """
        Method that tests the implementation of the method `sectorize`
        of the `RectCell` class.
        """
        # Add a circle to the cell
        radius = 0.2
        self.cell.add_circle(radius)
        # Test an exception is raised when passing an incorrect number of
        # sectorization options or the combination of sectors-angle is not
        # admitted
        with self.assertRaises(ValueError):
            self.cell.sectorize([4], [0])
        with self.assertRaises(ValueError):
            self.cell.sectorize([4, 5], [0, 10])
        # Verify the cell's sectorized face is not present
        self.assertIsNone(self.cell.sectorized_face)

        # Apply the sectorization without and with windmill and verify them
        self.__assess_sectorization(radius, ([4, 8], [0, 22.5]), False)
        self.__assess_sectorization(radius, ([4, 8], [0, 22.5]), True)

    def __assess_sectorization(
            self, radius, sect_opts: Tuple[List], windmill: bool) -> None:
        """
        Method that applies the sectorization to the instance of the
        `RectCell` class and verifies it has been applied correctly.

        Parameters
        ----------
        radius : float
            The radius of the only circle added to the cell.
        sect_opts : Tuple[List]
            A tuple containing two lists, the first indicating the number
            of sectors for each zone, the second the starting angle for
            the sectorization in each zone.
        windmill : bool
            Indicating whether the sectorization includes a windmill.
        """
        self.cell.sectorize(sect_opts[0], sect_opts[1], windmill=windmill)

        # Verify the sectorization happened successfully
        self.assertIsNotNone(self.cell.sectorized_face)
        n_zones = sum(sect_opts[0]) + (
            4 if self.cell.is_windmill_applied else 0)
        self.assertEqual(
            len(
                extract_sub_shapes(self.cell.sectorized_face,
                                   ShapeType.FACE)),
            n_zones
        )
        edges = extract_sub_shapes(self.cell.sectorized_face, ShapeType.EDGE)
        sect_edges = {}
        for e in edges:
            if str(get_kind_of_shape(e)[0]) != "SEGMENT":
                continue
            length = round(get_basic_properties(e)[0], 6)
            if (math.isclose(length, self.cell.figure.lx) or
                math.isclose(length, self.cell.figure.ly)):
                continue
            if length in sect_edges:
                sect_edges[length] += 1
            else:
                sect_edges[length] = 1

        # Assess the correct number of sectorization edges for the central
        # area, the outer one and the windmill, if applied
        l_outer = self.cell.figure.lx / math.cos(math.radians(22.5)) - radius
        l_wndml = (self.cell.figure.ly - self.cell.figure.lx / math.sin(
            math.radians(22.5))) / math.cos(math.pi/4)
        for e in sect_edges:
            if math.isclose(e, radius):
                self.assertEqual(sect_edges[e], sect_opts[0][0])
                continue
            elif math.isclose(e, l_outer):
                self.assertEqual(sect_edges[e], sect_opts[0][1])
                continue
            if self.cell.is_windmill_applied:
                if math.isclose(e, l_wndml):
                    self.assertEqual(sect_edges[e], 4)
                    continue


class TestHexCell(TestCell):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `HexCell` class.

    This test suite provides common setup and a set of tests to ensure that
    implementations of `HexCell` correctly handle the operations performed
    on a hexagonal cell, based on a `Hexagon` surface.
    Tests dealing with the cell's operations are declared in the `TestCell`
    class this class inherits from. They are run here, as this class declares
    the `cell` attribute.

    In addition, there are tests specific for a hexagonal-type cell.

    Attributes
    ----------
    In addition to the attributes declared in the `TestCell` superclass,
    there are the following ones:

    name : str
        The name of the cell when displayed in the SALOME viewer.
    surf : Hexagon
        The `Surface` subclass providing the characteristic shape for a
        hexagonal-type cell.
    cell : HexCell
        The `Cell` subclass that describes a hexagonal-type cell.
    sect_opts : Tuple[List[int], List[float]]
        Providing the sectorization options as the number of sectors and the
        starting angle for each cell-centered region.
    edge_length : float
        The length of the hexagon's edge.
    apothem : float
        The length of the hexagon's apothem.
    """
    def setUp(self):
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.name = "Hexagonal Cell"
        self.edge_length = 1.0
        self.apothem = self.edge_length * math.sin(math.pi/3)
        self.surf: Hexagon = Hexagon(
            center=get_point_coordinates(self.o),
            edge_length=1.0
        )
        self.cell: HexCell = HexCell(
            center=get_point_coordinates(self.o),
            edge_length=1.0,
            name=self.name
        )
        self.sect_opts: Tuple[List[int], List[float]] = (
            [1, 6, 6], [0, 0, 0])

    def test_check_radius_vs_cell_dim(self) -> None:
        """
        Method for testing if a circle can be added to the cell by verifying
        that the circle's radius does not exceed the cell's dimensions.
        """
        # Test the exception is raised when the radius is not valid
        with self.assertRaises(RuntimeError):
            self.cell._check_radius_vs_cell_dim(1.0)
        # Test the check passes with a valid radius
        self.cell._check_radius_vs_cell_dim(0.1)

    def test_init(self) -> None:
        """
        Method that tests the correct initialization of the `HexCell` class.
        """
        self.assertEqual(self.cell.edge_length, self.edge_length)
        self.assertEqual(self.cell.apothem, self.apothem)
        self.assertEqual(self.cell.figure.lx, self.edge_length)
        self.assertEqual(self.cell.figure.ly, self.apothem)
        self.assertEqual(self.cell.cell_type, CellType.HEX)

    def test_sectorize(self) -> None:
        """
        Method that tests the implementation of the method `sectorize`
        of the `HexCell` class.
        """
        # Add a circle to the cell
        radius = 0.2
        self.cell.add_circle(radius)
        # Test an exception is raised when passing an incorrect number of
        # sectorization options or the combination of sectors-angle is not
        # admitted
        with self.assertRaises(ValueError):
            self.cell.sectorize([6], [0])
        with self.assertRaises(ValueError):
            self.cell.sectorize([4, 5], [0, 10])
        # Verify the cell's sectorized face is not present
        self.assertIsNone(self.cell.sectorized_face)

        # Apply the sectorization and verify it
        self.__assess_sectorization(
            radius, (self.sect_opts[0][:2], self.sect_opts[1][:2]))

    def __assess_sectorization(self, radius, sect_opts: Tuple[List]) -> None:
        """
        Method that applies the sectorization to the instance of the
        `HexCell` class and verifies it has been applied correctly.

        Parameters
        ----------
        radius : float
            The radius of the only circle added to the cell.
        sect_opts : Tuple[List]
            A tuple containing two lists, the first indicating the number
            of sectors for each zone, the second the starting angle for
            the sectorization in each zone.
        """
        self.cell.sectorize(sect_opts[0], sect_opts[1])

        # Verify the sectorization happened successfully
        self.assertIsNotNone(self.cell.sectorized_face)
        n_zones = sum(sect_opts[0])
        self.assertEqual(
            len(
                extract_sub_shapes(self.cell.sectorized_face,
                                   ShapeType.FACE)),
            n_zones
        )
        edges = extract_sub_shapes(self.cell.sectorized_face, ShapeType.EDGE)
        sect_edges = {}
        for e in edges:
            if str(get_kind_of_shape(e)[0]) != "SEGMENT":
                continue
            length = round(get_basic_properties(e)[0], 6)
            if math.isclose(length, self.cell.figure.lx):
                continue
            if length in sect_edges:
                sect_edges[length] += 1
            else:
                sect_edges[length] = 1

        # Assess the correct number of sectorization edges for the central
        # area and the outer one
        l_outer = self.cell.figure.lx - radius
        for e in sect_edges:
            if math.isclose(e, radius):
                self.assertEqual(sect_edges[e], sect_opts[0][0])
                continue
            elif math.isclose(e, l_outer):
                self.assertEqual(sect_edges[e], sect_opts[0][1])
                continue


if __name__ == "__main__":
    unittest.main()