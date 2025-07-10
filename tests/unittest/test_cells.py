"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.cells` module have a valid implementation.
"""
from abc import ABC
from copy import deepcopy
import math
from typing import Any, Dict, Tuple, Union
import unittest

from glow.generator.support import PropertyType
from glow.geometry_layouts.cells import Cell, RectCell, Region
from glow.geometry_layouts.geometries import Rectangle, Surface
from glow.geometry_layouts.utility import are_same_shapes
from glow.interface.geom_interface import *
from support_funcs import make_ref_vectors


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
        self.cell.sectorize([4, 1], [0]*2)
        self.__assess_add_circle(
            None, 0.2, len(self.cell.inner_circles))

        # Restore the cell to the initial condition, so that other tests are
        # not affected by this one
        self.cell.restore()

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
        self.cell.sectorize([4, 1], [0]*2)
        self.__assess_remove_circle(
            None, radius, len(self.cell.inner_circles))
        self.__assess_remove_circle(
            position, radius, len(self.cell.inner_circles))
        # Verify the cell's face returned to the one of its figure
        self.assertTrue(
            are_same_shapes(
                self.cell.face, self.cell.figure.face, ShapeType.FACE)
        )

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
        ref_vectors1 = self.__build_cell_ref_vectors()
        self.cell.rotate(rot_angle)
        ref_vectors2 = self.__build_cell_ref_vectors()
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
        ref_vectors1 = self.__build_cell_ref_vectors()
        self.cell.rotate_from_axis(rot_angle, rot_axis)
        ref_vectors2 = self.__build_cell_ref_vectors()
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

    def __build_cell_ref_vectors(self) -> List[Any]:
        """
        Method that builds a list of vector objects on the first edge of
        each of the face objects belonging to a `Cell` instance, i.e. the
        whole cell's face and the faces stored in the dictionaries of
        properties and sectorization options.

        Returns
        -------
        List[Any]
            A list of vectors built of the first edge of the faces in the
            cell.
        """
        # List all the face objects comprising the cell's face and those
        # stored in the cell's dictionaries
        faces = [self.cell.face] + list(self.cell.tech_geom_props.keys())
        if self.cell.sectorized_face:
            faces += list(self.cell.tech_geom_sect_opts.keys())
        # Store the edge '0' for each face
        edges = [extract_sub_shapes(f, ShapeType.EDGE)[0] for f in faces]
        # Append the edge '0' for the cell's 'Surface' and reference circle
        edges.append(self.cell.figure.borders[0])
        edges.append(self.cell.figure.out_circle)
        # Build reference vectors for each edge
        return [
            make_vector_from_points(
                self.cell.figure.o,
                make_vertex_on_curve(e, 0.0)) for e in edges
        ]

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
            if are_same_shapes(face, f, ShapeType.FACE) and (
                get_basic_properties(face) == get_basic_properties(f)
            ):
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
    """
    def setUp(self):
        # Setup the common geometric elements
        super().setUp()
        # Setup the specific attributes for testing the `Surface` subclass
        self.name = "Cartesian Cell"
        self.surf: Rectangle = Rectangle(
            center=get_point_coordinates(self.o),
            height=1.0,
            width=1.0
        )
        self.cell: RectCell = RectCell(
            center=get_point_coordinates(self.o),
            height_x_width=(self.surf.ly, self.surf.lx),
            name=self.name
        )

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


if __name__ == "__main__":
    unittest.main()