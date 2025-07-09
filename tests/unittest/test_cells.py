"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.cells` module have a valid implementation.
"""
from typing import Any, Dict, Tuple, Union
import unittest

from glow.generator.support import PropertyType
from glow.geometry_layouts.cells import Region
from glow.geometry_layouts.utility import are_same_shapes
from glow.interface.geom_interface import *


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
