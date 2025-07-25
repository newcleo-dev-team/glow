import os
import sys
import unittest

from pathlib import Path

from glow.geometry_layouts.cells import RectCell
from glow.geometry_layouts.lattices import Lattice
from glow.main import analyse_and_generate_tdt
from glow.support.types import EDGE_NAME_VS_TYPE, EdgeType, GeometryType, \
    LatticeGeometryType, PropertyType


class TestMainFunction(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the function
    `analyse_and_generate_tdt` declared in the `main.py` module.
    The test verifies that the characteristics of the geometry layout
    are correctly extracted and the output TDT is generated.

    Attributes
    ----------
    lattice : Lattice
        The `Lattice` instance whose geometric characteristics are
        exported to file.
    file_name : str
        The path name of the output TDT file.
    geom_type : GeometryType
        The type of geometry for the lattice's cells, as element of the
        `GeometryType` enumeration.
    prop_type : PropertyType
        The type of property assigned to the lattice's regions, as element
        of the `PropertyType` enumeration.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the function
        `analyse_and_generate_tdt` of the 'main.py' module.
        It initializes the attributes common to all the tests.
        """
        cell = RectCell()
        cell.set_properties({PropertyType.MATERIAL: ['MAT']})
        self.lattice: Lattice = Lattice([cell])
        self.file_name: str = str(Path(__file__).parent / 'tdt_layout.dat')
        self.geom_type: GeometryType = GeometryType.TECHNOLOGICAL
        self.prop_type: PropertyType = PropertyType.MATERIAL
        self.lattice.build_regions(self.geom_type)

    def tearDown(self):
        """
        Method that is run after calling each of the tests. It removes the
        generated output TDT file, if exists.
        """
        # Remove the generated TDT file
        if os.path.exists(self.file_name):
            os.remove(self.file_name)

    def test_analyse_and_generate_tdt(self) -> None:
        """
        Method that tests the implementation of the function
        `analyse_and_generate_tdt` declared in the `main.py`
        module.
        """
        # Call the function to test
        analyse_and_generate_tdt(lattice=self.lattice,
                                 filename=self.file_name.split('.')[0],
                                 geom_type=self.geom_type,
                                 property_type=self.prop_type)

        # Verify that the output TDT file has been generated
        self.assertTrue(os.path.exists(self.file_name))
        # Read the content of the TDT file and verify the correct data is
        # present
        with open(self.file_name, "r") as f:
            content = f.read()

        # Verify the header section
        self.__assess_header_section(content)
        # Verify the regions section
        self.__assess_regions_section(content)
        # Verify the edges section (only segment-type edges are present)
        self.__assess_edges_section(content)
        # Verify the boundaries section
        self.__assess_boundaries_section(content)
        # Verify the properties section
        self.assertIn(f"# {1:2d} - {'MAT'}", content)
        self.assertIn(f"  {1}", content)

    def __assess_boundaries_section(self, content) -> None:
        """
        Method that tests that the boundaries section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        """
        # List storing the boundary axes for each edge
        boundaries = [
            (0, 1.0, 0.0),
            (-1.0, 0.0, 90.0),
            (0, -1.0, 0.0),
            (1.0, 0.0, 90.0),
        ]
        self.assertIn(f"  0, {len(boundaries)}, 0", content)
        self.assertIn("* albedo\n  1.0", content) # FIXME have attribute
        for i, bc in enumerate(boundaries):
            self.assertIn("* type  number of elements\n  2, 1", content)
            self.assertIn(f"*   elements\n{i+1}", content)
            self.assertIn("* tx, ty, angle", content)
            self.assertIn(f"{bc[0]:6f} {bc[1]:6f} {bc[2]:6f}", content)

    def __assess_edges_section(self, content) -> None:
        """
        Method that tests that the edges section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        """
        # List declaring the edge's XY coordinates of the starting point
        # and the dx, dy values
        edges = [
            (0.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0, 1.0),
            (1.0, 1.0, -1.0, 0.0),
            (0.0, 1.0, 0.0, -1.0)
        ]
        for i, edge in enumerate(edges):
            self.assertIn(
                f"* ELEM  {i+1}  {EDGE_NAME_VS_TYPE['SEGMENT'][1]}", content)
            self.assertIn(f" {EdgeType.SEGMENT.value}, 0, {1}", content)
            self.assertIn(
                f"  {edge[0]:6f}, {edge[1]:6f}, {edge[2]:6f}, {edge[3]:6f}",
                content)

    def __assess_header_section(self, content) -> None:
        """
        Method that tests that the header section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        """
        regions_nb = len(self.lattice.regions)
        edges_nb = 4
        self.assertIn(
            f"{LatticeGeometryType.RECTANGLE_TRAN.value:5d},{0:5d}, " + \
            f"{regions_nb:5d}, {edges_nb:5d}, " +\
            f"{1:5d},{regions_nb:5d}, {0:5d}, {1:5d}", content)
        # Verify the values for impressions and precisions
        self.assertIn(f"{0:5d}  {0:6d}  1", content)
        self.assertIn(f"{1e-5:7E}   {1e-5:7E}", content)

    def __assess_regions_section(self, content) -> None:
        """
        Method that tests that the regions section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        """
        regions_indices = ""
        for i in range(0, len(self.lattice.regions), 12):
            regions_indices += ",".join(f"{(j+1):4d}"
                for j, _ in enumerate(self.lattice.regions[i:i+12])) + ",\n"
        # Remove the last comma
        regions_indices = \
            "*   flux region number per geometry region (mesh)\n" + \
            regions_indices[0:-2] + "\n"
        self.assertIn(regions_indices, content)

