"""
Module containing unittest classes to assess that the functions of the
`glow.main` module have a valid implementation.
"""
import os
import unittest

from pathlib import Path
from typing import Any, List, Tuple

from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.geometries import Rectangle
from glow.geometry_layouts.lattices import CartesianLattice
from glow.main import TdtSetup, export_layout_to_tdt
from glow.support.types import EDGE_NAME_VS_TYPE, EdgeType, GeometryType, \
    LayoutGeometryType, PropertyType, SymmetryType
from tests.unittest.support_funcs import build_colorset, compute_hash


class TestTdtSetup(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `TdtSetup`
    dataclass. In particular, the logic that stores the value of the `albedo`
    and the `property_types` attributes is tested.
    """
    def test_init_default_values(self) -> None:
        """
        Method that tests the correct assignment of the default values.
        """
        tdt_setup = TdtSetup()
        self.assertEqual(tdt_setup.geom_type, GeometryType.TECHNOLOGICAL)
        self.assertEqual(tdt_setup.property_types, [PropertyType.MATERIAL])
        self.assertEqual(tdt_setup.albedo, None)
        self.assertEqual(tdt_setup.type_geo, LayoutGeometryType.ISOTROPIC)
        self.assertEqual(tdt_setup.symmetry_type, SymmetryType.FULL)
        # Verify the 'layout_type' attribute is not present at initialisation
        self.assertTrue(not hasattr(tdt_setup, "layout_type"))

    def test_init_valid_albedo(self) -> None:
        """
        Method that tests that a valid value for the albedo is accepted.
        """
        # Test lower/upper bounds
        self.assertEqual(TdtSetup(albedo=0.0).albedo, 0.0)
        self.assertEqual(TdtSetup(albedo=1.0).albedo, 1.0)
        # Test a value in the validity range
        self.assertEqual(TdtSetup(albedo=0.5).albedo, 0.5)
        # Test when None is provided or the default value is kept
        self.assertIsNone(TdtSetup(albedo=None).albedo)
        self.assertIsNone(TdtSetup().albedo)
        # Test values very close to bounds
        self.assertEqual(TdtSetup(albedo=1e-12).albedo, 1e-12)
        self.assertEqual(TdtSetup(albedo=1.0 - 1e-12).albedo, 1.0 - 1e-12)

    def test_init_invalid_albedo(self) -> None:
        """
        Method that tests that a invalid value for the albedo raises an
        exception.
        """
        with self.assertRaises(RuntimeError):
            TdtSetup(albedo=-0.1)
        with self.assertRaises(RuntimeError):
            TdtSetup(albedo=1.1)

    def test_init_properties(self) -> None:
        """
        Method that tests the assignment of the `property_types` attribute.
        """
        # Test that properties are stored as a list
        self.assertTrue(isinstance(TdtSetup().property_types, List))
        # Test that the default value is 'PropertyType.MATERIAL'
        self.assertEqual(TdtSetup().property_types, [PropertyType.MATERIAL])
        # Test the assignment of one property type
        self.assertEqual(
            TdtSetup(property_types=PropertyType.MACRO).property_types,
            [PropertyType.MACRO]
        )
        # Test the assignment of two property type
        self.assertEqual(
            TdtSetup(
                property_types=[PropertyType.MACRO, PropertyType.MATERIAL]
            ).property_types,
            [PropertyType.MACRO, PropertyType.MATERIAL]
        )

class TestMainFunction(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the function
    `export_layout_to_tdt` declared in the `main.py` module.
    The test verifies that the characteristics of the geometry layout
    are correctly extracted and the output TDT file is generated.

    Attributes
    ----------
    layout : CartesianCell
        The `CartesianCell` instance whose geometric characteristics are
        exported to file.
    file_name : str
        The path name of the output TDT file.
    geom_type : GeometryType
        The type of geometry for the layout, as element of the `GeometryType`
        enumeration.
    prop_type : PropertyType
        The type of property assigned to the layout's regions, as element
        of the `PropertyType` enumeration.
    format : str
        A string indicating the format with which the geometric data about
        the edges and the BCs is written to file.
    colorset : Lattice
        A `Lattice` instance made by assemblies each positioned appropriately
        to replicate a colorset.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the function
        `export_layout_to_tdt` of the `main.py` module.
        It initializes the attributes common to all the tests.
        """
        self.layout = CartesianCell(
            base_props={PropertyType.MATERIAL: 'MAT'}
        )
        self.file_name: str = str(Path(__file__).parent / 'tdt_layout.dat')
        self.geom_type: GeometryType = GeometryType.TECHNOLOGICAL
        self.prop_type: PropertyType = PropertyType.MATERIAL
        self.format: str = f"{{:.{7}E}}"
        self.colorset: CartesianLattice = build_colorset(self.layout)

    def tearDown(self) -> None:
        """
        Method that is run after calling each of the tests. It removes the
        generated output TDT file, if exists.
        """
        # Remove the generated TDT file
        if os.path.exists(self.file_name):
            os.remove(self.file_name)

    def test_export_layout_to_tdt_cell(self) -> None:
        """
        Method that tests the implementation of the `export_layout_to_tdt`
        function declared in the `main.py` module for a layout made by a
        single Cartesian cell.
        """
        # Call the function to test
        export_layout_to_tdt(
            layout=self.layout,
            filename=self.file_name.split('.')[0],
            tdt_setup=TdtSetup(
                self.geom_type,
                self.prop_type,
                0.0,
                LayoutGeometryType.RECTANGLE_TRAN
            )
        )

        # Verify that the output TDT file has been generated
        self.assertTrue(os.path.exists(self.file_name))
        # Read the content of the TDT file and verify the correct data is
        # present
        with open(self.file_name, "r") as f:
            content = f.read()

        # Declare the data to look for in the TDT file content
        # List storing the boundary axes for each edge
        boundaries = [
            (0, 1.0, 0.0),
            (-1.0, 0.0, 90.0),
            (0, -1.0, 0.0),
            (1.0, 0.0, 90.0),
        ]
        # List declaring the edges' XY coordinates of the starting point
        # and the dx, dy values
        edges = [
            (0.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0, 1.0),
            (1.0, 1.0, -1.0, 0.0),
            (0.0, 1.0, 0.0, -1.0)
        ]
        # Number of regions in the lattice
        no_regions = len(self.layout.get_regions())

        # Verify the header section
        self.__assess_header_section(
            content,
            no_regions,
            len(edges),
            LayoutGeometryType.RECTANGLE_TRAN
        )
        # Verify the regions section
        self.__assess_regions_section(content, no_regions)
        # Verify the edges section (only segment-type edges are present)
        self.__assess_edges_section(content, edges)
        # Verify the boundaries section
        self.__assess_boundaries_section(content, boundaries)
        # Verify the properties section
        self.assertIn(f"# {1:2d} - {'MAT'}", content)
        self.assertIn(f"  {1}", content)

    def test_export_layout_to_tdt_colorset(self) -> None:
        """
        Method that tests the implementation of the `export_layout_to_tdt`
        function declared in the `main.py` module when a colorset or its
        portion is provided.
        """
        # Test exceptions are raised
        # 1) colorset with inconsistency between layout-symmetry types and the
        #    'typgeo' value
        with self.assertRaises(RuntimeError):
            export_layout_to_tdt(
                layout=self.colorset,
                filename=self.file_name.split('.')[0],
                tdt_setup=TdtSetup(
                    self.geom_type,
                    self.prop_type,
                    0.0,
                    LayoutGeometryType.R120
                )
            )
        # 2) export a compound which is not a portion of the colorset
        portion = Rectangle((13.25, 3.25, 0.0), 3.25, 3.25)
        with self.assertRaises(RuntimeError):
            export_layout_to_tdt(
                layout=self.colorset,
                filename=self.file_name.split('.')[0],
                tdt_setup=TdtSetup(
                    self.geom_type, self.prop_type, 0.0
                ),
                compound_to_export=portion
            )

        # Apply a QUARTER symmetry so that the export can be done on a known
        # symmetry
        self.colorset.apply_symmetry(SymmetryType.QUARTER)
        # Call the function to test with a full colorset
        self.__assess_tdt_colorset(
            'e0190dde26720444e283d35e07b7c0e08959baf70137070b2fd7e0463192b48b'
        )
        # Call the function to test with a portion of the colorset
        self.__assess_tdt_colorset(
            '8a7c648b7e68c1b91d5b36bf5ae9d40fab4c58820652ee5c899f48cf65cd3a48',
            self.colorset * Rectangle((3.2, 3.2, 0.0), 3.2, 3.2)
        )

    def test_export_layout_to_tdt_macros(self) -> None:
        """
        Method that tests the implementation of the `export_layout_to_tdt`
        function declared in the `main.py` module when a layout with the
        `PropertyType.MACRO` is provided.
        """
        # Build a lattice of the same cell and apply different values for the
        # 'MACRO' to each cell
        lattice = CartesianLattice([self.layout])
        lattice.add_ring_of_cells(self.layout, 1)
        macro_index = 0
        for layer in lattice.layers:
            for cell in layer:
                # Increment the macro index for each cell
                macro_index += 1
                for cell_layer in cell.layers:
                    for cell_region in cell_layer:
                        cell_region.properties.update(
                            {PropertyType.MACRO: "MAC_00" + str(macro_index)}
                        )

        # Verify the exception is raised if not including 'MATERIAL' among
        # the properties to export
        with self.assertRaises(RuntimeError):
            export_layout_to_tdt(
                layout=lattice,
                filename=self.file_name.split('.')[0],
                tdt_setup=TdtSetup(
                    self.geom_type,
                    PropertyType.MACRO,
                )
            )

        # Export the lattice with both 'MATERIAL' and 'MACRO' properties
        export_layout_to_tdt(
            layout=lattice,
            filename=self.file_name.split('.')[0],
            tdt_setup=TdtSetup(
                self.geom_type,
                [PropertyType.MACRO, PropertyType.MATERIAL],
            )
        )
        # Read the content of the TDT file
        with open(self.file_name, "r") as f:
            content = f.read()
        # Verify the correctness of the macros section in the TDT file
        self.__assess_macro_regions_section(content, 9, "MAC_00")

    def __assess_boundaries_section(
            self, content: str, boundaries: List[Tuple[float]]) -> None:
        """
        Method that tests that the boundaries section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        boundaries : List[Tuple[float]]
            The XY coordinates and angle of the axes of the boundaries.
        """
        self.assertIn(f"  0, {len(boundaries)}, 0", content)
        self.assertIn("* albedo\n  0.0", content)
        for i, bc in enumerate(boundaries):
            self.assertIn("* type  number of elements\n  2, 1", content)
            self.assertIn(f"*   elements\n{i+1}", content)
            self.assertIn("* tx, ty, angle", content)
            self.assertIn(
                f"{self.format.format(bc[0])} {self.format.format(bc[1])}" + \
                f" {self.format.format(bc[2])}",
                content
            )

    def __assess_edges_section(
            self, content: str, edges: List[Tuple[float]]) -> None:
        """
        Method that tests that the edges section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        edges : List[Tuple[float]]
            The list of the edges' XY coordinates of the starting point
            and the dx, dy values
        """
        for i, edge in enumerate(edges):
            self.assertIn(
                f"* ELEM  {i+1}  {EDGE_NAME_VS_TYPE['SEGMENT'][1]}", content)
            self.assertIn(f" {EdgeType.SEGMENT.value}, 0, {1}", content)
            self.assertIn(
                f"  {self.format.format(edge[0])}, " + \
                f"{self.format.format(edge[1])}, " + \
                f"{self.format.format(edge[2])}, " + \
                f"{self.format.format(edge[3])}",
                content)

    def __assess_header_section(
            self,
            content: str,
            regions_nb: int,
            edges_nb: int,
            type_geo: LayoutGeometryType) -> None:
        """
        Method that tests that the header section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        regions_nb : int
            The number of regions.
        edges_nb : int
            The number of edges.
        type_geo : LayoutGeometryType
            The value of the `typgeo` parameter.
        """
        self.assertIn(
            f"{type_geo.value:5d},{0:5d}, " + \
            f"{regions_nb:5d}, {edges_nb:5d}, " +\
            f"{1:5d},{regions_nb:5d}, {0:5d}, {1:5d}", content)
        # Verify the values for impressions and precisions
        self.assertIn(f"{0:5d}  {0:6d}  1", content)
        self.assertIn(f"{1e-5:7E}   {1e-5:7E}", content)

    def __assess_macro_regions_section(
            self, content: str, no_regions: int, mac_name_prefix: str
        ) -> None:
        """
        Method that tests that the macro regions section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        no_regions : int
            The number of regions in the layout.
        """
        # Verify the presence of the names of the macros
        macro_names = ""
        for i in range(0, no_regions, 4):
            macro_names += " ".join(f"{mac_name_prefix}{j+1}"
                for j in range(i, min(i+4, no_regions))) + "\n "
        # Remove the last space
        macro_names = "*   names of macros\n " + macro_names[:-1]
        self.assertIn(macro_names, content)

        # Verify the presence of the indices of the macros ordered per flux
        # region
        macro_indices = ""
        for i in range(0, no_regions, 12):
            macro_indices += ",".join(f"{(j+1):3d}"
                for j in range(i, min(i+12, no_regions))) + ",\n"
        # Remove the last comma
        macro_indices = \
            "*   macro order number per flux region\n" + \
            macro_indices[0:-2] + "\n"
        self.assertIn(macro_indices, content)

    def __assess_regions_section(self, content: str, no_regions: int) -> None:
        """
        Method that tests that the regions section, contained in the
        input string, contains the expected data.

        Parameters
        ----------
        content : str
            The content string to search for the expected data in.
        no_regions : int
            The number of regions in the layout.
        """
        regions_indices = ""
        for i in range(0, no_regions, 12):
            regions_indices += ",".join(f"{(j+1):4d}"
                for j in range(i, min(i+12, no_regions))) + ",\n"
        # Remove the last comma
        regions_indices = \
            "*   flux region number per geometry region (mesh)\n" + \
            regions_indices[0:-2] + "\n"
        self.assertIn(regions_indices, content)

    def __assess_tdt_colorset(
            self, ref_hash: str, layout: Any | None = None) -> None:
        """
        Method that tests that the TDT file generated from a colorset or its
        portion is correct by comparing its SHA256 hash against the given
        pre-calculated hash of the reference file.
        """
        # Run the analysis and TDT file generation
        export_layout_to_tdt(
            layout=self.colorset,
            filename=self.file_name.split('.')[0],
            tdt_setup=TdtSetup(
                type_geo=LayoutGeometryType.SYMMETRIES_TWO,
                symmetry_type=SymmetryType.QUARTER
            ),
            compound_to_export=layout
        )

        # Verify that the output TDT file has been generated
        self.assertTrue(os.path.exists(self.file_name))
        # Read the content of the TDT file and verify the correct data is
        # present
        with open(self.file_name, "r") as f:
            _ = f.read()

        # Verify the SHA256 hash of the generated file is the same of the one
        # of the reference TDT file
        generated_tdt_hash = compute_hash(self.file_name)
        self.assertEqual(generated_tdt_hash, ref_hash)
