import os
import unittest

from typing import Any, List, Tuple
from pathlib import Path

from glow.geometry_layouts.cells import RectCell
from glow.geometry_layouts.geometries import Rectangle
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import ShapeType, add_to_study, extract_sub_shapes, get_bounding_box, get_kind_of_shape, make_common, make_compound, make_partition
from glow.main import TdtSetup, analyse_and_generate_tdt
from glow.support.types import EDGE_NAME_VS_TYPE, EdgeType, GeometryType, \
    LatticeGeometryType, PropertyType, SymmetryType
from support_funcs import build_colorset, compute_hash


class TestTdtSetup(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `TdtSetup`
    dataclass. In particular, the logic that stores the value of the `albedo`
    attribute is tested.
    """
    def test_init_valid_albedo(self):
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

    def test_init_invalid_albedo(self):
        """
        Method that tests that a invalid value for the albedo raises an
        exception.
        """
        with self.assertRaises(RuntimeError):
            TdtSetup(albedo=-0.1)
        with self.assertRaises(RuntimeError):
            TdtSetup(albedo=1.1)


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
    format : str
        A string indicating the format with which the geometric data about
        the edges and the BCs is written to file.
    colorset : List[Lattice]
        A list of ``Lattice`` instances each positioned appropriately to
        replicate a colorset.
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
        self.format: str = f"{{:.{7}E}}"
        self.colorset: List[Lattice] = build_colorset(self.lattice)

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
        analyse_and_generate_tdt(
            lattices=[self.lattice],
            filename=self.file_name.split('.')[0],
            tdt_config=TdtSetup(
                self.geom_type, self.prop_type, 0.0
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
        no_regions = len(self.lattice.regions)

        # Verify the header section
        self.__assess_header_section(
            content,
            no_regions,
            len(edges),
            LatticeGeometryType.RECTANGLE_TRAN
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

    def test_analyse_and_generate_tdt_colorset(self) -> None:
        """
        Method that tests the implementation of the function
        `analyse_and_generate_tdt` declared in the `main.py`
        module when a colorset or its portion is provided.
        """
        # Test exceptions are raised
        # 1) colorset with any lattice having a symmetry
        self.colorset[0].apply_symmetry(SymmetryType.EIGHTH)
        with self.assertRaises(RuntimeError):
            analyse_and_generate_tdt(
                lattices=self.colorset,
                filename=self.file_name.split('.')[0],
                tdt_config=TdtSetup(
                    self.geom_type, self.prop_type, 0.0
                )
            )
        # 2) colorset with any lattice having a symmetry and the analysis is
        #    done on a portion
        colorset_cmpd = make_partition(
            [l.lattice_cmpd for l in self.colorset], [], ShapeType.COMPOUND)
        colorset_portion = make_common(
                colorset_cmpd, Rectangle((3.2, 3.2, 0.0), 3.2, 3.2).face)
        with self.assertRaises(RuntimeError):
            analyse_and_generate_tdt(
                lattices=self.colorset,
                filename=self.file_name.split('.')[0],
                tdt_config=TdtSetup(
                    self.geom_type, self.prop_type, 0.0
                ),
                compound_to_export=colorset_portion
            )
        # 3) colorset with the analysis performed on a portion and
        #    inconsistency of values in TdTSetup
        with self.assertRaises(RuntimeError):
            analyse_and_generate_tdt(
                lattices=self.colorset,
                filename=self.file_name.split('.')[0],
                tdt_config=TdtSetup(
                    type_geo=LatticeGeometryType.RECTANGLE_TRAN,
                    symmetry_type=SymmetryType.QUARTER
                ),
                compound_to_export=colorset_portion
            )

        # Reset the symmetry application and set the typegeo of the reference
        # lattice (i.e. the first one in the list)
        self.colorset[0].apply_symmetry(SymmetryType.FULL)
        self.colorset[0].type_geo = LatticeGeometryType.RECTANGLE_TRAN
        # Call the function to test with a full colorset
        self.__assess_tdt_colorset(
            'de6955b64217cf08cb0875cf2da14d8a71e7583b3fc9c588fe778316bc545b09'
        )
        # Call the function to test with a portion of the colorset
        self.__assess_tdt_colorset(
            'b481a9fa3b4bf2d9af86ed9f2b674a9f1b9162d76fc09319c672c5b1bb3065d1',
            colorset_portion
        )

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
            type_geo: LatticeGeometryType) -> None:
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
        type_geo : LatticeGeometryType
            The value of the `typgeo` parameter.
        """
        self.assertIn(
            f"{type_geo.value:5d},{0:5d}, " + \
            f"{regions_nb:5d}, {edges_nb:5d}, " +\
            f"{1:5d},{regions_nb:5d}, {0:5d}, {1:5d}", content)
        # Verify the values for impressions and precisions
        self.assertIn(f"{0:5d}  {0:6d}  1", content)
        self.assertIn(f"{1e-5:7E}   {1e-5:7E}", content)

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
        analyse_and_generate_tdt(
            lattices=self.colorset,
            filename=self.file_name.split('.')[0],
            tdt_config=TdtSetup(
                type_geo=LatticeGeometryType.SYMMETRIES_TWO,
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

if __name__ == "__main__":
    cell = RectCell()
    cell.set_properties({PropertyType.MATERIAL: ['MAT']})
    lattice: Lattice = Lattice([cell])
    geom_type: GeometryType = GeometryType.TECHNOLOGICAL
    prop_type: PropertyType = PropertyType.MATERIAL
    lattice.build_regions(geom_type)
    format: str = f"{{:.{7}E}}"
    colorset: List[Lattice] = build_colorset(lattice)
    add_to_study(make_compound([l.lattice_cmpd for l in colorset]), "sa")

    analyse_and_generate_tdt([lattice], "pippo")

    Rectangle((1.6, 1.6, 0.0), 3.2, 3.2).show_face()

    colorset_cmpd = make_partition(
        [l.lattice_cmpd for l in colorset], [], ShapeType.COMPOUND)
    colorset_portion = make_partition(
        [make_common(
            colorset_cmpd, Rectangle((3.2, 3.2, 0.0), 3.2, 3.2).face)],
        [],
        ShapeType.COMPOUND
    )
    add_to_study(colorset_portion, "portion")

    analyse_and_generate_tdt(
        colorset,
        "test_colorset",
        TdtSetup(
            type_geo=LatticeGeometryType.SYMMETRIES_TWO,
            symmetry_type=SymmetryType.QUARTER
        ),
        colorset_portion
    )
