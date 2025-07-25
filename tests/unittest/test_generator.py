"""
Module containing unittest classes to assess that the classes and functions
of the `glow.generator.generator` module have a valid implementation.
"""
import io
import os
import unittest

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List

from glow.generator.generator import TdtData, _write_boundary_conditions, \
    _write_edges, _write_header, _write_properties, _write_regions, \
    write_tdt_file
from glow.generator.geom_extractor import Boundary, Edge, Face
from glow.geometry_layouts.geometries import Rectangle
from glow.interface.geom_interface import ShapeType, set_shape_name

from glow.support.types import EDGE_NAME_VS_TYPE, BoundaryType, EdgeType, \
    LatticeGeometryType, SymmetryType
from glow.support.utility import are_same_shapes


# Providing the 'Path' object for the parent folder of the unittests
UNITTEST_FLDR = Path(__file__).parent


class TestTdtData(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the class
    `TdtData` that stores  the geometrical data of need from a `Lattice`
    instance.

    Attributes
    ----------
    file_name : str
        The name of the TDT file to export.
    type_geo : LatticeGeometryType
        The type of geometry, as element of the `LatticeGeometryType`
        enumeration.
    symmetry_type : SymmetryType
        The type of symmetry, as element of the `SymmetryType` enumeration.
    edges : List[Edge]
        The list of `Edge` objects for each edge of the lattice.
    boundaries : List[Boundary]
        The list of `Boundary` objects for each border of the lattice.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the class `TdtData`.
        It initializes the attributes common to all the tests.
        """
        # Declare the data used for comparison purposes
        self.file_name: str = "layout"
        self.type_geo: LatticeGeometryType = \
            LatticeGeometryType.RECTANGLE_TRAN
        self.symmetry_type: SymmetryType = SymmetryType.FULL
        ref_shape = Rectangle()
        borders = ref_shape.borders
        set_shape_name(ref_shape.face, "FACE_1")
        self.ref_face: Face = Face(ref_shape.face, 'MAT')
        edge_names_vs_faces: Dict[str, List[Any | Face]]  = {
            'EDGE_1': [borders[0], self.ref_face],
            'EDGE_2': [borders[1], self.ref_face],
            'EDGE_3': [borders[2], self.ref_face],
            'EDGE_4': [borders[3], self.ref_face]
        }
        for b, n in zip(ref_shape.borders, edge_names_vs_faces.keys()):
            set_shape_name(b, n)
        self.edges: List[Edge] = [
            Edge(edge_face[0], *edge_face[1:])
                for edge_face in edge_names_vs_faces.values()
        ]
        self.boundaries: List[Boundary] = [
            Boundary(border, self.type_geo, ref_shape.o, (1.0, 1.0))
                for border in borders]

    def test_init(self) -> None:
        """
        Method that tests the initialization of the `TdtData` class.
        """
        # Instantiate the dataclass storing the needed lattice data
        tdt = TdtData(
            filename=self.file_name + ".dat",
            edges=self.edges,
            faces=[self.ref_face],
            boundaries=self.boundaries,
            type_geo=self.type_geo,
            type_sym=self.symmetry_type,
            impressions=(1, 1),
            precisions=(1e-6, 1e-6)
        )
        # Verify the correct assignment
        self.assertEqual(tdt.filename, self.file_name + '.dat')
        self.assertEqual(tdt.type_geo, self.type_geo)
        self.assertEqual(tdt.type_sym, self.symmetry_type)
        self.assertEqual(tdt.impressions, (1, 1))
        self.assertEqual(tdt.precisions, (1e-6, 1e-6))
        self.assertEqual(tdt.nb_folds, 0)
        for te, e in zip(tdt.edges, self.edges):
            self.assertEqual(te, e)
        for tb, b in zip(tdt.boundaries, self.boundaries):
            are_same_shapes(
                tb.border, b.border, ShapeType.EDGE
            )
        for tf, f in zip(tdt.faces, [self.ref_face]):
            are_same_shapes(
                tf.face, f.face, ShapeType.FACE
            )

        # Instantiate the dataclass with a geometry type used with a TISO
        # tracking type so to check that the 'nb_folds' attribute is
        # correctly set
        tdt = TdtData(
            filename=self.file_name + ".dat",
            edges=self.edges,
            faces=[self.ref_face],
            boundaries=self.boundaries,
            type_geo=LatticeGeometryType.SYMMETRIES_TWO,
            type_sym=SymmetryType.QUARTER,
            impressions=(1, 1),
            precisions=(1e-6, 1e-6)
        )
        self.assertEqual(tdt.nb_folds, 4)

    def test_build_properties_id(self) -> None:
        """
        Method that tests the implementation of the private method
        `__build_properties_id` of the `TdtData` class.
        """
        # Instantiate the 'TdtData' class without attributes and only
        # initialize the 'faces' one
        tdt = TdtData.__new__(TdtData)
        tdt.faces = []
        # Build four 'Face' objects each with a different property name but
        # the last one
        for i in range(4):
            face = deepcopy(self.ref_face.face)
            # Change the face's name
            set_shape_name(face, f'FACE_{i+1}')
            # Build and append a 'Face' object
            if i == 3:
                mat_name = "MAT_1"
            else:
                mat_name = f'MAT_{i+1}'
            tdt.faces.append(Face(face, mat_name))

        # Build the lists containing the names of the properties and the ID
        # of the property name associated to a face
        tdt._TdtData__build_properties_id()
        # Verify the data is correctly generated
        self.assertEqual(len(tdt.properties), 3)
        self.assertEqual(tdt.properties, ['MAT_1', 'MAT_2', 'MAT_3'])
        self.assertEqual(len(tdt.property_ids), 4)
        self.assertEqual(tdt.property_ids, [1, 2, 3, 1])


class TestGeneratorFunctions(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the functions
    declared in the `generator.py` module.

    Attributes
    ----------
    file_name : str
        The path to the TDT output file.
    type_geo : LatticeGeometryType
        The type of geometry, as element of the `LatticeGeometryType`
        enumeration.
    symmetry_type : SymmetryType
        The type of symmetry, as element of the `SymmetryType` enumeration.
    edges : List[Edge]
        The list of `Edge` objects for each edge of the lattice.
    boundaries : List[Boundary]
        The list of `Boundary` objects for each border of the lattice.
    face : Face
        A `Face` object representing a region of the lattice.
    tdt : TdtData
        The `TdtData` object collecting all the characteristics of the
        geometry layout.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the functions of the
        `generator.py` module.
        It initializes the attributes common to all the tests.
        """
        # Declare the data used for comparison purposes
        self.file_name: str = str(UNITTEST_FLDR / "test_layout.dat")
        self.type_geo: LatticeGeometryType = \
            LatticeGeometryType.RECTANGLE_TRAN
        self.symmetry_type: SymmetryType = SymmetryType.FULL
        ref_shape = Rectangle()
        borders = ref_shape.borders
        set_shape_name(ref_shape.face, "FACE_1")
        self.face: Face = Face(ref_shape.face, 'MAT')
        edge_names_vs_faces: Dict[str, List[Any | Face]]  = {
            'EDGE_1': [borders[0], self.face],
            'EDGE_2': [borders[1], self.face],
            'EDGE_3': [borders[2], self.face],
            'EDGE_4': [borders[3], self.face]
        }
        for b, n in zip(ref_shape.borders, edge_names_vs_faces.keys()):
            set_shape_name(b, n)
        self.edges: List[Edge] = [
            Edge(edge_face[0], *edge_face[1:])
                for edge_face in edge_names_vs_faces.values()
        ]
        self.boundaries: List[Boundary] = [
            Boundary(border, self.type_geo, ref_shape.o, (1.0, 1.0))
                for border in borders]
        self.tdt: TdtData = TdtData(
            filename=self.file_name,
            edges=self.edges,
            faces=[self.face],
            boundaries=self.boundaries,
            type_geo=self.type_geo,
            type_sym=SymmetryType.FULL,
            impressions=(1, 1),
            precisions=(1e-6, 1e-6)
        )

    def test_write_boundary_conditions(self) -> None:
        """
        Method that tests the implementation of the function
        `_write_boundary_conditions` declared in the `generator.py`
        module.
        """
        # Test function behavior for RECTANGLE_TRAN geometry
        tdt_tran = TdtData.__new__(TdtData)
        tdt_tran.boundaries = self.boundaries
        tdt_tran.type_geo = LatticeGeometryType.RECTANGLE_TRAN

        # Declare a buffer where the text is written to
        buffer = io.StringIO()
        # Call the function and retrieve the text from the buffer
        _write_boundary_conditions(buffer, tdt_tran)
        output = buffer.getvalue()
        # Verify the boundaries section contains the expected data for the
        # type of geometry
        self.__assess_boundaries(tdt_tran, output)

        # Test function behaviour for ISOTROPIC geometry: it should return
        # after writing the 'albedo' value
        tdt_iso = TdtData.__new__(TdtData)
        tdt_iso.boundaries = self.boundaries
        tdt_iso.type_geo = LatticeGeometryType.ISOTROPIC

        buffer_iso = io.StringIO()
        _write_boundary_conditions(buffer_iso, tdt_iso)
        output_iso = buffer_iso.getvalue()
        self.assertIn(
            "* boundaries conditions: defaul nbbcda allsur", output_iso)
        self.assertIn("  0, 0, 0", output_iso)
        self.assertIn("* albedo", output_iso)
        self.assertIn("  1.0", output_iso)
        # Should not contain boundary details
        self.assertNotIn("* type  number of elements", output_iso)

        # Test an exception is raised when boundaries have an unsupported
        # boundary type
        tdt_refl = TdtData.__new__(TdtData)
        refl_boundary = self.boundaries[0]
        refl_boundary.type = BoundaryType.REFL
        tdt_refl.boundaries = [refl_boundary]
        tdt_refl.type_geo = LatticeGeometryType.RECTANGLE_TRAN

        buffer_refl = io.StringIO()
        with self.assertRaises(RuntimeError):
            _write_boundary_conditions(buffer_refl, tdt_refl)

    def test_write_edges_segment(self) -> None:
        """
        Method that tests the implementation of the function `_write_edges`
        declared in the `generator.py` module for `EdgeType.SEGMENT`-type
        edges.
        """
        # Test the output data about edges having SEGMENT type
        # ----------------------------------------------------
        # Instantiate 'TdtData' with only the needed 'edges' attribute
        tdt_data = TdtData.__new__(TdtData)
        tdt_data.edges = []
        # Build the 'Edge' objects with only the needed attributes and
        # fixed values
        for i in range(4):
            edge = Edge.__new__(Edge)
            edge.kind = EdgeType.SEGMENT
            edge.no = i + 1
            edge.data = ['SEGMENT', 1.0, 2.0, 0.0, 3.0, 4.0, 0.0]
            edge.left = self.face
            edge.right = None
            tdt_data.edges.append(edge)

        # Declare a buffer where the text is written to
        buffer = io.StringIO()
        # Call the function and retrieve the text from the buffer
        _write_edges(buffer, tdt_data)
        output = buffer.getvalue()

        # Verify the edges section contains the expected data for
        # segment-type edges
        self.__assess_segment_edges(tdt_data, output)

    def test_write_edges_circle(self) -> None:
        """
        Method that tests the implementation of the function `_write_edges`
        declared in the `generator.py` module for `EdgeType.CIRCLE`-type
        edges.
        """
        # Test the output data about edges having CIRCLE type
        # ---------------------------------------------------
        # Instantiate 'TdtData' with only the needed 'edges' attribute
        tdt_data = TdtData.__new__(TdtData)
        tdt_data.edges = []
        # Build the 'Edge' objects with only the needed attributes and
        # fixed values
        for i in range(4):
            edge = Edge.__new__(Edge)
            edge.kind = EdgeType.CIRCLE
            edge.no = i + 1
            edge.data = ['CIRCLE', 5.0, 6.0, 0.0, 0.0, 0.0, 1.0, 2.5]
            edge.right = self.face
            edge.left = None
            tdt_data.edges.append(edge)

        # Declare a buffer where the text is written to
        buffer = io.StringIO()
        # Call the function and retrieve the text from the buffer
        _write_edges(buffer, tdt_data)
        output = buffer.getvalue()

        # Check header
        self.assertIn("   elements", output)
        # Check each edge is written with correct CIRCLE data
        for edge in sorted(tdt_data.edges):
            self.assertIn(
                f"* ELEM  {edge.no}  {EDGE_NAME_VS_TYPE['CIRCLE'][1]}",
                output)
            self.assertIn(
                f" {EdgeType.CIRCLE.value}, {edge.right.no}, 0", output)
            self.assertIn(
                f"  {edge.data[1]:6f}, {edge.data[2]:6f}, " + \
                f"{edge.data[7]:6f}, 0.0",
                output)

    def test_write_edges_arc_circle(self) -> None:
        """
        Method that tests the implementation of the function `_write_edges`
        declared in the `generator.py` module for `EdgeType.ARC_CIRCLE`-type
        edges.
        """
        # Test the output data about edges having ARC_CIRCLE type
        # -------------------------------------------------------
        # Instantiate 'TdtData' with only the needed 'edges' attribute
        tdt_data = TdtData.__new__(TdtData)
        tdt_data.edges = []
        # Build the 'Edge' objects with only the needed attributes and
        # fixed values
        for i in range(4):
            edge = Edge.__new__(Edge)
            edge.kind = EdgeType.ARC_CIRCLE
            edge.no = i + 1
            edge.data = ['ARC_CIRCLE', 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 3.0,
                         4.0, 1.0, 0.0, 1.0, 4.0, 0.0]
            edge.right = self.face
            edge.left = None
            tdt_data.edges.append(edge)

        # Declare a buffer where the text is written to
        buffer = io.StringIO()
        # Call the function and retrieve the text from the buffer
        _write_edges(buffer, tdt_data)
        output = buffer.getvalue()

        # Check header
        self.assertIn("   elements", output)
        # Check each edge is written with correct ARC_CIRCLE data
        for edge in sorted(tdt_data.edges):
            self.assertIn(
                f"* ELEM  {edge.no}  {EDGE_NAME_VS_TYPE['ARC_CIRCLE'][1]}",
                output)
            self.assertIn(
                f" {EdgeType.ARC_CIRCLE.value}, {edge.right.no}, 0", output)
            self.assertIn(
                f"  {edge.data[1]:6f}, {edge.data[2]:6f}, " + \
                f"{edge.data[7]:6f}, {0.0:6f}, {90.0:6f}",
                output)

    def test_write_header(self) -> None:
        """
        Method that tests the implementation of the function `_write_header`
        declared in the `generator.py` module.
        """
        # Declare a buffer where the text is written to
        buffer = io.StringIO()
        # Call the function and retrieve the text from the buffer
        _write_header(buffer, self.tdt)
        output = buffer.getvalue()

        # Verify that the header section contains the expected values
        self.__assess_header(self.tdt, output)

    def test_write_properties(self) -> None:
        """
        Method that tests the implementation of the function
        `_write_properties` declared in the `generator.py` module.
        """
        # Instantiate 'TdtData' with only the needed attributes
        tdt_data = TdtData.__new__(TdtData)
        tdt_data.properties = ['MAT_1', 'MAT_2', 'MAT_3']
        tdt_data.property_ids = [1, 2, 3, 1, 2, 3, 1]

        # Declare a buffer where the text is written to
        buffer = io.StringIO()
        # Call the function and retrieve the text from the buffer
        _write_properties(buffer, tdt_data)
        output = buffer.getvalue()

        # Verify the presence of the properties names and indices
        self.__assess_properties(tdt_data, output)

    def test_write_regions(self) -> None:
        """
        Method that tests the implementation of the function
        `_write_regions` declared in the `generator.py` module.
        """
        # Instantiate 'TdtData' with only the needed 'faces' attribute
        tdt_data = TdtData.__new__(TdtData)
        tdt_data.faces = []
        # Build the 'Face' objects with only the 'no' attribute
        for i in range(16):
            face = Face.__new__(Face)
            face.no = i + 1
            tdt_data.faces.append(face)

        # Declare a buffer where the text is written to
        buffer = io.StringIO()
        # Call the function and retrieve the text from the buffer
        _write_regions(buffer, tdt_data)
        output = buffer.getvalue()

        # Verify the presence of the regions' indices
        self.__assess_regions(tdt_data, output)

    def test_write_tdt_file(self) -> None:
        """
        Method that tests that the function `write_tdt_file` creates the
        output TDT file with the expected content.
        """
        # Write the TDT file in the same folder of the unittests
        write_tdt_file(self.tdt)
        # Verify that the TDT file has been created
        self.assertTrue(os.path.exists(self.file_name))
        # Read the content of the file
        with open(self.file_name, "r") as f:
            content = f.read()

        # Verify the presence of header, regions, edges, boundaries, and
        # properties information
        self.__assess_boundaries(self.tdt, content)
        self.__assess_header(self.tdt, content)
        self.__assess_properties(self.tdt, content)
        self.__assess_regions(self.tdt, content)
        self.__assess_segment_edges(self.tdt, content)

        # Remove the generated TDT file
        os.remove(self.file_name)

    def __assess_boundaries(self, tdt_data: TdtData, content: str) -> None:
        """
        Method that tests that the boundaries section, provided by the input
        string, contains the expected data indicating the characteristics
        of the boundaries of the geometry layout (FULL cartesian type).

        Parameters
        ----------
        tdt_data : TdtData
            Instance of the `TdtData` dataclass to get the characteristics
            of the boundaries of the geometry layout from.
        content : str
            The content string to search for the expected data in.
        """
        self.assertIn("* boundaries conditions: defaul nbbcda allsur", content)
        self.assertIn(f"  0, {len(tdt_data.boundaries)}, 0", content)
        self.assertIn("* albedo", content)
        self.assertIn("  1.0", content) # FIXME have attribute
        # Verify that each boundary is written with the correct data
        for bc in tdt_data.boundaries:
            self.assertIn("* type  number of elements", content)
            self.assertIn(
                f"  {bc.get_bc_type_number()}, {len(bc.edge_indxs)}", content)
            self.assertIn("*   elements", content)
            for edge_no in bc.edge_indxs:
                self.assertIn(f"{edge_no}", content)
                self.assertIn("* tx, ty, angle", content)
                self.assertIn(f"{bc.tx:6f} {bc.ty:6f} {bc.angle:6f}", content)

    def __assess_header(self, tdt_data: TdtData, content: str) -> None:
        """
        Method that tests that the header section, provided by the input
        string, contains the expected values indicating the characteristics
        of the geometry layout.

        Parameters
        ----------
        tdt_data : TdtData
            Instance of the `TdtData` dataclass to get the characteristics
            of the geometry layout from.
        content : str
            The content string to search for the expected values in.
        """
        # Verify the presence of the key header lines
        self.assertIn(
            "* typge, nbfo, node, elem, macr, nreg,    z, mac2", content)
        self.assertIn("* index  kindex", content)
        self.assertIn("*     eps    eps0", content)

        # Verify that the geometry type, nb_folds, nodes, elements, regions
        # are written as expected
        expected_typegeom = tdt_data.type_geo.value
        expected_nb_folds = tdt_data.nb_folds
        expected_nbnodes = len(tdt_data.faces)
        expected_nbelements = len(tdt_data.edges)
        expected_nbregions = len(tdt_data.faces)
        self.assertIn(
            f"{expected_typegeom:5d},{expected_nb_folds:5d}, " + \
            f"{expected_nbnodes:5d}, {expected_nbelements:5d}, " +\
            f"{1:5d},{expected_nbregions:5d}, {0:5d}, {1:5d}", content)
        # Verify the values for impressions and precisions
        self.assertIn(
            f"{tdt_data.impressions[0]:5d}  {tdt_data.impressions[1]:6d}  1",
            content
        )
        self.assertIn(
            f"{tdt_data.precisions[0]:7E}   {tdt_data.precisions[1]:7E}",
            content
        )

    def __assess_properties(self, tdt_data: TdtData, content: str) -> None:
        """
        Method that verifies that the names and the indices of the
        properties are correctly present in the given content.

        Parameters
        ----------
        tdt_data : TdtData
            Instance of the `TdtData` dataclass to get the properties from.
        content : str
            The content string to search for properties in.
        """
        # Verify the name of the properties, as well as their index,
        # are present
        for id, name in enumerate(tdt_data.properties):
            self.assertIn(f"# {id+1:2d} - {name}", content)
        for id in tdt_data.property_ids:
            self.assertIn(f"  {id}", content)

    def __assess_regions(self, tdt_data: TdtData, content: str) -> None:
        """
        Method that verifies that the indices of the regions are correctly
        present in the given content.

        Parameters
        ----------
        tdt_data : TdtData
            Instance of the `TdtData` dataclass to get the indices of the
            `Face` objects from.
        content : str
            The content string to search for regions' indices in.
        """
        regions_indices = ""
        regions_no = [f.no for f in tdt_data.faces]
        for i in range(0, len(regions_no), 12):
            regions_indices += ",".join(f"{n:4d}"
                                        for n in regions_no[i:i+12]) + ",\n"
        # Remove the last comma
        regions_indices = regions_indices[0:-2] + "\n"
        self.assertIn(regions_indices, content)

    def __assess_segment_edges(self, tdt_data: TdtData, content: str) -> None:
        """
        Method that tests that the edges section, provided by the input
        string, contains the expected data indicating the characteristics
        of the segment-type edges of the geometry layout.

        Parameters
        ----------
        tdt_data : TdtData
            Instance of the `TdtData` dataclass to get the characteristics
            of the segment-type edges of the geometry layout from.
        content : str
            The content string to search for the expected data in.
        """
        # Verify the header presence
        self.assertIn("   elements", content)
        # Verify that each edge is written with correct SEGMENT data
        for edge in sorted(tdt_data.edges):
            self.assertIn(
                f"* ELEM  {edge.no}  {EDGE_NAME_VS_TYPE['SEGMENT'][1]}",
                content)
            self.assertIn(
                f" {EdgeType.SEGMENT.value}, 0, {edge.left.no}", content)
            self.assertIn(
                f"  {edge.data[1]:6f}, {edge.data[2]:6f}, " + \
                f"{(edge.data[4]-edge.data[1]):6f}, " + \
                f"{(edge.data[5]-edge.data[2]):6f}",
                content)
