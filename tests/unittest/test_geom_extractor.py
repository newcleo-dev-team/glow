"""
Module containing unittest classes to assess that the classes and functions
of the `glow.generator.geom_extractor` module have a valid implementation.
"""
import io
import unittest

from contextlib import redirect_stdout
from copy import deepcopy
from math import cos, isclose, pi, sin, sqrt
from typing import Any, Dict, List

from glow.generator.export_data import EdgeData, FaceData, build_edge_id, \
    classify_layout_edges
from glow.generator.geom_extractor import LayoutDataExtractor, analyse_layout
from glow.geometry_layouts.cells import HexCell, CartesianCell
from glow.geometry_layouts.geometries import Circle, Hexagon, Rectangle
from glow.geometry_layouts.lattices import CartesianLattice, HexLattice, Lattice
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_entities import Face, Vertex, wrap_shape
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_shape_name, limit_tolerance, make_circle, make_common, make_compound, \
    make_edge, make_face, make_partition, make_translation, make_vector, \
    make_vertex, set_shape_name
from glow.main import TdtSetup
from glow.support.types import GeometryType, LayoutGeometryType, LayoutType, \
    PropertyType, SymmetryType
from glow.support.utility import are_same_shapes, build_compound_borders, \
    build_contiguous_edges
from tests.unittest.support_funcs import build_colorset


class TestLayoutDataExtractor(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the class
    `LayoutDataExtractor` that extracts the geometrical data of need
    from a `Lattice` instance.

    Attributes
    ----------
    lattice : CartesianLattice
        A `CartesianLattice` instance made by seven cartesian cells.
    hex_symm_vs_regions: Dict[SymmetryType, Dict[GeometryType, int]]
        Providing the number of regions for a boxed lattice with seven
        hexagonal cells according to the type of symmetry and the cell's
        geometry type.
    hex_symm_vs_shape : Dict[SymmetryType, Any]
        Providing the lattice's geometry layout for each symmetry type
        (seven hexagonal cells with a box)
    rect_symm_vs_regions: Dict[SymmetryType, Dict[GeometryType, int]]
        Providing the number of regions for a boxed Cartesian lattice with
        nine Cartesian cells according to the type of symmetry and the cell's
        geometry type.
    rect_symm_vs_shape : Dict[SymmetryType, Any]
        Providing the Cartesian lattice's geometry layout for each symmetry
        type (nine Cartesian cells with a box)
    colorset : CartesianLattice
        A `CartesianLattice` instance made by assemblies, each positioned
        appropriately to replicate a colorset.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the `LayoutDataExtractor`
        class. It initialises the attributes common to all the tests.
        """
        # Build a Cartesian lattice
        cell = CartesianCell()
        cell.add(Region(Circle(radius=0.25)))
        self.lattice: CartesianLattice = CartesianLattice([cell])
        self.lattice.add_ring_of_cells(cell, 1)
        self.lattice.update_hierarchical_structure()
        # Declare the number of regions for the technological and sectorised
        # geometry layout for a hexagonal assembly at different symmetry types
        self.hex_symm_vs_regions: Dict[
            SymmetryType, Dict[GeometryType, int]] = {
            SymmetryType.FULL: {GeometryType.TECHNOLOGICAL: 20,
                                GeometryType.SECTORIZED: 55},
            SymmetryType.THIRD: {GeometryType.TECHNOLOGICAL: 10,
                                 GeometryType.SECTORIZED: 23},
            SymmetryType.SIXTH: {GeometryType.TECHNOLOGICAL: 7,
                                 GeometryType.SECTORIZED: 14},
            SymmetryType.TWELFTH: {GeometryType.TECHNOLOGICAL: 5,
                                   GeometryType.SECTORIZED: 8},
        }
        # Declare the number of regions for the technological and sectorised
        # geometry layout for a Cartesian assembly at different symmetry types
        self.rect_symm_vs_regions: Dict[
            SymmetryType, Dict[GeometryType, int]] = {
            SymmetryType.FULL: {GeometryType.TECHNOLOGICAL: 10,
                                GeometryType.SECTORIZED: 37},
            SymmetryType.HALF: {GeometryType.TECHNOLOGICAL: 7,
                                GeometryType.SECTORIZED: 19},
            SymmetryType.QUARTER: {GeometryType.TECHNOLOGICAL: 5,
                                   GeometryType.SECTORIZED: 10},
            SymmetryType.EIGHTH: {GeometryType.TECHNOLOGICAL: 4,
                                  GeometryType.SECTORIZED: 7},
        }
        # Declare the shape of the symmetry for a hexagonal assembly at
        # different symmetry types
        box_ly = 5/2 + 0.1
        box_lx = box_ly / sin(pi/3)
        self.hex_symm_vs_shape: Dict[SymmetryType, Any] = {
            SymmetryType.FULL: Hexagon(edge_length=box_lx),
            SymmetryType.THIRD: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((box_lx, 0.0, 0.0)),
                make_vertex((3/2*box_lx, box_ly, 0.0)),
                make_vertex((1/2*box_lx, box_ly, 0.0)),
            ])),
            SymmetryType.SIXTH: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((box_lx, 0.0, 0.0)),
                make_vertex((box_lx/2, box_ly, 0.0))
            ])),
            SymmetryType.TWELFTH: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((box_lx, 0.0, 0.0)),
                make_vertex((sqrt(3)/2*box_ly, box_ly/2, 0.0))
            ]))
        }
        # Declare the shape of the symmetry for a Cartesian assembly at
        # different symmetry types
        box_lx = 3.2
        box_ly = 3.2
        self.rect_symm_vs_shape: Dict[SymmetryType, Any] = {
            SymmetryType.FULL: Rectangle(
                center=(1.6, 1.6, 0.0), height=box_ly, width=box_lx),
            SymmetryType.HALF: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((1/2*box_lx, 0.0, 0.0)),
                make_vertex((1/2*box_lx, box_ly, 0.0)),
                make_vertex((0.0, box_ly, 0.0)),
            ])),
            SymmetryType.QUARTER: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((1/2*box_lx, 0.0, 0.0)),
                make_vertex((1/2*box_lx, 1/2*box_ly, 0.0)),
                make_vertex((0.0, 1/2*box_ly, 0.0))
            ])),
            SymmetryType.EIGHTH: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((1/2*box_lx, 0.0, 0.0)),
                make_vertex((1/2*box_lx, 1/2*box_ly, 0.0))
            ]))
        }
        # Build a colorset as lattice of assemblies
        self.colorset: CartesianLattice = build_colorset(cell)

    def test_init(self) -> None:
        """
        Method that tests the initialisation of the `LayoutDataExtractor`
        class for a single lattice.
        """
        # Verify an exception is raised when the layout is empty
        with self.assertRaises(RuntimeError):
            _ = LayoutDataExtractor(Lattice(), TdtSetup())

        # Instantiate the 'LayoutDataExtractor' class for a full Cartesian
        # lattice
        tdt_setup = TdtSetup()
        tdt_setup.layout_type = LayoutType.RECT
        lde = LayoutDataExtractor(self.lattice, tdt_setup)
        # Verify the correct attributes assignment
        self.assertEqual(len(lde.regions), len(self.lattice.get_regions()))

        # Check if the lattice has been properly translated by verifying that
        # the shape built from the stored borders coincides with a rectangle
        # centered in (1.5, 1.5, 0.0) and that the extracted edges coincide
        # with those of the lattice compound
        center = (1.5, 1.5, 0.0)
        self.lattice.translate(center)
        self.assertTrue(
            are_same_shapes(
                make_face(lde.borders),
                Rectangle(center=center, height=3, width=3),
                ShapeType.FACE
            )
        )
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.layout_edges),
                make_partition([self.lattice], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )
        self.assertTrue(
            all(
                isclose(xyz1, xyz2, abs_tol=1e-6)
                for xyz1, xyz2 in zip(lde.layout_centre, center)
            )
        )
        self.assertEqual(lde.dimensions, (3.0, 3.0))
        self.assertEqual(lde.boundaries, [])
        self.assertEqual(lde.subfaces, [])
        self.assertEqual(lde.edges, [])
        self.assertEqual(len(lde.id_vs_edge.values()), len(lde.layout_edges))

    def test_init_colorset(self) -> None:
        """
        Method that tests the initialisation of the `LayoutDataExtractor`
        class for a full colorset.
        """
        # Instantiate the 'LayoutDataExtractor' class for a full Cartesian
        # colorset
        tdt_setup = TdtSetup()
        tdt_setup.layout_type = LayoutType.RECT
        lde = LayoutDataExtractor(self.colorset, tdt_setup)
        # Verify the correct attributes assignment
        self.assertEqual(len(lde.regions), len(self.colorset.get_regions()))

        # Translate the colorset so that the entire layout has its lower-left
        # corner in the origin
        center = (4.8, 4.8, 0.0)
        layout_cmpd = make_translation(self.colorset, make_vector(center))
        # Check if the colorset has been properly translated by verifying that
        # the shape built from the stored borders coincides with a rectangle
        # centered in (4.8, 4.8, 0.0) and that the extracted edges
        # coincide with those of the layout compound
        self.assertTrue(
            are_same_shapes(
                make_face(lde.borders),
                Rectangle(center=center, height=9.6, width=9.6),
                ShapeType.FACE
            ),
            f"{Face(make_face(lde.borders))}"
        )
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.layout_edges),
                make_partition([layout_cmpd], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )
        # Check the other instance attribute initialisation
        self.assertTrue(
            all(
                isclose(xyz1, xyz2) for xyz1, xyz2 in zip(
                    lde.layout_centre, center)
            )
        )
        self.assertEqual(lde.dimensions, (9.6, 9.6))
        self.assertEqual(lde.boundaries, [])
        self.assertEqual(lde.subfaces, [])
        self.assertEqual(lde.edges, [])
        self.assertEqual(len(lde.id_vs_edge.values()), len(lde.layout_edges))

    def test_init_colorset_portion(self) -> None:
        """
        Method that tests the initialisation of the `LayoutDataExtractor`
        class for a portion of a colorset.
        """
        # Extract a portion of the colorset
        shape = Rectangle((2.4, 2.4, 0.0), 4.8, 4.8)
        colorset_portion = make_common(self.colorset, shape)
        # Instantiate the 'LayoutDataExtractor' class for a portion of the
        # Cartesian colorset
        tdt_setup = TdtSetup()
        tdt_setup.layout_type = LayoutType.RECT
        lde = LayoutDataExtractor(self.colorset, tdt_setup, colorset_portion)

        # Verify the correct attributes assignment
        self.assertEqual(
            len(lde.regions),
            len(extract_sub_shapes(colorset_portion, ShapeType.FACE))
        )

        # Verify that the shape built from the stored borders coincides with
        # a rectangle centered in (2.4, 2.4, 0.0) and that the extracted edges
        # coincide with those of the colorset portion
        self.assertTrue(
            are_same_shapes(make_face(lde.borders), shape, ShapeType.FACE)
        )
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.layout_edges),
                make_partition([colorset_portion], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )
        # Check the other instance attribute initialisation
        self.assertTrue(
            all(
                isclose(xyz1, xyz2, abs_tol=1e-6) for xyz1, xyz2 in zip(
                    lde.layout_centre, (0.0, 0.0, 0.0))
            )
        )
        self.assertEqual(lde.dimensions, (4.8, 4.8))
        self.assertEqual(lde.boundaries, [])
        self.assertEqual(lde.subfaces, [])
        self.assertEqual(lde.edges, [])
        self.assertEqual(len(lde.id_vs_edge.values()), len(lde.layout_edges))

    def test_build_boundaries(self) -> None:
        """
        Method that tests the implementation of the method `build_boundaries`
        of the `LayoutDataExtractor` class.
        """
        # Verify no 'BoundaryData' objects are built when the layout type is
        # 'ISOTROPIC'
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        lde.boundaries = []
        lde.build_boundaries(LayoutGeometryType.ISOTROPIC)
        self.assertEqual(len(lde.boundaries), 0)

        # Verify 'BoundaryData' objects have been created for the borders of
        # the lattice when using a valid geometry layout type
        symm_type = SymmetryType.EIGHTH
        type_geo = LayoutGeometryType.RECTANGLE_EIGHT
        self.lattice.apply_symmetry(symm_type)
        tdt_setup = TdtSetup(
            symmetry_type=symm_type, type_geo=type_geo
        )
        tdt_setup.layout_type = LayoutType.RECT
        lde = LayoutDataExtractor(self.lattice, tdt_setup)
        lde.build_boundaries(type_geo)

        self.assertEqual(len(lde.boundaries), 3)
        self.assertTrue(
            are_same_shapes(
                make_face(build_contiguous_edges([
                    make_vertex((0.0, 0.0, 0.0)),
                    make_vertex((1.5, 0.0, 0.0)),
                    make_vertex((1.5, 1.5, 0.0)),
                ])),
                make_face([b.border for b in lde.boundaries]),
                ShapeType.FACE
            )
        )

    def test_build_edges(self) -> None:
        """
        Method that tests the implementation of the method `build_edges`
        of the `LayoutDataExtractor` class.
        """
        # Build the association between edge names and 'FaceData' objects for
        # the case of a rectangular region and its edges
        shared_shape = Rectangle()
        borders = shared_shape.borders
        region = Region(shared_shape, "FACE", {PropertyType.MATERIAL: "MAT"})
        face = FaceData(region, 1, [PropertyType.MATERIAL])
        edge_names_vs_faces: Dict[str, List[Any | FaceData]]  = {
            'EDGE_1': [borders[0], face],
            'EDGE_2': [borders[1], face],
            'EDGE_3': [borders[2], face],
            'EDGE_4': [borders[3], face]
        }
        for b, n in zip(shared_shape.borders, edge_names_vs_faces.keys()):
            set_shape_name(b, n)

        # Instantiate the 'LayoutDataExtractor' class without attributes
        # and only initialise the 'edges' one
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        lde.edges = []
        lde.build_edges(edge_names_vs_faces)

        # Verify the list of 'EdgeData' objects has been built and that the
        # left face has been assigned to each edge
        self.assertEqual(len(lde.edges), 4)
        self.assertTrue(all(e.left == face for e in lde.edges))

    def test_build_edges_and_faces_association(self) -> None:
        """
        Method that tests the implementation of the method
        `build_edges_and_faces_association` of the `LayoutDataExtractor`
        class on a rectangular shape.
        It verifies the resulting data structure by checking:

        - that has four entries;
        - if the edges' names (the keys) have a corresponding name among the
          layout's edges;
        - if the first element of the values of the resulting dictionary is
          equal to any of the layout's edges;
        - if the second element of the values of the resulting dictionary is
          equal to the reference face.
        """
        # Build a reference shape that mimics the layout
        ref_shape = Rectangle()
        set_shape_name(ref_shape, "FACE_1")
        ref_edges = ref_shape.borders
        ref_id_vs_edges = classify_layout_edges(ref_edges)
        ref_face = FaceData(
            Region(ref_shape, 'MAT', {PropertyType.MATERIAL: "MAT"}),
            1,
            [PropertyType.MATERIAL]
        )
        # Instantiate the 'LayoutDataExtractor' class without attributes
        # and only initialise the needed ones
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        lde.subfaces = [ref_face]
        lde.layout_edges = ref_edges
        lde.id_vs_edge = ref_id_vs_edges
        result = lde.build_edges_and_faces_association()

        # Verify the resulting data structure
        self.assertEqual(len(result), 4)
        for edge_name, edge_faces in result.items():
            self.assertTrue(
                any(get_shape_name(e) == edge_name
                    for e in ref_id_vs_edges.values())
            )
            self.assertTrue(
                any(
                    are_same_shapes(e, edge_faces[0], ShapeType.EDGE)
                    for e in ref_id_vs_edges.values()
                )
            )
            self.assertEqual(ref_face, edge_faces[1])

    def test_build_faces(self) -> None:
        """
        Method that tests the implementation of the method `build_faces`
        of the `LayoutDataExtractor` class.
        """
        # Instantiate the 'LayoutDataExtractor' class without attributes
        # and only initialise the needed ones
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        lde.subfaces = []
        lde.regions = [Region(Rectangle(), "FACE")]

        # Verify an exception is raised if building the 'FaceData' objects
        # from lattice's regions not having any 'PropertyType', or no value
        # is assigned for the indicated type.
        properties = [PropertyType.MATERIAL]
        with self.assertRaises(RuntimeError):
            lde.build_faces(properties)
        lde.regions[0].properties = {}
        with self.assertRaises(RuntimeError):
            lde.build_faces(properties)
        lde.regions[0].properties = {PropertyType.MATERIAL: "MAT"}
        with self.assertRaises(RuntimeError):
            lde.build_faces([PropertyType.MACRO])

        # Set the regions to those of the stored lattice with both the
        # 'MATERIAL' and 'MACRO' properties
        regions = self.lattice.get_regions()
        for region in regions:
            region.properties = {
                PropertyType.MACRO: "MAC", PropertyType.MATERIAL: "MAT"
            }
        lde.regions = regions
        # Call the method with only the 'MATERIAL' property
        lde.build_faces(properties)

        # Verify the 'FaceData' objects have been correctly generated
        self.assertEqual(len(lde.subfaces), len(regions))
        for i, face in enumerate(lde.subfaces):
            self.assertEqual(face.no, i+1)
            self.assertEqual(face.property_types, properties)

    def test_get_unique_edges(self) -> None:
        """
        Method that tests the implementation of the `__get_unique_edges`
        method of the `LayoutDataExtractor` class.

        Notes
        -----
        The test case is related to two adjacent Cartesian cells: one has
        the shared border subdivided into two edges, whereas the other no.
        When it comes to storing the edges, the ones of the cell having
        two subedges are saved.
        When the method being tested is called, an edge of one of the two
        regions is provided. For the second cell, the complete edge is
        looked for, but this is not present among the stored edges.
        The present test case verifies that the two edges of the first cell
        are correctly identified as part of the one of the second cell.
        """
        # Declare the edges and the IDs dictionary to use in the test
        layout_edges = [
            make_edge(v1, v2) for v1, v2 in [
                (make_vertex((0.0, 0.0, 0.0)), make_vertex((1.0, 0.0, 0.0))),
                (make_vertex((1.0, 0.0, 0.0)), make_vertex((1.0, 0.5, 0.0))),
                (make_vertex((1.0, 0.5, 0.0)), make_vertex((1.0, 1.0, 0.0))),
                (make_vertex((1.0, 1.0, 0.0)), make_vertex((0.0, 1.0, 0.0))),
                (make_vertex((0.0, 1.0, 0.0)), make_vertex((0.0, 0.0, 0.0))),
                (make_vertex((1.0, 0.0, 0.0)), make_vertex((2.0, 0.0, 0.0))),
                (make_vertex((2.0, 0.0, 0.0)), make_vertex((2.0, 1.0, 0.0))),
                (make_vertex((2.0, 1.0, 0.0)), make_vertex((1.0, 1.0, 0.0)))
            ]
        ]
        id_vs_edge = classify_layout_edges(layout_edges)

        # Instantiate the 'LayoutDataExtractor' class without attributes
        # and only initialise the needed attributes
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        lde.layout_edges = layout_edges
        lde.id_vs_edge = id_vs_edge

        # Call the method providing an edge that is present among the stored
        # ones
        edge_to_test = list(id_vs_edge.values())[0]
        edge_id = list(id_vs_edge.keys())[0]
        found_edges = lde._get_unique_edges(edge_to_test, edge_id)
        # Verify only one edge is returned and that is present among the
        # stored ones
        self.assertEqual(len(found_edges), 1)
        self.assertTrue(
            are_same_shapes(found_edges[0], edge_to_test, ShapeType.EDGE)
        )

        # Call the method providing the edge of the second cell that is not
        # present among the stored ones
        edge_to_test = make_edge(
            make_vertex((1.0, 1.0, 0.0)), make_vertex((1.0, 0.0, 0.0))
        )
        edge_id = build_edge_id(edge_to_test)
        found_edges = lde._get_unique_edges(edge_to_test, edge_id)
        # Verify two edges are returned and that are the ones belonging to
        # the first cell
        self.assertEqual(len(found_edges), 2)
        self.assertTrue(
            are_same_shapes(found_edges[0], layout_edges[1], ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(found_edges[1], layout_edges[2], ShapeType.EDGE)
        )

        # Verify an exception is raised when providing an edge not belonging
        # to any of the two cells
        invalid_edge = make_circle(make_vertex((0.0, 0.0, 0.0)), None, 1.0)
        with self.assertRaises(RuntimeError):
            _ = lde._get_unique_edges(
                invalid_edge, build_edge_id(invalid_edge)
            )

    def test_print_log_analysis(self) -> None:
        """
        Method that tests the implementation of the `print_log_analysis`
        method of the `LayoutDataExtractor` class.
        """
        # Build the association between edge names and 'Face' objects for
        # the case of a rectangular shape and its edges
        ref_shape = Rectangle()
        set_shape_name(ref_shape, "FACE_1")
        ref_edges = ref_shape.borders
        ref_face = FaceData(
            Region(ref_shape, properties={PropertyType.MATERIAL: "MAT"}),
            1,
            [PropertyType.MATERIAL]
        )
        edge_names_vs_faces: Dict[str, List[Any | Face]]  = {
            'EDGE_1': [ref_edges[0], ref_face],
            'EDGE_2': [ref_edges[1], ref_face],
            'EDGE_3': [ref_edges[2], ref_face],
            'EDGE_4': [ref_edges[3], ref_face]
        }
        # Instantiate the 'LayoutDataExtractor' class without attributes
        # and only initialise the needed attributes
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        lde.edges = [EdgeData.__new__(EdgeData)]*4
        lde.subfaces = [ref_face]

        # Call the printer method and capture its output
        f = io.StringIO()
        with redirect_stdout(f):
            lde.print_log_analysis(edge_names_vs_faces)
        # Split the lines
        output = f.getvalue().strip().splitlines()

        # Verify the correct information is printed
        self.assertIn("1", output[0])
        self.assertIn("4", output[1])
        self.assertIn("4", output[2])
        self.assertIn("0", output[3])
        self.assertIn("0", output[4])

    def test_preprocess_cart_assembly(self) -> None:
        """
        Method that tests the implementation of the method `_preprocess` of
        the `LayoutDataExtractor` class.
        An assembly made by a Cartesian lattice is used for testing purposes.
        The following cases are tested for a `TECHNOLOGICAL` geometry type:
        - full layout;
        - a half of layout;
        - a half of layout translated so that the lower-left corner
          is already in the XYZ origin;
        - a quarter of layout whose center coincides with the XYZ origin;
        - a quarter of layout with center not placed in the XYZ origin.

        The following cases are tested for a `SECTORIZED` geometry type:
        - a eighth of layout;
        - a eighth of layout with center not placed in the XYZ origin.
        """
        # Build the Cartesian assembly
        cell = CartesianCell()
        cell.sectorize([4], [0])
        lattice = CartesianLattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        assembly = CartesianCell(
            width_height=tuple(i + 0.1*2 for i in lattice.dimensions)
        )
        assembly.add(lattice)
        assembly.update_hierarchical_structure()

        # Instantiate the 'LayoutDataExtractor' class without attributes
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        # Verify the preprocess activities are performed correctly with a
        # full symmetry
        self.__assess_preprocess_rect_symm(assembly, lde)

        # Verify the preprocess activities are performed correctly with a
        # half symmetry that must be translated so that its lower-left
        # corner is in the XYZ origin
        assembly.apply_symmetry(SymmetryType.HALF)
        self.__assess_preprocess_rect_symm(assembly, lde)
        # Verify the preprocess activities are performed correctly with a
        # half symmetry for an already translated assembly
        assembly.translate((0.0, 1/2*assembly.dimensions[1], 0.0))
        self.__assess_preprocess_rect_symm(assembly, lde)
        assembly.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # quarter symmetry (no translation needed)
        assembly.apply_symmetry(SymmetryType.QUARTER)
        self.__assess_preprocess_rect_symm(assembly, lde)
        # Verify the preprocess activities are performed correctly with a
        # quarter symmetry that has been translated elsewhere (a translation
        # is needed to have the lower-left corned in the origin)
        assembly.translate(
            (assembly.dimensions[0]/2, assembly.dimensions[1]/2, 0.0)
        )
        self.__assess_preprocess_rect_symm(assembly, lde)
        assembly.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # eighth symmetry and sectorized geometry (a translation is not
        # needed). The 'show()' method is called to trigger the construction
        # of the sectorised geometry
        assembly.apply_symmetry(SymmetryType.EIGHTH)
        sec_geo = GeometryType.SECTORIZED
        assembly.geometry_maps[sec_geo] = assembly.get_geometry_map(sec_geo)
        self.__assess_preprocess_rect_symm(assembly, lde, sec_geo)
        # Verify the preprocess activities are performed correctly with a
        # eighth symmetry that has been translated elsewhere (a translation
        # is needed to have the lower-left corned in the origin)
        assembly.translate(
            (assembly.dimensions[0]/2, assembly.dimensions[1]/2, 0.0)
        )
        self.__assess_preprocess_rect_symm(assembly, lde, sec_geo)

    def test_preprocess_colorset(self) -> None:
        """
        Method that tests the implementation of the method `_preprocess` of
        the `LayoutDataExtractor` class.
        A colorset made by lattices of Cartesian cells is used for testing
        purposes. The following cases are tested:

        - full colorset;
        - a eighth of the colorset;
        - a compound which is not part of the full colorset.
        """
        # Properly translate the comparison compound so to have its
        # lower-left corner in the origin
        center = (4.8, 4.8, 0.0)
        layout_cmpd_trnsl = make_translation(
            self.colorset, make_vector(center)
        )
        # Verify the preprocess activities are performed correctly with a
        # full colorset
        self.__assess_preprocess_colorset(
            None,
            len(self.colorset.get_regions()),
            Rectangle(center, 9.6, 9.6),
            layout_cmpd_trnsl
        )

        # Build the shape for extracting the colorset portion
        eighth_shape = make_face(
            build_contiguous_edges(
                [
                    make_vertex((0.0, 0.0, 0.0)),
                    make_vertex((4.8, 0.0, 0.0)),
                    make_vertex((4.8, 4.8, 0.0))
                ]
            )
        )
        layout_symm = make_common(self.colorset, eighth_shape)
        # Verify the preprocess activities are performed correctly with a
        # eighth symmetry of the colorset
        self.__assess_preprocess_colorset(
            layout_symm,
            len(extract_sub_shapes(layout_symm, ShapeType.FACE)),
            eighth_shape,
            layout_symm
        )

        # Verify an exception is raised when the compound does not overlap
        # to the whole colorset
        layout = make_translation(layout_symm, make_vector((10.0, 10.0, 0.0)))
        # Instantiate the 'LayoutDataExtractor' class without attributes
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        # Setup the needed attributes
        lde.regions = []
        lde.geometry_layout = deepcopy(self.colorset)
        tdt_setup = TdtSetup(GeometryType.TECHNOLOGICAL)
        tdt_setup.layout_type = LayoutType.RECT
        with self.assertRaises(RuntimeError):
            lde._preprocess(tdt_setup, layout)

    def test_preprocess_hex_assembly(self) -> None:
        """
        Method that tests the implementation of the method `_preprocess` of
        the `LayoutDataExtractor` class.

        An assembly made by a hexagonal lattice is used for testing purposes.
        The following cases are tested for a `TECHNOLOGICAL` geometry type:
        - full layout;
        - one sixth of layout;
        - one sixth of layout translated so that the lower-left corner
          is already in the XYZ origin;
        - one twelfth of layout whose center coincides with the XYZ origin;
        - one twelfth of layout with center not placed in the XYZ origin.

        The following cases are tested for a `SECTORIZED` geometry type:
        - one third of layout;
        - one third of layout with lower-left corner already in the XYZ
          origin.
        """
        # Build the hexagonal assembly
        cell = HexCell()
        cell.rotate(90)
        cell.sectorize([6], [0])
        lattice = HexLattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        assembly = HexCell(side=(lattice.dimensions[1] + 0.1)/cos(pi/6))
        assembly.add(
            Region(Hexagon(edge_length=lattice.dimensions[1]/cos(pi/6)))
        )
        assembly.add(lattice)
        assembly.update_hierarchical_structure()

        # Instantiate the 'LayoutDataExtractor' class without attributes
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        # Verify the preprocess activities are performed correctly with a
        # full symmetry
        self.__assess_preprocess_hex_symm(assembly, lde)

        # Verify the preprocess activities are performed correctly with a
        # sixth symmetry that must be translated so that its lower-left
        # corner is in the XYZ origin
        assembly.apply_symmetry(SymmetryType.SIXTH)
        self.__assess_preprocess_hex_symm(assembly, lde)
        # Verify the preprocess activities are performed correctly with a
        # sixth symmetry for an already translated lattice
        assembly.translate(
            (assembly.dimensions[0]/2, assembly.dimensions[1], 0.0)
        )
        self.__assess_preprocess_hex_symm(assembly, lde)
        assembly.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # twelfth symmetry (no translation needed)
        assembly.apply_symmetry(SymmetryType.TWELFTH)
        self.__assess_preprocess_hex_symm(assembly, lde)
        # Verify the preprocess activities are performed correctly with a
        # twelfth symmetry that has been translated elsewhere (a translation
        # is needed to have the lower-left corned in the origin)
        assembly.translate(
            (assembly.dimensions[0]/2, assembly.dimensions[1], 0.0)
        )
        self.__assess_preprocess_hex_symm(assembly, lde)
        assembly.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # third symmetry and sectorized geometry (a translation is needed)
        assembly.apply_symmetry(SymmetryType.THIRD)
        sec_geo = GeometryType.SECTORIZED
        assembly.geometry_maps[sec_geo] = assembly.get_geometry_map(sec_geo)
        self.__assess_preprocess_hex_symm(assembly, lde, sec_geo)
        # Verify the preprocess activities are performed correctly with a
        # third symmetry for an already translated assembly
        assembly.translate(
            (assembly.dimensions[0]/2, assembly.dimensions[1], 0.0)
        )
        self.__assess_preprocess_hex_symm(assembly, lde, sec_geo)

    def test_preprocess_hex_assembly_numerical_precision(self) -> None:
        """
        Method that tests the implementation of the method `_preprocess` of
        the `LayoutDataExtractor` class.

        An assembly made by a hexagonal lattice is used for testing purposes.
        This test checks the `_preprocess` method against a specific layout
        that caused issues in determining the correct borders and edges.
        The problem was related to the tolerances applied to the sub-shapes
        which resulted in the presence of edges with length equal to the
        tolerance.

        This test builds the same layout and checks that the implementation
        processes the layout correctly for a one sixth and a third of the
        full layout.
        """
        # Build the hexagonal assembly
        cell = HexCell(side=0.7852193995)
        cell.rotate(90)
        cell.sectorize([6], [0])
        lattice = HexLattice([cell])
        lattice.add_rings_of_cells(cell, 6)
        assembly = HexCell(side=lattice.dimensions[0] + 2*0.05*cos(pi/3))
        assembly.add(
            Region(
                Hexagon(edge_length=lattice.dimensions[0] + 0.05*cos(pi/3))
                - lattice.shape,
                properties={PropertyType.MATERIAL: "CLADDING"}
            )
        )
        assembly.add(lattice)
        assembly.update_hierarchical_structure()

        # Set the number of regions and the shape of the simmetries used for
        # comparison purposes (no regions for the sectorised type)
        self.hex_symm_vs_regions = {
            SymmetryType.THIRD: {
                GeometryType.TECHNOLOGICAL: 66, GeometryType.SECTORIZED: None
            },
            SymmetryType.SIXTH: {
                GeometryType.TECHNOLOGICAL: 38, GeometryType.SECTORIZED: None
            }
        }
        box_lx = lattice.dimensions[0] + 2*0.05*cos(pi/3)
        box_ly = box_lx * sin(pi/3)
        self.hex_symm_vs_shape: Dict[SymmetryType, Any] = {
            SymmetryType.THIRD: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((box_lx, 0.0, 0.0)),
                make_vertex((3/2*box_lx, box_ly, 0.0)),
                make_vertex((1/2*box_lx, box_ly, 0.0)),
            ])),
            SymmetryType.SIXTH: make_face(build_contiguous_edges([
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((box_lx, 0.0, 0.0)),
                make_vertex((box_lx/2, box_ly, 0.0))
            ]))
        }

        # Instantiate the 'LayoutDataExtractor' class without attributes
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        # Verify the preprocess activities are performed correctly with a
        # sixth symmetry that must be translated so that its lower-left
        # corner is in the XYZ origin
        assembly.apply_symmetry(SymmetryType.SIXTH)
        self.__assess_preprocess_hex_symm(assembly, lde)
        # Verify the preprocess activities are performed correctly with a
        # third symmetry (a translation is needed)
        assembly.apply_symmetry(SymmetryType.THIRD)
        self.__assess_preprocess_hex_symm(assembly, lde)

    def __assess_preprocess_colorset(
            self,
            portion: Any | None,
            no_regions: int,
            cmpr_shape: Any,
            layout_cmpd: Any
        ) -> None:
        """
        Method that assesses the correctness of the method `__preprocess`
        in the case of a colorset from which either the full or a portion
        of the colorset is considered.

        Parameters
        ----------
        portion : Any | None
            The compound object representing the portion of the colorset to
            considered. If `None`, the whole colorset is taken.
        no_regions : int
            The number of layout's regions that should be present in the
            layout to analyse.
        cmpr_shape : Any
            A shape with which comparing the face build over the saved
            borders.
        layout_cmpd : Any
            A compound object from which edges are extracted and compared with
            those stored during the preprocess activities.
        """
        # Instantiate the 'LayoutDataExtractor' class without attributes
        lde = LayoutDataExtractor.__new__(LayoutDataExtractor)
        # Setup the needed attributes
        lde.regions = []
        lde.geometry_layout = deepcopy(self.colorset)
        tdt_setup = TdtSetup(GeometryType.TECHNOLOGICAL)
        tdt_setup.layout_type = LayoutType.RECT
        # Call the method
        lde._preprocess(tdt_setup, portion)
        # Assess the correct execution of the private method
        self.__assess_preprocess(
            lde,
            no_regions,
            cmpr_shape,
            make_partition([layout_cmpd], [], ShapeType.FACE)
        )

    def __assess_preprocess_hex_symm(
            self,
            assembly: Lattice,
            lde: LayoutDataExtractor,
            geom_type: GeometryType = GeometryType.TECHNOLOGICAL
        ) -> None:
        """
        Method that assesses the correctness of the method `_preprocess` in
        the case of a hexagonal assembly with any applied symmetry.

        Parameters
        ----------
        assembly : HexCell
            The `HexCell` instance to set the corresponding attribute in the
            `LayoutDataExtractor` object.
        lde : LayoutDataExtractor
            Instance of the `LayoutDataExtractor` class to test.
        geom_type : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry for the layout.
        """
        symm_type = assembly.state.symmetry_type
        tdt_setup = TdtSetup(geom_type, symmetry_type=symm_type)
        tdt_setup.layout_type = LayoutType.HEX
        # Set the stored geometry layout
        lde.geometry_layout = deepcopy(assembly)
        # Call the private method
        lde._preprocess(tdt_setup, None)
        # Get the compound corresponding to the applied symmetry
        cmpd = (
            assembly
            if symm_type == SymmetryType.FULL
            else wrap_shape(
                limit_tolerance(assembly)
            ) * assembly.symmetry_map[symm_type]
        )
        # Get the compound + the refinment edges, if needed
        if geom_type == GeometryType.SECTORIZED:
            cmpd = cmpd / assembly.geometry_maps[GeometryType.SECTORIZED]
        # Verify the preprocess with the applied symmetry
        self.__assess_preprocess(
            lde,
            self.hex_symm_vs_regions[symm_type][geom_type],
            self.hex_symm_vs_shape[symm_type],
            cmpd
        )

    def __assess_preprocess_rect_symm(
            self,
            assembly: CartesianCell,
            lde: LayoutDataExtractor,
            geom_type: GeometryType = GeometryType.TECHNOLOGICAL
        ) -> None:
        """
        Method that assesses the correctness of the private method
        `__preprocess` in the case of a cartesian lattice with any
        applied symmetry.

        Parameters
        ----------
        assembly : CartesianCell
            The `CartesianCell` instance to set the corresponding attribute
            in the `LayoutDataExtractor` object.
        lde : LayoutDataExtractor
            Instance of the `LayoutDataExtractor` class to test.
        geom_type : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry for the lattice's cells.
        """
        symm_type = assembly.state.symmetry_type
        tdt_setup = TdtSetup(geom_type, symmetry_type=symm_type)
        tdt_setup.layout_type = LayoutType.RECT
        # Set the stored geometry layout
        lde.geometry_layout = deepcopy(assembly)
        # Call the method
        lde._preprocess(tdt_setup, None)
        # Get the compound corresponding to the applied symmetry
        cmpd = (
            assembly
            if symm_type == SymmetryType.FULL
            else assembly * assembly.symmetry_map[symm_type]
        )
        # Get the compound + the refinment edges, if needed
        if geom_type == GeometryType.SECTORIZED:
            cmpd = cmpd / assembly.geometry_maps[GeometryType.SECTORIZED]
        # Verify the preprocess with the applied symmetry
        self.__assess_preprocess(
            lde,
            self.rect_symm_vs_regions[symm_type][geom_type],
            self.rect_symm_vs_shape[symm_type],
            cmpd
        )

    def __assess_preprocess(
            self,
            lde: LayoutDataExtractor,
            no_regions: int,
            cmpr_shape: Any,
            cmpd: Any
        ) -> None:
        """
        Method that verifies the preprocess activities performed during the
        initialisation of an instance of the class `LayoutDataExtractor`.
        It is checked whether the number of layout's regions equals the
        indicated number, the borders are those of the layout's outline
        and the edges are those of the layout portion to analyse (either
        full or a symmetry part).

        Parameters
        ----------
        lde : LayoutDataExtractor
            Instance of the `LayoutDataExtractor` class to test.
        no_regions : int
            The number of layout's regions that should be present in the
            layout to analyse.
        cmpr_shape: Any
            A shape with which comparing the face build over the saved
            borders.
        cmpd : Any
            A compound object from which edges are extracted and compared
            with those stored during the preprocess activities.
        """
        self.assertEqual(len(lde.regions), no_regions)
        self.assertEqual(
            len(lde.borders), len(build_compound_borders(cmpr_shape))
        )
        self.assertTrue(
            are_same_shapes(
                make_face(lde.borders),
                cmpr_shape,
                ShapeType.FACE
            )
        )
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.layout_edges),
                make_partition(
                    extract_sub_shapes(cmpd, ShapeType.EDGE),
                    [],
                    ShapeType.EDGE
                ),
                ShapeType.COMPOUND
            )
        )


class TestGeomExtractorFunctions(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the functions
    declared in the `geom_extractor.py` module.
    """
    def test_analyse_layout(self) -> None:
        """
        Method that tests the implementation of the function `analyse_layout`
        declared in the `geom_extractor.py` module when providing a single
        lattice.
        """
        cell = CartesianCell(base_props={PropertyType.MATERIAL: "MAT2"})
        cell.add(
            Region(
                Circle(radius=0.25),
                properties={PropertyType.MATERIAL: "MAT1"}
            )
        )
        lattice = CartesianLattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        lattice.update_hierarchical_structure()

        # Call the function extracting the geometric data from the lattice
        tdt_setup = TdtSetup()
        tdt_setup.layout_type = LayoutType.RECT
        lde = analyse_layout(lattice, tdt_setup)

        # Verify the 'LayoutDataExtractor' contains the needed information
        self.assertIsInstance(lde, LayoutDataExtractor)
        self.assertEqual(len(lde.borders), 4)
        self.assertEqual(len(lde.boundaries), 0)
        self.assertEqual(len(lde.subfaces), 18)
        self.assertEqual(len(lde.edges), 24+9)
        self.assertEqual(len(lde.id_vs_edge), 24+9)
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.layout_edges),
                make_partition([lattice], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )

    def test_analyse_layout_colorset(self) -> None:
        """
        Method that tests the implementation of the function `analyse_layout`
        declared in the `geom_extractor.py` module when providing a colorset
        or its portion.
        """
        # Build the colorset and its portion
        cell = CartesianCell(base_props={PropertyType.MATERIAL: "MAT2"})
        cell.add(
            Region(
                Circle(radius=0.25),
                properties={PropertyType.MATERIAL: "MAT1"}
            )
        )
        colorset = build_colorset(cell)
        eighth_shape = make_face(
            build_contiguous_edges(
                [
                    make_vertex((0.0, 0.0, 0.0)),
                    make_vertex((4.8, 0.0, 0.0)),
                    make_vertex((4.8, 4.8, 0.0))
                ]
            )
        )
        colorset_portion = make_common(colorset, eighth_shape)

        # Call the function extracting the geometric data from the colorset
        tdt_setup = TdtSetup(type_geo=LayoutGeometryType.RECTANGLE_TRAN)
        tdt_setup.layout_type = LayoutType.RECT
        lde = analyse_layout(colorset, tdt_setup)
        # Verify the 'LayoutDataExtractor' contains the needed information
        self.assertIsInstance(lde, LayoutDataExtractor)
        self.assertEqual(len(lde.borders), 4)
        self.assertEqual(len(lde.boundaries), 4)
        self.assertEqual(len(lde.subfaces), 171)
        self.assertEqual(len(lde.edges), 321)
        self.assertEqual(len(lde.id_vs_edge), 321)
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.layout_edges),
                make_partition([colorset], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )
        # Call the function extracting the geometric data from the colorset
        # portion
        tdt_setup = TdtSetup(
            type_geo=LayoutGeometryType.RECTANGLE_EIGHT,
            symmetry_type=SymmetryType.EIGHTH
        )
        tdt_setup.layout_type = LayoutType.RECT
        lde = analyse_layout(colorset, tdt_setup, colorset_portion)
        # Verify the 'LayoutDataExtractor' contains the needed information
        self.assertIsInstance(lde, LayoutDataExtractor)
        self.assertEqual(len(lde.borders), 3)
        self.assertEqual(len(lde.boundaries), 3)
        self.assertEqual(len(lde.subfaces), 33)
        self.assertEqual(len(lde.edges), 87)
        self.assertEqual(len(lde.id_vs_edge), 87)
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.layout_edges),
                make_partition([colorset_portion], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )
