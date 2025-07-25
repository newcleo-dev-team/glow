"""
Module containing unittest classes to assess that the classes and functions
of the `glow.generator.geom_extractor` module have a valid implementation.
"""
import io
import unittest

from contextlib import redirect_stdout
from copy import deepcopy
from math import pi, sin, sqrt
from typing import Any, Dict, List

from glow.generator.geom_extractor import Edge, Face, LatticeDataExtractor, \
    analyse_lattice, build_edge_id, classify_lattice_edges
from glow.geometry_layouts.cells import HexCell, RectCell
from glow.geometry_layouts.geometries import Hexagon, Rectangle
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_shape_name, make_arc_edge, make_circle, make_compound, make_edge, \
    make_face, make_partition, make_vertex, set_shape_name
from glow.support.types import GeometryType, LatticeGeometryType, \
    PropertyType, SymmetryType
from glow.support.utility import are_same_shapes, build_contiguous_edges


class TestLatticeDataExtractor(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the class
    `LatticeDataExtractor` that extracts the geometrical data of need
    from a `Lattice` instance.

    Attributes
    ----------
    lattice : Lattice
        A `Lattice` instance made by seven cartesian cells.
    hex_symm_vs_regions: Dict[SymmetryType, Dict[GeometryType, int]]
        Providing the number of regions for a boxed lattice with seven
        hexagonal cells according to the type of symmetry and the cell's
        geometry type.
    hex_symm_vs_shape : Dict[SymmetryType, Any]
        Providing the lattice's geometry layout for each symmetry type
        (seven hexagonal cells with a box)
    rect_symm_vs_regions: Dict[SymmetryType, Dict[GeometryType, int]]
        Providing the number of regions for a boxed lattice with nine
        cartesian cells according to the type of symmetry and the cell's
        geometry type.
    rect_symm_vs_shape : Dict[SymmetryType, Any]
        Providing the lattice's geometry layout for each symmetry type
        (nine hexagonal cells with a box)
    """
    def setUp(self):
        """
        Method that sets up the test environment for the class
        `LatticeDataExtractor`.
        It initializes the attributes common to all the tests.
        """
        cell = RectCell()
        cell.add_circle(0.25)
        self.lattice: Lattice = Lattice([cell])
        self.lattice.add_ring_of_cells(cell, 1)
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
        box_ly = 5/2 + 0.1
        box_lx = box_ly / sin(pi/3)
        self.hex_symm_vs_shape: Dict[SymmetryType, Any] = {
            SymmetryType.FULL: Hexagon(edge_length=box_lx).face,
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
        box_lx = 3.2
        box_ly = 3.2
        self.rect_symm_vs_shape: Dict[SymmetryType, Any] = {
            SymmetryType.FULL: Rectangle(
                center=(1.6, 1.6, 0.0), height=box_ly, width=box_lx).face,
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

    def test_init(self) -> None:
        """
        Method that tests the initialization of the `LatticeDataExtractor`
        class.
        """
        # Verify an exception is raised when the lattice has no cells
        with self.assertRaises(RuntimeError):
            _ = LatticeDataExtractor(Lattice(), GeometryType.TECHNOLOGICAL)

        # Instantiate the 'LatticeDataExtractor' class
        lde = LatticeDataExtractor(self.lattice, GeometryType.TECHNOLOGICAL)
        # Verify the correct attributes assignment
        self.assertEqual(len(lde.lattice.lattice_cells),
                         len(self.lattice.lattice_cells))
        # Translate the lattice so that the lower-left corner is in the origin
        self.lattice.translate((1.5, 1.5, 0.0))
        self.assertIsNone(lde.lattice.lattice_box)
        self.assertTrue(
            are_same_shapes(
                make_face(lde.borders),
                Rectangle(center=(1.5, 1.5, 0.0), height=3, width=3).face,
                ShapeType.FACE
            )
        )
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.lattice_edges),
                make_partition(
                    [self.lattice.lattice_cmpd], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )
        self.assertEqual(lde.boundaries, [])
        self.assertEqual(lde.subfaces, [])
        self.assertEqual(lde.edges, [])
        self.assertEqual(len(lde.id_vs_edge.values()), len(lde.lattice_edges))

    def test_build_boundaries(self) -> None:
        """
        Method that tests the implementation of the method `build_boundaries`
        of the `LatticeDataExtractor` class.
        """
        # Verify no 'Boundary' objects are built when the lattice type is
        # 'ISOTROPIC'
        self.lattice.type_geo = LatticeGeometryType.ISOTROPIC
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        lde.lattice = deepcopy(self.lattice)
        lde.boundaries = []
        lde.build_boundaries()
        self.assertEqual(len(lde.boundaries), 0)

        # Verify 'Boundary' objects have been created for the borders of the
        # lattice
        self.lattice.apply_symmetry(SymmetryType.EIGHTH)
        self.lattice.type_geo = LatticeGeometryType.RECTANGLE_EIGHT
        lde = LatticeDataExtractor(self.lattice, GeometryType.TECHNOLOGICAL)
        lde.build_boundaries()
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
        of the `LatticeDataExtractor` class.
        """
        # Build the association between edge names and 'Face' objects for
        # the case of a rectangular shape and its edges
        shared_shape = Rectangle()
        borders = shared_shape.borders
        set_shape_name(shared_shape.face, "FACE_1")
        face = Face(shared_shape.face, 'MAT')
        edge_names_vs_faces: Dict[str, List[Any | Face]]  = {
            'EDGE_1': [borders[0], face],
            'EDGE_2': [borders[1], face],
            'EDGE_3': [borders[2], face],
            'EDGE_4': [borders[3], face]
        }
        for b, n in zip(shared_shape.borders, edge_names_vs_faces.keys()):
            set_shape_name(b, n)

        # Instantiate the 'LatticeDataExtractor' class without attributes
        # and only initialize the 'edges' one
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        lde.edges = []
        lde.build_edges(edge_names_vs_faces)

        # Verify the list of 'Edge' objects has been built
        self.assertEqual(len(lde.edges), 4)
        self.assertTrue(
            all(
                are_same_shapes(
                    e.left.face, face.face, ShapeType.FACE) for e in lde.edges
            )
        )

    def test_build_edges_and_faces_association(self) -> None:
        """
        Method that tests the implementation of the method
        `build_edges_and_faces_association` of the
        `LatticeDataExtractor` class.
        """
        # Build a reference shape that mimics the lattice
        ref_shape = Rectangle()
        set_shape_name(ref_shape.face, "FACE_1")
        ref_edges = ref_shape.borders
        ref_id_vs_edges = classify_lattice_edges(ref_edges)
        ref_face = Face(ref_shape.face, 'MAT')
        # Instantiate the 'LatticeDataExtractor' class without attributes
        # and only initialize the needed attributes
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        lde.subfaces = [ref_face]
        lde.lattice_edges = ref_edges
        lde.id_vs_edge = ref_id_vs_edges
        result = lde.build_edges_and_faces_association()

        # Verify the resulting data structure by checking:
        # 1. that has four entries;
        # 2. if the edges' names (the keys) have a corresponding name among
        #    the lattice's edges;
        # 3. if the first element of the values of the resulting dictionary
        #    is equal to any of the lattice's edges;
        # 4. if the second element of the values of the resulting dictionary
        #    is equal to the reference face.
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
            self.assertTrue(
                are_same_shapes(
                    ref_face.face, edge_faces[1].face, ShapeType.FACE)
            )

    def test_build_faces(self) -> None:
        """
        Method that tests the implementation of the method `build_faces`
        of the `LatticeDataExtractor` class.
        """
        # Instantiate the 'LatticeDataExtractor' class without attributes
        # and only initialize the needed attributes
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        lde.lattice = deepcopy(self.lattice)
        lde.subfaces = []

        # Verify an exception is raised if building the 'Face' objects from
        # lattice's regions not having any 'PropertyType', or no value is
        # assigned for the indicated type.
        # N.B. Since only the 'PropertyType.MATERIAL' is present, it is not
        # possible to test the case where calling the method with a different
        # 'PropertyType' element.
        lde.lattice.build_regions(GeometryType.TECHNOLOGICAL)
        with self.assertRaises(RuntimeError):
            lde.build_faces(PropertyType.MATERIAL)
        for region in lde.lattice.regions:
            region.properties = {PropertyType.MATERIAL: ""}
        with self.assertRaises(RuntimeError):
            lde.build_faces(PropertyType.MATERIAL)

        # Assign the values for the 'PropertyType.MATERIAL' and call
        # the method
        for region in lde.lattice.regions:
            region.properties = {PropertyType.MATERIAL: "MAT"}
        lde.build_faces(PropertyType.MATERIAL)

        # Verify the 'Face' objects have been correctly generated
        self.assertEqual(len(lde.subfaces), len(lde.lattice.regions))
        for i, face in enumerate(lde.subfaces):
            self.assertEqual(face.no, i+1)

    def test_get_unique_edges(self) -> None:
        """
        Method that tests the implementation of the private method
        `__get_unique_edges` of the `LatticeDataExtractor` class.

        Notes
        -----
        The test case is related to two adjacent cartesian cells: one has
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
        lattice_edges = [
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
        id_vs_edge = classify_lattice_edges(lattice_edges)

        # Instantiate the 'LatticeDataExtractor' class without attributes
        # and only initialize the needed attributes
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        lde.lattice_edges = lattice_edges
        lde.id_vs_edge = id_vs_edge

        # Call the method providing an edge that is present among the stored
        # ones
        edge_to_test = list(id_vs_edge.values())[0]
        edge_id = list(id_vs_edge.keys())[0]
        found_edges = lde._LatticeDataExtractor__get_unique_edges(
            edge_to_test, edge_id)
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
        found_edges = lde._LatticeDataExtractor__get_unique_edges(
            edge_to_test, edge_id)
        # Verify two edges are returned and that are the ones belonging to
        # the first cell
        self.assertEqual(len(found_edges), 2)
        self.assertTrue(
            are_same_shapes(found_edges[0], lattice_edges[1], ShapeType.EDGE)
        )
        self.assertTrue(
            are_same_shapes(found_edges[1], lattice_edges[2], ShapeType.EDGE)
        )

        # Verify an exception is raised when providing an edge not belonging
        # to any of the two cells
        invalid_edge = make_circle(make_vertex((0.0, 0.0, 0.0)), None, 1.0)
        with self.assertRaises(RuntimeError):
            _ = lde._LatticeDataExtractor__get_unique_edges(
                invalid_edge, build_edge_id(invalid_edge))

    def test_print_log_analysis(self) -> None:
        """
        Method that tests the implementation of the method
        `print_log_analysis` of the `LatticeDataExtractor` class.
        """
        # Build the association between edge names and 'Face' objects for
        # the case of a rectangular shape and its edges
        ref_shape = Rectangle()
        set_shape_name(ref_shape.face, "FACE_1")
        ref_edges = ref_shape.borders
        ref_face = Face(ref_shape.face, 'MAT')
        edge_names_vs_faces: Dict[str, List[Any | Face]]  = {
            'EDGE_1': [ref_edges[0], ref_face],
            'EDGE_2': [ref_edges[1], ref_face],
            'EDGE_3': [ref_edges[2], ref_face],
            'EDGE_4': [ref_edges[3], ref_face]
        }
        # Instantiate the 'LatticeDataExtractor' class without attributes
        # and only initialize the needed attributes
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        lde.edges = [Edge.__new__(Edge)]*4
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

    def test_preprocess_cart_lattice(self) -> None:
        """
        Method that tests the implementation of the private method
        `__preprocess` of the `LatticeDataExtractor` class.
        A lattice made by cartesian cells is used for testing purposes.
        The following cases are tested for a `TECHNOLOGICAL` geometry type:
        - full lattice;
        - a half of lattice;
        - a half of lattice translated so that the lower-left corner
          is already in the XYZ origin;
        - a quarter of lattice whose center coincides with the XYZ origin;
        - a quarter of lattice with center not placed in the XYZ origin.

        The following cases are tested for a `SECTORIZED` geometry type:
        - a eighth of lattice;
        - a eighth of lattice with center not placed in the XYZ origin.
        """
        # Instantiate a Lattice with cartesian cells
        cell = RectCell()
        cell.sectorize([4], [0])
        lattice = Lattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        lattice.build_lattice_box([0.1])

        # Instantiate the 'LatticeDataExtractor' class without attributes
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        # Verify the preprocess activities are performed correctly with a
        # full symmetry
        self.__assess_preprocess_rect_symm(lattice, lde)

        # Verify the preprocess activities are performed correctly with a
        # half symmetry that must be translated so that its lower-left
        # corner is in the XYZ origin
        lattice.apply_symmetry(SymmetryType.HALF)
        self.__assess_preprocess_rect_symm(lattice, lde)
        # Verify the preprocess activities are performed correctly with a
        # half symmetry for an already translated lattice
        lattice.translate((0.0, 1/2*lattice.ly, 0.0))
        self.__assess_preprocess_rect_symm(lattice, lde)
        lattice.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # quarter symmetry (no translation needed)
        lattice.apply_symmetry(SymmetryType.QUARTER)
        self.__assess_preprocess_rect_symm(lattice, lde)
        # Verify the preprocess activities are performed correctly with a
        # quarter symmetry that has been translated elsewhere (a translation
        # is needed to have the lower-left corned in the origin)
        lattice.translate((lattice.lx/2, lattice.ly/2, 0.0))
        self.__assess_preprocess_rect_symm(lattice, lde)
        lattice.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # eighth symmetry and sectorized geometry (a translation is not
        # needed)
        lattice.apply_symmetry(SymmetryType.EIGHTH)
        lattice.build_regions(GeometryType.SECTORIZED)
        self.__assess_preprocess_rect_symm(
            lattice, lde, GeometryType.SECTORIZED)
        # Verify the preprocess activities are performed correctly with a
        # eighth symmetry that has been translated elsewhere (a translation
        # is needed to have the lower-left corned in the origin)
        lattice.translate((lattice.lx/2, lattice.ly/2, 0.0))
        self.__assess_preprocess_rect_symm(
            lattice, lde, GeometryType.SECTORIZED)

    def test_preprocess_hex_lattice(self) -> None:
        """
        Method that tests the implementation of the private method
        `__preprocess` of the `LatticeDataExtractor` class.
        A lattice made by hexagonal cells is used for testing purposes.
        The following cases are tested for a `TECHNOLOGICAL` geometry type:
        - full lattice;
        - one sixth of lattice;
        - one sixth of lattice translated so that the lower-left corner
          is already in the XYZ origin;
        - one twelfth of lattice whose center coincides with the XYZ origin;
        - one twelfth of lattice with center not placed in the XYZ origin.

        The following cases are tested for a `SECTORIZED` geometry type:
        - one third of lattice;
        - one third of lattice with lower-left corner already in the XYZ
          origin.
        """
        # Instantiate a Lattice with hexagonal cells
        cell = HexCell()
        cell.rotate(90)
        cell.sectorize([6], [0])
        lattice = Lattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        lattice.build_lattice_box([0.1])

        # Instantiate the 'LatticeDataExtractor' class without attributes
        lde = LatticeDataExtractor.__new__(LatticeDataExtractor)
        # Verify the preprocess activities are performed correctly with a
        # full symmetry
        self.__assess_preprocess_hex_symm(lattice, lde)

        # Verify the preprocess activities are performed correctly with a
        # sixth symmetry that must be translated so that its lower-left
        # corner is in the XYZ origin
        lattice.apply_symmetry(SymmetryType.SIXTH)
        self.__assess_preprocess_hex_symm(lattice, lde)
        # Verify the preprocess activities are performed correctly with a
        # sixth symmetry for an already translated lattice
        lattice.translate((lattice.lx/2, lattice.ly, 0.0))
        self.__assess_preprocess_hex_symm(lattice, lde)
        lattice.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # twelfth symmetry (no translation needed)
        lattice.apply_symmetry(SymmetryType.TWELFTH)
        self.__assess_preprocess_hex_symm(lattice, lde)
        # Verify the preprocess activities are performed correctly with a
        # twelfth symmetry that has been translated elsewhere (a translation
        # is needed to have the lower-left corned in the origin)
        lattice.translate((lattice.lx/2, lattice.ly, 0.0))
        self.__assess_preprocess_hex_symm(lattice, lde)
        lattice.translate((0.0, 0.0, 0.0))

        # Verify the preprocess activities are performed correctly with a
        # third symmetry and sectorized geometry (a translation is needed)
        lattice.apply_symmetry(SymmetryType.THIRD)
        lattice.build_regions(GeometryType.SECTORIZED)
        self.__assess_preprocess_hex_symm(
            lattice, lde, GeometryType.SECTORIZED)
        # Verify the preprocess activities are performed correctly with a
        # third symmetry for an already translated lattice
        lattice.translate((lattice.lx/2, lattice.ly, 0.0))
        self.__assess_preprocess_hex_symm(
            lattice, lde, GeometryType.SECTORIZED)

    def __assess_preprocess_hex_symm(
            self,
            lattice: Lattice,
            lde: LatticeDataExtractor,
            geom_type: GeometryType = GeometryType.TECHNOLOGICAL) -> None:
        """
        Method that assesses the correctness of the private method
        `__preprocess` in the case of a hexagonal lattice with any
        applied symmetry.

        Parameters
        ----------
        lattice : Lattice
            The `Lattice` instance to set the corresponding attribute in the
            `LatticeDataExtractor` object.
        lde : LatticeDataExtractor
            Instance of the `LatticeDataExtractor` class to test.
        geom_type : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry for the lattice's cells.
        """
        lde.lattice = deepcopy(lattice)
        # Call the private method
        lde._LatticeDataExtractor__preprocess(geom_type)
        cmpd = (lattice.lattice_cmpd
                if lattice.symmetry_type == SymmetryType.FULL
                else lattice.lattice_symm)
        self.__assess_preprocess(
            lde,
            self.hex_symm_vs_regions[lattice.symmetry_type][geom_type],
            self.hex_symm_vs_shape[lattice.symmetry_type],
            cmpd
        )

    def __assess_preprocess_rect_symm(
            self,
            lattice: Lattice,
            lde: LatticeDataExtractor,
            geom_type: GeometryType = GeometryType.TECHNOLOGICAL) -> None:
        """
        Method that assesses the correctness of the private method
        `__preprocess` in the case of a cartesian lattice with any
        applied symmetry.

        Parameters
        ----------
        lattice : Lattice
            The `Lattice` instance to set the corresponding attribute in the
            `LatticeDataExtractor` object.
        lde : LatticeDataExtractor
            Instance of the `LatticeDataExtractor` class to test.
        geom_type : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry for the lattice's cells.
        """
        lde.lattice = deepcopy(lattice)
        # Call the private method
        lde._LatticeDataExtractor__preprocess(geom_type)
        cmpd = (lattice.lattice_cmpd
                if lattice.symmetry_type == SymmetryType.FULL
                else lattice.lattice_symm)
        self.__assess_preprocess(
            lde,
            self.rect_symm_vs_regions[lattice.symmetry_type][geom_type],
            self.rect_symm_vs_shape[lattice.symmetry_type],
            cmpd
        )

    def __assess_preprocess(
            self,
            lde: LatticeDataExtractor,
            no_regions: int,
            cmpr_shape: Any,
            cmpd: Any) -> None:
        """
        Method that verifies the preprocess activities performed during the
        initialization of an instanc of the class`LatticeDataExtractor`.
        It is checked whether the number of lattice's regions equals the
        indicated number, the borders are those  of the lattice's outline
        and the edges are those of the lattice portion to analyse (either
        full or a symmetry part).

        Parameters
        ----------
        lde : LatticeDataExtractor
            Instance of the `LatticeDataExtractor` class to test.
        no_regions : int
            The number of lattice's regions that should be present in the
            lattice to analyse.
        cmpr_shape: Any
            A shape with which comparing the face build over the saved
            borders.
        cmpd : Any
            A compound object from which edges are extracted and compared
            with those stored during the preprocess activities.
        """
        self.assertEqual(len(lde.lattice.regions), no_regions)
        self.assertTrue(
            are_same_shapes(
                make_face(lde.borders),
                cmpr_shape,
                ShapeType.FACE
            ),
        )
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.lattice_edges),
                make_compound(extract_sub_shapes(cmpd, ShapeType.EDGE)),
                ShapeType.COMPOUND
            )
        )


class TestGeomExtractorFunctions(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the functions
    declared in the `geom_extractor.py` module.
    """
    def test_analyse_lattice(self) -> None:
        """
        Method that tests the implementation of the function `analyse_lattice`
        declared in the `geom_extractor.py` module.
        """
        cell = RectCell()
        cell.add_circle(0.25)
        cell.set_properties({PropertyType.MATERIAL: ['MAT1', 'MAT2']})
        lattice: Lattice = Lattice([cell])
        lattice.add_ring_of_cells(cell, 1)
        lattice.build_regions(GeometryType.TECHNOLOGICAL)

        # Call the function extracting the geometric data from the lattice
        lde = analyse_lattice(
            lattice, GeometryType.SECTORIZED, PropertyType.MATERIAL)

        # Verify the 'LatticeDataExtractor' contains the needed information
        self.assertIsInstance(lde, LatticeDataExtractor)
        self.assertEqual(len(lde.borders), 4)
        self.assertEqual(len(lde.boundaries), 0)
        self.assertEqual(len(lde.subfaces), 18)
        self.assertEqual(len(lde.edges), 24+9)
        self.assertEqual(len(lde.id_vs_edge), 24+9)
        self.assertTrue(
            are_same_shapes(
                make_compound(lde.lattice_edges),
                make_partition(
                    [lattice.lattice_cmpd], [], ShapeType.EDGE),
                ShapeType.COMPOUND
            )
        )

    def test_build_edge_id(self) -> None:
        """
        Method that tests the implementation of the function `build_edge_id`
        declared in the `geom_extractor.py` module.
        """
        # Declare the edges to test, one for each type
        sgmnt = make_edge(
            make_vertex((0.0, 0.0, 0.0)), make_vertex((1.0, 0.0, 0.0))
        )
        arc_crcl = make_arc_edge(
            make_vertex((0.0, 1.0, 0.0)),
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((1.0, 1.0, 0.0))
        )
        circle = make_circle(make_vertex((0.0, 0.0, 0.0)), None, 1.0)
        # Strings to compare the resulting ID with
        ref_sgmnt_id = "EDGE_SEGMENT_0_0_0_1_0_0"
        ref_arc_crcl_id = "EDGE_ARC_CIRCLE_0_1_0_0_0_1_1_0_0_0_1_1_0"
        ref_circle_id = "EDGE_CIRCLE_0_0_0_0_0_1_1"

        # Get the edges' information and verify its correctness
        sgmnt_id = build_edge_id(sgmnt)
        arc_crcl_id = build_edge_id(arc_crcl)
        circle_id = build_edge_id(circle)
        self.assertIn(ref_sgmnt_id, sgmnt_id)
        self.assertIn(ref_arc_crcl_id, arc_crcl_id)
        self.assertIn(ref_circle_id, circle_id)

        # Verify an exception is raised when providing a shape that is not
        # an edge
        with self.assertRaises(RuntimeError):
            _ = build_edge_id(make_face([circle]))

    def test_classify_lattice_edges(self) -> None:
        """
        Method that tests the implementation of the function
        `classify_lattice_edges` declared in the `geom_extractor.py`
        module.
        """
        # Declare the edges to use in the test, one for each type
        sgmnt = make_edge(
            make_vertex((0.0, 0.0, 0.0)), make_vertex((1.0, 0.0, 0.0))
        )
        arc_crcl = make_arc_edge(
            make_vertex((0.0, 1.0, 0.0)),
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((1.0, 1.0, 0.0))
        )
        circle = make_circle(make_vertex((0.0, 0.0, 0.0)), None, 1.0)
        # Reference dictionary of strings VS edges with which compare the
        # dictionary resulting from the function call
        ref_ids_edges = {
            "EDGE_SEGMENT_0_0_0_1_0_0": sgmnt,
            "EDGE_ARC_CIRCLE_0_1_0_0_0_1_1_0_0_0_1_1_0": arc_crcl,
            "EDGE_CIRCLE_0_0_0_0_0_1_1": circle
        }

        # Build the classificaton of the edges
        ids_edges = classify_lattice_edges([sgmnt, arc_crcl, circle])

        # Verify the correct classification information is present and that
        # the edges' names are set
        for ref_id_edge, id_edge in zip(ref_ids_edges, ids_edges):
            self.assertEqual(ref_id_edge, id_edge)
            self.assertTrue(
                are_same_shapes(
                    ref_ids_edges[ref_id_edge],
                    ids_edges[id_edge],
                    ShapeType.EDGE
                )
            )

        # Verify an exception is raised if any of the elements in the list
        # is not among the allowed edges' types
        with self.assertRaises(RuntimeError):
            _ = classify_lattice_edges([make_face(circle)])
