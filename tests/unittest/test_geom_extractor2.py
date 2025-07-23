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
    classify_lattice_edges
from glow.geometry_layouts.cells import HexCell, RectCell
from glow.geometry_layouts.geometries import Hexagon, Rectangle
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_shape_name, make_compound, make_face, make_vertex, set_shape_name
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
                make_compound(extract_sub_shapes(self.lattice.lattice_cmpd,
                                                 ShapeType.EDGE)),
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
        ref_id_vs_edges = classify_lattice_edges(ref_edges)
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

