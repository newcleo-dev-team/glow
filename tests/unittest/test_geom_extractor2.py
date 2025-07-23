"""
Module containing unittest classes to assess that the classes and functions
of the `glow.generator.geom_extractor` module have a valid implementation.
"""
from copy import deepcopy
from math import pi, sin, sqrt
from typing import Any, Dict
import unittest

from glow.generator.geom_extractor import LatticeDataExtractor
from glow.geometry_layouts.cells import HexCell, RectCell
from glow.geometry_layouts.geometries import Hexagon, Rectangle
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_shape_type, make_compound, make_cut, make_face, make_vertex
from glow.support.types import GeometryType, SymmetryType
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

