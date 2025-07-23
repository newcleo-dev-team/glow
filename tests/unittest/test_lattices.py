"""
Module containing unittest classes to assess that the classes and functions
of the `glow.geometry_layouts.lattice` module have a valid implementation.
"""
from copy import deepcopy
import math
import unittest

from typing import List

from glow.support.types import *
from glow.geometry_layouts.cells import Cell, HexCell, RectCell, get_region_info
from glow.geometry_layouts.geometries import Hexagon, Rectangle, build_hexagon
from glow.geometry_layouts.lattices import Lattice, get_changed_cells, \
    get_compound_from_geometry
from glow.support.utility import are_same_shapes, build_compound_borders
from glow.interface.geom_interface import *
from support_funcs import build_lattice_ref_vectors, set_up_hex_cells, \
    set_up_rect_cells


class TestLattice(unittest.TestCase):
    """
    Test case for verifying the geometric operations and visualization
    capabilities of the `Lattice` class.

    This test suite provides common setup and a set of tests to ensure that
    implementations of `Lattice` correctly handle the operations performed
    on a lattice, whose type is either cartesian or hexagonal.

    Attributes
    ----------
    name : str
        The name of the lattice when displayed in the SALOME viewer.
    hex_cells : List[HexCell]
        The list of `HexCell` objects the hexagonal lattice to test is made
        of.
    rect_cells : List[RectCell]
        The list of `RectCell` objects the cartesian lattice to test is made
        of.
    o : Any
        A vertex object positioned at the origin of the XYZ space.
    lattice : Lattice
        The `Lattice` object to test, made either by cartesian or hexagonal
        cells.
    box_layers : List[float]
        List storing the thickness for each layer the lattice's box is made of
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `Lattice` class.
        It initializes the attributes common to all the tests.
        """
        self.name: str = "Lattice"
        self.lattice: Lattice | None = Lattice(name=self.name)
        self.o = make_vertex((0.0, 0.0, 0.0))
        hex_cell = HexCell(name="Hexagonal Cell")
        rect_cell = RectCell(name="Cartesian Cell")
        self.hex_cells: List[HexCell] = set_up_hex_cells(hex_cell)
        self.rect_cells: List[RectCell] = set_up_rect_cells(rect_cell)
        self.box_layers: List[float] = [0.075, 0.075]

    def test_init(self) -> None:
        """
        Method that tests the correct initialization of the `Lattice` class.
        """
        # No cells case
        self.__assess_lattice_init([], [], 0)
        # No cells case with center different from XYZ origin
        center = (1.0, 1.0, 0.0)
        self.lattice = Lattice(center=center)
        self.__assess_lattice_init([], [], 0, center)
        # Hexagonal cells case
        self.lattice = Lattice(self.hex_cells)
        self.__assess_lattice_init(self.hex_cells, [], 1)
        # Cartesian case (odd number of cells)
        self.rect_cells = set_up_rect_cells(self.rect_cells[0])
        self.lattice = Lattice(self.rect_cells)
        self.__assess_lattice_init(self.rect_cells, [], 1)
        # Cartesian case (even number of cells)
        self.rect_cells = set_up_rect_cells(self.rect_cells[0], True)
        self.lattice = Lattice(self.rect_cells)
        self.__assess_lattice_init(self.rect_cells, [], 1)

        # Tests the instantiation of a lattice with hexagonal cells and the
        # thicknesses of its box layers
        self.lattice = Lattice(self.hex_cells, boxes_thick=self.box_layers)
        # Build the comparison figure
        box_face = self.__build_hex_box()
        # Verify the correctness of the box surface
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_box.face,
                            box_face,
                            ShapeType.COMPOUND)
        )

    def test_add_cell(self) -> None:
        """
        Method that tests the correct implementation of the method
        `add_cell` of the `Lattice` class.
        """
        # Instantiate the lattice without any cell
        self.lattice = Lattice()
        # Add an hexagonal cell
        cell = self.hex_cells[0]
        position = ()
        self.__assess_add_cell(
            cell, position, LatticeGeometryType.HEXAGON_TRAN)
        position = (2*cell.apothem, 0.0, 0.0)
        self.__assess_add_cell(
            cell, position, LatticeGeometryType.ISOTROPIC, 1, 2)
        # Verify correct exception handling
        self.__assess_add_cell(self.rect_cells[0], position)
        rotated_cell = deepcopy(cell)
        rotated_cell.rotate(30)
        self.__assess_add_cell(rotated_cell, position)

    def test_add_ring_of_cells(self) -> None:
        """
        Method that tests the correct implementation of the method
        `add_ring_of_cells` of the `Lattice` class.
        """
        # Instantiate the lattice with only a central hexagonal cell
        cell = self.hex_cells[0]
        self.lattice = Lattice([cell])
        # Add 1 ring of cells
        self.__assess_add_ring_of_cells(
            cell,
            ring_no=1,
            cells_no=6,
            ring_dist=[(6, 2*cell.figure.ly)])
        self.__assess_add_ring_of_cells(
            cell,
            ring_no=2,
            cells_no=12,
            ring_dist=[
                (6, 2*cell.figure.ly),
                (6, 4*cell.figure.ly),
                (6, 3*cell.figure.lx)
            ])

    def test_add_rings_of_cells(self) -> None:
        """
        Method that tests the correct implementation of the method
        `add_rings_of_cells` of the `Lattice` class.
        """
        # Instantiate the lattice with only a central hexagonal cell
        cell = self.hex_cells[0]
        self.lattice = Lattice([cell])
        # Add 2 rings of cells at once
        self.__assess_add_rings_of_cells(
            cell,
            ring_no=2,
            cells_no=18,
            ring_dist=[
                (6, 2*cell.figure.ly),
                (6, 4*cell.figure.ly),
                (6, 3*cell.figure.lx)
            ]
        )
        # Test a cell-centered lattice with cartesian-type cells
        cell = self.rect_cells[0]
        self.lattice = Lattice([cell])
        self.__assess_add_rings_of_cells(
            cell,
            ring_no=2,
            cells_no=24,
            ring_dist=[
                (4, cell.figure.ly),
                (4, math.sqrt(cell.figure.lx*cell.figure.lx +
                              cell.figure.ly*cell.figure.ly)),
                (4, 2*cell.figure.lx),
                (4, math.sqrt(4*cell.figure.lx*cell.figure.lx +
                              4*cell.figure.ly*cell.figure.ly)),
                (8, math.sqrt(4*cell.figure.lx*cell.figure.lx +
                              cell.figure.ly*cell.figure.ly))
            ]
        )
        # Test a lattice without any central cell
        cell = self.rect_cells[0]
        self.lattice = Lattice([])
        self.__assess_add_rings_of_cells(
            cell,
            ring_no=2,
            cells_no=16,
            ring_dist=[
                (4, math.sqrt(1/4*cell.figure.lx*cell.figure.lx +
                              1/4*cell.figure.ly*cell.figure.ly)),
                (4, math.sqrt(9/4*cell.figure.lx*cell.figure.lx +
                              9/4*cell.figure.ly*cell.figure.ly)),
                (8, math.sqrt(9/4*cell.figure.lx*cell.figure.lx +
                              1/4*cell.figure.ly*cell.figure.ly))
            ]
        )

    def test_apply_symmetry(self) -> None:
        """
        Method that tests the correct implementation of the method
        `apply_symmetry` of the `Lattice` class.
        """
        # Instantiate the lattice with a list of hexagonal cells
        self.lattice = Lattice(self.hex_cells)
        # Verify the lattice compound has been assembled
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_cmpd,
                            make_compound([c.face for c in self.hex_cells]),
                            ShapeType.COMPOUND)
        )
        # Verify an exception is raised when applying a symmetry to the
        # lattice of hexagonal cells without enclosing it in a box
        with self.assertRaises(AssertionError):
            self.lattice.apply_symmetry(SymmetryType.SIXTH)
        # Enclose the lattice in a box
        self.lattice.build_lattice_box([0.075])
        # Build the vertices of the shape of the symmetry
        symm_vs_vertices = {
            SymmetryType.SIXTH: [
                self.lattice.lattice_center,
                make_vertex((-self.lattice.lx/2, -self.lattice.ly, 0.0)),
                make_vertex((self.lattice.lx/2, -self.lattice.ly, 0.0))],
            SymmetryType.THIRD: [
                self.lattice.lattice_center,
                make_vertex((-self.lattice.lx/2, -self.lattice.ly, 0.0)),
                make_vertex((self.lattice.lx/2, -self.lattice.ly, 0.0)),
                make_vertex((self.lattice.lx, 0.0, 0.0))],
            SymmetryType.TWELFTH: [
                self.lattice.lattice_center,
                make_vertex((self.lattice.lx, 0.0, 0.0)),
                make_vertex((3/4*self.lattice.lx,
                             self.lattice.lx*math.sqrt(3)/4,
                             0.0))
            ]
        }
        # Verify the correctness of the symmetry application
        self.__assess_symmetry(SymmetryType.SIXTH,
                               symm_vs_vertices[SymmetryType.SIXTH])
        self.__assess_symmetry(SymmetryType.THIRD,
                               symm_vs_vertices[SymmetryType.THIRD])
        self.__assess_symmetry(SymmetryType.TWELFTH,
                               symm_vs_vertices[SymmetryType.TWELFTH])
        self.__assess_symmetry(
            SymmetryType.FULL,
            extract_sub_shapes(
                make_face(
                    build_compound_borders(self.lattice.lattice_box.face)),
                ShapeType.VERTEX))

        # Test the application of symmetries for a cartesian lattice with
        # and without a central cell
        self.__assess_rect_symmetry()
        self.__assess_rect_symmetry(True)

    def test_build_lattice_box(self) -> None:
        """
        Method that tests the correct implementation of the method
        `build_lattice_box` of the `Lattice` class.
        """
        # Instantiate the lattice with the hexagonal cells
        self.lattice = Lattice(self.hex_cells)
        # Build the comparison figure
        box_face = self.__build_hex_box()
        # Build the lattice box
        self.lattice.build_lattice_box(self.box_layers)
        # Verify the correctness of the box surface
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_box.face,
                            box_face,
                            ShapeType.COMPOUND)
        )
        # Instantiate the lattice with the cartesian cells
        self.lattice = Lattice(self.rect_cells)
        # Build the comparison figure
        box_face = self.__build_rect_box()
        # Build the lattice box
        self.lattice.build_lattice_box(self.box_layers)
        # Verify the correctness of the box surface
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_box.face,
                            box_face,
                            ShapeType.COMPOUND)
        )
        # Instantiate the lattice with the cartesian cells (even cells number)
        self.lattice = Lattice(
            set_up_rect_cells(self.rect_cells[0], True))
        # Build the comparison figure
        box_face = self.__build_rect_box()
        # Build the lattice box
        self.lattice.build_lattice_box(self.box_layers)
        # Verify the correctness of the box surface
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_box.face,
                            box_face,
                            ShapeType.COMPOUND)
        )

    def test_build_regions(self) -> None:
        """
        Method that tests the implementation of the method `build_regions`
        of the `Lattice` class.
        """
        # Instantiate the lattice without any cell, then add hexagonal cells
        self.lattice = Lattice()
        cell = self.hex_cells[0]
        self.lattice.add_cell(cell, ())
        self.lattice.add_ring_of_cells(cell, 1)
        # Verify no regions are present after adding the cells
        self.assertEqual(len(self.lattice.regions), 0)
        # Build the regions, using either the cells' technological geometry
        # or the sectorized one and verify the regions construction
        self.lattice.build_regions(GeometryType.TECHNOLOGICAL)
        self.__assess_build_regions(GeometryType.TECHNOLOGICAL)
        self.lattice.build_regions(GeometryType.SECTORIZED)
        self.__assess_build_regions(GeometryType.SECTORIZED)
        self.assertFalse(self.lattice.is_update_needed)

    def test_get_regions_info(self) -> None:
        """
        Method that tests the implementation of the method `get_regions_info`
        of the `Lattice` class.
        """
        # Instantiate the lattice with a list of hexagonal cells and
        # a box and build the regions
        self.lattice = Lattice(self.hex_cells)
        self.lattice.build_regions()

        # Verify an exception is raised when calling the method without
        # selecting any region
        with self.assertRaises(RuntimeError):
            self.lattice.get_regions_info()
        # Test the inner implementation by verifying an exception is raised
        # when providing a region not beloging to the lattice
        with self.assertRaises(RuntimeError):
            get_region_info(make_face([make_circle(self.o, None, 0.1)]),
                            self.lattice.regions)
        region_face = self.hex_cells[0].face
        # Check the message returns the name of the region and the associated
        # properties, if any are assigned
        self.__assess_get_regions_info(region_face,
                                       "No associated properties.")

        # Apply values for the PropertyType.MATERIAL to lattice's regions
        for layer in self.lattice.layers:
            for cell in layer:
                cell.set_properties({PropertyType.MATERIAL: ['MAT']})
        # Rebuilt the regions to apply the modifications
        self.lattice.build_regions()
        # Check the correct assignment
        self.__assess_get_regions_info(
            region_face, f"{PropertyType.MATERIAL.name}: MAT")

    def test_restore_cells(self) -> None:
        """
        Method that tests the implementation of the method `restore_cells`
        of the `Lattice` class.
        """
        # Set up a cartesian cell with circles and properties
        rect_cell = self.rect_cells[0]
        radii = [0.2, 0.3, 0.4]
        properties = {PropertyType.MATERIAL: ['MAT1', 'MAT2', 'MAT3', 'MAT4']}
        for r in radii:
            rect_cell.add_circle(r)
        rect_cell.set_properties(properties)
        # Instantiate the lattice without any cell
        self.lattice = Lattice()
        self.lattice.add_rings_of_cells(rect_cell, 2)
        # Add a cell so that it overlaps the four in the upper-right section
        self.lattice.add_cell(rect_cell, (1.0, 1.0, 0.0))
        self.lattice.build_regions()

        # Restore only the cells cut by the overlapping one
        cut_cells = [
            self.lattice.lattice_cells[3],
            self.lattice.lattice_cells[10],
            self.lattice.lattice_cells[11],
            self.lattice.lattice_cells[12],
        ]
        self.lattice.restore_cells(cut_cells, {PropertyType.MATERIAL: 'MAT1'})
        # Verify the cut cells have been restored
        self.__assess_cells_state(
            cut_cells, 0, {PropertyType.MATERIAL: ['MAT1']})
        # Verify no restore is performed if the cells have not been cut
        no_cut_cells = [self.lattice.lattice_cells[1],
                        self.lattice.lattice_cells[2]]
        self.lattice.restore_cells(
            no_cut_cells,
            {PropertyType.MATERIAL: 'MAT1'}
        )
        self.__assess_cells_state(no_cut_cells, 3, properties)
        # Verify cells not cut are restored, as indicated explicitly
        self.lattice.restore_cells(
            no_cut_cells,
            {PropertyType.MATERIAL: 'MAT1'},
            False
        )
        self.__assess_cells_state(
            no_cut_cells, 0, {PropertyType.MATERIAL: ['MAT1']})
        # Verify the lattice needs to be updated
        self.assertTrue(self.lattice.is_update_needed)

    def test_rotate(self) -> None:
        """
        Method that tests the implementation of the method `rotate` of the
        `Lattice` class.
        """
        # Instantiate the lattice with a list of hexagonal cells
        for cell in self.hex_cells:
            cell.sectorize([6], [0])
        self.lattice = Lattice(self.hex_cells)
        self.lattice.build_lattice_box([0.075])
        self.lattice.apply_symmetry(SymmetryType.SIXTH)
        # Store the initial rotation angle
        rotation_0 = self.lattice.cells_rot
        rot_angle = 90.0
        # Build reference vectors for each lattice element before and
        # after the rotation
        ref_vectors1 = build_lattice_ref_vectors(self.lattice)
        self.lattice.rotate(rot_angle)
        ref_vectors2 = build_lattice_ref_vectors(self.lattice)
        # Check the correct rotation happened by assessing that the angle
        # between the reference vectors is the rotation angle
        for v1, v2 in zip(ref_vectors1, ref_vectors2):
            self.assertTrue(
                math.isclose(
                    get_angle_between_shapes(v1, v2),
                    rot_angle),
                f"{get_angle_between_shapes(v1, v2)}, {rot_angle}"
            )
        # Check the 'rotation' attribute has updated
        self.assertTrue(
            math.isclose(
                self.lattice.cells_rot, rotation_0 + rot_angle)
        )

    def test_set_lattice_box_properties(self) -> None:
        """
        Method that tests the implementation of the method
        `set_lattice_box_properties` of the `Lattice` class.
        """
        # Instantiate the lattice with a list of hexagonal cells and
        # a box
        self.lattice = Lattice(self.hex_cells, boxes_thick=self.box_layers)

        # Verify exceptions are raised when incorrectly setting the box
        # properties (no properties or number of elements not matching
        # the box regions)
        with self.assertRaises(RuntimeError):
            self.lattice.set_lattice_box_properties({})
        with self.assertRaises(RuntimeError):
            self.lattice.set_lattice_box_properties(
                {PropertyType.MATERIAL: ['MAT']})

        # Assign the MATERIAL property to the regions of the box
        materials = {PropertyType.MATERIAL: ['MAT1', 'MAT2', 'MAT3']}
        self.lattice.set_lattice_box_properties(materials)
        # Verify the correct assignment of the materials to the regions
        # according to their distance from the lattice's center
        regions = sorted(
            extract_sub_shapes(self.lattice.lattice_box.face, ShapeType.FACE),
            key=lambda subface: get_min_distance(self.lattice.lattice_center,
                                                 subface)
        )
        for i, (r1, r2) in enumerate(
                zip(regions, self.lattice.lattice_box.tech_geom_props)):
            if are_same_shapes(r1, r2, ShapeType.FACE):
                self.assertEqual(
                    self.lattice.lattice_box.tech_geom_props[r2][
                        PropertyType.MATERIAL],
                    materials[PropertyType.MATERIAL][i]
                )

    def test_set_region_property(self) -> None:
        """
        Method that tests the implementation of the method
        `set_region_property` of the `Lattice` class.
        """
        # Instantiate the lattice with a list of hexagonal cells
        self.lattice = Lattice(self.hex_cells)
        property_value = (PropertyType.MATERIAL, "MAT")
        # Verify an exception is raised when calling the method without
        # selecting any region, with providing a non-existent one or when
        # the `SECTORIZED` geometry is displayed
        with self.assertRaises(RuntimeError):
            self.lattice.set_region_property(*property_value)
        with self.assertRaises(RuntimeError):
            self.lattice.set_region_property(
                *property_value,
                make_face([make_circle(self.o, None, 0.1)])
            )
        with self.assertRaises(RuntimeError):
            # Mock the displayed geometry is the sectorized one
            self.lattice.displayed_geom = GeometryType.SECTORIZED
            self.lattice.set_region_property(*property_value)

        # Call the method providing the region to update
        self.lattice.displayed_geom = GeometryType.TECHNOLOGICAL
        region_face = extract_sub_shapes(
            make_compound([self.lattice.lattice_cells[-1].face]),
            ShapeType.FACE)[0]
        region_pnt = make_vertex_inside_face(region_face)
        self.lattice.set_region_property(*property_value, region_face)
        # Check the correct assignment
        for layer in self.lattice.layers:
            for cell in layer:
                if not is_point_inside_shape(region_pnt, cell.face):
                    continue
                for region in cell.tech_geom_props:
                    if are_same_shapes(region, region_face, ShapeType.FACE):
                        self.assertTrue(
                            property_value[0] in cell.tech_geom_props[region]
                        )
                        self.assertEqual(
                            cell.tech_geom_props[region][property_value[0]],
                            property_value[1]
                        )
                        break
        self.assertTrue(self.lattice.is_update_needed)

    def test_set_type_geo(self) -> None:
        """
        Method that tests the implementation of the method `set_type_geo`
        of the `Lattice` class.
        """
        # Instantiate the lattice with cartesian cells
        self.lattice = Lattice(self.rect_cells)
        # Verify exceptions are raised if trying to set an invalid 'type_geo'
        # TRAN when the lattice is not enclosed in a box, any other type
        # indicating a symmetry when no symmetry is applied
        with self.assertRaises(RuntimeError):
            self.lattice.set_type_geo(LatticeGeometryType.RECTANGLE_TRAN)
        with self.assertRaises(RuntimeError):
            self.lattice.set_type_geo(LatticeGeometryType.RECTANGLE_SYM)
        # Apply a symmetry
        self.lattice.apply_symmetry(SymmetryType.EIGHTH)
        typ_geo = self.lattice.type_geo
        new_type_geo = LatticeGeometryType.SYMMETRIES_TWO
        self.lattice.set_type_geo(new_type_geo)
        # Verify the type_geo has changed
        self.assertTrue(
            typ_geo != self.lattice.type_geo and
            self.lattice.type_geo == new_type_geo
        )
        # Instantiate the lattice with hexagonal cells
        self.lattice = Lattice(self.hex_cells)
        # Verify exceptions are raised if trying to set an invalid 'type_geo'
        # TRAN when the lattice is not enclosed in a box, any other type
        # indicating a symmetry when no symmetry is applied
        with self.assertRaises(RuntimeError):
            self.lattice.set_type_geo(LatticeGeometryType.HEXAGON_TRAN)
        with self.assertRaises(RuntimeError):
            self.lattice.set_type_geo(LatticeGeometryType.R120)
        # Apply a symmetry
        self.lattice.build_lattice_box(self.box_layers)
        self.lattice.apply_symmetry(SymmetryType.SIXTH)
        typ_geo = self.lattice.type_geo
        new_type_geo = LatticeGeometryType.RA60
        self.lattice.set_type_geo(new_type_geo)
        # Verify the type_geo has changed
        self.assertTrue(
            typ_geo != self.lattice.type_geo and
            self.lattice.type_geo == new_type_geo
        )

    def test_show(self) -> None:
        """
        Method that tests the implementation of the method `show` of the
        `Lattice` class.
        """
        # Instantiate the lattice with a list of hexagonal cells
        for cell in self.hex_cells:
            cell.sectorize([6], [0])
        self.lattice = Lattice(self.hex_cells)
        # Check the cell's face has not been displayed yet
        self.assertIsNone(self.lattice.lattice_entry_id)
        # Verify lattice's regions are correctly shown
        self.__assess_show(None, GeometryType.TECHNOLOGICAL)
        self.__assess_show(None, GeometryType.SECTORIZED)
        # Verify an exception is raised when trying to show regions according
        # to a property type without assigning values to regions
        with self.assertRaises(RuntimeError):
            self.__assess_show(PropertyType.MATERIAL,
                               GeometryType.TECHNOLOGICAL)
        # Apply values for the PropertyType.MATERIAL to cells' regions
        for layer in self.lattice.layers:
            for cell in layer:
                cell.set_properties({PropertyType.MATERIAL: ['MAT1']})
        # Verify lattice's regions are shown according to the property values
        self.__assess_show(PropertyType.MATERIAL, GeometryType.TECHNOLOGICAL)
        self.__assess_show(PropertyType.MATERIAL, GeometryType.SECTORIZED)

    def test_translate(self) -> None:
        """
        Method that tests the implementation of the method `translate` of the
        `Lattice` class.
        """
        # Instantiate the lattice with a list of hexagonal cells
        for cell in self.hex_cells:
            cell.sectorize([6], [0])
        self.lattice = Lattice(self.hex_cells)
        self.lattice.build_lattice_box([0.075])
        self.lattice.apply_symmetry(SymmetryType.SIXTH)
        # Store the elements for assessing the correct translation
        new_cntr = (10.0, 10.0, 0.0)
        new_cntr_vrtx = make_vertex(new_cntr)
        cdg_pre = self.lattice.lattice_center
        distance = get_min_distance(self.lattice.lattice_center,
                                    new_cntr_vrtx)
        lattice_org = deepcopy(self.lattice)

        # Translate the lattice
        self.lattice.translate(new_cntr)
        # Check the correct translation happened
        self.assertTrue(
            are_same_shapes(
                self.lattice.lattice_center, new_cntr_vrtx, ShapeType.VERTEX)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(
                    make_cdg(lattice_org.lattice_cmpd),
                    make_cdg(self.lattice.lattice_cmpd)),
                distance)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(
                    make_cdg(lattice_org.lattice_box.face),
                    make_cdg(self.lattice.lattice_box.face)),
                distance)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(
                    make_cdg(lattice_org.lattice_symm),
                    make_cdg(self.lattice.lattice_symm)),
                distance)
        )
        self.assertTrue(
            math.isclose(
                get_min_distance(cdg_pre, self.lattice.lattice_center),
                distance)
        )
        for c1, c2 in zip(lattice_org.lattice_cells,
                          self.lattice.lattice_cells):
            self.assertTrue(
                math.isclose(
                    get_min_distance(c1.figure.o, c2.figure.o), distance)
            )
            self.assertTrue(
                math.isclose(
                    get_min_distance(
                        make_cdg(c1.sectorized_face),
                        make_cdg(c2.sectorized_face)),
                    distance)
            )
        for lyr1, lyr2 in zip(lattice_org.layers, self.lattice.layers):
            for c1, c2 in zip(lyr1, lyr2):
                self.assertTrue(
                    math.isclose(
                        get_min_distance(c1.figure.o, c2.figure.o), distance)
                )
        self.assertFalse(self.lattice.is_update_needed)

    def __assess_add_cell(
            self,
            cell: Cell,
            position: Tuple[float, float, float],
            type_geo: LatticeGeometryType = LatticeGeometryType.ISOTROPIC,
            ring_no: int = 0,
            cells_no: int = 1
        ) -> None:
        """
        Method that assesses the correct addition of a `Cell` object to the
        `Lattice` instance.

        Parameters
        ----------
        cell : Cell
            The cell to be added to the lattice.
        position : Tuple[float, float, float]
            The XYZ coordinates of the point where the cell is added to the
            lattice.
        type_geo : LatticeGeometryType = LatticeGeometryType.ISOTROPIC
            The lattice's type of geometry.
        ring_no : int = 0
            The number of lattice's rings.
        cells_no : int = 1
            The number of lattice's cells.
        """
        # Verify exceptions are raised when trying to add an invalid cell
        if (self.lattice.cells_type is not None and cell.cell_type !=
            self.lattice.cells_type):
            with self.assertRaises(RuntimeError):
                self.lattice.add_cell(cell, position)
            return
        if (self.lattice.cells_rot is not None and
            math.degrees(cell.rotation) != self.lattice.cells_rot):
            with self.assertRaises(RuntimeError):
                self.lattice.add_cell(cell, position)
            return
        # Add the valid cell to the lattice
        self.lattice.add_cell(cell, position)
        if not position:
            position = (0.0, 0.0, 0.0)
        self.assertEqual(self.lattice.cells_type, cell.cell_type)
        self.assertTrue(
            math.isclose(self.lattice.cells_rot, math.degrees(cell.rotation))
        )
        # Check the cell has been added in the correct position
        self.assertEqual(len(self.lattice.lattice_cells), cells_no)
        self.assertEqual(len(self.lattice.layers[cells_no]), 1)
        self.assertTrue(
            are_same_shapes(self.lattice.lattice_cells[0].face,
                            cell.face,
                            ShapeType.FACE)
        )
        self.assertTrue(
            get_min_distance(self.lattice.lattice_cells[cells_no-1].figure.o,
                             make_vertex(position)) < 1e-5
        )
        self.assertTrue(self.lattice.is_update_needed)
        self.assertEqual(self.lattice.type_geo, type_geo)
        self.assertEqual(self.lattice.rings_no, ring_no)

    def __assess_add_ring_of_cells(
            self,
            cell: Cell,
            type_geo: LatticeGeometryType = LatticeGeometryType.ISOTROPIC,
            ring_no: int = 1,
            cells_no: int = 1,
            ring_dist: List[Tuple[int, float]] = []
        ) -> None:
        """
        Method that assesses the correct addition of a ring of `Cell` objects
        to the `Lattice` instance.

        Parameters
        ----------
        cell : Cell
            The cell to be added to the lattice.
        type_geo : LatticeGeometryType = LatticeGeometryType.ISOTROPIC
            The lattice's type of geometry.
        ring_no : int = 1
            The number of lattice's rings.
        cells_no : int = 1
            The number of lattice's cells.
        ring_dist : List[Tuple[int, float]] = []
            Collecting the number of cells placed at a specific distance
            from the lattice center with the distance values themselves.
        """
        # Verify exceptions are raised when trying to add an invalid cell
        if (self.lattice.cells_type is not None and cell.cell_type !=
            self.lattice.cells_type):
            with self.assertRaises(RuntimeError):
                self.lattice.add_ring_of_cells(cell, ring_no)
            return
        if (self.lattice.cells_rot is not None and
            math.degrees(cell.rotation) != self.lattice.cells_rot):
            with self.assertRaises(RuntimeError):
                self.lattice.add_ring_of_cells(cell, ring_no)
            return
        if ring_no == 0:
            with self.assertRaises(RuntimeError):
                self.lattice.add_ring_of_cells(cell, ring_no)
            return
        # Store initial data
        n0 = len(self.lattice.lattice_cells)
        i0 = 1 if len(self.lattice.lattice_cells) > 0 else 0
        # Add a valid ring of cells to the lattice
        self.lattice.add_ring_of_cells(cell, ring_no)
        # Verify the correct addition of cells
        self.assertEqual(self.lattice.cells_type, cell.cell_type)
        self.assertTrue(
            math.isclose(self.lattice.cells_rot, math.degrees(cell.rotation))
        )
        # Check the cells have been added in the correct positions
        self.assertEqual(len(self.lattice.lattice_cells) - n0, cells_no)
        self.assertEqual(len(self.lattice.layers[-1]), cells_no)
        self.assertTrue(
            all(
                are_same_shapes(
                    self.lattice.lattice_cells[i].face,
                    cell.translate(
                        get_point_coordinates(
                            self.lattice.lattice_cells[i].figure.o)).face,
                    ShapeType.COMPOUND
                ) for i in range(n0, len(self.lattice.lattice_cells))
            )
        )
        self.assertTrue(self.lattice.is_update_needed)
        self.assertEqual(self.lattice.type_geo, type_geo)
        self.assertEqual(self.lattice.rings_no, ring_no)
        # Verify that the correct number of cells is placed at the right
        # distance from the lattice center
        self.__assess_cells_distances(i0, ring_dist)

    def __assess_add_rings_of_cells(
            self,
            cell: Cell,
            type_geo: LatticeGeometryType = LatticeGeometryType.ISOTROPIC,
            ring_no: int = 1,
            cells_no: int = 1,
            ring_dist: List[Tuple[int, float]] = []
        ) -> None:
        """
        Method that assesses the correct addition of several rings of `Cell`
        objects to the `Lattice` instance.

        Parameters
        ----------
        cell : Cell
            The cell to be added to the lattice.
        type_geo : LatticeGeometryType = LatticeGeometryType.ISOTROPIC
            The lattice's type of geometry.
        ring_no : int = 1
            The number of lattice's rings to add.
        cells_no : int = 1
            The number of lattice's cells.
        ring_dist : List[Tuple[int, float]] = []
            Collecting the number of cells placed at a specific distance
            from the lattice center with the distance value itself.
        """
        # Verify exceptions are raised when trying to add an invalid cell
        if (self.lattice.cells_type is not None and cell.cell_type !=
            self.lattice.cells_type):
            with self.assertRaises(RuntimeError):
                self.lattice.add_rings_of_cells(cell, ring_no)
            return
        if (self.lattice.cells_rot is not None and
            math.degrees(cell.rotation) != self.lattice.cells_rot):
            with self.assertRaises(RuntimeError):
                self.lattice.add_rings_of_cells(cell, ring_no)
            return
        if ring_no == 0:
            with self.assertRaises(RuntimeError):
                self.lattice.add_rings_of_cells(cell, ring_no)
            return
        # Store initial data
        n0 = len(self.lattice.lattice_cells)
        ring_0 = self.lattice.rings_no
        i0 = 1 if len(self.lattice.lattice_cells) > 0 else 0
        # Add a valid ring of cells to the lattice
        self.lattice.add_rings_of_cells(cell, ring_no)
        # Verify the correct addition of cells
        self.assertEqual(self.lattice.cells_type, cell.cell_type)
        self.assertTrue(
            math.isclose(self.lattice.cells_rot, math.degrees(cell.rotation))
        )
        # Check the cells have been added in the correct positions
        self.assertEqual(len(self.lattice.lattice_cells) - n0, cells_no)
        self.assertEqual(len(self.lattice.layers[-1]), cells_no)
        self.assertTrue(
            all(
                are_same_shapes(
                    self.lattice.lattice_cells[i].face,
                    cell.translate(
                        get_point_coordinates(
                            self.lattice.lattice_cells[i].figure.o)).face,
                    ShapeType.COMPOUND
                ) for i in range(n0, len(self.lattice.lattice_cells))
            )
        )
        self.assertTrue(self.lattice.is_update_needed)
        self.assertEqual(self.lattice.type_geo, type_geo)
        self.assertEqual(self.lattice.rings_no - ring_0, ring_no)
        # Verify that the correct number of cells is placed at the right
        # distance from the lattice center
        self.__assess_cells_distances(i0, ring_dist)

    def __assess_cells_distances(
            self,
            i0: int,
            ring_dist: List[Tuple[int, float]] = []) -> None:
        """
        Method that checks whether the cells have been added to the lattice
        correctly. This is done by grouping all the cells sharing the same
        distance from the lattice center and comparing the result with the
        given metric.

        Parameters
        ----------
        i0 : int
            Index indicating the starting element in the list of cells.
        ring_dist : List[Tuple[int, float]] = []
            Collecting the number of cells at a specific distance from the
            lattice center and the distance values themselves.
        """
        # Collect the distances and count the number of cells for each
        # distance
        distances = [
            round(
                get_min_distance(
                    cell.figure.o,
                    self.lattice.lattice_center),
                6) for cell in self.lattice.lattice_cells[i0:]
        ]
        unique_values = set(distances)
        result = []
        for value in unique_values:
            count = sum(
                1 for d in distances if math.isclose(d, value, abs_tol=1e-5))
            result.append((count, value))
        # Sort both collections in ascending order by the distance value
        result.sort(key=lambda d: d[1])
        ring_dist.sort(key=lambda d: d[1])
        # Check there is the correct number of cells for each distance
        for t1, t2 in zip(ring_dist, result):
            self.assertTrue(
                t1[0] == t2[0] and math.isclose(t1[1], t2[1], abs_tol=1e-5),
                f"{ring_dist}, {result}"
            )

    def __assess_cells_state(
            self,
            cells: List[Cell],
            no_circles: int,
            properties: Dict[PropertyType, List[str]]) -> None:
        """
        Method that verifies if all the provided cells have the indicated
        number of circles and properties.

        Parameters
        ----------
        cells : List[Cell]
            List of `Cell` objects to check.
        no_circles : int
            The number of circles all the given cells must have.
        properties : Dict[PropertyType, List[str]]
            The dictionary of properties all the given cells must have.
        """
        for cell in cells:
            self.assertEqual(len(cell.inner_circles), no_circles)
            for properties in cell.tech_geom_props.values():
                for p, v in properties.items():
                    self.assertEqual(v, properties[p])

    def __assess_rect_symmetry(self, is_even: bool = False) -> None:
        """
        Method that assesses the correct application of all the symmetries
        available for a lattice made by cartesian cells, with and without
        a central cell.

        Parameters
        ----------
        is_even : bool = False
            Flag indicating whether the cartesian lattice to test has an
            even number of cells (no central cell).
        """
        # Instantiate the lattice with a list of cartesian cells
        cells = set_up_rect_cells(self.rect_cells[0], is_even)
        self.lattice = Lattice(cells)
        # Verify the lattice compound has been assembled
        self.assertTrue(
            are_same_shapes(
                self.lattice.lattice_cmpd,
                make_compound([c.face for c in cells]),
                ShapeType.COMPOUND
            )
        )
        # Enclose the lattice in a box
        self.lattice.build_lattice_box([0.075])
        # Build the vertices of the shape of the symmetry
        symm_vs_vertices = {
            SymmetryType.HALF: [
                make_vertex((0.0, -self.lattice.ly/2, 0.0)),
                make_vertex((self.lattice.lx/2, -self.lattice.ly/2, 0.0)),
                make_vertex((self.lattice.lx/2, self.lattice.ly/2, 0.0)),
                make_vertex((0.0, self.lattice.ly/2, 0.0))],
            SymmetryType.QUARTER: [
                self.lattice.lattice_center,
                make_vertex((self.lattice.lx/2, 0.0, 0.0)),
                make_vertex((self.lattice.lx/2, self.lattice.ly/2, 0.0)),
                make_vertex((0.0, self.lattice.ly/2, 0.0))],
            SymmetryType.EIGHTH: [
                self.lattice.lattice_center,
                make_vertex((self.lattice.lx/2, 0.0, 0.0)),
                make_vertex((self.lattice.lx/2, self.lattice.ly/2, 0.0))]
        }
        # Verify the correctness of the symmetry application
        self.__assess_symmetry(SymmetryType.HALF,
                               symm_vs_vertices[SymmetryType.HALF])
        self.__assess_symmetry(SymmetryType.QUARTER,
                               symm_vs_vertices[SymmetryType.QUARTER])
        self.__assess_symmetry(SymmetryType.EIGHTH,
                               symm_vs_vertices[SymmetryType.EIGHTH])
        self.__assess_symmetry(
            SymmetryType.FULL,
            extract_sub_shapes(
                make_face(
                    build_compound_borders(self.lattice.lattice_box.face)),
                ShapeType.VERTEX))

    def __assess_show(
            self,
            property: PropertyType | None,
            geometry: GeometryType) -> None:
        """
        Method that assesses the correct display of lattice's regions
        according to whether regions should be colored by property values
        and to the geometry type, either `TECHNOLOGICAL` or `SECTORIZED`.

        Parameters
        ----------
        property : PropertyType | None
            The `PropertyType` according to which regions should be colored;
            if `None`, a default color must be used.
        geometry : GeometryType
            The `GeometryType` indicating whether regions are built from the
            technological or the sectorized geometry.
        """
        # Display the cell's technological geometry in the SALOME viewer
        self.lattice.show(property, geometry)
        # Verify an entry ID has been associated to the lattice's face
        self.assertIsNotNone(self.lattice.lattice_entry_id)
        # Verify the regions corresponding to the geometry have been built
        self.assertTrue(len(self.lattice.regions) > 0)
        self.__assess_build_regions(geometry)
        # Verify if regions have the same default color or the ones having
        # the same property value share the same color
        if property is None:
            default_color = (167, 167, 167)
            self.assertTrue(
                all(region.color == default_color
                    for region in self.lattice.regions)
            )
        else:
            prop_val_vs_color = {}
            for region in self.lattice.regions:
                value = region.properties[property]
                if value in prop_val_vs_color:
                    self.assertEqual(region.color, prop_val_vs_color[value])
                else:
                    prop_val_vs_color[value] = region.color
        # Verify the regions have been displayed by checking their entry ID
        # has been defined so that they are children of the lattice's face
        self.assertTrue(
            all(
                region.face_entry_id and region.face_entry_id.startswith(
                    self.lattice.lattice_entry_id + ':')
                for region in self.lattice.regions),
            f"{[region.face_entry_id for region in self.lattice.regions]}"
        )
        # Verify the attribute indicating the displayed geometry
        self.assertEqual(self.lattice.displayed_geom, geometry)

    def __assess_build_regions(self, geometry: GeometryType) -> None:
        """
        Method that assesses that the correct regions have been built
        according to the given geometry.

        Parameters
        ----------
        geometry : GeometryType
            The type of geometry of the cells to build regions from.
        """
        if geometry == GeometryType.TECHNOLOGICAL:
            self.assertEqual(
                len(self.lattice.regions),
                sum(len(cell.tech_geom_props.keys())
                    for cell in self.lattice.lattice_cells)
            )
            for region in self.lattice.regions:
                found = False
                for cell in self.lattice.lattice_cells:
                    for tr in cell.tech_geom_props:
                        if are_same_shapes(region.face, tr, ShapeType.FACE):
                            found = True
                            break
                    else:
                        continue
                    self.assertTrue(found)
                    break
        elif geometry == GeometryType.SECTORIZED:
            for region in self.lattice.regions:
                found = False
                for cell in self.lattice.lattice_cells:
                    for tr in extract_sub_shapes(
                        cell.sectorized_face, ShapeType.FACE):
                        if are_same_shapes(region.face, tr, ShapeType.FACE):
                            found = True
                            break
                    else:
                        continue
                    self.assertTrue(found)
                    break

    def __assess_get_regions_info(
            self, region_face: Any, str_to_find: str) -> None:
        """
        Method that assesses whether the right message is shown when
        displaying information about a specified region of the lattice.

        Parameters
        ----------
        region_face : Any
            The lattice's region to retrieve information from.
        str_to_find: str
            The informative string to find for in the one being produced.
        """
        # Retrieve the informative message about a region
        mssg = get_region_info(region_face, self.lattice.regions)
        for region in self.lattice.regions:
            if are_same_shapes(region.face, region_face, ShapeType.FACE):
                self.assertTrue(region.name in mssg)
                self.assertTrue(str_to_find in mssg)
                return

    def __assess_symmetry(
            self, sym_type: SymmetryType, vertices: List[Any]) -> None:
        """
        Method that assesses the correct application of the given symmetry
        to the lattice.

        Parameters
        ----------
        sym_type : SymmetryType
            The type of symmetry to verify.
        vertices : List[Any]
            The list of vertices identifying the shape of the symmetry.
        """
        # Build the edges identifying the shape of the symmetry
        edges = [make_edge(vertices[i], vertices[(i+1) % len(vertices)])
                    for i in range(len(vertices))]
        # Apply the symmetry
        self.lattice.apply_symmetry(sym_type)
        # Verify the correct symmetry application
        self.assertTrue(
            are_same_shapes(
                make_face(build_compound_borders(self.lattice.lattice_symm)),
                make_face(edges),
                ShapeType.FACE
            )
        )

    def __assess_lattice_init(
            self,
            cells: List[Cell],
            box_layers: List[float],
            no_rings: int,
            center: Tuple[float] | None = None) -> None:
        """
        Method that assesses the correct initialization of a `Lattice`
        instance depending on whether any cells, any box layers, or center
        have been provided.

        Parameters
        ----------
        cells : List[Cell]
            List of cells in the lattice.
        box_layers : List[float]
            List of the thicknesses of the box enclosing the lattice.
        no_rings : int
            The number of rings of cells in the lattice.
        center : Tuple[float] | None = None
            The center of the lattice, if any.
        """
        if center is None:
            center = (0.0, 0.0, 0.0)
        self.assertTrue(
            are_same_shapes(
                self.lattice.lattice_center,
                make_vertex(center),
                ShapeType.VERTEX)
        )
        if cells:
            lx = cells[0].figure.lx
            ly = cells[0].figure.ly
            self.assertTrue(self.lattice.layers and all(
                layer for layer in self.lattice.layers))
            self.assertEqual(self.lattice.cells_type, cells[0].cell_type)
            self.assertEqual(
                self.lattice.cells_rot, math.degrees(cells[0].rotation))
            self.assertAlmostEqual(
                self.lattice.distance,
                max(math.dist(
                    get_point_coordinates(self.lattice.lattice_center),
                    get_point_coordinates(c.figure.o)) for c in cells)
            )
        else:
            lx = 0.0
            ly = 0.0
            self.assertFalse(self.lattice.layers and all(
                layer for layer in self.lattice.layers))
            self.assertIsNone(self.lattice.cells_type)
            self.assertIsNone(self.lattice.cells_rot)
            self.assertEqual(self.lattice.distance, 0.0)

        self.assertEqual(len(self.lattice.lattice_cells), len(cells))
        self.assertEqual(len(self.lattice.layers[0]), len(cells))
        self.assertEqual(self.lattice.name, self.name)
        self.assertIsNone(self.lattice.lattice_entry_id)
        self.assertEqual(self.lattice.rings_no, no_rings)
        self.assertEqual(self.lattice.type_geo, LatticeGeometryType.ISOTROPIC)
        self.assertEqual(self.lattice.symmetry_type, SymmetryType.FULL)
        self.assertEqual(self.lattice.lx, lx)
        self.assertEqual(self.lattice.ly, ly)
        self.assertEqual(len(self.lattice.box_layers), len(box_layers))
        self.assertIsNone(self.lattice.lattice_symm)
        self.assertEqual(len(self.lattice.regions), 0)
        self.assertEqual(self.lattice.displayed_geom,
                         GeometryType.TECHNOLOGICAL)
        self.assertFalse(self.lattice.is_update_needed)

        if box_layers:
            self.assertIsNotNone(self.lattice.lattice_box)
        else:
            self.assertIsNone(self.lattice.lattice_box)

    def __build_hex_box(self) -> Any:
        """
        Method that builds a face object made by several encapsulated
        hexagon surfaces in number equal to the box layers.

        Returns
        -------
        Any
            A face object made by assembling several hexagons together.
        """
        box_surfaces: List[Hexagon] = []
        box_apothem = self.hex_cells[0].edge_length * (
            2 + math.sin(math.pi/6))
        if self.box_layers[0] > 0:
            box_layers = [0.0] + self.box_layers
        for thick in box_layers:
            box_apothem += thick
            box_surfaces.append(
                build_hexagon(
                    box_apothem,
                    get_point_coordinates(self.lattice.lattice_center)))
        # Perform a 'partition' operation to assemble the lattice box
        return make_partition(
            [hex.face for hex in box_surfaces], [], ShapeType.FACE)

    def __build_rect_box(self) -> Any:
        """
        Method that builds a face object made by several encapsulated
        rectangle surfaces in number equal to the box layers.

        Returns
        -------
        Any
            A face object made by assembling several rectangles together.
        """
        # Declare the starting dimensions of the box
        x_min, x_max, y_min, y_max = get_bounding_box(
            make_compound([cell.face for cell in self.lattice.lattice_cells]))
        height = y_max - y_min
        width = x_max - x_min
        box_surfaces: List[Rectangle] = []
        if self.box_layers[0] > 0:
            box_layers = [0.0] + self.box_layers
        # Build a rectangle for each layer
        for thick in box_layers:
            height += 2*thick
            width += 2*thick
            box_surfaces.append(
                Rectangle(get_point_coordinates(self.lattice.lattice_center),
                          height,
                          width))
        # Perform a 'partition' operation to assemble the lattice box
        return make_partition(
            [rect.face for rect in box_surfaces], [], ShapeType.FACE)


class TestLatticesFunctions(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the functions
    declared in the `lattices.py` module.

    Attributes
    ----------
    cells : List[Cell]
        A list of `Cell` objects the lattice is made of.
    lattice : Lattice
        The `Lattice` object to test, made either by cartesian or hexagonal
        cells.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the functions declared
        in the `lattices.py` module.
        It initializes the attributes common to all the tests.
        """
        # Declare an hexagonal cell with a circular region and a sectorization
        cell = HexCell()
        cell.add_circle(0.25)
        cell.sectorize([1, 6], [0, 0])
        # Initialize the attributes
        self.cells: List[Cell] = set_up_hex_cells(cell)
        self.lattice: Lattice = Lattice(self.cells)

    def test_get_compound_from_geometry(self) -> None:
        """
        Method that tests the implementation of the function
        `get_compound_from_geometry` declared in the `lattices.py`
        module.
        """
        # Extract and verify the correct compound is returned depending on
        # the type of geometry
        cmpd = get_compound_from_geometry(GeometryType.TECHNOLOGICAL,
                                          self.lattice.lattice_cells)
        no_edges = 7
        self.assertEqual(
            len(extract_sub_shapes(cmpd, ShapeType.EDGE)),
            no_edges*len(self.lattice.lattice_cells)
        )
        self.assertEqual(
            get_basic_properties(cmpd)[1],
            len(self.lattice.lattice_cells) * get_basic_properties(
                self.lattice.lattice_cells[0].face)[1]
        )
        cmpd = get_compound_from_geometry(GeometryType.SECTORIZED,
                                          self.lattice.lattice_cells)
        no_edges = 18
        self.assertEqual(
            len(extract_sub_shapes(cmpd, ShapeType.EDGE)),
            no_edges*len(self.lattice.lattice_cells)
        )
        self.assertEqual(
            get_basic_properties(cmpd)[1],
            len(self.lattice.lattice_cells) * get_basic_properties(
                self.lattice.lattice_cells[0].face)[1]
        )

    def test_get_changed_cells(self) -> None:
        """
        Method that tests the implementation of the function
        `get_changed_cells` declared in the `lattices.py`
        module.
        """
        # Set up a cartesian cell with circles and properties
        rect_cell = RectCell()
        radii = [0.2, 0.3, 0.4]
        for r in radii:
            rect_cell.add_circle(r)
        # Instantiate the lattice without any cell
        self.lattice = Lattice()
        self.lattice.add_rings_of_cells(rect_cell, 2)
        # Add a cell so that it overlaps the four in the upper-right section
        self.lattice.add_cell(rect_cell, (1.0, 1.0, 0.0))
        self.lattice.build_regions()

        # Get the cells whose geometry layout has changed
        changed_cells = get_changed_cells(self.lattice)
        # Verify the correct cells are returned
        self.assertEqual(len(changed_cells), 4)
        indices = ["3", "10", "11", "12"]
        for cell in changed_cells:
            self.assertFalse(
                are_same_shapes(cell.face, rect_cell.face, ShapeType.COMPOUND)
            )
            for indx in indices:
                self.assertTrue(indx in cell.name)
                # Remove the index as already checked
                indices.pop(indices.index(indx))
                break


if __name__ == "__main__":
    unittest.main()
