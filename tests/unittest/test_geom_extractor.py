"""
Module containing unittest classes to assess that the classes and functions
of the `glow.generator.geom_extractor` module have a valid implementation.
"""
import math
import unittest

from typing import Any, Dict, List, Tuple

from glow.generator.geom_extractor import Boundary, Edge, Face, build_edge_id
from glow.geometry_layouts.geometries import Rectangle
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_kind_of_shape, is_point_inside_shape, make_arc_edge, make_circle, \
    make_compound, make_edge, make_face, make_partition, make_vertex, \
    set_shape_name
from glow.support.types import NAME_EDGE_TYPE, CellType, EdgeType, \
    LatticeGeometryType, SymmetryType
from glow.support.utility import are_same_shapes
from support_funcs import BoundaryData, build_boundary_data


class TestBoundary(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `Boundary`
    class that provides the data structure for a GEOM edge object being a
    border of the lattice.

    Attributes
    ----------
    border : Any
        An edge object representing a border of a lattice.
    lattice_center : Any
        A vertex object representing the center of the lattice.
    dimensions : Dict[CellType, Tuple[float, float]]
        Storing the lattice characteristic dimensions used in the tests
        according to the type of cells.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the `Boundary` class.
        It initializes the attributes common to all the tests.
        """
        self.border: Any = make_edge(
            make_vertex((0.0, 0.0, 0.0)), make_vertex((2.0, 0.0, 0.0))
        )
        self.lattice_center: Any = make_vertex((0.0, 0.0, 0.0))
        self.dimensions: Dict[CellType, Tuple[float, float]] = {
            CellType.HEX: (1.0, math.sqrt(3)/2),
            CellType.RECT: (2.0, 1.0)
        }

    def test_init_exceptions(self) -> None:
        """
        Method that tests that exceptions are raised when initializing the
        `Boundary` class wrongly, i.e. with an edge having an angle different
        from 0° or 90° for a `RECTAGLE_TRAN` type of lattice, or with an edge
        belonging to a Y-oriented hexagon.
        """
        # Cartesian case with type_geo = 'RECTAGLE_TRAN' (angle != 0°, 90°)
        with self.assertRaises(RuntimeError):
            _ = Boundary(
                make_edge(
                    make_vertex((0.5, 0.0, 0.0)), make_vertex((1.0, 0.5, 0.0))
                ),
                LatticeGeometryType.RECTANGLE_TRAN,
                self.lattice_center,
                self.dimensions[CellType.RECT]
            )
        # Hexagonal case with type_geo = 'HEXAGON_TRAN' (angle = 90°)
        with self.assertRaises(RuntimeError):
            _ = Boundary(
                make_edge(
                    make_vertex((0.0, 0.0, 0.0)), make_vertex((0.0, 1.0, 0.0))
                ),
                LatticeGeometryType.HEXAGON_TRAN,
                self.lattice_center,
                self.dimensions[CellType.HEX]
            )
        # Initialization with a shape whose type is not 'EDGE'
        with self.assertRaises(RuntimeError):
            _ = Boundary(
                make_face(make_circle(self.lattice_center, None, 1.0)),
                LatticeGeometryType.RECTANGLE_TRAN,
                self.lattice_center,
                self.dimensions[CellType.RECT]
            )

    def test_init_hex(self) -> None:
        """
        Method that tests the initialization of the `Boundary` class for
        a hexagonal-type lattice with different symmetries and type of
        geometries.
        """
        # Full hexagon case
        self.__assess_init_full_hex()

        # Sixth of symmetry cases
        self.__assess_init_boundary(
            CellType.HEX, SymmetryType.SIXTH, LatticeGeometryType.SA60)
        self.__assess_init_boundary(
            CellType.HEX, SymmetryType.SIXTH, LatticeGeometryType.RA60)
        self.__assess_init_boundary(
            CellType.HEX,
            SymmetryType.SIXTH,
            LatticeGeometryType.SYMMETRIES_TWO)
        self.__assess_init_boundary(
            CellType.HEX, SymmetryType.SIXTH, LatticeGeometryType.ROTATION)

        # Third of symmetry cases (center redeclared to fit the case)
        x, y = self.dimensions[CellType.HEX]
        self.lattice_center = make_vertex((x/2, y, 0.0))
        self.__assess_init_boundary(
            CellType.HEX, SymmetryType.THIRD, LatticeGeometryType.R120)
        self.__assess_init_boundary(
            CellType.HEX, SymmetryType.THIRD, LatticeGeometryType.ROTATION)

        # Twelfth of symmetry cases
        self.__assess_init_boundary(
            CellType.HEX, SymmetryType.TWELFTH, LatticeGeometryType.S30)
        self.__assess_init_boundary(
            CellType.HEX,
            SymmetryType.TWELFTH,
            LatticeGeometryType.SYMMETRIES_TWO)

    def test_init_rect(self) -> None:
        """
        Method that tests the initialization of the `Boundary` class for
        a cartesian-type lattice with different symmetries and type of
        geometries.
        """
        # Verify an exception is raised when instantiating with an ISOTROPIC
        # type of geometry
        with self.assertRaises(RuntimeError):
            _ = Boundary(
                self.border,
                LatticeGeometryType.ISOTROPIC,
                self.lattice_center,
                self.dimensions
            )
        # Assess the initialization for a full hexagon
        self.__assess_init_boundary(
            CellType.RECT,
            SymmetryType.FULL,
            LatticeGeometryType.RECTANGLE_TRAN)

        # Half of symmetry cases
        self.__assess_init_boundary(
            CellType.RECT,
            SymmetryType.HALF,
            LatticeGeometryType.RECTANGLE_SYM)
        self.__assess_init_boundary(
            CellType.RECT,
            SymmetryType.HALF,
            LatticeGeometryType.SYMMETRIES_TWO)

        # Quarter of symmetry cases
        self.__assess_init_boundary(
            CellType.RECT,
            SymmetryType.QUARTER,
            LatticeGeometryType.RECTANGLE_SYM)
        self.__assess_init_boundary(
            CellType.RECT,
            SymmetryType.QUARTER,
            LatticeGeometryType.SYMMETRIES_TWO)

        # Eighth of symmetry cases
        self.__assess_init_boundary(
            CellType.RECT,
            SymmetryType.EIGHTH,
            LatticeGeometryType.RECTANGLE_EIGHT)
        self.__assess_init_boundary(
            CellType.RECT,
            SymmetryType.EIGHTH,
            LatticeGeometryType.SYMMETRIES_TWO)

    def test_find_edges_on_border(self) -> None:
        """
        Method that tests the correct implementation of the method
        `find_edges_on_border` of the `Boundary` class.
        """
        # Initialize a 'Boundary' instance
        boundary = Boundary(
            self.border,
            LatticeGeometryType.RECTANGLE_TRAN,
            self.lattice_center,
            self.dimensions[CellType.RECT])
        # Declare edges laying on the borders
        edges_on_border = [
            make_edge(
                make_vertex((0.0, 0.0, 0.0)), make_vertex((0.1, 0.0, 0.0))),
            make_edge(
                make_vertex((0.1, 0.0, 0.0)), make_vertex((0.5, 0.0, 0.0))),
            make_edge(
                make_vertex((0.5, 0.0, 0.0)), make_vertex((1.5, 0.0, 0.0))),
            make_edge(
                make_vertex((1.5, 0.0, 0.0)), make_vertex((2.0, 0.0, 0.0))),
            make_edge(
                make_vertex((2.0, 0.0, 0.0)), make_vertex((2.0, 0.5, 0.0))),
            make_edge(
                make_vertex((2.0, 0.5, 0.0)), make_vertex((2.0, 1.0, 0.0))),
            make_edge(
                make_vertex((2.0, 1.0, 0.0)), make_vertex((1.0, 1.0, 0.0))),
            make_edge(
                make_vertex((1.0, 1.0, 0.0)), make_vertex((0.0, 1.0, 0.0))),
            make_edge(
                make_vertex((0.0, 1.0, 0.0)), make_vertex((0.0, 0.0, 0.0))),
        ]
        # Build the dictionary ID VS edge object
        id_vs_edges = {}
        for i, e in enumerate(edges_on_border):
            # Set the edge's name and add an entry to the dictionary
            set_shape_name(e, f"EDGE_{i + 1}")
            id_vs_edges[build_edge_id(e)] = e

        # Identify which edges lay on the border the 'Boundary' instance
        # refers to
        boundary.find_edges_on_border(
            make_compound(edges_on_border), id_vs_edges)
        # Verify the correct edges have been identified
        indexes = set(boundary.edge_indxs)
        self.assertEqual(len(indexes), 4)
        for indx in indexes:
            self.assertIn(indx, [1, 2, 3, 4])

        # Verify an exception is raised when the edges contains an arc of
        # circle (the boundary itself must be an arc of circle)
        boundary = Boundary(
            make_arc_edge(make_vertex((0.1, 0.1, 0.0)),
                          make_vertex((0.0, 0.1, 0.0)),
                          make_vertex((0.1, 0.0, 0.0))),
            LatticeGeometryType.RECTANGLE_SYM,
            self.lattice_center,
            self.dimensions[CellType.RECT])
        edges_on_border = [
            make_arc_edge(make_vertex((0.1, 0.1, 0.0)),
                          make_vertex((0.0, 0.1, 0.0)),
                          make_vertex((0.05, 0.0, 0.0)))]
        with self.assertRaises(RuntimeError):
            boundary.find_edges_on_border(
                make_compound(edges_on_border), id_vs_edges)

    def test_get_bc_type_number(self) -> None:
        """
        Method that tests the correct implementation of the method
        `get_bc_type_number` of the `Boundary` class.
        """
        # Initialize a 'Boundary' instance
        boundary = Boundary(
            self.border,
            LatticeGeometryType.RECTANGLE_TRAN,
            self.lattice_center,
            self.dimensions[CellType.RECT])

        # Verify the correct 'BoundaryType' index is returned, i.e. 2
        self.assertEqual(boundary.get_bc_type_number(), 2)

    def __assess_boundary(
            self,
            bd: BoundaryData,
            i: int,
            border: Any,
            type_geo: LatticeGeometryType) -> None:
        """
        Method that verifies if the attributes of a `Boundary` instance have
        been correctly assigned by comparing them with the ones stored in
        a `BoundaryData` data structure.

        Parameters
        ----------
        bd : BoundaryData
            Data structure collecting all the assess values.
        i : int
            Index to retrieve the values corresponding to the `Boundary`
            instance in the `BoundaryData` object.
        border : Any
            The GEOM edge object representing one of the lattice's borders.
        type_geo : LatticeGeometryType
            Element of the `LatticeGeometryType` enumeration indicating the
            lattice type of geometry.
        """
        # Initialize a 'Boundary' instance
        boundary = Boundary(
                border,
                type_geo,
                self.lattice_center,
                bd.dimensions)
        # Verify the correct assignment
        self.assertTrue(
            math.isclose(boundary.angle, bd.angles[i]),
            f"{i}) {boundary.angle}, {bd.angles[i]}"
        )
        self.assertEqual(boundary.type, bd.bd_type[i])
        self.assertEqual(boundary.tx, bd.axis[i][0])
        self.assertEqual(boundary.ty, bd.axis[i][1])
        self.assertEqual(len(boundary.edge_indxs), 0)

    def __assess_init_full_hex(self) -> None:
        """
        Method that assesses the initialization of a `Boundary` instance for
        a full hexagonal geometry.
        It verifies that the initialization with an `ISOTROPIC` geometry type
        raises a `RuntimeError`.
        It also assesses that the `Boundary` initialization for each edge of
        the full hexagonal geometry with `HEXAGON_TRAN` type of geometry
        correctly assignes the values to its attributes.
        """
        # Verify an exception is raised when instantiating with an ISOTROPIC
        # type of geometry
        with self.assertRaises(RuntimeError):
            _ = Boundary(
                self.border,
                LatticeGeometryType.ISOTROPIC,
                self.lattice_center,
                self.dimensions
            )
        # Assess the initialization for a full hexagon
        self.__assess_init_boundary(
            CellType.HEX,
            SymmetryType.FULL,
            LatticeGeometryType.HEXAGON_TRAN)

    def __assess_init_boundary(
            self,
            cells_type: CellType,
            symm_type: SymmetryType,
            type_geo: LatticeGeometryType) -> None:
        """
        Method that assesses the initialization of a `Boundary` instance for
        a generic lattice, having either hexagonal or cartesian cells.
        It creates a reference `BoundaryData` object by calling a builder
        function with the lattice characteristic dimensions, the type of
        cells, the applied symmetry and the lattice type of geometry.

        The correct assignment of the `Boundary` attributes is verified by
        comparing them with the data stored in the built `BoundaryData`
        object.

        Parameters
        ----------
        cells_type : CellType
            An element of the `CellType` enumeration.
        symm_type : SymmetryType
            The type of symmetry applied to the lattice.
        type_geo : LatticeGeometryType
            The type of lattice geometry.
        """
        # Instantiate a 'BoundaryData' object for comparison purposes
        bd: BoundaryData = build_boundary_data(
            self.dimensions[cells_type],
            cells_type,
            symm_type,
            type_geo)
        for i, e in enumerate(bd.edges):
            # Initialize and verify each 'Boundary' instance
            self.__assess_boundary(bd, i, e, type_geo)


class TestEdge(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `Edge`
    class that provides the data structure for any GEOM edge object
    contained in the lattice.

    Attributes
    ----------
    edge_faces : List[Tuple[Any, Face]]
        List of tuples, each having an edge and the associated `Face`
        objects
    shape : Any
        A shape built by combining a square and a circle surfaces.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `Edge` class.
        It initializes the attributes common to all the tests.
        """
        # Build two adjacent rectangles sharing one edge
        r1 = Rectangle((0.5, 0.5, 0.0))
        r2 = Rectangle((1.5, 0.5, 0.0))
        set_shape_name(r1.face, "FACE_1")
        set_shape_name(r2.face, "FACE_2")
        prop_val = 'MAT'
        # Build a list of the edges with the face(s) they belong to
        self.edge_faces: List[Tuple[Any, Face]] = [
            (r1.borders[0], Face(r1.face, prop_val)),
            (r1.borders[1], Face(r1.face, prop_val), Face(r2.face, prop_val)),
            (r1.borders[2], Face(r1.face, prop_val)),
            (r1.borders[3], Face(r1.face, prop_val)),
            (r2.borders[0], Face(r2.face, prop_val)),
            (r2.borders[1], Face(r2.face, prop_val)),
            (r2.borders[2], Face(r2.face, prop_val))
        ]
        # Assign a name to the edges so that they have a global index
        for i in range(len(self.edge_faces)):
            set_shape_name(self.edge_faces[i][0], f"EDGE_{i+1}")

        # Build a surface made by a circle inside a square
        self.shape: Any = make_partition(
            [r1.face],
            [make_circle(make_vertex((0.5, 0.5, 0.0)), None, 0.25)],
            ShapeType.FACE
        )

    def test_init(self) -> None:
        """
        Method that tests the initialization of the `Edge` class.
        """
        # Verify an exception is raised when the edge does not have a name
        # or no index can be retrieved from its name or the provided shape
        # is not an edge
        edge_faces = (
            make_edge(make_vertex((0.0, 0.0, 0.0)),
                      make_vertex((1.0, 0.0, 0.0))),
            [])
        set_shape_name(edge_faces[0], "")
        with self.assertRaises(RuntimeError):
            _ = Edge(*edge_faces)
        set_shape_name(edge_faces[0], "EDGE1")
        with self.assertRaises(RuntimeError):
            _ = Edge(*edge_faces)
        with self.assertRaises(RuntimeError):
            _ = Edge(
                make_face(
                    make_circle(make_vertex((0.0, 0.0, 0.0)), None, 1.0)),
                None
            )

        # Instantiate an 'Edge' object and verify the correct attributes
        # assignment
        for i, ef in enumerate(self.edge_faces):
            edge = Edge(*ef)
            self.assertTrue(
                are_same_shapes(edge.edge, ef[0], ShapeType.EDGE)
            )
            self.assertEqual(edge.data, get_kind_of_shape(ef[0]))
            self.assertEqual(edge.no, i+1)
            self.assertEqual(edge.kind, EdgeType.SEGMENT)
            self.assertEqual(edge.left, ef[1])
            if i == 1:
                self.assertEqual(edge.right, ef[2])
            else:
                self.assertEqual(edge.right, None)

    def test_add_face(self) -> None:
        """
        Method that tests the implementation of the private method
        `__add_face` of the `Edge` class.
        """
        # Get a segment and a circle edge from the saved shape
        edges = extract_sub_shapes(self.shape, ShapeType.EDGE)
        sgmnt = [e for e in edges
                 if str(get_kind_of_shape(e)[0]) == 'SEGMENT'][0]
        circle = [e for e in edges
                  if str(get_kind_of_shape(e)[0]) == 'CIRCLE'][0]
        # Get the corresponding faces
        faces = extract_sub_shapes(self.shape, ShapeType.FACE)
        sgmnt_face = Face.__new__(Face)
        sgmnt_face.face = [
            f for f in faces
            if str(get_kind_of_shape(f)[0]) != 'DISK_CIRCLE'][0]
        circle_face = Face.__new__(Face)
        circle_face.face = [
            f for f in faces
            if str(get_kind_of_shape(f)[0]) == 'DISK_CIRCLE'][0]

        # Instantiate an `Edge` object without any parameter
        edge = Edge.__new__(Edge)
        # Set the attributes so to replicate an edge segment
        edge.data = get_kind_of_shape(sgmnt)
        edge.kind = NAME_EDGE_TYPE[str(edge.data[0])][0]
        edge.edge = sgmnt
        edge.right = None
        edge.left = None
        # Verify the correct assignment of the segment face to the 'left'
        # attribute
        edge._Edge__add_face(sgmnt_face)
        self.assertTrue(
            are_same_shapes(edge.left.face, sgmnt_face.face, ShapeType.FACE)
        )

        # Instantiate an `Edge` object without any parameter
        edge = Edge.__new__(Edge)
        # Set the attributes so to replicate an edge circle
        edge.data = get_kind_of_shape(circle)
        edge.kind = NAME_EDGE_TYPE[str(edge.data[0])][0]
        edge.edge = circle
        edge.right = None
        edge.left = None
        # Verify the correct assignment of the 'left' and 'right' attributes
        edge._Edge__add_face(sgmnt_face)
        self.assertTrue(
            are_same_shapes(edge.left.face, sgmnt_face.face, ShapeType.FACE)
        )
        edge._Edge__add_face(circle_face)
        self.assertTrue(
            are_same_shapes(edge.right.face, circle_face.face, ShapeType.FACE)
        )

    def test_origin(self) -> None:
        """
        Method that tests the implementation of the private method `__origin`
        of the `Edge` class.
        """
        # Get a segment and a circle edge from the saved shape
        edges = extract_sub_shapes(self.shape, ShapeType.EDGE)
        sgmnt = [e for e in edges
                 if str(get_kind_of_shape(e)[0]) == 'SEGMENT'][0]
        circle = [e for e in edges
                  if str(get_kind_of_shape(e)[0]) == 'CIRCLE'][0]
        faces = extract_sub_shapes(self.shape, ShapeType.FACE)
        sgmnt_face = [f for f in faces
                      if str(get_kind_of_shape(f)[0]) != 'DISK_CIRCLE'][0]
        circle_face = [f for f in faces
                       if str(get_kind_of_shape(f)[0]) == 'DISK_CIRCLE'][0]
        # Assign names to edges and faces
        set_shape_name(sgmnt, "EDGE_1")
        set_shape_name(circle, "EDGE_2")
        set_shape_name(sgmnt_face, "FACE_1")
        set_shape_name(circle_face, "FACE_2")
        # Instantiate an 'Edge' object for both
        e1 = Edge(sgmnt, Face(sgmnt_face, 'MAT'))
        e2 = Edge(circle, Face(circle_face, 'MAT'))
        # Verify the correct value for the edge's origin is returned
        self.assertEqual(e1._Edge__origin(), (0.0, 0.0, 0.0))
        self.assertEqual(e2._Edge__origin(), (0.75, 0.5, 0.0))


class TestFace(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `Face`
    class that provides the data structure for any GEOM face object
    contained in the lattice.
    """
    def test_init(self) -> None:
        """
        Method that tests the initialization of the `Face` class.
        """
        # Build a rectangular surface and set the face's name
        rect_face = Rectangle((0.5, 0.5, 0.0))
        prop_val = 'MAT'

        # Verify an exception is raised when instantiating a 'Face' object
        # from a shape without any name or with an invalid format or if the
        # provided shape is not a face
        set_shape_name(rect_face.face, "FACE")
        with self.assertRaises(RuntimeError):
            face = Face(rect_face.face, prop_val)
        set_shape_name(rect_face.face, "FACE1")
        with self.assertRaises(RuntimeError):
            face = Face(rect_face.face, prop_val)
        with self.assertRaises(RuntimeError):
            _ = Face(rect_face.borders[0], prop_val)

        # Instantiate the 'Face' object with the valid format for the name
        set_shape_name(rect_face.face, "FACE_1")
        face = Face(rect_face.face, prop_val)

        # Verify the correct initialization of the attributes
        self.assertTrue(
            are_same_shapes(face.face, rect_face.face, ShapeType.FACE)
        )
        self.assertEqual(face.property, prop_val)
        self.assertEqual(face.no, 1)
        self.assertEqual(face.sort_index, 1)
        self.assertTrue(
            is_point_inside_shape(
                make_vertex(face.inner_point), rect_face.face)
        )
        self.assertTrue(
            all(build_edge_id(b) in face.edge_vs_id.values()
                for b in rect_face.borders)
        )
