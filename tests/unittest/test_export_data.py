"""
Module containing unittest classes to assess that the classes and functions
of the `glow.generator.geom_extractor` module have a valid implementation.
"""
import math
import unittest

from typing import Any, Dict, List, Tuple

from glow.generator.geom_extractor import classify_layout_edges
from glow.generator.geom_extractor import BoundaryData, EdgeData, FaceData, \
    build_edge_id
from glow.geometry_layouts.geometries import Circle, Rectangle
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_kind_of_shape, is_point_inside_shape, make_arc_center, \
    make_circle, make_compound, make_edge, make_face, make_vertex, \
    set_shape_name
from glow.support.types import EDGE_NAME_VS_TYPE, BoundaryType, EdgeType, \
    LayoutGeometryType, LayoutType, PropertyType, SymmetryType
from glow.support.utility import are_same_shapes
from tests.unittest.support_funcs import BoundaryInfo, build_boundary_info


class TestBoundaryData(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `BoundaryData`
    class that provides the data structure for a GEOM edge object being a
    border of the geometry layout.

    Attributes
    ----------
    border : Any
        An edge object representing a border of the geometry layout.
    centre : Any
        A vertex object representing the centre of the geometry layout.
    dimensions : Dict[LayoutType, Tuple[float, float]]
        Storing the characteristic dimensions of the geometry layout used in
        the tests according to the type of layout.
    """
    def setUp(self):
        """
        Method that sets up the test environment for the `BoundaryData` class.
        It initialises the attributes common to all the tests.
        """
        self.border: Any = make_edge(
            make_vertex((0.0, 0.0, 0.0)), make_vertex((2.0, 0.0, 0.0))
        )
        self.centre: Any = make_vertex((0.5, math.sqrt(3)/2, 0.0))
        self.dimensions: Dict[LayoutType, Tuple[float, float]] = {
            LayoutType.HEX: (1.0, math.sqrt(3)/2),
            LayoutType.RECT: (1.0, 1.0)
        }

    def test_init_exceptions(self) -> None:
        """
        Method that tests that exceptions are raised when initialising the
        `BoundaryData` class wrongly, i.e. with an edge having an angle
        different from 0° or 90° for a `RECTAGLE_TRAN` type of layout, or
        with an edge belonging to a Y-oriented hexagon.
        """
        # Cartesian case with type_geo = 'RECTAGLE_TRAN' (angle != 0°, 90°)
        with self.assertRaises(RuntimeError):
            _ = BoundaryData(
                make_edge(
                    make_vertex((0.5, 0.0, 0.0)), make_vertex((1.0, 0.5, 0.0))
                ),
                LayoutGeometryType.RECTANGLE_TRAN,
                self.centre,
                self.dimensions[LayoutType.RECT]
            )
        # Hexagonal case with type_geo = 'HEXAGON_TRAN' (angle = 90°)
        with self.assertRaises(RuntimeError):
            _ = BoundaryData(
                make_edge(
                    make_vertex((0.0, 0.0, 0.0)), make_vertex((0.0, 1.0, 0.0))
                ),
                LayoutGeometryType.HEXAGON_TRAN,
                self.centre,
                self.dimensions[LayoutType.HEX]
            )
        # Initialisation with a shape whose type is not 'EDGE'
        with self.assertRaises(RuntimeError):
            _ = BoundaryData(
                make_face(make_circle(self.centre, None, 1.0)),
                LayoutGeometryType.RECTANGLE_TRAN,
                self.centre,
                self.dimensions[LayoutType.RECT]
            )

    def test_init_hex(self) -> None:
        """
        Method that tests the initialisation of the `BoundaryData` class for
        a hexagonal-type layout with different symmetries and type of
        geometries.
        """
        # Full hexagon case
        self.__assess_init_full_hex()

        # Sixth of symmetry cases
        self.__assess_init_boundary(
            LayoutType.HEX, SymmetryType.SIXTH, LayoutGeometryType.SA60
        )
        self.__assess_init_boundary(
            LayoutType.HEX, SymmetryType.SIXTH, LayoutGeometryType.RA60
        )
        self.__assess_init_boundary(
            LayoutType.HEX,
            SymmetryType.SIXTH,
            LayoutGeometryType.SYMMETRIES_TWO
        )
        self.__assess_init_boundary(
            LayoutType.HEX, SymmetryType.SIXTH, LayoutGeometryType.ROTATION
        )

        # Third of symmetry cases
        self.__assess_init_boundary(
            LayoutType.HEX, SymmetryType.THIRD, LayoutGeometryType.R120
        )
        self.__assess_init_boundary(
            LayoutType.HEX, SymmetryType.THIRD, LayoutGeometryType.ROTATION
        )

        # Twelfth of symmetry cases
        self.__assess_init_boundary(
            LayoutType.HEX, SymmetryType.TWELFTH, LayoutGeometryType.S30
        )
        self.__assess_init_boundary(
            LayoutType.HEX,
            SymmetryType.TWELFTH,
            LayoutGeometryType.SYMMETRIES_TWO
        )

    def test_init_hex_precision(self) -> None:
        """
        Method that tests the initialisation of the `BoundaryData` class for
        a hexagonal-type layout with `THIRD` and `SIXTH` symmetry types,
        adopting a rotational BC on the internal borders and a translation BC
        on the external one. The geometry layout exhibit numerical floating
        point precision noise on the vertices. The test aims to verify that
        the code correctly assigns the BC type `BoundaryType.ROTATION` to the
        internal borders and the `BoundaryType.TRANSLATION` to the external
        one.
        """
        # Declare the XY dimensions of the hexagonal case
        lx = 20.3214824999
        ly = lx/2 * math.tan(math.pi/3)
        self.centre = make_vertex((lx/2.0, ly-1e-5, 0.0))
        # R120 case
        bd = BoundaryInfo(
            vertices=[
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((lx, 0.0, 0.0)),
                make_vertex((3/2*lx, ly, 0.0)),
                make_vertex((lx/2, ly, 0.0))
            ],
            axis=[(0.0, 0.0), (lx, 0.0), (lx/2, ly), (0.0, 0.0)],
            angles=[0.0, 60.0, 0.0, 60.0],
            dimensions=(lx, ly),
            bd_type=[
                BoundaryType.TRANSLATION,
                BoundaryType.TRANSLATION,
                BoundaryType.ROTATION,
                BoundaryType.ROTATION
            ]
        )
        # Initialise and verify the 'BoundaryData' for each edge
        for i, e in enumerate(bd.edges):
            self.__assess_boundary(bd, i, e, LayoutGeometryType.R120)

        # RA60 case
        bd = BoundaryInfo(
            vertices=[
                make_vertex((0.0, 0.0, 0.0)),
                make_vertex((lx, 0.0, 0.0)),
                make_vertex((lx/2, ly, 0.0))
            ],
            axis=[(0.0, 0.0), (lx, 0.0), (0.0, 0.0)],
            angles=[0.0, 120.0, 60.0],
            dimensions=(lx, ly),
            bd_type=[
                BoundaryType.TRANSLATION,
                BoundaryType.ROTATION,
                BoundaryType.ROTATION
            ]
        )
        # Initialise and verify the 'BoundaryData' for each edge
        for i, e in enumerate(bd.edges):
            self.__assess_boundary(bd, i, e, LayoutGeometryType.RA60)

    def test_init_rect(self) -> None:
        """
        Method that tests the initialisation of the `BoundaryData` class for
        a Cartesian-type layout with different symmetries and type of
        geometries.
        """
        # Assess the initialisation for a full Cartesian layout
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.FULL,
            LayoutGeometryType.RECTANGLE_TRAN)

        # Half of symmetry cases
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.HALF,
            LayoutGeometryType.RECTANGLE_SYM
        )
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.HALF,
            LayoutGeometryType.SYMMETRIES_TWO
        )

        # Half of symmetry, along the diagonal, cases
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.DIAG,
            LayoutGeometryType.RECTANGLE_SYM
        )
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.DIAG,
            LayoutGeometryType.SYMMETRIES_TWO
        )

        # Quarter of symmetry cases
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.QUARTER,
            LayoutGeometryType.RECTANGLE_SYM
        )
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.QUARTER,
            LayoutGeometryType.SYMMETRIES_TWO
        )

        # Eighth of symmetry cases
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.EIGHTH,
            LayoutGeometryType.RECTANGLE_EIGHT
        )
        self.__assess_init_boundary(
            LayoutType.RECT,
            SymmetryType.EIGHTH,
            LayoutGeometryType.SYMMETRIES_TWO
        )

    def test_find_edges_on_border(self) -> None:
        """
        Method that tests the correct implementation of the method
        `find_edges_on_border` of the `BoundaryData` class.
        """
        # Initialise a 'BoundaryData' instance for the border of a full
        # rectangular layout
        boundary = BoundaryData(
            self.border,
            LayoutGeometryType.RECTANGLE_TRAN,
            self.centre,
            self.dimensions[LayoutType.RECT]
        )
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

        # Identify which edges lay on the border the 'BoundaryData' instance
        # refers to
        boundary.find_edges_on_border(
            make_compound(edges_on_border), id_vs_edges
        )
        # Verify the correct edges have been identified
        indexes = set(boundary.edge_indxs)
        self.assertEqual(len(indexes), 4)
        for indx in indexes:
            self.assertIn(indx, [1, 2, 3, 4])

    def test_find_edges_on_border_with_arcs(self) -> None:
        """
        Method that tests the `find_edges_on_border` method to verify that
        an exception is raised when looking for an arc of circle edge on a
        border which is an arc of circle itself.
        """
        # Declare a 'BoundaryData' for an arc of circle border
        boundary = BoundaryData(
            make_arc_center(
                make_vertex((0.1, 0.1, 0.0)),
                make_vertex((0.0, 0.1, 0.0)),
                make_vertex((0.1, 0.0, 0.0))
            ),
            LayoutGeometryType.RECTANGLE_SYM,
            self.centre,
            self.dimensions[LayoutType.RECT]
        )
        # Declare an arc of circle sub-edge on the border
        edges_on_border = [
            make_arc_center(
                make_vertex((0.1, 0.1, 0.0)),
                make_vertex((0.0, 0.1, 0.0)),
                make_vertex((0.05, 0.0, 0.0))
            )
        ]
        # An empty dictionary of edges IDs vs edge object is passed as not
        # needed for the test
        with self.assertRaises(RuntimeError):
            boundary.find_edges_on_border(make_compound(edges_on_border), {})

    def test_get_bc_type_number(self) -> None:
        """
        Method that tests the correct implementation of the method
        `get_bc_type_number` of the `BoundaryData` class.
        """
        # Initialise a 'BoundaryData' instance
        boundary = BoundaryData(
            self.border,
            LayoutGeometryType.RECTANGLE_TRAN,
            self.centre,
            self.dimensions[LayoutType.RECT]
        )

        # Verify the correct 'BoundaryType' index is returned, i.e. 2
        self.assertEqual(boundary.get_bc_type_number(), 2)

    def __assess_boundary(
            self,
            boundary_info: BoundaryInfo,
            i: int,
            border: Any,
            type_geo: LayoutGeometryType
        ) -> None:
        """
        Method that verifies if the attributes of a `BoundaryData` instance
        have been correctly assigned by comparing them with the ones stored
        in a `BoundaryInfo` data structure.

        Parameters
        ----------
        boundary_info : BoundaryInfo
            Data structure collecting all the comparison data.
        i : int
            Index to retrieve the values corresponding to the `BoundaryInfo`
            instance in the `BoundaryData` object.
        border : Any
            The GEOM edge object representing one of the layouts's borders.
        type_geo : LayoutGeometryType
            Element of the `LayoutGeometryType` enumeration indicating the
            layouts type of geometry.
        """
        # Initialise a 'BoundaryData' instance
        boundary = BoundaryData(
            border,
            type_geo,
            self.centre,
            boundary_info.dimensions
        )
        # Verify the correct assignment
        self.assertTrue(
            math.isclose(
                boundary.angle, boundary_info.angles[i], abs_tol=1e-5
            )
        )
        self.assertEqual(
            boundary.type, boundary_info.bd_type[i]
        )
        self.assertTrue(
            math.isclose(boundary.tx, boundary_info.axis[i][0])
        )
        self.assertTrue(
            math.isclose(boundary.ty, boundary_info.axis[i][1])
        )
        self.assertEqual(len(boundary.edge_indxs), 0)

    def __assess_init_boundary(
            self,
            layout_type: LayoutType,
            symm_type: SymmetryType,
            type_geo: LayoutGeometryType
        ) -> None:
        """
        Method that assesses the initialisation of a `BoundaryData` instance
        for a generic layout.
        It creates a reference `BoundaryData` object by calling a builder
        function receiving the layouts's characteristic dimensions, the type
        of layout, the applied symmetry and the type of geometry.

        The correct assignment of the `BoundaryData` attributes is verified by
        comparing them with the data stored in the built `BoundaryInfo`
        object.

        Parameters
        ----------
        layout_type : LayoutType
            An element of the `LayoutType` enumeration.
        symm_type : SymmetryType
            The type of symmetry applied to the layout.
        type_geo : LayoutGeometryType
            The type of geometry of the layout according to DRAGON5.
        """
        # Instantiate a 'BoundaryInfo' object for comparison purposes
        bd = build_boundary_info(
            self.dimensions[layout_type],
            layout_type,
            symm_type,
            type_geo
        )
        for i, e in enumerate(bd.edges):
            # Initialise and verify the 'BoundaryData' for each edge
            self.__assess_boundary(bd, i, e, type_geo)

    def __assess_init_full_hex(self) -> None:
        """
        Method that assesses the initialisation of a `BoundaryData` instance
        for a full hexagonal geometry.
        It verifies that the initialisation with an `ISOTROPIC` geometry type
        raises a `RuntimeError`.
        It also assesses that the `BoundaryData` initialisation for each edge
        of the full hexagonal geometry with `HEXAGON_TRAN` type of geometry
        correctly assignes the values to its attributes.
        """
        # Verify an exception is raised when instantiating with an ISOTROPIC
        # type of geometry
        with self.assertRaises(RuntimeError):
            _ = BoundaryData(
                self.border,
                LayoutGeometryType.ISOTROPIC,
                self.centre,
                self.dimensions
            )
        # Assess the initialisation for a full hexagon
        self.__assess_init_boundary(
            LayoutType.HEX,
            SymmetryType.FULL,
            LayoutGeometryType.HEXAGON_TRAN
        )


class TestEdgeData(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `EdgeData`
    class that provides the data structure for any GEOM edge object
    contained in the layout.

    Attributes
    ----------
    edge_faces : List[Tuple[Any, FaceData]]
        List of tuples, each having an edge and the associated `FaceData`
        objects.
    shape : Any
        A shape built by combining a square and a circle surfaces.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `EdgeData` class.
        It initialises the attributes common to all the tests.
        """
        # Build two adjacent rectangular regions sharing one edge
        r1 = Rectangle((0.5, 0.5, 0.0))
        r2 = Rectangle((1.5, 0.5, 0.0))
        reg_1 = Region(r1, "FACE_1", {PropertyType.MATERIAL: "MAT_1"})
        reg_2 = Region(r2, "FACE_2", {PropertyType.MATERIAL: "MAT_2"})
        prop_types = [PropertyType.MATERIAL]
        # Build a list of the edges with the face(s) they belong to
        self.edge_faces: List[Tuple[Any, FaceData]] = [
            (r1.borders[0], FaceData(reg_1, 1, prop_types)),
            (
                r1.borders[1],
                FaceData(reg_1, 1, prop_types),
                FaceData(reg_2, 2, prop_types)
            ),
            (r1.borders[2], FaceData(reg_1, 1, prop_types)),
            (r1.borders[3], FaceData(reg_1, 1, prop_types)),
            (r2.borders[0], FaceData(reg_2, 2, prop_types)),
            (r2.borders[1], FaceData(reg_2, 2, prop_types)),
            (r2.borders[2], FaceData(reg_2, 2, prop_types))
        ]
        # Assign a name to the edges so that they have a global index
        for i in range(len(self.edge_faces)):
            set_shape_name(self.edge_faces[i][0], f"EDGE_{i+1}")

        # Build a surface made by a circle inside a square
        self.shape: Any = r1 // Circle((0.5, 0.5, 0.0), 0.25)

    def test_init(self) -> None:
        """
        Method that tests the initialisation of the `EdgeData` class.
        """
        # Instantiate an 'EdgeData' object and verify the correct attributes
        # assignment
        for i, ef in enumerate(self.edge_faces):
            edge = EdgeData(*ef)
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

    def test_init_exceptions(self) -> None:
        """
        Method that tests the initialisation of the `EdgeData` class raises
        the expected `RuntimeError` exception when:

        - the edge does not have a name;
        - no index can be retrieved from the edge's name;.
        - the provided shape is not an edge.
        """
        # Build an edge object
        edge = make_edge(
            make_vertex((0.0, 0.0, 0.0)), make_vertex((1.0, 0.0, 0.0))
        )
        # Verify exceptions are raised when invalid names are assigned to the
        # edge
        with self.assertRaises(RuntimeError):
            _ = EdgeData(edge, None)
        set_shape_name(edge, "")
        with self.assertRaises(RuntimeError):
            _ = EdgeData(edge, None)
        set_shape_name(edge, "EDGE1")
        with self.assertRaises(RuntimeError):
            _ = EdgeData(edge, None)
        with self.assertRaises(RuntimeError):
            _ = EdgeData(self.shape, None)

    def test_determine_face_position(self) -> None:
        """
        Method that tests the implementation of the `_determine_face_position`
        method of the `EdgeData` class.
        """
        # Get a segment and a circle edge from the saved shape
        edges = extract_sub_shapes(self.shape, ShapeType.EDGE)
        sgmnt = [
            e for e in edges if str(get_kind_of_shape(e)[0]) == 'SEGMENT'
        ][0]
        circle = [
            e for e in edges if str(get_kind_of_shape(e)[0]) == 'CIRCLE'
        ][0]
        # Get the corresponding faces and build the 'FaceData' objects
        faces = extract_sub_shapes(self.shape, ShapeType.FACE)
        sgmnt_face = FaceData.__new__(FaceData)
        sgmnt_face.region = [
            Region(f) for f in faces
            if str(get_kind_of_shape(f)[0]) != 'DISK_CIRCLE'
        ][0]
        circle_face = FaceData.__new__(FaceData)
        circle_face.region = [
            Region(f) for f in faces
            if str(get_kind_of_shape(f)[0]) == 'DISK_CIRCLE'
        ][0]

        # Instantiate an `EdgeData` object without any parameter
        edge = EdgeData.__new__(EdgeData)
        # Set the attributes so to replicate an edge segment
        edge.data = get_kind_of_shape(sgmnt)
        edge.kind = EDGE_NAME_VS_TYPE[str(edge.data[0])][0]
        edge.edge = sgmnt
        edge.right = None
        edge.left = None
        # Verify the correct assignment of the segment face to the 'left'
        # attribute
        edge._determine_face_position(sgmnt_face)
        self.assertTrue(
            are_same_shapes(
                edge.left.region, sgmnt_face.region, ShapeType.FACE
            )
        )

        # Instantiate an `EdgeData` object without any parameter
        edge = EdgeData.__new__(EdgeData)
        # Set the attributes so to replicate an edge circle
        edge.data = get_kind_of_shape(circle)
        edge.kind = EDGE_NAME_VS_TYPE[str(edge.data[0])][0]
        edge.edge = circle
        edge.right = None
        edge.left = None
        # Verify the correct assignment of the 'left' and 'right' attributes
        edge._determine_face_position(sgmnt_face)
        self.assertTrue(
            are_same_shapes(
                edge.left.region, sgmnt_face.region, ShapeType.FACE
            )
        )
        edge._determine_face_position(circle_face)
        self.assertTrue(
            are_same_shapes(
                edge.right.region, circle_face.region, ShapeType.FACE
            )
        )

    def test_origin(self) -> None:
        """
        Method that tests the implementation of the method `_get_origin`
        of the `EdgeData` class.
        """
        # Get a segment and a circle edge from the saved shape
        edges = extract_sub_shapes(self.shape, ShapeType.EDGE)
        sgmnt = [
            e for e in edges if str(get_kind_of_shape(e)[0]) == 'SEGMENT'
        ][0]
        circle = [
            e for e in edges if str(get_kind_of_shape(e)[0]) == 'CIRCLE'
        ][0]
        # Get the corresponding faces and build the 'FaceData' objects
        faces = extract_sub_shapes(self.shape, ShapeType.FACE)
        sgmnt_face = FaceData.__new__(FaceData)
        sgmnt_face.region = [
            Region(f) for f in faces
            if str(get_kind_of_shape(f)[0]) != 'DISK_CIRCLE'
        ][0]
        circle_face = FaceData.__new__(FaceData)
        circle_face.region = [
            Region(f) for f in faces
            if str(get_kind_of_shape(f)[0]) == 'DISK_CIRCLE'
        ][0]
        # Assign names to edges
        set_shape_name(sgmnt, "EDGE_1")
        set_shape_name(circle, "EDGE_2")
        # Instantiate an 'Edge' object for both
        e1 = EdgeData(sgmnt, sgmnt_face)
        e2 = EdgeData(circle, circle_face)
        # Verify the correct value for the edge's origin is returned
        self.assertEqual(e1._get_origin(), (0.0, 0.0, 0.0))
        self.assertEqual(e2._get_origin(), (0.75, 0.5, 0.0))


class TestFaceData(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the `FaceData`
    class that provides the data structure for any GEOM face object
    contained in the geometry layout.
    """
    def test_init(self) -> None:
        """
        Method that tests the initialisation of the `FaceData` class.
        """
        # Build a rectangular region
        rect_face = Rectangle((0.5, 0.5, 0.0))
        region = Region(rect_face, "FACE", {PropertyType.MATERIAL: "MAT"})

        # Instantiate the 'FaceData' object
        prop_types = [PropertyType.MATERIAL]
        face_data = FaceData(region, 1, prop_types)

        # Verify the correct initialisation of the attributes
        self.assertTrue(
            are_same_shapes(face_data.region, rect_face, ShapeType.FACE)
        )
        self.assertEqual(face_data.property_types, prop_types)
        self.assertEqual(face_data.no, 1)
        self.assertEqual(face_data.sort_index, 1)
        self.assertTrue(
            is_point_inside_shape(
                make_vertex(face_data.inner_point), rect_face
            )
        )
        self.assertTrue(
            all(build_edge_id(b) in face_data.edge_vs_id.values()
                for b in rect_face.borders)
        )

    def test_init_exceptions(self) -> None:
        """
        Method that tests the initialisation of the `FaceData` class raises
        the expected `RuntimeError` exception when:

        - no properties are assigned to the region of the `FaceData` instance;
        - any of the requested properties is not assigned to the region of
          the `FaceData` instance;
        - no value for any of the requested properties is assigned to the
          region of the `FaceData` instance.
        """
        # Build a rectangular region
        region = Region(
            Rectangle((0.5, 0.5, 0.0)),
            "FACE"
        )

        # Verify exceptions are raised when invalid properties are assigned
        with self.assertRaises(RuntimeError):
            _ = FaceData(region, 1, [PropertyType.MATERIAL])

        region.properties = {}
        with self.assertRaises(RuntimeError):
            _ = FaceData(region, 1, [PropertyType.MATERIAL])

        region.properties = {PropertyType.MATERIAL: "MAT"}
        with self.assertRaises(RuntimeError):
            _ = FaceData(region, 1, [PropertyType.MACRO])

        region.properties.update({PropertyType.MACRO: ""})
        with self.assertRaises(RuntimeError):
            _ = FaceData(
                region, 1, [PropertyType.MACRO, PropertyType.MATERIAL]
            )


class TestExportDataFunctions(unittest.TestCase):
    """
    Test case for verifying the correct implementation of the functions
    declared in the `export_data.py` module.

    Attributes
    ----------
    segment : Any
        A segment-type edge object.
    arc_circle : Any
        A arc of circle-type edge object.
    circle : Any
        A circle-type edge object.
    """
    def setUp(self) -> None:
        """
        Method that sets up the test environment for the `EdgeData` class.
        It initialises the attributes common to all the tests.
        """
        # Declare the edges to test, one for each type
        self.segment: Any = make_edge(
            make_vertex((0.0, 0.0, 0.0)), make_vertex((1.0, 0.0, 0.0))
        )
        self.arc_circle: Any = make_arc_center(
            make_vertex((0.0, 1.0, 0.0)),
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((1.0, 1.0, 0.0))
        )
        self.circle: Any = make_circle(
            make_vertex((0.0, 0.0, 0.0)), None, 1.0
        )

    def test_build_edge_id(self) -> None:
        """
        Method that tests the implementation of the function `build_edge_id`.
        """
        # Strings to compare the resulting IDs with
        ref_sgmnt_id = "EDGE_SEGMENT_0_0_0_1_0_0"
        ref_arc_crcl_id = "EDGE_ARC_CIRCLE_0_1_0_0_0_1_1_0_0_0_1_1_0"
        ref_circle_id = "EDGE_CIRCLE_0_0_0_0_0_1_1"

        # Get the edges' information and verify its correctness
        sgmnt_id = build_edge_id(self.segment)
        arc_crcl_id = build_edge_id(self.arc_circle)
        circle_id = build_edge_id(self.circle)
        self.assertIn(ref_sgmnt_id, sgmnt_id)
        self.assertIn(ref_arc_crcl_id, arc_crcl_id)
        self.assertIn(ref_circle_id, circle_id)

        # Verify an exception is raised when providing a shape that is not
        # an edge
        with self.assertRaises(RuntimeError):
            _ = build_edge_id(make_face([self.circle]))

    def test_classify_layout_edges(self) -> None:
        """
        Method that tests the implementation of the `classify_layout_edges`
        function.
        """
        # Reference dictionary of strings VS edges with which compare the
        # dictionary resulting from the function call
        ref_ids_edges = {
            "EDGE_SEGMENT_0_0_0_1_0_0": self.segment,
            "EDGE_ARC_CIRCLE_0_1_0_0_0_1_1_0_0_0_1_1_0": self.arc_circle,
            "EDGE_CIRCLE_0_0_0_0_0_1_1": self.circle
        }

        # Build the classificaton of the edges
        ids_edges = classify_layout_edges(
            [self.segment, self.arc_circle, self.circle]
        )

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
            _ = classify_layout_edges([make_face(self.circle)])
