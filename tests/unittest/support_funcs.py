"""
Module declaring functions to support the execution of the unit tests.
"""
import hashlib

from copy import deepcopy
from dataclasses import dataclass, field
from math import degrees, sqrt, atan2
from pathlib import Path
from typing import Any, List, Tuple

from glow.geometry_layouts.cells import Cell, HexCell, RectCell
from glow.geometry_layouts.geometries import Surface
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    make_circle, make_edge, make_vector_from_points, make_vertex,\
    make_vertex_on_curve
from glow.support.types import BoundaryType, CellType, LatticeGeometryType, \
    PropertyType, SymmetryType
from glow.support.utility import build_contiguous_edges


@dataclass
class BoundaryData():
    """
    Dataclass storing geometric information about the boundaries of a lattice
    in terms of vertices, contiguous edges, axes, angles and characteristic
    dimensions of the lattice.

    Attributes
    ----------
    vertices : List[Any]
        List of vertex objects defining the lattice's borders.
    edges : List[Any]
        List of edge objects automatically constructed from the vertices.
    axis : List[Tuple[float, float]]
        Providing the XY directions representing the axes for each border.
    angles : List[float]
        List of angles (in degrees) associated with each axis.
    dimensions : tuple of float
        Characteristic dimensions of the lattice.
    bd_type : List[BoundaryType]
        Type associated to each lattice's border provided as element of the
        `BoundaryType` enumeration.
    """
    vertices: List[Any]
    edges: List[Any] = field(init=False)
    axis: List[Tuple[float, float]]
    angles: List[float]
    dimensions: Tuple[float, float]
    bd_type: List[BoundaryType]

    def __post_init__(self) -> None:
        """
        Method run after the dataclass initialization for building the
        contiguous edges related to the stored vertices.
        """
        self.edges = build_contiguous_edges(self.vertices)


def build_bd_full_hex(lx: float, ly: float) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing a regular
    hexagon centered at the origin. The information related to the borders,
    that is contained in the `BoundaryData` object, derives from the
    characteristic dimensions of the hexagon.
    The geometric data of the hexagon is provided in counter-clockwise order
    starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the hexagon (i.e. its edge length).
    ly : float
        Characteristic Y-dimension of the hexagon (i.e. its apothem length).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the six vertices of the full
        hexagon, the XY directions of the borders' axes and the corresponding
        angles (as needed by DRAGON), and the type of boundary, assigned as
        `TRANSLATION` for all the borders.
    """
    return BoundaryData(
        vertices=[
            make_vertex((-lx/2, -ly, 0.0)),
            make_vertex((lx/2, -ly, 0.0)),
            make_vertex((lx, 0.0, 0.0)),
            make_vertex((lx/2, ly, 0.0)),
            make_vertex((-lx/2, ly, 0.0)),
            make_vertex((-lx, 0.0, 0.0))
        ],
        axis=[
            (0.0, 2*ly),
            (-3/2*lx, ly),
            (-3/2*lx, -ly),
            (0.0, -2*ly),
            (3/2*lx, -ly),
            (3/2*lx, ly),
        ],
        angles=[0.0, 60.0, 120.0, 0.0, 60.0, 120.0],
        dimensions=(lx, ly),
        bd_type=[BoundaryType.TRANSLATION]*6
    )


def build_bd_sixth_hex(
        lx: float, ly: float, type_geo: LatticeGeometryType) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing a sixth
    symmetry of a regular hexagon. The information related to the borders of
    the triangular portion of the hexagon, that is contained in the
    `BoundaryData` object, derives from the characteristic dimensions of the
    full hexagon.
    The geometric data of the triangular shape is provided in
    counter-clockwise order starting from its bottom left corner that
    coincides with the XYZ origin.
    The type of boundary condition assigned to each edge depends on the
    `type_geo` parameter.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the hexagon (i.e. its edge length).
    ly : float
        Characteristic Y-dimension of the hexagon (i.e. its apothem length).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the three vertices of a sixth of
        a full hexagon, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned to each border depending on the `type_geo` parameter.
    """
    if type_geo in [LatticeGeometryType.SA60,
                    LatticeGeometryType.SYMMETRIES_TWO]:
        bd_type = [BoundaryType.AXIAL_SYMMETRY]*3
    elif type_geo in [LatticeGeometryType.RA60, LatticeGeometryType.ROTATION]:
        bd_type = [
            BoundaryType.ROTATION,
            BoundaryType.TRANSLATION,
            BoundaryType.ROTATION
        ]
    return BoundaryData(
        vertices=[
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((lx, 0.0, 0.0)),
            make_vertex((lx/2, ly, 0.0))
        ],
        axis=[(0.0, 0.0), (lx, 0.0), (0.0, 0.0)],
        angles=[0.0, 120.0, 60.0],
        dimensions=(lx, ly),
        bd_type=bd_type
    )


def build_bd_third_hex(lx: float, ly: float) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing a third
    symmetry of a regular hexagon. The information related to the borders of
    the quadrilateral portion of the hexagon, that is contained in the
    `BoundaryData` object, derives from the characteristic dimensions of the
    full hexagon.
    The geometric data of the quadrilateral shape is provided in
    counter-clockwise order starting from its bottom left corner that
    coincides with the XYZ origin.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the hexagon (i.e. its edge length).
    ly : float
        Characteristic Y-dimension of the hexagon (i.e. its apothem length).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the four vertices of a third of
        a full hexagon, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        being either `TRANSLATION` or `ROTATION`.
    """
    return BoundaryData(
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


def build_bd_twelfth_hex(lx: float, ly: float) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing a twelfth
    symmetry of a regular hexagon. The information related to the borders of
    the triangular portion of the hexagon, that is contained in the
    `BoundaryData` object, derives from the characteristic dimensions of the
    full hexagon.
    The geometric data of the triangular shape is provided in
    counter-clockwise order starting from its bottom left corner that
    coincides with the XYZ origin.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the hexagon (i.e. its edge length).
    ly : float
        Characteristic Y-dimension of the hexagon (i.e. its apothem length).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the three vertices of a twelfth of
        a full hexagon, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        all being `AXIAL_SYMMETRY`.
    """
    return BoundaryData(
        vertices=[
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((lx, 0.0, 0.0)),
            make_vertex((sqrt(3)/2*ly, ly/2, 0.0))
        ],
        axis=[(0.0, 0.0), (lx, 0.0), (0.0, 0.0)],
        angles=[0.0, 120.0, 30.0],
        dimensions=(lx, ly),
        bd_type=[BoundaryType.AXIAL_SYMMETRY]*6
    )


def build_bd_full_rect(lx: float, ly: float) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing a full
    rectangle placed so that its bottom left corner coincides with the
    XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryData` object, derives from the characteristic dimensions
    of the rectangle.
    The geometric data of the rectangle is provided in counter-clockwise
    order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the four vertices of the rectangle,
        the XY directions of the borders' axes and the corresponding angles
        (as needed by DRAGON), and the type of boundary, assigned as
        `TRANSLATION` for all the borders.
    """
    return BoundaryData(
        vertices=[
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((lx, 0.0, 0.0)),
            make_vertex((lx, ly, 0.0)),
            make_vertex((0.0, ly, 0.0))
        ],
        axis=[(0.0, ly), (-lx, 0.0), (0.0, -ly), (lx, 0.0)],
        angles=[0.0, 90.0, 0.0, 90.0],
        dimensions=(lx, ly),
        bd_type=[BoundaryType.TRANSLATION]*4
    )


def build_bd_half_rect(lx: float, ly: float) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing a half
    symmetry of a rectangle. The vertices are defined so that the shape's
    bottom left corner coincides with the XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryData` object, derives from the characteristic dimensions of
    the rectangle.
    The geometric data of the rectangular shape is provided in
    counter-clockwise order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the four vertices of a half of
        the rectangle, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned as `AXIAL_SYMMETRY` for all the borders.
    """
    return BoundaryData(
        vertices=[
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((lx/2, 0.0, 0.0)),
            make_vertex((lx/2, ly, 0.0)),
            make_vertex((0.0, ly, 0.0))
        ],
        axis=[(0.0, 0.0), (lx/2, 0.0), (0.0, ly), (0.0, 0.0)],
        angles=[0.0, 90.0, 0.0, 90.0],
        dimensions=(lx, ly),
        bd_type=[BoundaryType.AXIAL_SYMMETRY]*4
    )


def build_bd_quarter_rect(lx: float, ly: float) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing a quarter
    symmetry of a rectangle. The vertices are defined so that the shape's
    bottom left corner coincides with the XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryData` object, derives from the characteristic dimensions of
    the rectangle.
    The geometric data of the rectangular shape is provided in
    counter-clockwise order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the four vertices of a quarter of
        the rectangle, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned as `AXIAL_SYMMETRY` for all the borders.
    """
    return BoundaryData(
        vertices=[
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((lx/2, 0.0, 0.0)),
            make_vertex((lx/2, ly/2, 0.0)),
            make_vertex((0.0, ly/2, 0.0))
        ],
        axis=[(0.0, 0.0), (lx/2, 0.0), (0.0, ly/2), (0.0, 0.0)],
        angles=[0.0, 90.0, 0.0, 90.0],
        dimensions=(lx, ly),
        bd_type=[BoundaryType.AXIAL_SYMMETRY]*4
    )


def build_bd_eighth_rect(lx: float, ly: float) -> BoundaryData:
    """
    Function that creates a `BoundaryData` object representing an eighth
    symmetry of a rectangle. The vertices are defined so that the shape's
    bottom left corner coincides with the XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryData` object, derives from the characteristic dimensions of
    the rectangle.
    The geometric data of the triangular shape is provided in
    counter-clockwise order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryData
        A `BoundaryData` object built from the three vertices of an eighth of
        the rectangle, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned as `AXIAL_SYMMETRY` for all the borders.
    """
    return BoundaryData(
        vertices=[
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((lx/2, 0.0, 0.0)),
            make_vertex((lx/2, ly/2, 0.0))
        ],
        axis=[(0.0, 0.0), (lx/2, 0.0), (0.0, 0.0)],
        angles=[0.0, 90.0, degrees(atan2(ly, lx))],
        dimensions=(lx, ly),
        bd_type=[BoundaryType.AXIAL_SYMMETRY]*3
    )


def build_boundary_data(
        dimensions: Tuple[float, float],
        cell_type: CellType,
        symm_type: SymmetryType,
        type_geo: LatticeGeometryType) -> BoundaryData:
    """
    Function that constructs a `BoundaryData` object providing the boundary
    characteristics to use as a reference for test purposes.
    The instance is built depending on the type of cells, the applied symmetry
    and the corresponding lattice type of geometry.


    Parameters
    ----------
    dimensions : Tuple[float, float]
        The lattice X-Y characteristic dimensions.
    cell_type : CellType
        The type of cells in the lattice.
    symm_type : SymmetryType
        The type of symmetry applied to the lattice; it drives the selection
        of the builder function to generate the `BoundaryData` instance.
    type_geo : LatticeGeometryType
        The lattice geometry type, used to further specialize the
        `BoundaryData` instance in the case of a sixth symmetry.

    Returns
    -------
    BoundaryData
        A boundary representation including vertex positions, edges, axis
        directions, corresponding angles, and associated BC types.
    """
    if cell_type == CellType.HEX:
        match symm_type:
            case SymmetryType.FULL:
                return build_bd_full_hex(*dimensions)
            case SymmetryType.SIXTH:
                return build_bd_sixth_hex(*dimensions, type_geo)
            case SymmetryType.THIRD:
                return build_bd_third_hex(*dimensions)
            case SymmetryType.TWELFTH:
                return build_bd_twelfth_hex(*dimensions)
    else:
        match symm_type:
            case SymmetryType.FULL:
                return build_bd_full_rect(*dimensions)
            case SymmetryType.HALF:
                return build_bd_half_rect(*dimensions)
            case SymmetryType.QUARTER:
                return build_bd_quarter_rect(*dimensions)
            case SymmetryType.EIGHTH:
                return build_bd_eighth_rect(*dimensions)


def build_cell_ref_vectors(cell: Cell) -> List[Any]:
    """
    Function that builds a list of vector objects on the first edge of
    each of the face objects belonging to a `Cell` instance, i.e. the
    whole cell's face and the faces stored in the dictionaries of
    properties and sectorization options.

    Parameters
    ----------
    cell : Cell
        The cell to build reference vectors for.

    Returns
    -------
    List[Any]
        A list of vectors built of the first edge of the faces in the
        cell.
    """
    # List all the face objects comprising the cell's face and those
    # stored in the cell's dictionaries
    faces = [cell.face] + list(cell.tech_geom_props.keys())
    if cell.sectorized_face:
        faces += list(cell.tech_geom_sect_opts.keys())
    # Store the edge '0' for each face
    edges = [extract_sub_shapes(f, ShapeType.EDGE)[0] for f in faces]
    # Append the edge '0' for the cell's 'Surface' and reference circle
    edges.append(cell.figure.borders[0])
    edges.append(cell.figure.out_circle)
    # Build reference vectors for each edge
    return [
        make_vector_from_points(
            cell.figure.o,
            make_vertex_on_curve(e, 0.0)) for e in edges
    ]


def build_colorset(lattice: Lattice) -> List[Lattice]:
    """
    Function that facilitates the construction of a colorset as a list of
    ``Lattice`` instances.

    Parameters
    ----------
    lattice : Lattice
        The reference lattice.

    Returns
    -------
    List[Lattice]
        A list of ``Lattice`` instances each positioned appropriately to
        replicate a colorset.
    """
    # Deepcopy the given lattice
    lattice = deepcopy(lattice)
    # Build a ring of cells around the only cell present in the given
    # lattice and add a box
    lattice.add_ring_of_cells(lattice.lattice_cells[0], 1, 0)
    lattice.build_lattice_box([0.1])
    lattice.set_lattice_box_properties(
        {PropertyType.MATERIAL: ["MAT1"]})
    # Build the regions
    lattice.build_regions()
    # Build the positions of the lattices in the colorset
    lattices: List[Lattice] = [lattice]
    pos = [
        (lattice.lattice_box.figure.lx, 0.0, 0.0),
        (lattice.lattice_box.figure.lx,
            lattice.lattice_box.figure.ly,
            0.0),
        (0.0, lattice.lattice_box.figure.ly, 0.0),
        (-lattice.lattice_box.figure.lx,
            lattice.lattice_box.figure.ly,
            0.0),
        (-lattice.lattice_box.figure.lx, 0.0, 0.0),
        (-lattice.lattice_box.figure.lx,
            -lattice.lattice_box.figure.ly,
            0.0),
        (0.0, -lattice.lattice_box.figure.ly, 0.0),
        (lattice.lattice_box.figure.lx,
            -lattice.lattice_box.figure.ly,
            0.0)
    ]
    # Copy and translate the given lattice to the calculated positions
    for xyz in pos:
        l = deepcopy(lattice)
        l.translate(xyz)
        lattices.append(l)
    return lattices


def build_hex_geom_elements(
        center: Any, edge_length: float) -> Tuple[List[Any], List[Any]]:
    """
    Function that, given the center and the length of the hexagon edge,
    builds the vertex and edge objects that represent a hexagon.

    Parameters
    ----------
    center : Any
        A vertex object representing the center of the resulting hexagonal
        shape.
    edge_length : float
        The length of the hexagon's edge.

    Returns
    -------
    Tuple[List[Any], List[Any]]
        A tuple with the list of vertices and edges of the hexagon.
    """
    vertices = [
        make_vertex_on_curve(
            make_circle(center, None, edge_length), i/6) for i in range(6)
    ]
    edges = [
            make_edge(vertices[i], vertices[(i+1) % 6]) for i in range(6)
    ]
    return vertices, edges


def build_lattice_ref_vectors(lattice: Lattice) -> List[Any]:
    """
    Function that builds a list of vector objects on the first edge of
    each of the face objects belonging to a `Lattice` instance, i.e. the
    whole lattice compound and the one representing the symmetry. The same
    goes for the cells in the lattice.

    Parameters
    ----------
    lattice : Lattice
        The lattice to build reference vectors for.

    Returns
    -------
    List[Any]
        A list of vectors built of the first edge of the elements of the
        lattice and its cells.
    """
    # List all the face objects comprising the lattice's compound and its
    # cells' faces
    faces = [lattice.lattice_cmpd]
    if lattice.lattice_box:
        faces.append(lattice.lattice_box.face)
    if lattice.lattice_symm:
        faces.append(lattice.lattice_symm)
    for cell in lattice.lattice_cells:
        faces.append(cell.face)
    # Store the edge '0' for each face
    edges = [extract_sub_shapes(f, ShapeType.EDGE)[0] for f in faces]
    # Build reference vectors for each edge
    return [
        make_vector_from_points(
            lattice.lattice_center,
            make_vertex_on_curve(e, 0.0)) for e in edges
    ]


def compute_hash(file_path: Path) -> str:
    """
    Function to compute the SHA256 hash of a file.

    Parameters
    ----------
    file_path : Path
        The ``Path`` object of the file to process.

    Returns
    -------
    str
        The SHA256 hash of a given file.
    """
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256.update(block)
    return sha256.hexdigest()


def make_ref_vectors(surf: Surface) -> Tuple[Any, Any]:
    """
    Function that returns two vector objects built on the `Surface` first
    border element and on its corresponding circle the surface is inscribed
    into.

    Parameters
    ----------
    surf : Surface
        The `Surface` object used to build the reference vectors.

    Returns
    -------
    Tuple[Any, Any]
        Tuple providing two reference vector objects, the first built
        on the surface border, the second on the circle the surface is
        inscribed into.
    """
    face_ref_vect = make_vector_from_points(
        surf.o, make_vertex_on_curve(surf.borders[0], 0.0))
    out_circle_ref_vect = make_vector_from_points(
        surf.o, make_vertex_on_curve(surf.out_circle, 0.0))
    return face_ref_vect, out_circle_ref_vect


def set_up_hex_cells(hex_cell: HexCell) -> List[HexCell]:
    """
    Function that builds a list of hexagonal cells with a central cell
    surrounded by six other cells.

    Parameters
    ----------
    hex_cell : HexCell
        The `HexCell` object representing an hexagonal cell.

    Returns
    -------
    List[HexCell]
        A list of hexagonal cells with a central one surrounded by six
        cells.
    """
    dx = hex_cell.apothem
    dy = 3/2*hex_cell.edge_length
    hex_cell.rotate(90.0)
    return [
        hex_cell,
        hex_cell.translate((dx, dy, 0)),
        hex_cell.translate((-dx, dy, 0)),
        hex_cell.translate((-dx, -dy, 0)),
        hex_cell.translate((dx, -dy, 0)),
        hex_cell.translate((2*dx, 0, 0)),
        hex_cell.translate((-2*dx, 0, 0))
    ]

def set_up_rect_cells(rect_cell: RectCell,
                      is_even: bool = False) -> List[RectCell]:
    """
    Function that builds a list of cartesian cells. Depending on the
    boolean flag `is_even`, the resulting list is made by a central
    cell surrounded by eight other cells, if `False`, or by cells
    without any central one replicating a pattern with an even number
    of cells.

    Parameters
    ----------
    rect_cell : RectCell
        The `RectCell` object representing a cartesian cell.
    is_even : bool
        Boolean flag indicating the kind of pattern of cells (either with
        an odd or even number of cells).

    Returns
    -------
    List[RectCell]
        A list of cartesian cells with a specific pattern.
    """
    dx = rect_cell.width
    dy = rect_cell.height
    if is_even:
        dx /= 2
        dy /= 2
        return [
            rect_cell.translate((dx, dy, 0)),
            rect_cell.translate((-dx, dy, 0)),
            rect_cell.translate((-dx, -dy, 0)),
            rect_cell.translate((dx, -dy, 0)),
            rect_cell.translate((2*dx, 0, 0)),
            rect_cell.translate((2*dx, dy, 0)),
            rect_cell.translate((2*dx, 2*dy, 0)),
            rect_cell.translate((dx, 2*dy, 0)),
            rect_cell.translate((-dx, 2*dy, 0)),
            rect_cell.translate((-2*dx, 2*dy, 0)),
            rect_cell.translate((-2*dx, dy, 0)),
            rect_cell.translate((-2*dx, 0, 0)),
            rect_cell.translate((-2*dx, -dy, 0)),
            rect_cell.translate((-2*dx, -2*dy, 0)),
            rect_cell.translate((-dx, -2*dy, 0)),
            rect_cell.translate((dx, -2*dy, 0)),
            rect_cell.translate((2*dx, -2*dy, 0)),
            rect_cell.translate((2*dx, -dy, 0)),
        ]
    return [
        rect_cell,
        rect_cell.translate((dx, 0, 0)),
        rect_cell.translate((dx, dy, 0)),
        rect_cell.translate((0, dy, 0)),
        rect_cell.translate((-dx, dy, 0)),
        rect_cell.translate((-dx, 0, 0)),
        rect_cell.translate((-dx, -dy, 0)),
        rect_cell.translate((0, -dy, 0)),
        rect_cell.translate((dx, -dy, 0))
    ]
