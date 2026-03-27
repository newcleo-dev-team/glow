"""
Module declaring functions to support the execution of the unit tests.
"""
import io
import hashlib
import sys

from dataclasses import dataclass, field
from math import degrees, sqrt, atan2
from pathlib import Path
from typing import Any, Callable, List, Tuple

from glow.geometry_layouts.cells import Cell, HexCell, CartesianCell
from glow.geometry_layouts.geometries import Surface
from glow.geometry_layouts.lattices import CartesianLattice, Lattice
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_point_coordinates,  make_circle, make_edge, make_vector_from_points, \
    make_vertex, make_vertex_on_curve
from glow.support.types import BoundaryType, LayoutType, LayoutGeometryType, \
    PropertyType, SymmetryType
from glow.support.utility import build_contiguous_edges


@dataclass
class BoundaryInfo():
    """
    Dataclass storing geometric information about the boundaries of a
    geometry layout in terms of vertices, contiguous edges, axes, angles
    and characteristic dimensions of the layout.

    Attributes
    ----------
    vertices : List[Any]
        List of vertex objects defining the layout's borders.
    edges : List[Any]
        List of edge objects automatically constructed from the vertices.
    axis : List[Tuple[float, float]]
        Providing the XY directions representing the axes for each border.
    angles : List[float]
        List of angles (in degrees) associated with each axis.
    dimensions : tuple of float
        Characteristic dimensions of the layout.
    bd_type : List[BoundaryType]
        Type associated to each layout's border provided as element of the
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


def build_bd_full_hex(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a regular
    hexagon centered at the origin. The information related to the borders,
    that is contained in the `BoundaryInfo` object, derives from the
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
    BoundaryInfo
        A `BoundaryInfo` object built from the six vertices of the full
        hexagon, the XY directions of the borders' axes and the corresponding
        angles (as needed by DRAGON), and the type of boundary, assigned as
        `TRANSLATION` for all the borders.
    """
    return BoundaryInfo(
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
        lx: float, ly: float, type_geo: LayoutGeometryType
    ) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a sixth
    symmetry of a regular hexagon. The information related to the borders of
    the triangular portion of the hexagon, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions of the
    full hexagon.
    The geometric data of the triangular shape is provided in counterclockwise
    order starting from its bottom left corner that coincides with the XYZ
    origin.
    The type of boundary condition assigned to each edge depends on the
    `type_geo` parameter.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the hexagon (i.e. its edge length).
    ly : float
        Characteristic Y-dimension of the hexagon (i.e. its apothem length).
    type_geo : LayoutGeometryType
        The type of geometry of the layout, according to DRAGON5.

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the three vertices of a sixth of
        a full hexagon, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned to each border depending on the `type_geo` parameter.
    """
    if type_geo in [LayoutGeometryType.SA60,
                    LayoutGeometryType.SYMMETRIES_TWO]:
        bd_type = [BoundaryType.AXIAL_SYMMETRY]*3
    elif type_geo in [LayoutGeometryType.RA60, LayoutGeometryType.ROTATION]:
        bd_type = [
            BoundaryType.ROTATION,
            BoundaryType.TRANSLATION,
            BoundaryType.ROTATION
        ]
    return BoundaryInfo(
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


def build_bd_third_hex(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a third
    symmetry of a regular hexagon. The information related to the borders of
    the quadrilateral portion of the hexagon, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions of the
    full hexagon.
    The geometric data of the quadrilateral shape is provided in
    counterclockwise order starting from its bottom left corner that coincides
    with the XYZ origin.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the hexagon (i.e. its edge length).
    ly : float
        Characteristic Y-dimension of the hexagon (i.e. its apothem length).

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the four vertices of a third of
        a full hexagon, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        being either `TRANSLATION` or `ROTATION`.
    """
    return BoundaryInfo(
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


def build_bd_twelfth_hex(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a twelfth
    symmetry of a regular hexagon. The information related to the borders of
    the triangular portion of the hexagon, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions of the
    full hexagon.
    The geometric data of the triangular shape is provided in counterclockwise
    order starting from its bottom left corner that coincides with the XYZ
    origin.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the hexagon (i.e. its edge length).
    ly : float
        Characteristic Y-dimension of the hexagon (i.e. its apothem length).

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the three vertices of a twelfth of
        a full hexagon, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        all being `AXIAL_SYMMETRY`.
    """
    return BoundaryInfo(
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


def build_bd_full_rect(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a full
    rectangle placed so that its bottom left corner coincides with the
    XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions
    of the rectangle.
    The geometric data of the rectangle is provided in counterclockwise
    order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the four vertices of the rectangle,
        the XY directions of the borders' axes and the corresponding angles
        (as needed by DRAGON), and the type of boundary, assigned as
        `TRANSLATION` for all the borders.
    """
    return BoundaryInfo(
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


def build_bd_diag_rect(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a half
    symmetry of a rectangle along its diagonal. The vertices are defined so
    that the shape's bottom left corner coincides with the XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions of
    the rectangle.
    The geometric data of the triangular shape is provided in counterclockwise
    order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the three vertices of a half of
        the rectangle along its diagonal, the XY directions of the borders'
        axes and the corresponding angles (as needed by DRAGON), and the type
        of boundary, assigned as `AXIAL_SYMMETRY` for all the borders.
    """
    return BoundaryInfo(
        vertices=[
            make_vertex((0.0, 0.0, 0.0)),
            make_vertex((lx, 0.0, 0.0)),
            make_vertex((lx, ly, 0.0))
        ],
        axis=[(0.0, 0.0), (lx, 0.0), (0.0, 0.0)],
        angles=[0.0, 90.0, degrees(atan2(ly, lx))],
        dimensions=(lx, ly),
        bd_type=[BoundaryType.AXIAL_SYMMETRY]*3
    )


def build_bd_half_rect(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a half
    symmetry of a rectangle. The vertices are defined so that the shape's
    bottom left corner coincides with the XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions of
    the rectangle.
    The geometric data of the rectangular shape is provided in
    counterclockwise order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the four vertices of a half of
        the rectangle, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned as `AXIAL_SYMMETRY` for all the borders.
    """
    return BoundaryInfo(
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


def build_bd_quarter_rect(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing a quarter
    symmetry of a rectangle. The vertices are defined so that the shape's
    bottom left corner coincides with the XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions of
    the rectangle.
    The geometric data of the rectangular shape is provided in
    counterclockwise order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the four vertices of a quarter of
        the rectangle, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned as `AXIAL_SYMMETRY` for all the borders.
    """
    return BoundaryInfo(
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


def build_bd_eighth_rect(lx: float, ly: float) -> BoundaryInfo:
    """
    Function that creates a `BoundaryInfo` object representing an eighth
    symmetry of a rectangle. The vertices are defined so that the shape's
    bottom left corner coincides with the XYZ space origin.
    The information related to the borders, that is contained in the
    `BoundaryInfo` object, derives from the characteristic dimensions of
    the rectangle.
    The geometric data of the triangular shape is provided in counterclockwise
    order starting from its bottom left corner.

    Parameters
    ----------
    lx : float
        Characteristic X-dimension of the rectangle (i.e. its width).
    ly : float
        Characteristic Y-dimension of the rectangle (i.e. its height).

    Returns
    -------
    BoundaryInfo
        A `BoundaryInfo` object built from the three vertices of an eighth of
        the rectangle, the XY directions of the borders' axes and the
        corresponding angles (as needed by DRAGON), and the type of boundary,
        assigned as `AXIAL_SYMMETRY` for all the borders.
    """
    return BoundaryInfo(
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


def build_boundary_info(
        dimensions: Tuple[float, float],
        layout_type: LayoutType,
        symm_type: SymmetryType,
        type_geo: LayoutGeometryType
    ) -> BoundaryInfo:
    """
    Function that constructs a `BoundaryInfo` object providing the boundary
    characteristics to use as a reference for test purposes.
    The instance is built depending on the type of layout, the applied
    symmetry and the corresponding type of geometry.

    Parameters
    ----------
    dimensions : Tuple[float, float]
        The layout X-Y characteristic dimensions.
    layout_type : LayoutType
        The type of layout, depending on its characteristic shape.
    symm_type : SymmetryType
        The type of symmetry applied to the layout; it drives the selection
        of the builder function to generate the `BoundaryInfo` instance.
    type_geo : LayoutGeometryType
        The layout geometry type, used to further specialize the
        `BoundaryInfo` instance in the case of a sixth symmetry.

    Returns
    -------
    BoundaryInfo
        A boundary representation including vertex positions, edges, axis
        directions, corresponding angles, and associated BC types.
    """
    if layout_type == LayoutType.HEX:
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
            case SymmetryType.DIAG:
                return build_bd_diag_rect(*dimensions)
            case SymmetryType.QUARTER:
                return build_bd_quarter_rect(*dimensions)
            case SymmetryType.EIGHTH:
                return build_bd_eighth_rect(*dimensions)


def build_colorset(cell: Cell) -> CartesianLattice:
    """
    Function that facilitates the construction of a colorset as a
    `CartesianLattice` instance made of the given `Cell` objects.

    Parameters
    ----------
    cell : Cell
        The reference cell.

    Returns
    -------
    CartesianLattice
        A `CartesianLattice` instance where cells are positioned appropriately
        to replicate a colorset.
    """
    # Clone the given cell
    cell = cell.clone()
    # Build a lattice made by a ring of cells around the central one and
    # include it in another cell to replicate an assembly
    lattice = CartesianLattice([cell])
    lattice.add_ring_of_cells(cell, 1, 0)
    assembly = CartesianCell(
        width_height=tuple(xy + 0.25 for xy in lattice.dimensions),
        base_props={PropertyType.MATERIAL: "MAT1"}
    )
    assembly.add(lattice)

    # Assemble the colorset as a lattice of assemblies
    colorset = CartesianLattice([assembly])
    colorset.add_ring_of_cells(assembly, 1)
    colorset.update_hierarchical_structure()
    return colorset


def build_hex_geom_elements(center: Any, edge_length: float) -> List[Any]:
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
    List[Any]
        The list of edges of the hexagon.
    """
    vertices = [
        make_vertex_on_curve(
            make_circle(center, None, edge_length), i/6) for i in range(6)
    ]
    edges = [
            make_edge(vertices[i], vertices[(i+1) % 6]) for i in range(6)
    ]
    return edges


def capture_output(
        func: Callable[..., Any], *args: Any, **kwargs: Any
    ) -> str:
    """
    Function that captures and returns all standard output produced by a
    given `Callable` object.

    This function temporarily redirects `sys.stdout` to an in-memory
    buffer, executes the provided callable with the given positional and
    keyword arguments, restores the original stdout stream, and returns
    the captured text.

    Parameters
    ----------
    func : Callable[..., Any]
        The function or callable object whose printed output should be
        captured.
    *args : Any
        Positional arguments forwarded to `func`.
    **kwargs : Any
        Keyword arguments forwarded to `func`.

    Returns
    -------
    str
        A string containing everything written to standard output during
        the execution of `func`.
    """
    old_stdout = sys.stdout
    buffer = io.StringIO()
    sys.stdout = buffer
    try:
        func(*args, **kwargs)
    finally:
        sys.stdout = old_stdout
    return buffer.getvalue()


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


def make_ref_vectors(surf: Surface) -> Any:
    """
    Function that returns a vector objects built on the `Surface` first
    border element.

    Parameters
    ----------
    surf : Surface
        The `Surface` object used to build the reference vector.

    Returns
    -------
    Any
        The reference vector object built on the surface border.
    """
    face_ref_vect = make_vector_from_points(
        surf.o, make_vertex_on_curve(surf.borders[0], 0.0)
    )
    return face_ref_vect


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
    # Retrieve the cell's centre coordinates and calculate the X-Y shifts
    o_x, o_y, _ = get_point_coordinates(hex_cell.o)
    dx = hex_cell.dimensions[1]
    dy = 3/2*hex_cell.dimensions[0]
    hex_cell.rotate(90.0)
    # Build the list of positions of the cells
    centres = [
        (o_x + dx, o_y + dy, 0),
        (o_x - dx, o_y + dy, 0),
        (o_x - dx, o_y - dy, 0),
        (o_x + dx, o_y - dy, 0),
        (o_x + 2*dx, o_y, 0),
        (o_x - 2*dx, o_y, 0)
    ]
    # Populate the list of translated cells
    cells = [hex_cell]
    for centre in centres:
        cell = hex_cell.clone()
        cell.translate(centre)
        cells.append(cell)
    # Return the list of translated cells
    return cells


def set_up_rect_cells(
        rect_cell: CartesianCell, is_even: bool = False
    ) -> List[CartesianCell]:
    """
    Function that builds a list of Cartesian cells. Depending on the
    boolean flag `is_even`, the resulting list is made by a central
    cell surrounded by eight other cells, if `False`, or by cells
    without any central one replicating a pattern with an even number
    of cells.

    Parameters
    ----------
    rect_cell : CartesianCell
        The `CartesianCell` object representing a Cartesian cell.
    is_even : bool
        Boolean flag indicating the kind of pattern of cells (either with
        an odd or even number of cells).

    Returns
    -------
    List[CartesianCell]
        A list of Cartesian cells with a specific pattern.
    """
    # Retrieve the cell's centre coordinates and calculate the X-Y shifts
    o_x, o_y, _ = get_point_coordinates(rect_cell.o)
    dx, dy = rect_cell.dimensions
    cells = []
    if is_even:
        dx /= 2
        dy /= 2
        # Build the list of positions of the cells
        centres = [
            (o_x + dx, o_y + dy, 0),
            (o_x - dx, o_y + dy, 0),
            (o_x - dx, o_y - dy, 0),
            (o_x + dx, o_y - dy, 0),
            (o_x + 3*dx, o_y + dy, 0),
            (o_x + 3*dx, o_y + 3*dy, 0),
            (o_x + dx, o_y + 3*dy, 0),
            (o_x - dx, o_y + 3*dy, 0),
            (o_x - 3*dx, o_y + 3*dy, 0),
            (o_x - 3*dx, o_y + dy, 0),
            (o_x - 3*dx, o_y - dy, 0),
            (o_x - 3*dx, o_y - 3*dy, 0),
            (o_x - dx, o_y - 3*dy, 0),
            (o_x + dx, o_y - 3*dy, 0),
            (o_x + 3*dx, o_y - 3*dy, 0),
            (o_x + 3*dx, o_y - dy, 0)
        ]
    else:
        centres = [
            (o_x + dx, o_y, 0),
            (o_x + dx, o_y + dy, 0),
            (o_x, dy, o_y),
            (o_x - dx, o_y + dy, 0),
            (o_x - dx, o_y, 0),
            (o_x - dx, o_y - dy, 0),
            (o_x, o_y - dy, 0),
            (o_x + dx, o_y - dy, 0)
        ]
        cells.append(rect_cell)

    # Populate the list of translated cells
    for centre in centres:
        cell = rect_cell.clone()
        cell.translate(centre)
        cells.append(cell)
    # Return the list of translated cells
    return cells
