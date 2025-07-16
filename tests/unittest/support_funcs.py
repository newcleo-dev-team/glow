"""
Module declaring functions to support the execution of the unit tests.
"""
from typing import Any, List, Tuple

from glow.geometry_layouts.cells import Cell
from glow.geometry_layouts.geometries import Surface
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, make_circle, make_edge, \
    make_vector_from_points, make_vertex_on_curve


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
