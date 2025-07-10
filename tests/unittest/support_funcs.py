"""
Module declaring functions to support the execution of the unit tests.
"""
from typing import Any, List, Tuple

from glow.geometry_layouts.geometries import Surface
from glow.interface.geom_interface import make_circle, make_edge, \
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