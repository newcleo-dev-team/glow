"""
Module containing classes providing an interface towards the topological
entities of the GEOM module of SALOME.
"""
from abc import ABC
from typing import Any, List, Self, Sequence, Type

# Import SALOME element
from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_point_coordinates, get_shape_name, get_shape_type, make_common, \
    make_cut, make_fuse, make_partition
from glow.support.utility import check_shape_expected_types


class GeomWrapper(ABC):
    """
    Abstract class acting as a wrapper for low-level GEOM topological shape
    objects.
    It provides:
        - type validation of the provided GEOM object;
        - attribute delegation, so that this class behaves like the provided
          GEOM object when using any GEOM function;
        - geometric operator overloads. Math operators are overloaded to
          provide boolean operations between two ``GeomWrapper`` objects.

    Parameters
    ----------
    geom_obj : Any | None
        The underlying GEOM_Object (or None). If provided, its shape type is
        checked against ``expected_types``.
    expected_types : List[ShapeType]
        List of allowed shape types to check the provided GEOM_Object against.
        If the actual type of ``geom_obj`` is not in this list, initialization
        raises a ``ValueError`` exception.

    Notes
    -----
    - Attribute access is delegated to the underlying GEOM object via
      ``__getattr__``. During initialization this delegation is guarded to
      avoid recursion if the GEOM object has not yet been stored.
    - The wrapper implements a few binary operator overloads that produce
      new wrapped shapes:
      - ``+`` (``__add__``) performs a fuse (union) of two shapes.
      - ``-`` (``__sub__``) performs a cut (difference).
      - ``*`` (``__mul__``) produces the common (intersection).
      - ``/`` (``__truediv__``) perform a partition operation where the
        GEOM object is subdivided in sub shapes by given shape objects.
      - ``//`` (``__floordiv__``) perform a partition operation where the
        resulting shape is the combination of all the shapes after being
        partitioned.
      These operators call helper functions that interface with the
      corresponding GEOM ones, and then wrap the resulting GEOM_Object in an
      object of the subclasses of ``GeomWrapper``.

    Raises
    ------
    ValueError
        If ``geom_obj`` is not None and its shape type is not contained among
        the expected ones.
    """
    def __init__(
            self, geom_obj: Any | None, expected_types: List[ShapeType]
        ) -> None:
        # Check if the type of the GEOM object is valid
        if geom_obj is not None:
            check_shape_expected_types(geom_obj, expected_types)
        # Store the types for successive validation
        self.expected_types = expected_types
        # Directly modify the internal __dict__ to avoid triggering
        # __setattr__ which would result in an infinite recursion
        self.__dict__["_geom_obj"] = geom_obj

    @property
    def geom_obj(self) -> Any:
        """
        Access or set the corresponding GEOM_Object.

        Parameters
        ----------
        geom_obj : Any
            The GEOM_Object this instance should be set to.
        expected_types : List[ShapeType]
            List of allowed shape types to check the provided GEOM_Object
            against.

        Returns
        -------
        Any
            The GEOM_Object this instance is based on.

        Raises
        ------
        ValueError
            If the type of the provided GEOM_Object does not match with the
            expected one.
        """
        return self._geom_obj

    @geom_obj.setter
    def geom_obj(self, geom_obj: Any) -> None:
        # Check if the type of the GEOM object is valid
        if geom_obj is not None:
            check_shape_expected_types(geom_obj, self.expected_types)
        self._geom_obj = geom_obj

    def __getattr__(self, name: str) -> Any:
        # Avoid recursion if _geom_obj is not yet set
        if "_geom_obj" not in self.__dict__:
            raise AttributeError(
                f"Error while checking for {name} attribute in "
                f"'<{self.__class__.__name__}>' object: the attribute "
                "'_geom_obj' has not been set yet.")
        return getattr(self.__dict__["_geom_obj"], name)

    def __repr__(self) -> str:
        """
        Outputs the class name followed by the ``name`` attribute of the
        wrapped GEOM object.

        Returns
        -------
        str
            Descriptive string containing the class name followed by the
            ``name`` attribute of the wrapped GEOM object.
        """
        try:
            name = get_shape_name(self.geom_obj)
        except Exception:
            name = "Unnamed"
        return f"<{self.__class__.__name__}: {name}>"

    def __add__(self, other: Self | Sequence[Self]) -> Self:
        """
        Returns an instance of the ``GeomWrapper`` subclasses resulting from
        the fusion of this shape and the given ones.

        Parameters
        ----------
        other : Self | Sequence[Self]
            One or more ``GeomWrapper`` objects to fuse the current object
            with.

        Returns
        -------
        Self
            An instance of the ``GeomWrapper`` subclasses representing the
            fused shape.
        """
        return wrap_shape(
            make_fuse(
                [self._geom_obj] + [
                    o for o in (
                        other if isinstance(other, Sequence) else [other])
                ]
            )
        )

    def __sub__(self, other: Self) -> Self:
        """
        Returns an instance of the ``GeomWrapper`` subclasses resulting from
        cutting this shape with the given one.

        Parameters
        ----------
        other : Self
            A ``GeomWrapper`` object to cut this one with.

        Returns
        -------
        Self
            An instance of the ``GeomWrapper`` subclasses representing the
            cut shape.
        """
        return wrap_shape(make_cut(self._geom_obj, other._geom_obj))

    def __mul__(self, other: Self) -> Self:
        """
        Returns an instance of the ``GeomWrapper`` subclasses resulting as
        the common part between this shape and the given one.

        Parameters
        ----------
        other : Self
            A ``GeomWrapper`` object to extract the common part.

        Returns
        -------
        Self
            An instance of the ``GeomWrapper`` subclasses representing the
            common part between the current shape and the given one.
        """
        return wrap_shape(make_common(self._geom_obj, other._geom_obj))

    def __truediv__(self, other: Self | Sequence[Self]) -> Self:
        """
        Returns an instance of the ``GeomWrapper`` subclasses resulting from
        partitioning this shape with the given ones. The resulting shape
        includes the current one subdivided by the given shapes.

        Parameters
        ----------
        other : Self | Sequence[Self]
            One or more ``GeomWrapper`` objects to partition the current
            object with.

        Returns
        -------
        Self
            An instance of the ``GeomWrapper`` subclasses representing the
            current partitioned shape.
        """
        return wrap_shape(
            make_partition(
                [self._geom_obj],
                [o for o in (
                    other if isinstance(other, Sequence) else [other])],
                ShapeType.COMPOUND
            )
        )

    def __floordiv__(self, other: Self | Sequence[Self]) -> Self:
        """
        Returns an instance of the ``GeomWrapper`` subclasses resulting from
        partitioning this shape with the given ones. The resulting shape
        includes every shape involved in the operation.

        Parameters
        ----------
        other : Self | Sequence[Self]
            One or more ``GeomWrapper`` objects to partition the current
            object with.

        Returns
        -------
        Self
            An instance of the ``GeomWrapper`` subclasses representing the
            combination of all the shapes after being partitioned.
        """
        return wrap_shape(
            make_partition(
                [self._geom_obj] + \
                [o for o in (
                    other if isinstance(other, Sequence) else [other])],
                [],
                ShapeType.COMPOUND
            )
        )


class Compound(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM compound object.
    """
    def __init__(self, geom_obj: Any):
        super().__init__(geom_obj, [ShapeType.COMPOUND])

    def __repr__(self):
        faces = [
            wrap_shape(f) for f in extract_sub_shapes(self, ShapeType.FACE)]
        edges = [
            wrap_shape(e) for e in extract_sub_shapes(self, ShapeType.EDGE)]
        # Collect the information about its faces and edges
        faces_str = "\n   ".join(repr(f) for f in faces)
        edge_str = "\n   ".join(repr(e) for e in edges)
        return super().__repr__() + "\n   " + faces_str + "\n   " + edge_str


class Vertex(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM vertex object.
    """
    def __init__(self, geom_obj: Any):
        super().__init__(geom_obj, [ShapeType.VERTEX])

    def __repr__(self):
        return (
            super().__repr__()
            + f" - (x, y, z) = {get_point_coordinates(self)}"
        )

class Edge(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM edge object.
    """
    def __init__(self, geom_obj: Any):
        super().__init__(geom_obj, [ShapeType.EDGE])

    def __repr__(self):
        vertices = [
            wrap_shape(v) for v in extract_sub_shapes(self, ShapeType.VERTEX)]
        # Collect the information about its vertices
        vertices_str = "\n   ".join(repr(v) for v in vertices)
        return super().__repr__() + "\n   " + vertices_str

class Face(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM face object.
    """
    def __init__(self, geom_obj: Any):
        super().__init__(geom_obj, [ShapeType.FACE])

    def __repr__(self):
        edges = [
            wrap_shape(e) for e in extract_sub_shapes(self, ShapeType.EDGE)]
        # Collect the information about its edges
        edge_str = "\n   ".join(repr(e) for e in edges)
        return super().__repr__() + "\n   " + edge_str


def wrap_shape(obj: Any) -> GeomWrapper:
    """
    Function that wraps a generic GEOM object in the appropriate wrapper
    class.

    Parameters
    ----------
    geom_obj : Any
        The generic GEOM object to wrap.

    Returns
    -------
    GeomWrapper
        An instance of a ``GeomWrapper`` subclass built from the given
        GEOM object.

    Raises
    ------
    TypeError
        If the given GEOM object is not a valid GEOM shape.
    ValueError
        If the given GEOM object does not have a corresponding wrapper class.
    """
    # Check if the given object is a valid GEOM shape
    try:
        shape_type = get_shape_type(obj)
    except Exception as e:
        raise TypeError("Provided object is not a valid GEOM shape.") from e
    # Map shape types to wrapper classes
    wrapper_map: dict[ShapeType, Type[GeomWrapper]] = {
        ShapeType.VERTEX: Vertex,
        ShapeType.EDGE: Edge,
        ShapeType.FACE: Face,
        ShapeType.COMPOUND: Compound,
    }
    # Get the wrapper class, if any
    wrapper_cls = wrapper_map.get(shape_type)
    if wrapper_cls is None:
        raise ValueError(f"No wrapper defined for shape type '{shape_type}'")
    # Return the instantiated wrapper class
    return wrapper_cls(obj)
