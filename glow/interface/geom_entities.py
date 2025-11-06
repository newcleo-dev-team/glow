"""
Module containing classes providing an interface towards the topological
entities of the GEOM module of SALOME.
"""
from abc import ABC
import math
from typing import Any, List, Self, Sequence, Type

from glow.interface.geom_interface import ShapeType, extract_sub_shapes, \
    get_point_coordinates, get_shape_name, get_shape_type, make_common, \
    make_cut, make_face, make_fuse, make_partition, set_shape_name
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
        The GEOM object to wrap (or ``None``). If provided, its shape type
        is checked against the expected types.
    expected_types : List[ShapeType]
        List of allowed shape types to check the provided GEOM object against.

    Attributes
    ----------
    geom_obj : Any | None
        The wrapped GEOM object (or ``None``).
    name : str
        The name of the wrapped GEOM object.

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
        super().__init__()
        # Check if the type of the GEOM object is valid
        if geom_obj is not None:
            check_shape_expected_types(geom_obj, expected_types)
        # Store the types for successive validation
        self.__expected_types = expected_types
        # Directly modify the internal __dict__ to avoid triggering
        # __setattr__ which would result in an infinite recursion
        self.__dict__["_geom_obj"] = geom_obj
        self.__name: str = ""

    @property
    def geom_obj(self) -> Any:
        """
        Access or set the corresponding GEOM object.

        Parameters
        ----------
        geom_obj : Any
            The GEOM object this instance should be set to.
        expected_types : List[ShapeType]
            List of allowed shape types to check the provided GEOM object
            against.

        Returns
        -------
        Any
            The GEOM object this instance is based on.

        Raises
        ------
        ValueError
            If the type of the provided GEOM object does not match with the
            expected one.
        """
        return self._geom_obj

    @geom_obj.setter
    def geom_obj(self, geom_obj: Any | None) -> None:
        # Check if the type of the GEOM object is valid
        if geom_obj is not None:
            check_shape_expected_types(geom_obj, self.__expected_types)
        self._geom_obj = geom_obj

    @property
    def name(self) -> str:
        """
        Get or set the name of the wrapped GEOM object.

        Parameters
        ----------
        new_name : str
            The name the wrapped GEOM object should be set to.

        Returns
        -------
        str
            The name assigned to the wrapped GEOM object.
        """
        return self.__name

    @name.setter
    def name(self, new_name: str) -> None:
        self.__name = new_name
        # Set the GEOM_Object name, if the GEOM_Object is defined
        if self.geom_obj:
            set_shape_name(self.geom_obj, new_name)

    def __add__(self, other: Self | Sequence[Self]) -> Self:
        """
        Return an instance of the ``GeomWrapper`` subclasses resulting from
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

    def __floordiv__(self, other: Self | Sequence[Self]) -> Self:
        """
        Return an instance of the ``GeomWrapper`` subclasses resulting from
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

    def __getattr__(self, name: str) -> Any:
        """
        Delegate attribute access to the wrapped GEOM object.

        This method is called called when an attribute is missing on the
        instance. It delegates lookup to the internal `_geom_obj` if
        initialized; otherwise raises an exception.

        Parameters
        ----------
        name : str
            The name of the attribute to access.

        Returns
        -------
        Any
            The value of the requested attribute from the wrapped GEOM object.

        Raises
        ------
        AttributeError
            If no GEOM object has still be wrapped or the indicated attribute
            does not exist on the GEOM object.
        """
        if "_geom_obj" not in self.__dict__:
            raise AttributeError(
                f"Error while checking for {name} attribute in "
                f"'<{self.__class__.__name__}>' object: the attribute "
                "'_geom_obj' has not been set yet.")
        return getattr(self.__dict__["_geom_obj"], name)

    def __mul__(self, other: Self) -> Self:
        """
        Return an instance of the ``GeomWrapper`` subclasses resulting as
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

    def __repr__(self) -> str:
        """
        Return a descriptive string of the ``GeomWrapper`` instance containing
        the class name followed by the ``name`` attribute of the wrapped GEOM
        object.

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

    def __sub__(self, other: Self) -> Self:
        """
        Return an instance of the ``GeomWrapper`` subclasses resulting from
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

    def __truediv__(self, other: Self | Sequence[Self]) -> Self:
        """
        Return an instance of the ``GeomWrapper`` subclasses resulting from
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


class Compound(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM compound object.

    Parameters
    ----------
    geom_obj : Any | None
        The GEOM compound object to wrap (or ``None``). If provided, its
        shape type is checked against the expected ``ShapeType.COMPOUND``.

    Attributes
    ----------
    geom_obj : Any | None
        The wrapped GEOM compound object (or ``None``).
    name : str
        The name of the wrapped GEOM compound object.
    """
    def __init__(self, geom_obj: Any | None) -> None:
        super().__init__(geom_obj, [ShapeType.COMPOUND])

    def __add__(self, other: Self | Sequence[Self]) -> Self:
        """
        Return a new ``Compound`` object resulting from fusing the current
        GEOM compound with the given ones.

        Parameters
        ----------
        other : Self | Sequence[Self]
            One or more ``Compound`` objects to fuse the current object with.

        Returns
        -------
        Self
            A new ``Compound`` instance resulting from fusing the current
            GEOM compound with the given ones.
        """
        return wrap_shape(
            make_face([
                make_fuse(
                    [self.geom_obj] + [
                        o for o in (
                            other if isinstance(other, Sequence) else [other])
                    ]
                )
            ])
        )

    def __repr__(self) -> str:
        """
        Return a descriptive string of the ``Compound`` instance containing
        a descriptive string for all the ``Face`` and ``Edge`` objects
        included in the compound.

        Returns
        -------
        str
            Descriptive string of the current instance.
        """
        faces = [
            wrap_shape(f) for f in extract_sub_shapes(self, ShapeType.FACE)]
        edges = [
            wrap_shape(e) for e in extract_sub_shapes(self, ShapeType.EDGE)]
        # Collect the information about its faces and edges
        faces_str = "\n   ".join(repr(f) for f in faces)
        edge_str = "\n   ".join(repr(e) for e in edges)
        return super().__repr__() + "\n   " + faces_str + "\n   " + edge_str


class Edge(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM edge object.

    Parameters
    ----------
    geom_obj : Any | None
        The GEOM edge object to wrap (or ``None``). If provided, its
        shape type is checked against the expected ``ShapeType.EDGE``.

    Attributes
    ----------
    geom_obj : Any | None
        The wrapped GEOM edge object (or ``None``).
    name : str
        The name of the wrapped GEOM edge object.
    """
    def __init__(self, geom_obj: Any) -> None:
        super().__init__(geom_obj, [ShapeType.EDGE])

    def __repr__(self) -> str:
        """
        Return a descriptive string of the ``Edge`` instance containing the
        XYZ coordinates of the GEOM vertex objects the edge is made of.

        Returns
        -------
        str
            Descriptive string of the current instance.
        """
        vertices = [
            wrap_shape(v) for v in extract_sub_shapes(self, ShapeType.VERTEX)]
        # Collect the information about its vertices
        vertices_str = "\n   ".join(repr(v) for v in vertices)
        return super().__repr__() + "\n   " + vertices_str


class Face(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM face object.

    Parameters
    ----------
    geom_obj : Any | None
        The GEOM face object to wrap (or ``None``). If provided, its
        shape type is checked against the expected ``ShapeType.FACE``.

    Attributes
    ----------
    geom_obj : Any | None
        The wrapped GEOM face object (or ``None``).
    name : str
        The name of the wrapped GEOM face object.
    """
    def __init__(self, geom_obj: Any) -> None:
        super().__init__(geom_obj, [ShapeType.FACE])

    def __add__(self, other: Self | Sequence[Self]) -> Self:
        """
        Return a new ``Face`` object resulting from fusing the current
        GEOM face with the given ones.

        Parameters
        ----------
        other : Self | Sequence[Self]
            One or more ``Face`` objects to fuse the current object with.

        Returns
        -------
        Self
            A new ``Face`` instance resulting from fusing the current GEOM
            face with the given ones.
        """
        return wrap_shape(
            make_face([
                make_fuse(
                    [self.geom_obj] + [
                        o for o in (
                            other if isinstance(other, Sequence) else [other])
                    ]
                )
            ])
        )

    def __repr__(self) -> str:
        """
        Return a descriptive string of the ``Face`` instance containing a
        descriptive string for each of the GEOM edge objects the face is
        made of.

        Returns
        -------
        str
            Descriptive string of the current instance.
        """
        edges = [
            wrap_shape(e) for e in extract_sub_shapes(self, ShapeType.EDGE)]
        # Collect the information about its edges
        edge_str = "\n   ".join(repr(e) for e in edges)
        return super().__repr__() + "\n   " + edge_str


class Vertex(GeomWrapper):
    """
    Class acting as a wrapper for the GEOM vertex object.

    Parameters
    ----------
    geom_obj : Any | None
        The GEOM vertex object to wrap (or ``None``). If provided, its
        shape type is checked against the expected ``ShapeType.VERTEX``.

    Attributes
    ----------
    geom_obj : Any | None
        The wrapped GEOM vertex object (or ``None``).
    name : str
        The name of the wrapped GEOM vertex object.
    """
    def __init__(self, geom_obj: Any) -> None:
        super().__init__(geom_obj, [ShapeType.VERTEX])

    def __eq__(self, other: Self) -> bool:
        """
        Return whether the current ``Vertex`` object is equal to the given
        one.

        Parameters
        ----------
        other : Self
            The ``Vertex`` object to check the current one for equality.

        Returns
        -------
        bool
            ``True`` if the current and given ``Vertex`` objects are the
            same (in terms of their coordinates); ``False`` otherwise.
        """
        if get_shape_type(other) != ShapeType.VERTEX:
            raise RuntimeError("The given shape must be a vertex object.")
        return all(
            math.isclose(i, j, abs_tol=1e-6) for i, j in zip(
                get_point_coordinates(self), get_point_coordinates(other)))

    def __repr__(self) -> str:
        """
        Return a descriptive string of the ``Vertex`` instance containing
        the XYZ coordinates of the GEOM vertex object.

        Returns
        -------
        str
            Descriptive string of the current instance.
        """
        return (
            super().__repr__()
            + f" - (x, y, z) = {get_point_coordinates(self)}"
        )


def wrap_shape(geom_obj: Any) -> GeomWrapper:
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
        shape_type = get_shape_type(geom_obj)
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
    return wrapper_cls(geom_obj)
