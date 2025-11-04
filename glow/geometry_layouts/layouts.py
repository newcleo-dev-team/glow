"""
Module containing the classes enabling the creation of the visualisation of
geometry layouts built in GLOW.
"""
from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Any, Dict, Self, Sequence, Tuple

from glow.interface.geom_entities import Compound, Edge, Face, Vertex, \
    wrap_shape
from glow.interface.geom_interface import make_cdg, make_vertex
from glow.support.types import PropertyType


DEFAULT_REGION_COLOR: Tuple[int, int, int] = (167, 167, 167)
"""
Default RGB color for a ``Region`` object when displayed in the SALOME viewer.
The code corresponds to light grey.
"""


class Layout(ABC):
    """
    Abstract class that provides characteristic data and behaviour a generic
    layout object should have.
    Attributes are common to all subclasses of ``Layout``; they provide
    geometric data and the ID used by SALOME when displaying the GEOM object
    the layout refers to. The GEOM object of the layout is provided by
    subclasses of ``Layout``.
    Methods herein declared are abstract, meaning only their signature is
    provided. Subclasses must provide dedicated implementations for each of
    them. These methods enable the application of Euclidean transformations
    to generic layouts, as well as the capability to display and update the
    layout.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the GEOM object.
    entry_id : str | None
        The ID attributed by SALOME when the GEOM object is added to the
        study.
    o : Vertex
        The ``Vertex`` object representing the centre of the GEOM object.
    rot_angle : float
        The rotation angle of the GEOM object wrt the X-axis.
    """
    def __init__(self) -> None:
        super().__init__()
        # Initialize the attributes
        self.dimensions: Tuple[float, float] = (0.0, 0.0)
        self.entry_id: str | None = None
        self.o: Vertex = Vertex(make_vertex((0.0, 0.0, 0.0)))
        self.rot_angle: float = 0.0

    @abstractmethod
    def rotate(self, angle: float) -> None:
        """
        Abstract method for rotating the layout by the given angle (in
        degrees) around the axis perpendicular to the layout and passing
        through its centre.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        """

    @abstractmethod
    def rotate_from_axis(self, angle: float, axis: Edge) -> None:
        """
        Abstract method for rotating the layout by the given angle (in
        degrees) around the given axis.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        axis : Edge
            An ``Edge`` object representing the rotation axis.
        """

    @abstractmethod
    def scale(self, factor: float) -> None:
        """
        Abstract method for scaling the layout by the given factor.

        Parameters
        ----------
        factor : float
            The scaling factor.
        """

    @abstractmethod
    def show(self, *args: Any, **kwargs: Any) -> None:
        """
        Abstract method for displaying the layout in the 3D viewer of SALOME
        according to the given settings.

        Parameters
        ----------
        *args : Any
            Positional arguments.
        **kwargs : Any
            Key arguments.
        """

    @abstractmethod
    def translate(self, new_cntr: Tuple[float, float, float]) -> None:
        """
        Abstract method for translating the layout so that its centre
        is positioned at the given coordinates.

        Parameters
        ----------
        new_cntr : Tuple[float, float, float]
            The XYZ coordinates the layout centre should be placed at.
        """

    @abstractmethod
    def update(self, layout: Compound | Face) -> None:
        """
        Abstract method for updating the GEOM object the layout refers to.

        Parameters
        ----------
        layout : Compound | Face
            The new layout to update the current with.
        """


class Region(Face):
    """
    Class that groups information regarding a generic region of the cell's
    or the lattice's geometry layout.
    A ``Region`` instance represents any 2D surface bounded by one or two
    edges that is filled by properties, e.g. the material.
    Properties are expressed as a dictionary of items of the ``PropertyType``
    enumeration vs the corresponding names.

    This class inherits from the ``Face`` wrapper class; this guarantees a
    ``Region`` instance behaves like the corresponding GEOM face object when
    used with GEOM functions.

    Parameters
    ----------
    geom_obj : Any
        The GEOM face object the region represent.
    name : str
        The name of the region used in the current SALOME study.
    properties: Dict[PropertyType, str] | None = None
        A dictionary of items of the ``PropertyType`` enumeration and the
        corresponding values. They indicate the properties associated to
        the region (e.g. the material).

    Attributes
    ----------
    color : Tuple[int, int, int]
        The tuple whose values represent the RGB color code associated to the
        region when displayed in the SALOME viewer according to a property
        map.
    entry_id : str | None
        The ID of the region used in the current SALOME study.
    geom_obj : Any | None
        The internal `GEOM_Object` representative of the region's surface.
    name : str | None = None
        The name associated to the region used in the current SALOME study.
    properties : Dict[PropertyType, str] | None = None
        The dictionary of property types and value names associated to the
        region.
    region_id : int
        The ID used to globally identify a ``Region`` object.
    """
    def __init__(
            self,
            geom_obj: Any,
            name: str | None = None,
            properties: Dict[PropertyType, str] | None = None
        ) -> None:
        super().__init__(geom_obj)
        # Initialize attributes
        self.color: Tuple[int, int, int] = DEFAULT_REGION_COLOR
        self.entry_id: str | None = None
        self.properties: Dict[PropertyType, str] | None = properties
        self.region_id: int = id(self._geom_obj)
        # Set the underlying GEOM face object name
        self.name = name if name else f"Region {self.region_id}"

    def clone(self) -> Self:
        """
        Method that returns a copy of the current ``Region`` instance with
        the same ID, properties and color.

        Returns
        -------
        Self
            A copy of the current instance sharing the same ID, properties
            and color.
        """
        # Store ID and color
        prev_r_id = self.region_id
        prev_col = self.color
        # Create a new Region with the same GEOM_Object, name and properties
        cloned_region = Region(
            self.geom_obj, self.name, deepcopy(self.properties))
        # Restore the ID and color
        cloned_region.region_id = prev_r_id
        cloned_region.color = prev_col
        return cloned_region

    def reset_region_color(self) -> None:
        """
        Method for resetting the region color to the (167, 167, 167) RGB
        value, which corresponds to light grey.
        """
        self.color = DEFAULT_REGION_COLOR

    def set_region_color(self, color: Tuple[int, int, int]) -> None:
        """
        Method for associating a specific RGB color to the region.

        Parameters
        ----------
        color : Tuple[int, int, int]
            RGB values identifying the color to associate the region with.
        """
        # Check the given color
        if (not color
            or len(color) != 3
            or any(value < 0 or value > 255 for value in color)
        ):
            raise ValueError("The 'color' parameter must be provided as a "
                             "tuple of 3 integers in the 0:255 range.")
        # Store the RGB color
        self.color = color

    def __add__(self, other: Self | Sequence[Self]) -> Self:
        """
        Return a new ``Region`` instance resulting from the fusion of the
        GEOM face associated with the current region and the given ones.
        The fuse operation cannot be performed if the involved regions have
        different values for the same properties.

        Parameters
        ----------
        other : Self | Sequence[Self]
            One or more ``Region`` objects to fuse the current object with.

        Returns
        -------
        Self
            A new ``Region`` instance resulting from fusing the current region
            with the given ones.

        Raises
        ------
        RuntimeError
            When all the involved regions provide the same property type, but
            at least one of them uses a different value.
        """
        # Extract a list of the given regions
        regions = [
            o for o in (other if isinstance(other, Sequence) else [other])]
        # Properties belonging to the fused region
        fused_properties = deepcopy(self.properties)
        # If the regions have different values for the same properties,
        # the fuse operation cannot be performed
        for region in regions:
            for prop, val in region.properties.items():
                if prop not in fused_properties:
                    fused_properties[prop] = val
                    continue
                if val != fused_properties[prop]:
                    raise RuntimeError(
                        f"No fuse between the {self.name} and {region.name} "
                        "is possible as they do not share the same value for "
                        f"the {prop}.")
        # Fuse all the GEOM_Objects and return a new Region out of it
        return Region(
            geom_obj=(
                wrap_shape(self.geom_obj)
                + [wrap_shape(region) for region in regions]
            ).geom_obj,
            properties=fused_properties
        )

    def __repr__(self) -> str:
        """
        Return a descriptive string about the current ``Region`` instance
        indicating its characteristics, i.e.:
        - the region name;
        - the associated properties;
        - the color associated to the GEOM face the region corresponds to;
        - the XYZ coordinates of the region's CDG.

        Returns
        -------
        str
            A descriptive string containing information about the current
            ``Region`` instance.
        """
        return f"{self.name}, {self.properties}, {self.color}, " + \
            f"{wrap_shape(make_cdg(self.geom_obj))}"

    def __sub__(self, other: Self) -> Self:
        """
        Return a new ``Region`` instance resulting from cutting the current
        region with the given one.

        Parameters
        ----------
        other : Self
            A ``Region`` object to cut the current one with.

        Returns
        -------
        Self
            A new ``Region`` instance resulting from cutting the current
            region with the given one.
        """
        return Region(
            wrap_shape(self.geom_obj) - wrap_shape(other).geom_obj,
            properties=deepcopy(self.properties)
        )
