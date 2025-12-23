"""
Module containing the classes enabling the creation of the visualisation of
geometry layouts built in GLOW.
"""
import math

from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Any, Dict, List, Self, Sequence, Set, Tuple

from glow.interface.geom_entities import Compound, Edge, Face, Vertex, \
    wrap_shape
from glow.interface.geom_interface import ShapeType, add_to_study, add_to_study_in_father, clear_view, \
    display_shape, extract_sub_shapes, get_basic_properties, get_bounding_box, \
    get_closed_free_boundary, get_min_distance, get_object_from_id, \
    get_point_coordinates, get_shape_type, is_point_inside_shape, make_cdg, make_common, \
    make_compound, make_cut, make_face, make_partition, make_rotation, \
    make_scale, make_translation, make_vector_from_points, make_vertex, \
    make_vertex_inside_face, remove_from_study, set_color_face, update_salome_study
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.support.utility import generate_unique_random_colors


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
        The rotation angle (in degrees) of the GEOM object wrt the X-axis.
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
    def _rotate_from_axis(self, angle: float, axis: Edge) -> None:
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
    def show(self, *args: Any) -> None:
        """
        Abstract method for displaying the layout in the 3D viewer of SALOME
        according to the given settings.

        Parameters
        ----------
        *args : Any
            Positional arguments providing the display settings.
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


class Region(Face, Layout):
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
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the GEOM object.
    entry_id : str | None
        The ID of the region used in the current SALOME study.
    geom_obj : Any | None
        The internal `GEOM_Object` representative of the region's surface.
    name : str | None = None
        The name associated to the region used in the current SALOME study.
    o : Vertex
        The ``Vertex`` object representing the centre of the GEOM object.
    properties : Dict[PropertyType, str] | None = None
        The dictionary of property types and value names associated to the
        region.
    region_id : int
        The ID used to globally identify a ``Region`` object.
    rot_angle : float
        The rotation angle (in degrees) of the GEOM object wrt the X-axis.
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
        self.properties: Dict[PropertyType, str] | None = properties
        self.region_id: int = id(self._geom_obj)
        self.name = name if name else f"Region {self.region_id}"
        # Initialize superclass attributes
        self.entry_id = None
        self.o = wrap_shape(make_cdg(self.geom_obj))

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

    def rotate(self, angle: float) -> None:
        """
        Method for rotating the layout by the given angle (in degrees)
        around the axis perpendicular to the layout and passing through
        its centre.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        """
        # Get the figure center coordinates
        center = get_point_coordinates(self.o)
        # Build the Z-axis of rotation positioned in the figure center
        z_axis = wrap_shape(
            make_vector_from_points(
                self.o, make_vertex((center[0], center[1], 1))
            )
        )
        # Rotate the surface elements
        self._rotate_from_axis(angle, z_axis)

    def _rotate_from_axis(self, angle: float, axis: Edge) -> None:
        """
        Method for rotating the region by the given angle (in degrees)
        around the given axis.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        axis : Edge
            An ``Edge`` object representing the rotation axis.
        """
        # Convert the rotation angle in radians
        self.rot_angle = math.radians(angle)
        # Rotate the geometric elements of the surface
        self.geom_obj = make_rotation(self, axis, self.rot_angle)

    def scale(self, factor: float) -> None:
        """
        Method for scaling the region by the given factor.

        Parameters
        ----------
        factor : float
            The scaling factor.
        """
        self.geom_obj = make_scale(self.geom_obj, self.o, factor)

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

    def show(self, *args: Any) -> None:
        """
        Method for displaying the region in the 3D viewer of SALOME
        according to the given settings.

        Parameters
        ----------
        *args : Any
            Positional arguments providing the display settings. Must not
            be provided.
        """
        # Check that no args are provided
        if args:
            raise ValueError(f"No arguments are accepted for 'show()'.")
        # Erase all objects from the current view
        clear_view()
        # Delete the region from the study if already present
        if self.entry_id and get_object_from_id(self.entry_id):
            remove_from_study(self.entry_id)
        # Add the region's GEOM face to the current SALOME study
        self.entry_id = add_to_study(self, self.name)
        display_shape(self.entry_id)
        # Update the SALOME view
        update_salome_study()

    def translate(self, new_cntr: Tuple[float, float, float]) -> None:
        """
        Method for translating the region so that its centre is positioned
        at the given coordinates.

        Parameters
        ----------
        new_cntr : Tuple[float, float, float]
            The XYZ coordinates the region centre should be placed at.
        """
        # Build a vector from the current center to the new one
        transl_vect = make_vector_from_points(self.o, make_vertex(new_cntr))
        # Translate the wrapped GEOM face object
        self.geom_obj = make_translation(self.geom_obj, transl_vect)
        self.o = wrap_shape(make_translation(self.o, transl_vect))

    def update(self, layout: Face) -> None:
        """
        Method for updating the GEOM object the region refers to.

        Parameters
        ----------
        layout : Face
            The new layout to update the current region with.
        """
        self.geom_obj = layout.geom_obj
        self.o = wrap_shape(make_cdg(layout))

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

    def __mul__(self, other: Self) -> Self:
        """
        Return a new ``Region`` instance resulting as the common part between
        the current region and the given one.
        The result of this operation shares the same properties of the given
        region.

        Parameters
        ----------
        other : Self
            A ``Region`` object to extract the common part with the current
            one.

        Returns
        -------
        Self
            A new ``Region`` instance resulting as the common part between
            the current region and the given one.
        """
        return Region(
            wrap_shape(self.geom_obj) * wrap_shape(other).geom_obj,
            properties=deepcopy(other.properties)
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


class Fillable(Compound, Layout):
    """
    Class to represent any geometry layout that can be filled with regions
    each associated to properties, such as the material.

    This class inherits from the ``Compound`` wrapper class; this guarantees
    that subclasses of ``Fillable`` behave like the corresponding GEOM
    compound object when used with GEOM functions.
    In addition, this class inherits from the ``Layout`` class, meaning that
    a proper implementation of each abstract method of ``Layout`` is included
    herein or in subclasses of ``Fillable``.

    This class offers to its subclasses the capability to represent their
    geometry layout as the superimposition of multiple layers made either by
    ``Region`` or ``Fillable`` objects. By relying on this layer concept,
    the whole layout can be represented in a tree-like structure where nodes
    are represented by ``Fillable`` objects, while leaves by ``Region`` ones.
    In addition, it enables a mapping between different geometry views,
    according to the ``GeometryType`` enumeration.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the GEOM object.
    displayed_geom : GeometryType
        The currently ``GeometryType`` displayed in the SALOME 3D viewer.
    entry_id : str | None
        The ID attributed by SALOME when the GEOM object is added to the
        study.
    geometry_maps : Dict[GeometryType, Compound]
        A mapping from ``GeometryType`` values to ``Compound`` objects. Each
        entry provides a different representation for the geometry layout this
        instance refers. It is used to switch between different visualisation
        types (e.g., technological, sectorized).
    is_update_needed : bool
        Flag that indicates whether an update is required (e.g., geometry
        layout rebuilding).
    layers : List[List[Region | Self]]
        A list of layers, each layer itself being a list of ``Region`` objects
        or nested ``Fillable`` instances. Layers represent the hierarchical
        structure of the geometry layout.
    o : Vertex
        The ``Vertex`` object representing the centre of the GEOM object.
    regions : List[Region]
        A flat list of Region objects contained by the current instance.
        Maintained in addition to the layered `layers` structure to allow
        quick iteration or lookups over all regions.
    rot_angle : float
        The rotation angle (in degrees) of the GEOM object wrt the X-axis.
    """
    def __init__(self) -> None:
        super().__init__(None)
        # Initialize attributes
        self.layers: List[List[Region | Self]] = []
        self.geometry_maps: Dict[GeometryType, Compound] = {}
        self.is_update_needed: bool = False
        self.displayed_geom: GeometryType = GeometryType.TECHNOLOGICAL
        self.regions: List[Region] = []

    def add(self,
            layout: Region | Self,
            position: Tuple[float, float, float] | None = None,
            layer_index: int | None = None
        ) -> None:
        """
        Method that adds a generic layout, i.e. either a ``Region`` or a
        `Fillable` object, to the technological geometry layout of this
        instance at the indicated layer.

        The given layout object is first translated, if needed, so that the
        coordinates of its CDG match the indicated position. If no position
        is provided, the layout object is placed in the CDG of this instance.
        The given layout is stored at the end of the sublist of the ``layers``
        attribute that is specified by the value of the indicated parameter
        ``layer_index``, if any. Otherwise, a new sublist (i.e. a new layer)
        is created with the given layout.
        If the layer index value is not valid, an exception is raised.
        Please note that, when providing layouts on the same layer, they must
        not overlap.

        Parameters
        ----------
        layout : Region | Self
            The generic layout object to add.
        position : Tuple[float, float, float] | None
            The XYZ coordinates of the layout's centre, if any. It defaults
            to ``None``, meaning that the layout is added at the current
            instance centre.
        layer_index : int | None
            The index identifying the layer at which the given layout should
            be added. It defaults to ``None``, meaning that the layout is
            stored in a new layer.

        Raises
        ------
        ValueError
            If the indicated layer index is not valid.

        Notes
        -----
        This method simply updates the list of layers of the technological
        geometry layout with the given layout without collapsing the layers
        and updating the entire GEOM compound object this instance refers to.
        To collapse the layers and build the regions of this instance without
        displaying the geometry in the 3D viewer of SALOME, call the method
        ``build_regions``. To build and display the geometry layout, call the
        method ``show``.
        """
        # Set the given layout name to include the one of the current compound
        layout.name = f"{self.name}_{layout.name}"
        # Set the given layout position to the current compound centre, if no
        # position is provided.
        if not position:
            position = get_point_coordinates(self.o)
        # Translate the given layout if its position differs from the compound
        # centre
        if not all(math.isclose(i, 0.0) for i in position):
            layout = layout.clone()
            layout.translate(position)

        # Include the given layout at the end of the sublist specified by the
        # indicated index, if any; otherwise, either create a new sublist or
        # raise an exception
        if layer_index is None or layer_index == len(self.layers):
            # Create a new layer
            self.layers.append([layout])
        elif layer_index >= 0 or layer_index < len(self.layers):
            # Add to existing layer
            self.layers[layer_index].append(layout)
        else:
            raise ValueError(f"Invalid layer index {layer_index}.")

        # Indicate the need to update the layout by building its regions
        self.is_update_needed = True

    def build_regions(self) -> None:
        """
        Method that processes the layers of the current instance object to
        construct a list of ``Region`` objects that are representative of
        the technological geometry layout.
        These regions are the result of collapsing all the layers by
        traversing the entire hierarchical structure in reverse order, and
        cutting the regions of the layers that are overlapped. This operation
        is skipped and the method exits without changes if there is no need
        to update the technological geometry layout (i.e. the value of the
        attribute ``is_update_needed`` is ``False``).

        The operation of assembling all the layers requires that each layout
        is up-to-date. This means that each ``Fillable`` object in the list
        of layers is further processed by recursively calling this method to
        collapse its own layers.

        Given the updated layout objects for each layer, the corresponding
        ``Region`` objects are retrieved and a flat list is populated.
        Finally, the GEOM compound object representative of this instance is
        updated by performing a partition with all the built ``Region``
        objects.
        """
        # Return immediately if there is no need to update the layout
        if not self.is_update_needed:
            return
        # Reverse the layers
        reversed_layers = self.layers[::-1]
        # Ensure the layout is up-to-date for each 'Fillable' object by
        # recursively calling this method
        for layer in reversed_layers:
            for layout in layer:
                if isinstance(layout, Fillable):
                    layout.build_regions()
        # Collapse all the layers and cut out the overlapping regions
        self._collapse_layers()
        # Clear and rebuild the list of regions from those of each layer
        self.regions.clear()
        for layer in reversed_layers:
            for layout in layer:
                if isinstance(layout, Region):
                    self.regions.append(layout)
                elif isinstance(layout, Fillable):
                    self.regions.extend(layout.regions)

        # Update the GEOM compound object of this instance
        self.geom_obj = make_partition(self.regions, [], ShapeType.COMPOUND)
        # Update the flag stating there is no need to rebuild the regions
        self.is_update_needed = False

    def clone(self) -> Self:
        """
        Method that returns a copy of the current ``Fillable`` instance.

        Returns
        -------
        Self
            A copy of the current instance.
        """
        return deepcopy(self)

    def show(self, *args: Any) -> None:
        """
        Method for displaying the geometry layout of the cell in the 3D viewer
        of SALOME according to the ``PropertyType`` and ``GeometryType``
        settings provided as input regardless of the order.
        Regions are displayed using a color map according to the values of
        ``PropertyType`` associated to them. If no ``PropertyType`` is
        specified, the regions are displayed without any color.
        By default this method builds (if needed) and displays the regions of
        the technological geometry. If a different ``GeometryType`` is given,
        the method also displays the GEOM compound of the edges decribing the
        cell's refined geometry.

        Parameters
        ----------
        *args : Any
            Positional arguments providing the display settings.

        Notes
        -----
        Allowed arguments to this method are:

        - ``PropertyType``, used to derive the color map associated to the
          regions of the technological geometry;
        - ``GeometryType``, used to indicate the type of geometry of the
          layouts of the cell.
        """
        # Initialize the display settings with default values
        prop_type = None
        geom_type = GeometryType.TECHNOLOGICAL
        # Extract the display settings from the parameters of the method
        settings = list(args)
        for setting in settings:
            if isinstance(setting, PropertyType):
                prop_type = setting
            elif isinstance(setting, GeometryType):
                geom_type = setting
            else:
                raise RuntimeError(f"Unknown '{setting}'.")

        # Erase all objects from the current view
        clear_view()
        # If the GEOM compound is already present in the current SALOME study,
        # remove it
        if self.entry_id and get_object_from_id(self.entry_id):
            remove_from_study(self.entry_id)
        # Re-build the regions, if needed
        self.build_regions()
        # Add the updated GEOM compound of the cell to the study
        add_to_study(self.geom_obj, self.name)

        # Assign the same color to all the regions having the same property
        # value, if any has been specified to show
        try:
            associate_colors_to_regions(prop_type, self.regions)
        except Exception as e:
            # Add all the regions of the layout to the study, with the faulty
            # ones (i.e. those without the property) colored in red
            self._show_regions()
            # Re-raise the exception
            raise RuntimeError from e
        # Add all the regions of the layout to the study, each with its color,
        # if any
        self._show_regions()
        # Update the displayed geometry type
        self.displayed_geom = geom_type
        # Set update flag to False
        self.is_update_needed = False

    def update(self, layout: Compound | Face) -> None:
        """
        Method for updating the GEOM object this instance refers to.

        Parameters
        ----------
        layout : Compound | Face
            The new layout to update the current instance with.
        """
        self.geom_obj = layout.geom_obj
        self.o = wrap_shape(make_cdg(layout))

    def _apply_cut_to_layouts_in_layer(
            self, layer: List[Region | Self], cutting_tool: Face) -> None:
        """
        Method that applies a cutting operation to all layout objects within
        a given layer using the given ``Face`` object as cutting tool.

        A loop through all the layout objects of the given layer is performed
        to handle the overlap operations:
        - if the distance between the cutting tool and the layout object is
          greater than 0 (with a tolerance), no overlapping is expected and a
          new layout is considered;
        - the cut operation is performed and, if the resulting shape does not
          contain any GEOM face object, it means the layout object is
          completely overlapped and its index in the layer list is stored;
        - if the layout object is a ``Region`` instance, its GEOM face is
          updated with the result of the cut, if it is a GEOM face, otherwise
          if the cut results in multiple faces, each is inserted as a new
          ``Region`` and the original instance is removed;
        - if the layout object is a ``Fillable`` instance, its GEOM compound
          is updated with the result of the cut and the same operation is
          recursively applied to its ``Region`` objects;
        - lastly, the layout objects completely overlapped are removed by the
          list of objects of the given sub-layer.

        Parameters
        ----------
        layer : List[Region | Self]
            The list of layout objects in the layer to be cut.
        cutting_tool : Face
            The geometric shape used to cut the layouts in the layer.
        """
        # List storing the indices of the layouts to remove from the
        # sub-layer, as completely overlapped
        layouts_to_remove: List[int] = []
        # Loop through all the layout objects of the sub-layer and perform
        # the cuts
        for i, layout in enumerate(layer):
            # Continue with the next layout in the sub-layer if the current
            # layout is not close to the layer
            if get_min_distance(layout, cutting_tool) > 1e-6:
                continue
            # Cut the layout with the layer shape and check the result
            cut_layout = make_cut(layout, cutting_tool)
            if not extract_sub_shapes(
                make_compound([cut_layout]), ShapeType.FACE
            ):
                layouts_to_remove.append(i)
                continue
            # Update the layout object according to its type and the result
            # of the cut operation
            is_cut_a_cmpd = get_shape_type(cut_layout) == ShapeType.COMPOUND
            if isinstance(layout, Region):
                # Substitute the Region with those resulting from the cut
                if is_cut_a_cmpd:
                    layer.pop(i)
                    for j, face in enumerate(
                        extract_sub_shapes(cut_layout, ShapeType.FACE)
                    ):
                        # Create a new region for each subface
                        layer.insert(
                            i+j, Region(face, layout.name, layout.properties))
                    continue
                layer[i].update(wrap_shape(cut_layout))
            elif isinstance(layout, Fillable):
                layer[i].update(wrap_shape(cut_layout))
                # Recursively apply the cut
                self._apply_cut_to_layouts_in_layer(
                    layer[i].regions, cutting_tool)

        # Remove the layouts completely overlapped by the superior layer
        for index in sorted(layouts_to_remove, reverse=True):
            layer.pop(index)

    def _collapse_layers(self) -> None:
        """
        Method that collapses all the layers (a list containing lists of
        ``Region`` or ``Fillable`` objects) of the current instance by
        traversing them in reverse order.
        If any compound object of a layer overlaps any layer below it, its
        ``Region`` objects are cut by the higher layer compound.

        Notes
        -----
        Layout objects belonging to the same layer should not overlap each
        other.
        """
        # Reverse the list of layers
        reversed_layers = self.layers[::-1]
        for i, layer in enumerate(reversed_layers):
            # Skip the layer, if empty
            if not layer:
                continue
            # Loop through all the layers below the current one to cut all the
            # layout objects of each sublayer, if any is overlapped.
            for sub_layer in reversed_layers[i + 1:]:
                # Skip the layer, if empty
                if not sub_layer:
                    continue
                # Overlap the current layer onto the layers below
                self._overlap_layer_to(reversed_layers[i], sub_layer)

    def _overlap_layer_to(
            self, layer: List[Region | Self], sub_layer: List[Region | Self]
        ) -> None:
        """
        Method that overlaps a layer, whose shape is derived from its layout
        objects, onto a layer below it, which is given as a list of ``Region``
        and/or ``Fillable`` objects.
        The shape of the layer, which acts as the cutting tool for the layout
        objects of the given sub-layer is determined on the basis of the
        number of closed boundaries that can be extracted by the compound of
        the layer. When more than one are found, a partition operation is
        performed to reduce the number of closed boundaries to the minimum.
        In any case, a GEOM face is built from the found boundaries. If
        multiple boundaries are found, the resulting shape may have holes.
        The cut operation is then performed on the sublayer by using the built
        shape as cutting tool.

        Parameters
        ----------
        layer : List[Region | Self]
            The superior layer representing the cutting tool.
        sub_layer : List[Region | Self]
            The inferior layer whose layout objects are cut by the superior
            layer.
        """
        # Build the shape of the layer compound
        layer_cmpd = make_compound(layer)
        boundaries = get_closed_free_boundary(layer_cmpd)
        if len(boundaries) > 1:
            boundaries = get_closed_free_boundary(
                make_partition(layer, [], ShapeType.FACE)
            )
        layer_shape = make_face(boundaries)
        # Apply the cut on the layout objects of the sub-layer
        self._apply_cut_to_layouts_in_layer(sub_layer, layer_shape)

    def _show_regions(self) -> None:
        """
        Method that adds all the regions of the layout to the current SALOME
        study. In the Object Browser they are available as children of the
        GEOM compound representative of the technological geometry of the
        layout.
        Each region is displayed with a colour defined beforehands.
        """
        for region in self.regions:
            # Set the region color in the viewer
            set_color_face(region.geom_obj, region.color)
            # Add the cell region to the study
            region.entry_id = add_to_study_in_father(
                self.geom_obj, region, region.name)
            # Display the region in the current view, if needed
            display_shape(region.entry_id)

        # Show everything on the SALOME application
        update_salome_study()


def is_layout_contained(
        container: Compound | Face,
        candidate: Compound | Face,
        tolerance: float = 1e-6) -> bool:
    """
    Check if `candidate` layout is entirely contained within `container`
    layout. The containment check consists of:

    - Area comparison: ``False`` is returned if the candidate area is larger
      than container area.
    - Bounding box check: ``False`` is returned if the candidate extends
      beyond container bounds.
    - Precise containment: The common area is computed and compared with the
      candidate area. ``True`` is returned if the two areas are the same.

    Parameters
    ----------
    container : Compound | Face
        The geometry layout acting as the container.
    candidate : Compound | Face
        The geometry layout to check for containment.
    tolerance : float
        Numerical tolerance for area and bounding box comparisons. Default
        valus is 1e-6.

    Returns
    -------
    bool
        ``True`` if `candidate` layout is fully contained within the layout of
        the `container` by the given tolerance; ``False`` otherwise.
    """
    # Area check
    area_container = get_basic_properties(container)[1]
    area_candidate = get_basic_properties(candidate)[1]
    if area_candidate > area_container + tolerance * area_container:
        return False

    # Bounding box check
    bbox_container = get_bounding_box(container)
    bbox_candidate = get_bounding_box(candidate)
    if any([
        bbox_candidate[i] - bbox_container[i] <= -tolerance for i in [0, 2]
    ]) or any([
        bbox_candidate[i] - bbox_container[i] >= tolerance for i in [1, 3]
    ]):
        return False

    # Precise geometric intersection
    common = make_common(container, candidate)
    area_common = get_basic_properties(common)[1]
    return abs(area_common - area_candidate) < tolerance * area_candidate


def associate_colors_to_regions(
        property_type: PropertyType | None, regions: List[Region]) -> None:
    """
    Method that assigns the same color to all the regions having the same
    value for the given property type.

    Parameters
    ----------
    property_type : PropertyType | None
        The type of property for which colors must be assigned to regions
        having the same property type value. If ``None``, the regions
        color is reset to its default value.
    regions : List[Region]
        The list of ``Region`` objects to colour according to the value of
        the associated property type.
    """
    # If no colorset to display, reset the region colors
    if not property_type:
        for region in regions:
            # Set the region color to its default value
            region.reset_region_color()
        return
    # Extract the unique values of the given property type associated to
    # each region
    values = get_unique_values_for_property(property_type, regions)
    # Generate a specific amount of colors as the number of different
    # values for the same given property type
    colors = generate_unique_random_colors(len(values))
    # Build a dictionary of values for the given property type VS color
    property_vs_color = dict(zip(values, colors))
    # Loop through all the regions and assign a color corresponding to
    # the value of the given property type
    for region in regions:
        # Get the value of the given property type associated to the
        # region
        value = region.properties[property_type]
        # Set the region color
        region.set_region_color(property_vs_color[value])


def get_unique_values_for_property(
        property_type: PropertyType, regions: List[Region]) -> List[str]:
    """
    Method that gets the unique values of the given property type for the
    given regions. If any ``Region`` object does not have any property
    or the given property type is missing, a reference point for the
    region is stored for logging purposes.
    An exception showing the coordinates of the points of the problematic
    regions is raised.

    Parameters
    ----------
    property_type : PropertyType
        The type of property whose unique values to collect.
    regions : List[Region]
        The list of ``Region`` objects to select the unique values associated
        to the indicated ``PropertyType``.

    Returns
    -------
    List[str]
        A list of the unique names for the given property type that have
        been associated to the given regions.

    Raises
    ------
    RuntimeError
        Showing the coordinates of the points of the regions having any
        issue.
    """
    values = set()
    missing_regions_points = []
    for region in regions:
        if not region.properties or property_type not in region.properties:
            missing_regions_points.append(
                get_point_coordinates(
                    make_vertex_inside_face(region)))
            region.set_region_color((255, 0, 0))
            continue
        values.add(region.properties[property_type])
    # Raise an exception if there are regions with missing property
    if missing_regions_points:
        message = (
            f"No {property_type.name} property type has been found "
            "for the regions identified by the inner points with "
            f"coordinates: ")
        for point in missing_regions_points[:-1]:
            message += f"'{point}', "
        message += (f"'{missing_regions_points[-1]}'. Please, call the "
            + "'show()' method to show the regions, select the one with "
            + "a missing property and call the 'set_region_property()' "
            + "method to assign the property to.")
        raise RuntimeError(message)
    return list(values)
