"""
Module containing the class enabling the creation and the visualisation of
geometry layouts that can be described according to a hierarchical structure.
"""
from copy import deepcopy
import math
from typing import Any, Dict, List, Self, Tuple
from glow.geometry_layouts.layouts import Layout, Region, \
    associate_colors_to_regions
from glow.interface.geom_entities import Compound, Face, wrap_shape
from glow.interface.geom_interface import ShapeType, add_to_study, \
    add_to_study_in_father, clear_view, display_shape, extract_sub_shapes, \
    get_closed_free_boundary, get_min_distance, get_object_from_id, \
    get_point_coordinates, get_shape_type, make_cdg, make_compound, make_cut, \
    make_face, make_partition, remove_from_study, set_color_face, \
    update_salome_study
from glow.support.types import GeometryType, PropertyType


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
