"""
Module containing the class enabling the creation and the visualisation of
geometry layouts that can be described according to a hierarchical structure.
"""
import logging
import math

from abc import abstractmethod
from copy import deepcopy
from typing import Any, Dict, Iterator, List, Self, Tuple

from glow.geometry_layouts.geometries import Surface
from glow.geometry_layouts.layouts import Layout, LayoutState, Region, \
    associate_colors_to_regions, is_layout_contained
from glow.geometry_layouts.symmetry_management import SymmetryDomain, \
    build_common_symmetry_shape, use_symmetry_logic
from glow.interface.geom_entities import Compound, Edge, Face, Vertex, \
    wrap_shape
from glow.interface.geom_interface import ShapeType, add_to_study, \
    add_to_study_in_father, clear_view, display_shape, extract_sub_shapes, \
    get_closed_free_boundary, get_min_distance, get_object_from_id, \
    get_point_coordinates, get_shape_name, get_shape_type, \
    is_point_inside_shape, make_cdg, make_common, make_compound, make_cut, \
    make_face, make_partition, make_rotation, make_scale, make_translation, \
    make_vector_from_points, make_vertex, make_vertex_inside_face, \
    remove_from_study, set_color_face, update_salome_study
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.support.utility import are_same_shapes, build_z_axis_from_vertex, \
    compute_point_by_reference, flatten_list, retrieve_selected_object


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
    geometry layout as the superimposition of multiple layers made either of
    ``Region`` or ``Fillable`` objects. By relying on this layer concept,
    the whole layout can be represented in a tree-like structure where nodes
    are represented by ``Fillable`` objects, while leaves by ``Region`` ones.
    In addition, it enables a mapping between different geometry views,
    according to the ``GeometryType`` enumeration.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the GEOM object.
    entry_id : str | None
        The ID attributed by SALOME when the GEOM object is added to the
        study.
    geom_obj : Any | None
        The internal `GEOM_Object` representative of the layout this instance
        refers to.
    geometry_maps : Dict[GeometryType, Compound]
        A mapping from ``GeometryType`` values to ``Compound`` objects. Each
        entry provides a different representation for the geometry layout this
        instance refers. It is used to switch between different visualisation
        types (e.g., technological, sectorized).
    layers : List[List[Region | Self]]
        A list of layers, each layer itself being a list of ``Region`` objects
        or nested ``Fillable`` instances. Layers represent the hierarchical
        structure of the geometry layout.
    name : str | None = None
        The name of the GEOM object identifying the layout this instance
        refers to. It is used when the layout is added to the current SALOME
        study.
    o : Vertex
        The ``Vertex`` object representing the centre of the GEOM object.
    regions : List[Region]
        A flat list of Region objects contained by the current instance.
        Maintained in addition to the layered `layers` structure to allow
        quick iteration or lookups over all regions.
    rot_angle : float
        The rotation angle (in degrees) of the GEOM object wrt the X-axis.
    shape : Surface | None = None
        The ``Surface`` object representing the characteristic shape of the
        layout.
    state : LayoutState
        Providing the state of the layout in the SALOME study.
    symmetry_map : Dict[SymmetryType, Surface]
        A mapping from ``SymmetryType`` values to ``Surface`` objects. Each
        entry provides the characteristic shape of the corresponding symmetry
        type.
    """
    def __init__(self) -> None:
        super().__init__(None)
        # Initialize attributes
        self.layers: List[List[Region | Self]] = []
        self.geometry_maps: Dict[GeometryType, Compound] = {}
        self.state = LayoutState()
        self.regions: List[Region] = []
        self.symmetry_map: Dict[SymmetryType, Surface] = {}
        self.shape: Surface | None = None

    def add(self,
            layout: Region | Self,
            position: Tuple[float, float, float] | None = None,
            layer_index: int | None = None
        ) -> None:
        """
        Method that adds a generic layout, i.e. either a ``Region`` or a
        ``Fillable`` object, to the technological geometry layout of this
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
        To collapse the layers and update all the contained layouts of this
        instance without displaying the geometry in the 3D viewer of SALOME,
        call the method ``update_hierarchical_structure``.
        To update and display the geometry layout with all the regions of the
        technological geometry, call the method ``show``.
        """
        # Set the given layout position to the current compound centre, if no
        # position is provided.
        if not position:
            position = get_point_coordinates(self.o)
        # Update the hierarchical structure of the given layout, if needed,
        # and clone it
        if isinstance(layout, Fillable):
            layout.update_hierarchical_structure()
        layout = layout.clone()
        # Set the given layout name to include the one of the current compound
        layout.name = f"{self.name}_{layout.name}"
        # Translate the given layout if its position differs from the layout
        # centre
        if not all(math.isclose(i, 0.0, abs_tol=1e-6) for i in position):
            layout.translate(position)

        # Include the given layout at the end of the sublist specified by the
        # indicated index, if any; otherwise, either create a new sublist or
        # raise an exception
        if layer_index is None or layer_index == len(self.layers):
            # Create a new layer
            self.layers.append([layout])
        elif layer_index >= 0 and layer_index < len(self.layers):
            # Add to existing layer
            self.layers[layer_index].append(layout)
        else:
            raise ValueError(f"Invalid layer index {layer_index}.")

        # Indicate the need to update the layout by building its regions
        self.state.is_update_needed = True

    def apply_symmetry(self, symmetry: SymmetryType) -> None:
        """
        Method for deriving the portion of the geometry layout that represents
        the given type of symmetry.
        According to the type, the shape that identifies the symmetry is
        derived and stored in the corresponding instance attribute mapping
        types and shapes.
        The ``state`` attribute, which indicates the state of the layout, is
        modified by properly setting the active symmetry type.

        Parameters
        ----------
        symmetry : SymmetryType
            The type of symmetry to handle.

        Raises
        ------
        RuntimeError
            If the indicated symmetry type is not supported for the current
            instance.
        """
        # Update the layout characteristic shape from its regions, if needed
        if self.geom_obj is None:
            self.update_hierarchical_structure()
        # Build the shape of the symmetry
        try:
            symm_shape = self._build_symmetry_shape(
                symmetry,
                SymmetryDomain(get_point_coordinates(self.o), self.shape)
            )
        except RuntimeError as e:
            raise RuntimeError(
                "Error while building the symmetry shape for a "
                f"'{self.__class__.__name__}' instance."
            ) from e
        # Rotate the shape of the symmetry, if needed
        symm_shape.rotate(self.rot_angle, build_z_axis_from_vertex(self.o))
        # Store the shape of the symmetry in the mapping
        self.symmetry_map[symmetry] = symm_shape
        # Update the state of the layout
        self.state.symmetry_type = symmetry

    def clone(self) -> Self:
        """
        Method that returns a copy of the current ``Fillable`` instance.

        Returns
        -------
        Self
            A copy of the current instance.
        """
        return deepcopy(self)

    def get_geometry_map(self, geom_type: GeometryType) -> Compound:
        """
        Method that gets the ``Compound`` object associated to the given
        ``GeometryType`` for the current instance. It defaults to the
        compound made by collecting all the mappings between the given
        ``GeometryType`` and the corresponding ``Compound`` object for every
        ``Fillable`` in the hierarchy tree of this instance.
        If an empty object is found, an exception is raised.

        Parameters
        ----------
        geom_type : GeometryType
            The type of geometry to look for.

        Returns
        -------
        Compound
            Grouping all the mappings between the given ``GeometryType`` and
            the corresponding ``Compound`` object for every ``Fillable`` in
            the hierarchy tree of this instance.

        Raises
        ------
        RuntimeError
            If no ``Compound`` object is associated to the indicated
            ``GeometryType`` element.
        """
        # Get the compound associated to the geometry type, it defaults to the
        # compound collecting those of its layouts in the layers
        geom_map = self.geometry_maps.get(
            geom_type,
            wrap_shape(
                make_compound(
                    list(self._iterate_over_geom_mappings(geom_type))
                )
            )
        )
        # Raise an exception if the 'Compound' has no edges
        if not extract_sub_shapes(geom_map, ShapeType.EDGE):
            raise RuntimeError(
                f"No mapping from '{geom_type}' to geometry layout is "
                "currently present."
            )
        return geom_map

    def get_regions(self) -> List[Region]:
        """
        Method that collects and returns all the ``Region`` objects from the
        layers of the current instance.
        For each layer, it iterates through its layouts and collect any
        ``Region`` instances directly, or recursively calls this method for
        each ``Fillable`` object found along the hierarchy.
        The result of the collection is a flat list of all the ``Region``
        objects that are representative of the technological geometry layout.

        Returns
        ----------
        List[Region]
            A list containing the ``Region`` objects representative of the
            technological geometry layout.
        """
        return [
            r
            for layer in self.layers
            for layout in layer
            for r in (
                [layout] if isinstance(layout, Region)
                else layout.get_regions() if isinstance(layout, Fillable)
                else []
            )
        ]

    def get_regions_with_symmetry(
            self, symmetry: SymmetryType
        ) -> List[Region]:
        """
        Method that collects and returns all the ``Region`` objects from the
        layers of the current instance that are representative of the
        technological geometry layout.
        If any symmetry type other than ``SymmetryType.FULL`` is given, only
        the regions (or their portions) in common with the shape identifying
        the symmetry type are returned.
        If no shape is stored for the indicated symmetry type, an exception
        is raised.

        Parameters
        ----------
        symmetry : SymmetryType
            The type of symmetry for which the corresponding ``Region``
            objects of the layout are returned.

        Returns
        ----------
        List[Region]
            A list containing the ``Region`` objects representative of the
            technological geometry layout according to the given type of
            symmetry.

        Raises
        ------
        RuntimeError
            If no shape is stored in the current instance for the indicated
            symmetry type.
        """
        # Get the regions of the full technological layout
        regions = self.get_regions()
        # Return immediately if the provided symmetry is FULL
        if symmetry == SymmetryType.FULL:
            return regions
        # Get the shape of the symmetry. Raise an exception if not present
        if symmetry not in self.symmetry_map:
            raise RuntimeError(
                f"The indicated '{symmetry.name}' type of symmetry is "
                f"currently not available for '{self.name}'. To enable it, "
                "please call the method 'apply_symmetry()' first with the "
                "desired symmetry type."
            )
        symm_shape = self.symmetry_map[symmetry]
        # Return the regions in common with the shape of the symmetry
        common_regions = []
        for region in regions:
            common = wrap_shape(make_common(region, symm_shape))
            # Proceed only if the common part is still a face object
            if isinstance(common, Face):
                r = region.clone()
                r.update(common)
                common_regions.append(r)
        return common_regions

    def print_region_info(self, shape: Any | None = None) -> None:
        """
        Method that retrieves and prints descriptive information about the
        region of the layout this instance refers to according to the given
        GEOM face, if any, otherwise the one currently selected in the SALOME
        study.
        The printed information is related to:

        - the tree representation of the hierarchical structure up to the
          target ``Region`` object;
        - the values for each of the ``PropertyType`` items associated to the
          target ``Region`` object.

        Parameters
        ----------
        shape : Any | None = None
            The GEOM face indicating the region belonging to the technological
            geometry of the layout this instance refers to and whose
            information have to be retrieved and printed. When not provided,
            the GEOM face currently selected in the SALOME study is
            considered.

        Raises
        ------
        RuntimeError
            When no GEOM face, or more than one, is selected.
        RuntimeError
            If the given GEOM face does not correspond to any of the
            ``Region`` objects belonging to the hierarichal tree of this
            instance.
        """
        # Extract the geometrical object currently selected in the current
        # SALOME study, if any
        if shape is None:
            shape = retrieve_selected_object(
                "Please, select a single region whose data to show."
            )
        # Get the path to the region that corresponds to the given shape, if
        # any can be found
        path = find_region_path_in_tree(self.layers, shape)
        if not path:
            raise RuntimeError(
                f"The indicated region, named '{get_shape_name(shape)}', "
                "could not be found among any of the regions along the "
                f"hierarchical structure of the layout named '{self.name}'."
            )
        # Get the region by following the path up to it
        region = follow_path_to_node(self, path)

        # Print information about the subtree from the root to the target
        print_subtree(self, region, path)

        # Print the information about the properties of the region
        print(f"Properties of '{region.name}':")
        if not region.properties:
            print("   No associated properties.")
            return
        for prop_type, value in region.properties.items():
            print(f"   {prop_type.name}: {value}\n")

    def restore(self) -> None:
        """
        Method that restores the geometry layout of this instance.
        It removes all the layouts stored in the layers while keeping the
        characteristic shape intact.
        Any symmetry and geometry types mappings are completely cleared.
        The result is that the layout is set to a GEOM face built from the
        borders of the current layout, while the ``layers`` and ``regions``
        lists are filled only by the corresponding ``Region`` without any
        associated property.
        """
        # Clear all the mappings and the lists
        self.geometry_maps.clear()
        self.symmetry_map.clear()
        self.layers.clear()
        self.regions.clear()
        # Update the GEOM object this instance refers to by building a face
        # from the borders of the current layout
        boundaries = get_closed_free_boundary(self.geom_obj)
        if len(boundaries) > 1:
            boundaries = get_closed_free_boundary(
                make_partition([self.geom_obj], [], ShapeType.FACE)
            )
        region = Region(make_face(boundaries))
        self.update(region)
        # Update the list of regions
        self.regions = [region]
        # Update the layers with only one containing the region based on the
        # restored layout
        self.layers = [[region]]

    def rotate(self, angle: float, axis: Edge | None = None) -> None:
        """
        Method for rotating the layout by the given angle (in degrees) around
        the given axis, if any is provided, otherwise around the axis
        perpendicular to the layout and passing through its centre.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        axis : Edge | None = None
            The ``Edge`` object representing the rotation axis, if any.
        """
        # Return immediately if the angle is zero
        if math.isclose(angle, 0.0, abs_tol=1e-6):
            return
        if self.geom_obj is None:
            self.update_hierarchical_structure()
        # Build a Z-axis, if none is provided
        if not axis:
            # Build the Z-axis of rotation positioned in the figure center
            axis = wrap_shape(build_z_axis_from_vertex(self.o))
        # Rotate the surface elements
        self._rotate_from_axis(angle, axis)

    def scale(self, factor: float, origin: Vertex | None = None) -> None:
        """
        Method for scaling the layout by the given factor wrt the ``Vertex``
        object, if any. If no reference vertex is provided, the scaling is
        performed wrt the centre of the region.
        The GEOM object this instance refers to, all the layout objects in
        each layer, and the mappings, are scaled by the same factor.

        Parameters
        ----------
        factor : float
            The scaling factor.
        origin : Vertex | None = None
            Identifying the point wrt the scaling is performed. If ``None``,
            the reference point is the region's centre.

        Raises
        ------
        ValueError
            If the scaling factor is less than or equal to zero.
        """
        # Check the validity of the scaling factor
        if factor <= 0.0:
            raise ValueError(
                f"The indicated scaling factor of {factor} is not valid. "
                "Please, provide a value greater than zero."
            )
        if self.geom_obj is None:
            self.update_hierarchical_structure()
        # Perform the scaling wrt to the given origin, otherwise the layout's
        # centre
        if origin is None:
            origin = self.o
        # Scale the GEOM object this instance refers to together with the
        # layouts in the layers
        self.geom_obj = make_scale(self, origin, factor)
        for layer in self.layers:
            for layout in layer:
                layout.scale(factor, origin)

        # Scale the Compound objects in the mappings
        self.geometry_maps.update({
            geom_type: wrap_shape(make_scale(layout, origin, factor))
            for geom_type, layout in self.geometry_maps.items()
        })
        self.symmetry_map.update({
            symm_type: wrap_shape(make_scale(layout, origin, factor))
            for symm_type, layout in self.symmetry_map.items()
        })

    def set_region_properties(
            self,
            properties: Dict[PropertyType, str],
            region: Any | None = None
        ) -> None:
        """
        Method that enables setting the values of the given property types
        for the region of the layout whose GEOM face is passed as input.
        Alternatively, if no one is provided, the method gets the GEOM face
        currently selected in the SALOME study, if any.
        The whole hierarchical tree of the layout is traversed to find the
        ``Region`` object that matches the GEOM face and its properties
        updated with the given dictionary.

        Parameters
        ----------
        properties : Dict[PropertyType, str]
            The dictionary associating the ``PropertyType`` items to the
            corresponding values with which the properties of the region
            have to be updated.
        region : Any | None = None
            The GEOM face indicating the region belonging to the technological
            geometry of the layout this instance refers to and whose
            properties have to be updated. When not provided, the GEOM face
            currently selected in the SALOME study is considered.

        Raises
        ------
        RuntimeError
            When no GEOM face, or more than one, is selected.
        RuntimeError
            If the given GEOM face does not correspond to any of the
            ``Region`` objects belonging to the hierarichal tree of this
            instance.
        """
        # If no geometric object is passed, get the one currently selected in
        # the SALOME study, if any
        if region is None:
            region = retrieve_selected_object(
                "Please, select a single region."
            )

        # Get the 'Region' object that corresponds to the GEOM face by
        # traversing the tree up to the region, if any
        found_region = get_region_in_tree(self.layers, region)
        if not found_region:
            raise RuntimeError(
                f"The indicated region, named '{get_shape_name(region)}', "
                "could not be found among any of the regions along the "
                f"hierarchical structure of the layout named '{self.name}'."
            )
        # Set the properties of the found 'Region' object
        if found_region.properties is None:
            found_region.properties = properties
        else:
            found_region.properties.update(properties)

    def show(self, *args: Any) -> None:
        """
        Method for displaying the geometry layout in the 3D viewer of SALOME
        according to the ``PropertyType`` and ``GeometryType`` settings,
        provided as input regardless of their order.
        Regions representative of the technological geometry layout are
        displayed using a color map according to the values of the indicated
        ``PropertyType`` associated to them. If no ``PropertyType`` is
        specified, the regions are displayed without any color.

        In addition, the displayed regions correspond to those of the last
        applied ``SymmetryType``. Given the regions of the entire layout, a
        common operation with the shape of the symmetry provides the regions
        to display.

        By default this method displays the regions of the technological
        geometry. If a different ``GeometryType`` is given, the method also
        displays the GEOM compound of the edges describing the layout's
        refined geometry.
        If specified differently, this method only adds the GEOM face objects
        of the regions of technological geometry in the Object Browser of
        SALOME without triggering their visualization. This can be used to
        avoid overheads when layouts characterised by lots of face objects
        are shown in SALOME.

        Parameters
        ----------
        *args : Any
            Positional arguments providing the display settings.

        Raises
        ------
        RuntimeError
            If any parameter other than the expected ones in provided.
        RuntimeError
            If any displayed ``Region`` object does not declare a value for
            the indicated ``PropertyType``, if any.
        RuntimeError
            If this instance does not have any ``Compound`` object associated
            to the given ``GeometryType``.

        Notes
        -----
        Allowed arguments to this method are:

        - ``PropertyType``, used to derive the colour map associated to the
          regions of the technological geometry;
        - ``GeometryType``, used to indicate the type of geometry of this
          layout.
        - boolean flag, used to avoid triggering the automatic visualisation
          of the regions and the sectorised edges, if ``False``.
        """
        # Initialize the display settings with default values
        prop_type = None
        geom_type = GeometryType.TECHNOLOGICAL
        display_regions = True
        # Extract the display settings from the parameters of the method
        for setting in args:
            if isinstance(setting, PropertyType):
                prop_type = setting
            elif isinstance(setting, GeometryType):
                geom_type = setting
            elif isinstance(setting, bool):
                display_regions = setting
            else:
                raise RuntimeError(f"Unknown '{setting}'.")

        # Erase all objects from the current view
        clear_view()
        # If the GEOM compound is already present in the current SALOME study,
        # remove it
        if self.entry_id and get_object_from_id(self.entry_id):
            remove_from_study(self.entry_id)
        # Update the layout's tree, if needed
        self.update_hierarchical_structure()
        # Add the updated GEOM compound of the layout to the study
        self.entry_id = add_to_study(self.geom_obj, self.name)
        logging.info(
            f"{self.__class__.__name__} - ({self.name}): "
            + "Added GEOM object to study"
        )

        # Collect the regions according to the applied symmetry, if any, and
        # display them with colour map associated to the property type
        self.regions = self.get_regions_with_symmetry(
            self.state.symmetry_type
        )
        logging.info(
            f"{self.__class__.__name__} - ({self.name}): "
            + f"Got {len(self.regions)} regions of the layout"
        )
        try:
            # Assign the same color to all the regions having the same
            # property value, if any has been specified to show
            associate_colors_to_regions(prop_type, self.regions)
            logging.info(
                f"{self.__class__.__name__} - ({self.name}): "
                + f"Associated {prop_type.name if prop_type else 'None'} "
                + "color map to regions"
            )
            # Add all the regions of the layout to the study
            self._show_regions(display_regions)
        except RuntimeError as e:
            # Add all the regions of the layout to the study, with the faulty
            # ones (i.e. those without the property) coloured in red
            self._show_regions(display_regions)
            # Re-raise the exception
            raise RuntimeError(
                "Error while displaying the regions of the geometry layout."
            ) from e

        # Handle the visualization of the edges of the selected geometry type,
        # if different from the technological one
        try:
            self._show_geometry_type_edges(geom_type)
        except RuntimeError as e:
            # Re-raise the exception
            raise RuntimeError(
                "Error while displaying the edges associated to the "
                f"indicated '{geom_type.name}' geometry type."
            ) from e

        # Update the state of the layout
        self.state.displayed_geom = geom_type
        self.state.is_update_needed = False

    def translate(self, new_cntr: Tuple[float, float, float]) -> None:
        """
        Method for translating the geometric elements of the current layout,
        i.e. the corresponding GEOM compound, the centre and the layouts in
        each layer. If any of the layouts of a layer is an instance of the
        subclasses of ``Fillable``, this method is applied recursively.

        Parameters
        ----------
        new_cntr : Tuple[float, float, float]
            The XYZ coordinates of the new center of the layout.
        """
        # Return immediately if the new centre coincides with the current one
        new_cntr_vrtx = make_vertex(new_cntr)
        if are_same_shapes(self.o, new_cntr_vrtx, ShapeType.VERTEX):
            return
        if self.geom_obj is None:
            self.update_hierarchical_structure()
        # Build a vector from the current center to the new one
        transl_vect = make_vector_from_points(self.o, new_cntr_vrtx)
        prev_o = self.o
        # Translate all the characteristic geometric elements
        self.o = wrap_shape(make_translation(self.o, transl_vect))
        self._translate_layout_specific_elems(new_cntr)
        # Translate the layouts in each layer so that its relative position
        # wrt the translated centre is kept
        for layer in self.layers:
            for layout in layer:
                layout_center = compute_point_by_reference(
                    layout.o, prev_o, new_cntr
                )
                layout.translate(layout_center)

        # Update the GEOM compounds representing the different geometry
        # layout types
        self.geom_obj = make_translation(self, transl_vect)
        self.geometry_maps.update({
            geom_type: wrap_shape(make_translation(layout, transl_vect))
            for geom_type, layout in self.geometry_maps.items()
        })
        self.symmetry_map.update({
            symm_type: wrap_shape(make_translation(layout, transl_vect))
            for symm_type, layout in self.symmetry_map.items()
        })

    def update(self, layout: Compound | Face) -> None:
        """
        Method for updating the current instance by updating the GEOM object
        this instance refers to, the centre of the layout and any mapping
        between ``GeometryType`` values and ``Compound`` objects.

        Parameters
        ----------
        layout : Compound | Face
            The new layout to update the current instance with.

        Notes
        -----
        This method does not update the layouts stored in the ``layers``
        instance attribute accordingly with the GEOM object of this instance.
        """
        self.geom_obj = make_compound(layout.geom_obj)
        self.o = wrap_shape(make_cdg(layout))
        # Update the mappings
        self.geometry_maps = {
            geom_type: make_common(
                self.geometry_maps[geom_type], self.geom_obj
            ) for geom_type in self.geometry_maps
        }

    def update_hierarchical_structure(
            self, collapse_layers: bool = False
        ) -> None:
        """
        Method that processes the layers of the current instance object to
        update the entire geometric hierarchical structure which represents
        the technological geometry layout.
        This method collapses all the layers by traversing the entire tree
        in reverse order, and cutting the ``Region`` or the ``Fillable``
        objects that are overlapped.
        This operation is skipped and the method exits without changes if
        there is no need to update the technological geometry layout (i.e.
        the ``is_update_needed`` value of the ``state`` attribute is
        ``False``).

        The operation of assembling all the layers requires that each layout
        is up-to-date. This means that each ``Fillable`` object in the list
        of layers is further processed by recursively calling this method to
        collapse its own layers.

        Finally, the GEOM compound object representative of this instance is
        updated by collecting in a GEOM compound all the ``Region`` objects
        retrieved from the layouts in each layer.

        To simplify the hierarchical tree, a ``True`` value can be provided.
        If so, the hierarchical tree is collapsed and the layers reduced to a
        sigle layer collecting all the ``Region`` objects in the tree.

        Parameters
        ----------
        collapse_layers : bool = False
            Flag stating whether the current hierarchical tree should be
            collapsed and reduced to a list of ``Region`` objects only.
        """
        # Return immediately if there is no need to update the layout
        if not self.state.is_update_needed:
            logging.info(
                f"{self.__class__.__name__} - ({self.name}): "
                "No need to update...")
            return
        # Reverse the layers
        reversed_layers = self.layers[::-1]
        # Ensure the layout is up-to-date for each 'Fillable' object by
        # recursively calling this method
        for layer in reversed_layers:
            for layout in layer:
                if isinstance(layout, Fillable):
                    layout.update_hierarchical_structure()
        # Collapse all the layers and cut out the overlapping layouts
        if len(reversed_layers) > 1:
            self._collapse_layers(reversed_layers)
        # Collapse the hierarchical tree to only one layer made of regions
        if collapse_layers:
            # Store any geometry map present in the 'Fillable' objects along
            # the hierarchical tree as collapsing it will lose any reference
            # to them. If no mapping is present for a geometry type, the
            # exception is caught but not re-raised as not an error.
            for geom_type in GeometryType:
                try:
                    if geom_type == GeometryType.TECHNOLOGICAL:
                        continue
                    # Update the mapping of the current instance
                    self.geometry_maps[geom_type] = self.get_geometry_map(
                        geom_type
                    )
                except RuntimeError:
                    continue
            self.layers = [self.get_regions()]

        # Update the GEOM compound object of this instance
        self.geom_obj = make_compound(self.get_regions())
        # Update the state of the layout
        self.state.is_update_needed = False

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
          recursively applied to its layers;
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
                layout.update(wrap_shape(cut_layout))
            elif isinstance(layout, Fillable):
                layout.update(wrap_shape(cut_layout))
                # Recursively apply the cut to each layer of the layout
                for l in layout.layers:
                    self._apply_cut_to_layouts_in_layer(l, cutting_tool)

        # Remove the layouts completely overlapped by the superior layer
        for index in sorted(layouts_to_remove, reverse=True):
            layer.pop(index)

    @use_symmetry_logic(build_common_symmetry_shape)
    def _build_symmetry_shape(
            self,
            symmetry_type: SymmetryType,
            context: SymmetryDomain
        ) -> Surface:
        """
        Method that builds the geometric shape that corresponds to the given
        symmetry type and domain for a generic layout.
        This method supports only the symmetry types that are common to all
        layouts, indipendently from their characteristic shape.

        Supported symmetries are:

        - FULL: returns the characteristic shape of the layout.
        - HALF: builds a rectangle representing half of the layout.
        - QUARTER: builds a rectangle representing one quarter of the layout.

        Parameters
        ----------
        symmetry : SymmetryType
            The symmetry type for which the characteristic shape is built.
        domain : SymmetryDomain
            Instance providing the domain of the full layout in terms of
            XY-bounding extents, full shape and its centre.

        Returns
        -------
        Surface
            A geometric shape representing the requested symmetry.

        Raises
        ------
        RuntimeError
            If the indicated symmetry type is not supported for a generic
            layout.

        Notes
        -----
        The method is decorated so that it calls the function handling the
        construction of the symmetry shape for generic-type layouts. For this
        reason, no implementation is included here.
        """
        pass

    def _collapse_layers(self, layers: List[List[Region | Self]]) -> None:
        """
        Method that collapses all the given layers (a list containing lists
        of ``Region`` or ``Fillable`` objects).
        If any compound object of a layer overlaps any layer below it, the
        layouts contained in the overlapped layer are cut by the layer
        compound at a higher level.

        Notes
        -----
        Layout objects belonging to the same layer should not overlap each
        other.
        """
        logging.info(
            f"{self.__class__.__name__} - ({self.name}): "
            "Collapsing its layers...")
        # Loop through all the layers
        for i, layer in enumerate(layers):
            # Skip the layer, if empty
            if not layer:
                continue
            # Loop through all the layers below the current one to cut all the
            # layout objects of each sublayer, if any is overlapped.
            for sub_layer in layers[i + 1:]:
                # Skip the layer, if empty
                if not sub_layer:
                    continue
                # Overlap the current layer onto the layers below
                self._overlap_layer_to(layer, sub_layer)

    def _iterate_over_geom_mappings(
            self, geom_type: GeometryType) -> Iterator[Compound]:
        """
        Method that traverses the hierarchical structure of this instance to
        collect all the geometry mappings for the given geometry type.

        This method recursively descends into all nested ``Fillable`` objects
        contained in this instance's ``layers`` attribute. For each nested
        ``Fillable``, the needed geometry mapping is yielded before yielding
        the one of the current object (i.e. a post-order traversal is
        adopted).

        Parameters
        ----------
        geom_type : GeometryType
            The geometry type whose compounds should be returned.

        Yields
        ------
        Compound
            The flattened sequence of ``Compound`` objects associated with
            the indicated ``GeometryType`` for all descendant ``Fillable``
            objects, followed by the one of the current instance.
        """
        # Recurse into 'Fillable' objects of each layer
        for layer in self.layers:
            for layout in layer:
                if isinstance(layout, Fillable):
                    # Depth-first into nested 'Fillable' objects
                    yield from layout._iterate_over_geom_mappings(geom_type)

        # Yield this Fillable's compound after children
        yield from flatten_list(self.geometry_maps.get(geom_type, []))

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
        # Build the boundaries of the layer compound
        layer_cmpd = make_compound(layer)
        boundaries = get_closed_free_boundary(layer_cmpd)
        if len(boundaries) > 1:
            boundaries = get_closed_free_boundary(
                make_partition(layer, [], ShapeType.FACE)
            )
        # Try to build a planar face from the extracted boundaries; if the
        # operation raises an exception, fall back to the given compound
        try:
            layer_shape = make_face(boundaries)
        except:
            layer_shape = layer_cmpd
        # Apply the cut on the layout objects of the sub-layer
        self._apply_cut_to_layouts_in_layer(sub_layer, layer_shape)

    def _rotate_from_axis(self, angle: float, axis: Edge) -> None:
        """
        Method for rotating the layout and its geometric elements by the
        given angle (in degrees) around the given axis.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        axis : Edge
            The ``Edge`` object representing the rotation axis.
        """
        self.rot_angle += angle
        # Convert the rotation angle in radians
        rot_angle = math.radians(angle)
        # Rotate the layouts in each layer araound the given axis
        for layer in self.layers:
            for layout in layer:
                layout.rotate(angle, axis)

        # Update the GEOM compounds representing the different geometry
        # layout types
        self.geom_obj = make_rotation(self, axis, rot_angle)
        self.o = make_rotation(self.o, axis, rot_angle)
        self.geometry_maps.update({
            geom_type: wrap_shape(make_rotation(layout, axis, rot_angle))
            for geom_type, layout in self.geometry_maps.items()
        })
        self.symmetry_map.update({
            symm_type: wrap_shape(make_rotation(layout, axis, rot_angle))
            for symm_type, layout in self.symmetry_map.items()
        })

    def _show_geometry_type_edges(self, geom_type: GeometryType) -> None:
        """
        Method that handles the visualization of the edges of the given
        geometry type, if different from the technological one.
        It gets the ``Compound`` object associated to the ``GeometryType``,
        extracting the common part with the shape of the symmetry currently
        applied, and it adds it to the SALOME study.
        In the Object Browser, the added compound is available as child of
        the GEOM compound representative of the technological geometry of the
        layout this instance refers to.

        Parameters
        ----------
        geom_type : GeometryType
            The geometry type whose ``Compound`` object to display.

        Raises
        ------
        RuntimeError
            If this instance does not have any ``Compound`` object associated
            to the given ``GeometryType``.
        """
        # Return if the given geometry type is the technological one
        if geom_type == GeometryType.TECHNOLOGICAL:
            return
        # Get the 'Compound' object associated to the indicated geometry type
        geom_map = self.get_geometry_map(geom_type)
        # Update the mapping of the current instance
        self.geometry_maps[geom_type] = geom_map
        # Get only the portion in common with the symmetry shape
        if self.state.symmetry_type != SymmetryType.FULL:
            geom_map = make_common(
                geom_map, self.symmetry_map[self.state.symmetry_type]
            )
        # Display the compound by colouring the edges with black color
        set_color_face(geom_map, (0,0,0))
        display_shape(
            add_to_study_in_father(
                self.geom_obj, geom_map, f"{geom_type.name} edges"
            )
        )
        # Update the SALOME's study to display all the GEOM compound object
        update_salome_study()

    def _show_regions(self, display: bool = True) -> None:
        """
        Method that adds all the regions of the layout to the current SALOME
        study. In the Object Browser they are available as children of the
        GEOM compound representative of the technological geometry of the
        layout.
        Each region is displayed with a colour defined beforehands.

        Parameters
        ----------
        display : bool = True
            If ``False``, the GEOM face objects of the regions are not
            automatically displayed in the viewer of SALOME. Defaults to
            ``True``.
        """
        for i, region in enumerate(self.regions):
            # Set the region color in the viewer
            set_color_face(region.geom_obj, region.color)
            # Add the layout region to the study
            region.entry_id = add_to_study_in_father(
                self.geom_obj, region, region.name)
            logging.info(
                f"{self.__class__.__name__} - ({self.name}): "
                + f"Added {region.name} (no. {i}) to the study"
            )
        # Display the GEOM face objects of the regions, if needed
        if display:
            for region in self.regions:
                display_shape(region.entry_id)
        # Update the SALOME's study to display all the GEOM objects
        update_salome_study()

    @abstractmethod
    def _translate_layout_specific_elems(
        self, new_cntr: Tuple[float, float, float]
    ) -> None:
        """
        Abstract method for translating layout-specific elements, considering
        the centre of the translated layout is positioned at the given XYZ
        coordinates.

        Parameters
        ----------
        new_cntr : Tuple[float, float, float]
            The XYZ coordinates of the centre of the translated layout.
        """


# -------------------------------------------------------------------------- #
#                                FUNCTIONS                                   #
# -------------------------------------------------------------------------- #

def find_region_path_in_tree(
        layers: List[List[Region | Fillable]],
        shape: Any,
        path: List[Tuple[int, int]] | None = None
    ) -> List[Tuple[int, int]] | None:
    """
    Function that traverses the given hierarchical tree as a list of possibly
    nested lists to find the path to the ``Region`` object that matches the
    given shape.
    A positive match occurs either if their shapes are identical or if a
    vertex built within the input shape is within a ``Region`` object in the
    tree. If no match is found, ``None`` is returned.

    For each ``Fillable`` object that contains the given GEOM face, the
    indices of the list of lists of its parent ``Fillable`` in which the
    current one is found are stored.
    The final result is a list of tuples of indices that helps to reconstruct
    the entire path from the root to the matching ``Region``.

    Parameters
    ----------
    layers : List[List[Region | Fillable]]
        The hierarchical tree as list of nested lists of ``Region`` (leaves)
        and ``Fillable`` objects (nodes).
    shape : Any
        The GEOM face whose corresponding ``Region`` has to be found.
    path : List[Tuple[int, int]] | None = None
        The path storing the indices to progressively reach the parent of the
        current ``Fillable`` object. It defaults to ``None``, indicating the
        beginning of the search.

    Returns
    -------
    List[Tuple[int, int]] | None
        The list of tuples of indices used in each of the nested lists from
        the root ``Fillable`` to the matching ``Region``. If no match is found
        ``None`` is returned.
    """
    # Initialize the list of indices if not passed
    if path is None:
        path = []
    # Traverse the list of lists of the current node
    for layer_i, layer in enumerate(layers):
        for elem_i, layout in enumerate(layer):
            # Skip layouts that do not contain the shape
            if not is_layout_contained(layout, shape):
                continue
            # Update the path with the indices to the layout containing the
            # given shape
            path += [(layer_i, elem_i)]
            # Recursively call this function to descend deeper in the tree, if
            # the container is a 'Fillable'
            if isinstance(layout, Fillable):
                return find_region_path_in_tree(
                    layout.layers,
                    shape,
                    path
                )
            else:
                # Check if a match with any of the regions of the node occurs
                if (are_same_shapes(shape, layout, ShapeType.FACE)
                    or is_point_inside_shape(
                        make_vertex_inside_face(shape), layout
                    )
                ):
                    # Return the path to the found region
                    return path
    # Return None if no matches are found in this tree
    return None


def follow_path_to_node(
        root: Fillable, path: List[Tuple[int, int]]
    ) -> Region | Fillable:
    """
    Function that traverses the nested hierarchical trees starting from the
    root ``Fillable`` object up to the node that results from following the
    provided path.

    Parameters
    ----------
    root : Fillable
        Root node of the hierarchical tree.
    path : List[Tuple[int, int]]
        Sequence of pairs of indices describing the path up to the node.

    Returns
    -------
    Region | Fillable
        The layout (either a ``Region`` or a ``Fillable`` object) at the end
        of the path.

    Raises
    ------
    IndexError
        If a layer or element index is out of range.
    """
    # Set the root as the starting node
    node = root
    # Traverse the tree from the root object
    for depth, (i, j) in enumerate(path):
        # Return if the current node is a 'Region' object
        if isinstance(node, Region):
            return node
        # Try accessing the node at the current indices
        try:
            node = node.layers[i][j]
        except IndexError as e:
            raise IndexError(
                f"Invalid path at depth {depth}: "
                f"layers[{i}][{j}] is out of range."
            ) from e

    return node


def get_region_in_tree(
        tree: List[List[Region | Fillable]],
        shape: Any,
    ) -> Region | None:
    """
    Function that traverses the given hierarchical tree as a list of possibly
    nested lists to find the ``Region`` object that matches the given shape.
    A positive match occurs either if their shapes are identical or if a
    vertex built within the input shape is within a ``Region`` object in the
    tree. If no match is found, ``None`` is returned.

    Parameters
    ----------
    tree : List[List[Region | Fillable]]
        The hierarchical tree as list of nested lists of ``Region`` (leaves)
        and ``Fillable`` objects (nodes).
    shape : Any
        The GEOM face whose corresponding ``Region`` has to be found.

    Returns
    -------
    Region | None
        The ``Region`` contained in the hierarchical tree matching the given
        shape, or ``None`` if no one is found.
    """
    # Get the region that corresponds to the given shape
    for subtree in tree:
        for layout in subtree:
            # Continue with the next layout if the current one does not
            # contain the shape
            if not is_layout_contained(layout, shape):
                continue
            if isinstance(layout, Fillable):
                # Return the found region, if any, by traversing the tree of
                # the current layout
                return get_region_in_tree(layout.layers, shape)
            else:
                # Return the region if the matching conditions are met
                if (are_same_shapes(shape, layout, ShapeType.FACE)
                    or is_point_inside_shape(
                        make_vertex_inside_face(shape), layout)
                ):
                    return layout
    # No region has been found, hence None is returned
    return None


def print_subtree(
        root: Fillable,
        target: Region | Fillable,
        path: List[Tuple[int, int]]
    ) -> None:
    """
    Function that prints to console the portion of the hierarchical tree from
    the given root ``Fillable`` object to the target (either a ``Region`` or
    a ``Fillable`` object) according to the path containing the list of tuples
    of indices allowing to traverse the tree up to the target.

    Parameters
    ----------
    root : Fillable
        The root ``Fillable`` object.
    target: Region | Fillable
        The target (either a ``Region`` or a ``Fillable`` object) the path
        points to.
    path: List[Tuple[int, int]]
        The sequence of tuples of indices that allows to traverse the tree
        from root to target.
    """
    print(
        f"\n--- TREE PATH FROM {root.__class__.__name__} '{root.name}' "
        f"TO {target.__class__.__name__} '{target.name}' ---\n"
    )
    # Print visual indented tree and the access command to target
    print(f"{root.name} (root)")
    access_code = "<root-instance-name>"
    for depth, (l, i) in enumerate(path, 1):
        indent = "  " * depth
        print(
            f"{indent}└── layers[{l}][{i}] -> "
            f"{follow_path_to_node(root, path[:depth]).name}"
        )
        access_code += f".layers[{l}][{i}]"
    print(f"\nAccess Code: {access_code}\n")
