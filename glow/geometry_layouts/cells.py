"""
Module containing classes providing the means for creating cells from basic
geometries in SALOME.
"""
import math

from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Self, Tuple, Union

from glow.geometry_layouts.geometries import GenericSurface, Circle, \
    Rectangle, Hexagon, update_relative_pos
from glow.interface.geom_interface import ShapeType, add_to_study, \
    add_to_study_in_father, clear_view, display_shape, extract_sub_shapes, \
    get_kind_of_shape, get_min_distance, get_object_from_id, \
    get_point_coordinates, get_selected_object, get_shape_name, \
    is_point_inside_shape, make_cdg, make_common, make_cut, make_fuse, \
    make_line, make_partition, make_rotation, make_translation, \
    make_vector_from_points, make_vertex, make_vertex_inside_face, \
    make_vertex_on_curve, make_vertex_on_lines_intersection, \
    remove_from_study, set_color_face, update_salome_study
from glow.generator.support import CellType, GeometryType, PropertyType, \
    generate_unique_random_colors


@dataclass
class Region():
    """
    Grouping information regarding a generic region of the cell
    technological geometry.

    Attributes
    ----------
    - cdg           : Union[Any, None] = field(init=False)
                      A GEOM vertex object representing the cell region
                      CDG, if any has been set
    - face          : Any = None
                      A GEOM face object of the cell region
    - face_entry_id : Union[str, None] = None
                      The ID of the cell region used in the current SALOME
                      study
    - inner_point   : Union[Any, None] = None
                      A GEOM vertex object representing a point within
                      the cell region
    - name          : Union[str, None] = None
                      The name associated to the cell region used in the
                      current SALOME study
    - properties    : Union[Dict[PropertyType, str], None] = None
                      A dictionary of property types and names associated
                      to the cell region
    - colors        : Union[Dict[PropertyType, Tuple[int, int, int]], None]
                      A dictionary associating an RGB color to each property
                      type of the cell region
    """
    cdg: Union[Any, None] = field(init=False)
    face: Any = None
    face_entry_id: Union[str, None] = None
    inner_point: Union[Any, None] = None
    name: Union[str, None] = None
    properties: Union[Dict[PropertyType, str], None] = None
    color: Tuple[int, int, int] = field(init=False)

    def __post_init__(self) -> None:
        """
        Dataclasses specific method called automatically after the class
        initialization for setting attributes not required at class
        instantiation.
        It calculates the region CDG, given the face object, if any; in
        addition, it associates a fixed color for all the property types
        set at instantiation.
        """
        self.cdg = self.__calculate_cdg()
        self.color = (167, 167, 167)
        # # Associate a color to each property type, if any has been set
        # if self.properties:
        #     rgb_colors = [(255, 0, 0)] * len(self.properties.keys())
        #     # Associate the same color for each property type
        #     self.colors = dict(zip(self.properties.keys(), rgb_colors))

    def __calculate_cdg(self) -> Union[Any, None]:
        """
        Method for calculating the CDG of the face this sector is referring
        to. If no face has been set yet, the CDG, as a GEOM vertex object,
        is left to None.

        Returns
        -------
        A GEOM point object representing the face CDG, if any face has been
        set; None otherwise.
        """
        if not self.face:
            return None
        return make_cdg(self.face)

    def reset_color_to_default(self):
        """
        Method for resetting the region color to the (167, 167, 167) RGB
        value, which corresponds to light grey.
        """
        self.color = (167, 167, 167)

    def set_property_color(self,
                           color: Tuple[int, int, int]) -> None:
        """
        Method for associating a specific RGB color to the region.

        Parameters
        ----------
        color     : color: Tuple[int, int, int]
                    RGB values identifying the color to associate to the
                    region
        """
        # Check the given color
        if not color or not isinstance(color, Tuple) or len(color) != 3 or \
            any(value < 0 or value > 255 for value in color):
            raise ValueError("The 'color' parameter must be provided as a "
                             "tuple of 3 integers in the 0:255 range.")
        # Store the RGB color
        self.color = color


class GenericCell(ABC):
    """
    Abstract class representative of a generic cell.

    Parameters
    ----------
    center  : Tuple[float, float, float]
        The X-Y-Z coordinates of the cell center
    figure  : GenericSurface
        The figure representing the cell shape
    name    : str
        The name of the cell to be used in the current study

    Attributes
    ----------
    cell_type       : Union[CellType, None]
        The type of cell expressed as a value of the 'CellType' enumeration
    face            : Any
        The face object providing the cell technological geometry
        representation
    face_entry_id   : Union[str, None]
        The ID of the cell face object used in the current study
    figure          : GenericSurface
        The figure representing the cell shape
    inner_circles   : List[Circle]
        The list of 'Circle' objects representing the cell concentric circles
    is_windmill_applied : bool
        Boolean flag stating if the windmill sectorization is applied
    name            : str
        The name of the cell to be used in the current study
    regions         : List[Region]
        The list of 'Region' objects storing information about either the
        regions of the cell technological geometry or the ones resulting from
        its sectorization
    rotation        : float
        The rotation angle of the cell expressed in radians
    sectorized_face : Any
        The face object providing the cell sectorized geometry representation
        (used for calculations)
    tech_geom_props : Dict[Any, Dict[PropertyType, str]]
        A dictionary between the regions of the cell technological geometry
        VS the dictionary of property types VS value associated to each zone
    tech_geom_sect_opts : Dict[Any, Tuple[int, float]]
        A dictionary between the regions of the cell technological geometry
        VS the tuple identifying the associated sectorization options

    VALID_SECTOR_NO_VS_ANGLE  : Dict[int, List[float]]
        Dictionary providing the number of valid sectors each cell zone
        could be subdivided VS the accepted angles the sectorization can
        start from.
    """
    # Number of valid sectors each cell zone could be subdivided VS the
    # accepted angles the sectorization can start from. Values are
    # specific to the cell type
    VALID_SECTOR_NO_VS_ANGLE: Dict[int, List[float]] = {}

    def __init__(self, figure: GenericSurface, name: str) -> None:
        super().__init__()
        # Check that the class attribute providing the number of accepted
        # sectors and angles for the sectorization has been defined
        print("VALID_SECTOR_NO_VS_ANGLE = ", self.VALID_SECTOR_NO_VS_ANGLE)
        if not self.VALID_SECTOR_NO_VS_ANGLE:
            raise ValueError(f"{self.__name__} must define the "
                             "'VALID_SECTOR_NO_VS_ANGLE' attribute.'")
        # Store the geometrical figure this cell represents
        self.figure: GenericSurface = figure
        # Build the cell face from its borders
        self.figure.build_face()
        # Initialize the cell-related instance attributes
        self.sectorized_face: Union[Any, None] = None
        self._initialize_specific_cell()
        self.face: Any = self.figure.face
        self.inner_circles: List[Circle] = []
        self.name: str = name
        # Extract the cell subfaces, if any; otherwise use the cell face for
        # initialization purposes
        subfaces = self.extract_subfaces()
        if not subfaces:
            subfaces = [self.face]
        self.tech_geom_props: Dict[Any, Dict[PropertyType, str]] = {
            region: {} for region in subfaces
        }
        self.tech_geom_sect_opts: Dict[Any, Tuple[int, float]] = {
            region: (1, 0) for region in subfaces
        }
        self.face_entry_id: Union[str, None] = None
        self.rotation: float = 0.0
        self.regions: List[Region] = []
        self.cell_type: Union[CellType, None] = None
        self.is_windmill_applied: bool = False
        self.displayed_geom: GeometryType = GeometryType.TECHNOLOGICAL
        self.added_edges: List[Any] = []

    @abstractmethod
    def _initialize_specific_cell(self) -> None:
        """
        Abstract method for initializing instance attributes that are
        specific for the type of cell subclasses of this class represent.
        """

    @abstractmethod
    def _check_radius_vs_cell_dim(self, radius: float) -> None:
        """
        Abstract method for assessing if the circle, whose radius is given
        as input, can be added to the cell. The radius is checked against
        the characteristic dimensions of the cell.
        Since specific to the cell geometry, each subclass must provide its
        own implementation to this method.

        Parameters
        ----------
        radius  : float
                  The radius of the circle to check against the cell
                  dimensions
        """

    def extract_subfaces(self) -> List[Any]:
        """
        Method that extracts the subfaces the cell face is made by. These
        subfaces are the result of the construction of the cell technological
        geometry.
        In addition, the list is ordered in terms of the distances of each
        subface from the cell center, in ascending order.

        Returns
        -------
        A list of all the subfaces that comprise the cell face, sorted in
        ascending order by means of their distance from the cell center.
        """
        return self.__extract_subfaces_from_face(self.face)

    def __extract_subfaces_from_face(self, face: Any) -> List[Any]:
        """
        Method that extracts the subfaces the given face is made by.
        In addition, it returns a list of these subfaces which is
        ordered in terms of the distances of each subface from the
        cell center, in ascending order.

        Parameters
        ----------
        face  : Any
                A face object the subfaces have to be extracted from

        Returns
        -------
        A list of all the subfaces that comprise the given cell face, sorted
        in ascending order by means of their distance from the cell center.
        """
        subfaces = extract_sub_shapes(face, ShapeType.FACE)
        return sorted(subfaces,
                      key=lambda subface: get_min_distance(self.figure.o,
                                                           subface))

    def add_circle(
            self,
            radius: float = 1.0,
            position: Union[Tuple[float, float, float], None] = None) -> None:
        """
        Method that allows to add a circle to the cell technological geometry
        only, if not already present.
        Given the radius, a corresponding face object is built in the given
        position within the cell, if any is specified; otherwise the circle
        is added in the cell center.
        In any case, a partition operation between the cell face and the new
        circle is performed. The result is a geometric surface that comprises
        both faces.
        In addition, the dictionaries associating the regions of the cell
        technological geometry with the properties and with the sectorization
        options respectively are updated.

        Parameters
        ----------
        radius    : float
                    The radius of the circle to add
        position  : Union[Tuple[float, float, float], None]
                    The X-Y-Z coordinates of the circle center to be added
                    within the cell
        """
        # Check if the given radius is acceptable, i.e. the resulting
        # circle is within the cell
        self._check_radius_vs_cell_dim(radius)
        # Check if there is already one circle in the cell with the same
        # radius (with a 1e-10 tolerance) and same center position
        if self.__check_circle_presence(radius):
            raise RuntimeError(
                f"Another circle with the same radius '{radius}' and in the "
                f"same posizion '{position if position else (0, 0, 0)}' is "
                "already present.")
        # Build a 'Circle' object
        circle = Circle(center=position, radius=radius)
        # Build the corresponding face
        circle.build_face()

        # Add the circle to the cell
        self.__add_circle_to_pos(circle)

    def __add_circle_to_pos(self, circle: Circle) -> None:
        """
        Method for adding the given 'Circle' object to the cell at any
        position (indicated by the 'Circle' attribute) within the cell.
        A partition operation between the cell face and the circle one
        is performed. The result is a geometric surface that comprises
        both faces.
        In addition, the dictionaries associating the regions of the cell
        technological geometry with the properties and with the sectorization
        options respectively are updated.

        Parameters
        ----------
        circle  : Circle
            The 'Circle' object to add to the cell
        """
        # Add the 'Circle' object to the list storing the circles within
        # the cell
        self.inner_circles.append(circle)
        # Sort the list of circles by the size of the radius (lower to greater)
        self.inner_circles = sorted(self.inner_circles,
                                    key=lambda item: item.radius)

        # Build a partition between the cell and the circle faces
        self.face = make_partition(
            [self.face, circle.face], [], ShapeType.FACE)

        # ------------------------------------------------------------------
        # Substitute the region where the circle has been added with the one
        # resulting by cutting it with the added circle. The common part
        # between the substituted subface and the added circle is included
        # as well with a default value.
        # ------------------------------------------------------------------
        # Build a point on the border of the added circle to identify the
        # region to cut
        ref_point = make_vertex_on_curve(circle.borders[0], 0)
        # Substitute the region with its cut version and add the common part
        # with default value in the properties dictionary
        for subface in self.tech_geom_props:
            if is_point_inside_shape(ref_point, subface):
                self.tech_geom_props[
                    make_cut(subface, circle.face)
                        ] = self.tech_geom_props.pop(subface)
                self.tech_geom_props[
                    make_common(subface, circle.face)
                        ] = {}
                break

        # ----------------------------------------------------------------
        # Update the sectorization options dictionary for circles centered
        # in the cell only. In any case, perform a partition between the
        # sectorized face and the added circle, if any sectorization has
        # been applied yet.
        # ----------------------------------------------------------------
        # If no sectorization has been applied yet, return without updating
        # the dictionary
        if not self.sectorized_face:
            return
        # Update the dictionary for cell-centered circles
        if get_min_distance(self.figure.o, circle.o) < 1e-5:
            # Substitute the region with its cut version and add the common
            # part with default value in the sectorization options dictionary,
            # if any has been performed
            for subface in self.tech_geom_sect_opts:
                if is_point_inside_shape(ref_point, subface):
                    self.tech_geom_sect_opts[
                        make_cut(subface, circle.face)
                            ] = self.tech_geom_sect_opts.pop(subface)
                    self.tech_geom_sect_opts[
                        make_common(subface, circle.face)
                            ] = (1, 0)
                    break

            # Extract the regions sorted according to their distance from
            # the cell center
            sorted_faces = sorted(self.tech_geom_sect_opts.keys(),
                key=lambda item: get_min_distance(
                  self.figure.o, make_vertex_inside_face(item)))
            # Extract the sectorization options for each region
            sectors_angles = []
            for region in sorted_faces:
                sectors_angles.append(self.tech_geom_sect_opts[region])
            # Apply the sectorization with the extracted options for each
            # region
            # FIXME simplify this passage by updating the dictionary of
            # sectorization options directly. The partition operation, since
            # valid both in the cell-centered circle and not can then be
            # performed outside of the if-clause.
            self._sectorize_cell([value[0] for value in sectors_angles],
                                 [value[1] for value in sectors_angles],
                                 self.is_windmill_applied)
            # Partition the sectorization result with the externally added
            # edges, if any
            # TODO check the necessity to include an attribute for edges
            # added to the sectorized face from within SALOME. When created
            # in SALOME, they could be directly added to the sectorized face
            if self.added_edges:
                self.sectorized_face = make_partition(
                    [self.sectorized_face],
                    self.added_edges,
                    ShapeType.FACE)
        else:
            # Perform a partition between the sectorized face and the added
            # cicle only, whose position is not in the cell center
            self.sectorized_face = make_partition(
                [self.sectorized_face], [circle.face], ShapeType.FACE)

    def remove_circle(self, radius: float = 1.0) -> None:
        """
        Method that allows to remove a cell-centered circle from the cell
        technological geometry, given the radius.
        All references to the circle in the cell are removed, which means
        that the corresponding 'Circle' object is removed from the
        'inner_circles' list attribute and the whole face object of the
        cell is re-built.
        The dictionaries associating the regions of the cell technological
        geometry with the properties and with the sectorization options
        respectively are updated as well.

        Parameters
        ----------
        radius  : float
                  The radius of the circle to remove
        """
        # TODO instead of passing the radius of the circle, which could result
        # in problems as two circles can have the same radius but different
        # positions, it should be better to retrieve the currently selected
        # object in the SALOME viewer. A check on the actual presence of the
        # circle among the stored ones in necessary, also because the selected
        # object could not be a full circle but an annulus, whose external
        # radius has to be retrieved.
        if not self.inner_circles:
            raise RuntimeError("No circles are present in the cell.")
        # Retrieve the original cell face from the figure it is based on
        self.face = self.figure.face
        # Retrieve the index of the circle in the list of cell circles, if
        # any, but only among the cell centered ones
        indx = self.__check_circle_presence(radius)
        if indx is None:
            raise RuntimeError("No circle could be found within the cell for "
                               f"the given radius '{radius}'")
        # Remove the 'Circle' object from the list
        circle = self.inner_circles.pop(indx)

        # Re-build the partition between the face and the remaining circles
        circle_faces = [circle.face for circle in self.inner_circles]
        self.face = make_partition(
            [self.face] + circle_faces, [], ShapeType.FACE)

        # ----------------------------------------------------------------
        # Substitute the region where the circle has been removed with the
        # one resulting by fusing it with the removed circle
        # ----------------------------------------------------------------
        # Build a point on the border of the removed circle to identify the
        # region to fuse
        ref_point = make_vertex_on_curve(circle.borders[0], 0)
        # Substitute the region with its fused version in the properties
        # dictionary. Its properties are assigned to the fused region.
        for subface in self.tech_geom_props:
            if is_point_inside_shape(ref_point, subface):
                self.tech_geom_props[
                    make_fuse([subface, circle.face])
                        ] = self.tech_geom_props.pop(subface)
                break
        # Substitute the region with its fused version in the sectorization
        # options dictionary, if any has been performed. Its options are
        # assigned to the fused region.
        if self.sectorized_face:
            for subface in self.tech_geom_sect_opts:
                if is_point_inside_shape(ref_point, subface):
                    self.tech_geom_sect_opts[
                        make_fuse([subface, circle.face])
                            ] = self.tech_geom_sect_opts.pop(subface)
                    break

        # Re-build the regions-properties/sectorization options association
        # and re-apply the sectorization operation, if needed
        self.__reapply_associations()

    def __extract_sectorization_option_faces(self) -> List[Any]:
        """
        Method that extracts a list of face objects corresponding to the
        cell regions of the technological geometry from which the circles
        not centered in the cell are excluded.
        This list is sorted according to the distance of each face from
        the cell center.
        """
        # Remove entries corresponding to regions not centered in the cell
        supp_dict = {}
        # Loop through all the regions of the cell technological geometry
        for region in self.tech_geom_props:
            # Build a point inside the region to be used to retrieve the
            # corresponding 'Circle' object, if the region is based on a
            # circle, if not, none is found
            point = make_vertex_inside_face(region)
            found_circle = None
            # Loop through all the cell 'Circle' objects to get the one that
            # corresponds to the built region point, if the region comes from
            # a circle
            for circle in self.inner_circles:
                # Use the min distance to get the circle the point lays on
                if get_min_distance(point, circle.face) < 1e-5:
                    found_circle = circle
                    break
            # Filter out any found 'Circle' object that is not centered in the
            # cell
            if found_circle:
                if get_min_distance(make_cdg(found_circle.face),
                                    self.figure.o) < 1e-5:
                    # Given the region that corresponds to the found 'Circle'
                    # object build an entry associating the region face to its
                    # properties
                    supp_dict[region] = deepcopy(
                        self.tech_geom_props[region])
            else:
                # Add the entry associating the region face to its properties;
                # here, the region does not have a corresponding 'Circle'
                # object
                supp_dict[region] = deepcopy(self.tech_geom_props[region])
        # Sort the regions according to the distance from the cell center
        return sorted(supp_dict.keys(),
            key=lambda item: get_min_distance(
                self.figure.o, make_vertex_inside_face(item)))

    def __reapply_associations(self, to_skip: bool = False):
        """
        Method that re-builds the association of regions with properties and
        sectorization options (i.e. number of sectors and starting angle) and
        re-apply the sectorization operation, if it has already been applied
        beforehands.
        If specified, the second update operation can be skipped.

        Parameters
        ----------
        to_skip : bool
            If True, skip the sectorization options dictionary update
        """
        # Extract the regions of the cell technological geometry
        regions = self.extract_subfaces()
        # Associate each region of the cell technological geometry to its
        # properties
        self.__update_regions_association(regions, self.tech_geom_props, {})
        # If specified, skip the association of regions to sectorization
        # options
        # FIXME check the actual necessity of this parameter
        if to_skip:
            return
        self.__update_regions_association(
            regions, self.tech_geom_sect_opts, (1, 0))

        # Re-apply the sectorization, if it was present
        if self.sectorized_face:
            # Extract a list of face objects for the regions of the cell
            # technological geometry sorted according to the distance from
            # the cell center
            sorted_faces = self.__extract_sectorization_option_faces()
            supp_dict = {}
            for region in sorted_faces:
                if region in self.tech_geom_sect_opts:
                    supp_dict[region] = deepcopy(
                        self.tech_geom_sect_opts[region])
            # Update the dictionary associating to each cell region its
            # sectorization options
            self.tech_geom_sect_opts.clear()
            self.tech_geom_sect_opts.update(supp_dict)
            # Extract a list of tuples containing the numbers of sectors and
            # the sectorization starting angles
            sectors_angles = list(self.tech_geom_sect_opts.values())
            # Re-apply the sectorization with the updated options values
            self._sectorize_cell([value[0] for value in sectors_angles],
                                 [value[1] for value in sectors_angles],
                                 self.is_windmill_applied)
            # Partition the sectorization result with the externally added
            # edges, if any
            if self.added_edges:
                self.sectorized_face = make_partition([
                    self.sectorized_face], self.added_edges, ShapeType.FACE)

    def __update_regions_association(self,
                                     tech_regions: List[Any],
                                     regions_dict: Dict,
                                     default_value: Any) -> None:
        """
        Method that updates the given dictionary providing the association
        of the regions of the cell technological geometry with a specific
        option type (e.g. the properties, the sectorization options) that
        depends on the dictionary itself.
        If no option has been assigned yet to a region, the given default
        value is considered.

        Parameters
        ----------
        tech_regions  : List[Any]
            The list of face objects extracted from the cell technological
            geometry
        regions_dict  : Dict
            The dictionary associating the cell regions to the specific
            option type
        default_value : Any
            A default value to use for defining an entry in the given
            dictionary
        """
        # Declare a support dictionary
        support_dict = {}
        # Loop through all the current regions of the cell technological
        # geometry
        for region in tech_regions:
            # Initialize the entry
            support_dict[region] = default_value
            # Loop through all the regions and associated values
            for zone, values in regions_dict.items():
                # Check if a point inside the previous zone belongs to the
                # new region; if so, add a new entry using the values once
                # belonging to the zone
                if get_min_distance(make_vertex_inside_face(zone),
                                    region) < 1e-7:
                    support_dict[region] = values
                    break
        # Update the dictionary of regions of the technological geometry VS
        # associated values
        regions_dict.clear()
        regions_dict.update(support_dict)

    def __check_circle_presence(self, radius: float) -> Union[int, None]:
        """
        Method that loops through all the circles in the cell and checks
        if there is one whose radius value is the same of the given input
        positioned in the cell center.

        Parameters
        ----------
        radius  : float
                  The radius of the circle to look for in the list of
                  circles of the cell

        Returns
        -------
        An integer representing the index of the found circle in the list
        or None, if none could be found.
        """
        for indx, circle in enumerate(self.inner_circles):
            if (math.isclose(circle.radius, radius) and
                (get_point_coordinates(circle.o) ==
                 get_point_coordinates(self.figure.o))
                ):
                return indx
        return None

    def rotate(self, angle: float) -> None:
        """
        Method that rotates all the geometrical elements the cell is made of
        by a given angle in degrees.
        The cell GEOM face is rebuild from the rotated elements and so both
        the subfaces of the cell technological geometry and those coming from
        the cell sectorization are derived as well.

        Parameters
        ----------
        angle : float
                The angle of rotation in degree
        """
        # Get the cell center coordinates
        center = get_point_coordinates(self.figure.o)
        # Build the Z-axis of rotation positioned in the cell center
        z_axis = make_vector_from_points(
            self.figure.o, make_vertex((center[0], center[1], 1)))
        rotation = math.radians(angle)

        # Rotate the main figure (face, borders, vertices) of the cell
        self.figure.rotate(angle)
        # Rotate the cell inner circles
        for circle in self.inner_circles:
            circle.rotate(angle)
        # Rotate the whole face object representing the cell technological
        # geometry
        self.face = make_rotation(self.face, z_axis, rotation)

        # Re-build the dictionary of cell technological regions VS properties
        # with rotated zones
        new_dict = {}
        for f, p in self.tech_geom_props.items():
            # Add a new entry with the rotated region
            new_dict[make_rotation(f, z_axis, rotation)] = p
        # Copy the just built dictionary back into the old one
        self.tech_geom_props = deepcopy(new_dict)

        # Rotate the sectorized cell face, if any, without re-applying the
        # sectorization
        if self.sectorized_face:
            self.sectorized_face = make_rotation(
                self.sectorized_face, z_axis, rotation)

        # Update the cell rotation
        self.rotation += math.radians(angle)

    def translate(self, new_pos: Tuple[float, float, float]) -> Self:
        """
        Method that builds a copy of this class translated with respect to
        the current cell center position.
        All the geometrical elements describing the cell geometry, both its
        technological and its sectorized one, are translated accordingly.

        Parameters
        ----------
        new_pos : Tuple[float, float, float]
                  The new position where the cell has to be moved

        Returns
        -------
        A copy of this class instance with all the geometrical elements
        describing the cell geometry have been moved to the new position.
        """
        # Deep copy of this class instance
        cell = deepcopy(self)
        # Build a vector from the current center to the new one
        transl_vect = make_vector_from_points(cell.figure.o,
                                              make_vertex(new_pos))

        # Move the cell figure and all the characteristic geometrical elements
        cell.figure.translate(new_pos)
        # Move the inner circles
        for circle in cell.inner_circles:
            circle.translate(new_pos)

        # Move the cell face
        cell.face = make_translation(cell.face, transl_vect)

        # Move the cell sectorized face, if any
        if cell.sectorized_face:
            cell.sectorized_face = make_translation(cell.sectorized_face,
                                                    transl_vect)

        # Re-build the dictionary of the cell technological geometry zones
        # VS properties with translated regions
        new_dict = {}
        for tech_face, prop in cell.tech_geom_props.items():
            # Calculate the region new position so that it keeps its relative
            # position wrt the translated cell center
            region_to_center = update_relative_pos(
                make_cdg(tech_face), self.figure.o, new_pos)
            # Build a translation vector between the region previous CDG and
            # its new position
            vector = make_vector_from_points(
                make_cdg(tech_face), make_vertex(region_to_center))
            # Add a new entry with the translated region
            new_dict[make_translation(tech_face, vector)] = prop
        # Copy the just built dictionary back into the old one
        cell.tech_geom_props = deepcopy(new_dict)

        # Return the copied class instance whose position has been translated
        return cell

    def __build_regions(self) -> None:
        """
        Method that builds the 'Region' objects that corresponds to the
        regions coming from the cell technological geometry.
        These objects store information about the face object and associated
        properties.
        """
        # Re-initialize the list of 'Region' objects
        self.regions.clear()
        # For each cell region build a new 'Region' object and add to the
        # corresponding list
        for i, (tech_face, props) in enumerate(self.tech_geom_props.items()):
            # Add a new region
            self.regions.append(Region(
                face=tech_face,
                inner_point=make_vertex_inside_face(tech_face),
                name=f"Region {i}",
                properties=deepcopy(props)
            ))

    def __build_lines_intersection(
            self, lines1: List[Any], lines2: List[Any]) -> List[Any]:
        """
        Method that builds a list of points resulting from the intersection
        of two lines.
        For each element of the two given list of GEOM line objects, the
        corresponding intersection point, as a GEOM vertex object is built,
        if any can be found.

        Parameters
        ----------
        lines1  : List[Any]
            The first list of line objects
        lines2  : List[Any]
            The second list of line objects

        Returns
        -------
        A list of GEOM vertex objects representing the intersection points
        found between all the two lists of lines.
        """
        # Initialize the list of intersection points
        intersections = []
        # Loop over all the lines of the first list
        for l1 in lines1:
            # Loop over all the lines of the second list
            for l2 in lines2:
                try:
                    # Build the intersection point (as GEOM object) over
                    # the current lines and add it to the list
                    intersections.append(
                        make_vertex_on_lines_intersection(l1, l2))
                    break
                except RuntimeError:
                    # In case of no intersection, an exception is raised,
                    # but skipped as not considered a problem
                    pass
        # Return the built list of intersection points
        return intersections

    @abstractmethod
    def sectorize(
        self, sectors_no: List[int], angles: List[float], **kwargs) -> None:
        """
        Abstract method that subdivides the cell into sectors. Given the
        number of sectors for each cell zone and the values of the angles
        to start the sectorization from, points are built on the geometric
        surfaces of the cell.
        Lines are drawn between those points to define the sectors and the
        corresponding subfaces are extracted.

        Parameters
        ----------
        sectors_no  : List[int]
                      List of integers representing the number of subdivision
                      for each cell zone coming from the technological geometry
        angles      : List[float]
                      List of angles (in degree) the sectorization should start
                      from for each cell zone coming from the technological
                      geometry
        kwargs      : Any
                      Additional parameters specific to the type of cell to be
                      sectorized
        """

    def _sectorize_cell(self,
                        sectors_no: List[int],
                        angles: List[float],
                        windmill: bool = False) -> None:
        """
        Method that performs the sectorization of the cell zones coming from
        the technological geometry.
        Given the number of sectors for each cell zone and the values of the
        angles to start the sectorization from, points are built on the
        geometric surfaces of the cell.
        Lines are drawn between those points to define the sectors and the
        corresponding subfaces are extracted.

        Parameters
        ----------
        sectors_no  : List[int]
                      List of integers representing the number of subdivision
                      for each cell zone coming from the technological geometry
        angles      : List[float]
                      List of angles (in degree) the sectorization should start
                      from for each cell zone coming from the technological
                      geometry
        windmill    : bool = False
                      Flag indicating if a windmill sectorization of the cell
                      needs to be considered
        """
        # Check the correctness of the lists storing the information for
        # performing the sectorization
        self.__check_sectorization_elements_len(
            len([c.borders[0] for c in self.inner_circles
                   if get_min_distance(c.o, self.figure.o) < 1e-5]) + 1,
            len(sectors_no),
            len(angles))
        # Check the correctness of the number of sectors, as specific values
        # are admitted only according to the cell type
        self.__check_sectors_no(sectors_no)
        # Check the correctness of the angle the sectorization starts from
        # for each cell zone
        self.__check_sectorization_angle(angles, sectors_no)
        # Re-initialize the cell face object and store the windmill option
        self.sectorized_face = None
        self.is_windmill_applied = windmill

        # Handle the sectorization operation: the resulting list holds all
        # the cell edges identifying:
        # . the cell borders;
        # . the sectors between consecutive circles;
        # . the circles themselves;
        # . the lines connecting consecutive intersection points (if the
        #   'windmill' option is active)
        edges = self.__build_sectorization_edges(sectors_no, angles, windmill)

        # Build the cell face over the edges defining the cell sectorization
        self.sectorized_face = make_partition([self.face],
                                              edges,
                                              ShapeType.FACE)
        # Update the sectorization options for each zone of the cell
        # technological geometry
        self.tech_geom_sect_opts.clear()
        for region, opts in zip(self.__extract_sectorization_option_faces(),
                                zip(sectors_no, angles)):
            self.tech_geom_sect_opts[region] = opts

    def __build_sectorization_edges(self,
                                    sectors_no: List[int],
                                    angles: List[float],
                                    windmill: bool) -> List[Any]:
        """
        Method that performs the sectorization operation on the cell and
        builds the edge objects resulting from the application of this
        operation. This list stores:
        - the cell borders;
        - the sectors between consecutive circles;
        - the circles themselves;
        - the lines connecting consecutive intersection points (if the
          'windmill' option is active)

        Parameters
        ----------
        sectors_no  : List[int]
            List of integers representing the number of subdivision for
            each cell zone coming from the technological geometry
        angles      : List[float]
            List of angles (in degree) the sectorization should start
            from for each cell zone coming from the technological geometry
        windmill    : bool = False
            Flag indicating if a windmill sectorization of the cell needs
            to be considered

        Returns
        -------
        A list of edge objects resulting from the sectorization operation.
        """
        # -------------------------------------------------------------
        # Initialization of support lists (containing the sectorization
        # construction elements)
        # -------------------------------------------------------------
        # Extract a list of circles (as GEOM edge objects) from the stored
        # 'Circle' objects representing the circles centered in the cell
        # geometry
        circles = [c.borders[0] for c in self.inner_circles
                   if get_min_distance(c.o, self.figure.o) < 1e-5]
        # Declare the list storing all the cell sectorization edges
        edges = self.figure.borders + circles

        # ---------------------------------
        # Construction of sectors and zones
        # ---------------------------------
        # Determine the parameter identifying where a point on the circles
        # will be added. It based on the rotation of the cell and ranges
        # from 0 to 1. In case the rotation is negative, the opposite value
        # is taken into account.
        u_rotation = self.rotation / (2 * math.pi)
        rotation = u_rotation if self.rotation >= 0 else 1 + u_rotation
        zones_no = len(circles) + 1
        # Loop over all the sector numbers and angles
        for zone_no in range(zones_no):
            # Extract the number of subdivision sectors and the angle to
            # start the subdivision from
            sector = sectors_no[zone_no]
            angle = angles[zone_no]

            if sector == 1:
                # No subdivision has to be performed on the current cell zone,
                # hence append an empty list to the one storing all the cell
                # edges resulting from the sectorization
                edges += []
                # Continue with the next zone
                continue
            # --------------------------------------------------------
            # Handle the sectorization if the number of subdivision is
            # greater than 1
            # --------------------------------------------------------
            # Build subdivision points on the 'previous' circle
            if zone_no > 0:
                # Build a list of GEOM points on the 'previous' circle
                # curve (wrt the zone number), each on a different
                # position given by the sector index and the angle
                sbdv_pnts_pre = self.__build_sectorization_points(
                    circles[zone_no - 1], angle, sector, rotation)
            else:
                # Build a list of GEOM points for the cell center
                sbdv_pnts_pre = [self.figure.o] * sector

            # Build subdivision points on the 'current' circle curve
            # (wrt the zone number)
            if zone_no < zones_no - 1:
                sbdv_pnts_cur = self.__build_sectorization_points(
                    circles[zone_no], angle, sector, rotation)
            else:
                # Build subdivision points on the construction circle the
                # cell is inscribed into
                sbdv_pnts_out = self.__build_sectorization_points(
                    self.figure.out_circle, angle, sector, 0.0)
                # Build a list of GEOM lines between the just built points
                # on the construction circle and the origin
                lines_o_out = [
                    make_line(self.figure.o, v0) for v0 in sbdv_pnts_out]
                # Build a list of points (as a GEOM object each) on the
                # intersections between the line connecting the
                # construction circle subdivision points with the cell
                # center and the cell borders
                sbdv_pnts_cur = self.__build_lines_intersection(
                    lines_o_out, self.figure.borders)

                # -------------
                # Windmill case
                # -------------
                # Build the wings of the sectorization, if required, and
                # add to the list holding all the sectorization edges
                if windmill:
                    edges += self.__handle_windmill_sectorization(
                        sector, angle, sbdv_pnts_cur)

            # ------------------------
            # Creation of sector edges
            # ------------------------
            if len(sbdv_pnts_pre) != len(sbdv_pnts_cur):
                raise RuntimeError(
                    f"ERROR while building sectors for {zone_no} sector: "
                    f"sector: no. sectors {sector}, starting angle: "
                    f"{angle}.")
            # Given the subdivision points of the 'previous' circle and
            # the ones of the 'current' circle (or the intersection
            # points), lines are built between corresponding points and
            # added to the list of edges resulting from the sectorization
            edges += [make_line(v1, v2)
                for v1, v2 in zip(sbdv_pnts_pre, sbdv_pnts_cur)]

        # Return the built list of edges
        return edges

    def __handle_windmill_sectorization(self,
                                        sector: int,
                                        angle: float,
                                        sbdv_pnts_cur: List[Any]) -> List[Any]:
        """
        Method that builds the wings of the windmill sectorization, if needed.
        The method connects two successive intersection points (on adjacent
        cell borders), if the number of sectors for the outmost region of the
        cartesian cell technological geometry is 8 or 16.

        Parameters
        ----------
        sector        : int
            The number of sector for the current region of the cartesian cell
            technological geometry
        angle         : float
            The starting angle of the sectorization
        sbdv_pnts_cur : List[Any]
            The list of vertex objects to connect

        Returns
        -------
        The list of line objects built on the subdivision point, if required,
        otherwise an empty list.
        """
        # Build a list of line object giving the windmill wings
        if sector == 8 and angle == 22.5:
            return [make_line(sbdv_pnts_cur[i], sbdv_pnts_cur[i+1])
                for i in (0, 2, 4, 6)]
        if sector == 16 and angle == 0:
            return [make_line(sbdv_pnts_cur[i], sbdv_pnts_cur[i+2])
                for i in (1, 5, 9, 13)]
        # If here, return an empty list
        return []

    def __build_sectorization_points(self,
                                     circle: Any,
                                     starting_angle: float,
                                     no_sectors: int,
                                     rotation: float) -> List[Any]:
        """
        Method that builds a list of vertex objects on the given circle
        object, each on a different position given by the sector index,
        the starting angle and the rotation parameter.

        Parameters
        ----------
        circle          : Any
            The circle object on which points have to be built
        starting_angle  : float
            The angle (in degrees) from which to build points
        no_sectors      : int
            The number of subdivision pointds to build
        rotation        : float
            A parameter in the [0-1] range identifying where to built
            a point. It based on the rotation of the cell

        Returns
        -------
        A list of vertex objects subdividing the given circle.
        """
        return [
            make_vertex_on_curve(
                circle, starting_angle/360.0 + float(i)/no_sectors + rotation)
            for i in range(no_sectors)
        ]

    def __check_sectorization_elements_len(self,
                                           no_regions: int,
                                           sectors_no_len: int,
                                           angles_len: int) -> None:
        """
        Method for checking the correctness of the lenght of the lists
        for the sectorization operation. It is verified that:
        - the lists length is the same;
        - their size coincides with the number of cell zones coming from
          the cell technological geometry.

        Parameters
        ----------
        no_regions      : int
                          The number of regions of the cell technological
                          geometry
        sectors_no_len  : int
                          The length of the list providing the number of
                          sectors each cell zone has to be subdivided into
        angles_len      : int
                          The length of the list providing the angle the
                          sectorization starts from
        """
        # Check that the number of sectors and angles is the same
        if not sectors_no_len == angles_len:
            raise ValueError(f"The numbers of sectors ({sectors_no_len}) and "
                             f"angles to start the sectorization from ("
                             f"{angles_len}) do not coincides.")
        if not sectors_no_len == no_regions:
            raise ValueError("The number of elements in the lists for the "
                             f"sectorization ({sectors_no_len}) does not "
                             f"coincide with the number of cell zones ("
                             f"{no_regions}).")

    def __check_sectors_no(self, sectors_no: List[int]) -> None:
        """
        Method that checks if the number of sectors each cell zone has
        to be subdivided is valid according to the specific cell type.
        An exception is raised if an invalid number is found.

        Parameters
        ----------
        sectors_no  : List[int]
                      List of integers identifying the number of sectors
                      each cell zone has to be subdivided
        """
        for i, sec_no in enumerate(sectors_no):
            if sec_no not in self.VALID_SECTOR_NO_VS_ANGLE:
                raise ValueError(f"The zone #{i} has an invalid "
                    "sectorization: the accepted values are "
                    f"{self.VALID_SECTOR_NO_VS_ANGLE.keys()}")

    def __check_sectorization_angle(self,
                                    angles: List[float],
                                    sectors_no: List[int]) -> None:
        """
        Method that checks if the angle the sectorization starts from for
        each cell zone is valid according to the specific cell type and
        to the given number of sectors.

        An exception is raised if an invalid value is found.

        Parameters
        ----------
        angles      : List[float]
                      List of angles the sectorization starts from for each
                      cell zone
        sectors_no  : List[int]
                      List of integers identifying the number of sectors
                      each cell zone has to be subdivided
        """
        for i, (angle, sec_no) in enumerate(zip(angles, sectors_no)):
            if not angle in self.VALID_SECTOR_NO_VS_ANGLE[sec_no]:
                raise ValueError(f"The zone #{i} has an invalid sectorization "
                                 "starting angle: the accepted value for the "
                                 f"given {sec_no} number of sectors are "
                                 f"{self.VALID_SECTOR_NO_VS_ANGLE[sec_no]}")

    def set_properties(
            self, properties: Dict[PropertyType, List[str]]) -> None:
        """
        Method that allows to set the properties for all the zones the cell
        is made of (technological geometry).
        The accepted convention considers that the elements in the list for
        each property represent the values associated to the regions starting
        from the inner circle to the last element (space between last circle
        and cell borders).

        Parameters
        ----------
        properties  : Dict[PropertyType, List[str]]
                      Dictionary collecting the properties for each cell
                      region; different types, provided by the 'PropertyType'
                      enumeration, can be provided
        """
        # Extract the subfaces identifying the cell technological geometry
        tech_geom_faces = self.extract_subfaces()
        # Check if the number of cell regions coincides with the number of
        # property elements for all the given property types
        no_regions = len(tech_geom_faces)
        for prop_type, values in properties.items():
            if not no_regions == len(values):
                raise ValueError(
                    "There is no correspondence between the number "
                    f"of cell zones ({no_regions}) and that of the "
                    f"properties with type '{prop_type.name}' (no. "
                    f"{len(values)})")

        # Clear any previously set entry in the dictionary associating the
        # regions to the properties, as this method sets the properties for
        # all.
        self.tech_geom_props.clear()

        # Re-build the dictionary of cell zones VS properties
        for i, tech_face in enumerate(tech_geom_faces):
            # Collect all properties for a region
            prop_types = {
                type: values[i] for type, values in properties.items()}
            # Associate a region with its properties
            self.tech_geom_props[tech_face] = prop_types

    def __build_sector_regions(self) -> None:
        """
        Method for building the 'Region' objects that corresponds to the
        cell regions after applying a sectorization operation.
        By comparing the distance between a point on each sector face and
        the regions coming from the techonological geometry, the associated
        properties are assigned to the newly built 'Region' instance.
        """
        if not self.sectorized_face:
            raise RuntimeError("No sectorization has been applied: regions "
                               "cannot be extracted.")
        # Clear the previously stored regions of the cell geometry
        self.regions.clear()
        # Extract the subfaces of the sectorized cell
        subfaces = self.__extract_subfaces_from_face(self.sectorized_face)
        indx = 0

        # Loop through all the subfaces resulting from the cell sectorization
        for subface in subfaces:
            # Get a vertex inside the subface
            subface_point = make_vertex_inside_face(subface)
            # Loop through all the regions of the cell technological geometry
            for region, props in self.tech_geom_props.items():
                # Check if the lattice subface overlaps the cell one
                if get_min_distance(region, subface_point) < 1e-7:
                    # Update the index
                    indx += 1
                    # Build a region for the lattice
                    self.regions.append(Region(
                        face=subface,
                        inner_point=subface_point,
                        name=f"Region {indx}",
                        properties=deepcopy(props)
                    ))
                    break

        # Check if the number of subfaces coincides with the one of regions
        if len(subfaces) != len(self.regions):
            raise AssertionError("Mismatch between the number of cells "
                                 f"subfaces ({len(subfaces)}) and that of "
                                 f"the found regions ({len(self.regions)})")

    def show(self,
             property_type_to_show: Union[PropertyType, None] = None,
             geometry_type_to_show: GeometryType = GeometryType.TECHNOLOGICAL,
             update_view: bool = True) -> None:
        """
        Method that shows the cell and its regions to the current SALOME
        study.
        Given the property type, the viewe shows the region with a color
        that has been associated to the property type of the region. If
        None is provided, the regions are shown without any color.

        Parameter
        ---------
        property_type_to_show : Union[PropertyType, None] = None
            The type of property to show for each cell region; the viewer
            shows them in terms of a spefic color that each region associates
            to the indicated type
        geometry_type_to_show : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry to show: regions are shown accordingly
            with the selected type
        update_view           : bool
            Flag stating if the viewer should be updated by displaying the
            cell regions
        """
        # Erase all objects from the current view
        clear_view()

        if self.name == "":
            self.name = get_shape_name(self.face)
        # Delete the face from the study if already present
        if self.face_entry_id:
            remove_from_study(self.face_entry_id)
        # Add the cell to the current SALOME study
        self.face_entry_id = add_to_study(self.face, self.name)

        try:
            # Show in the viewer the cell regions according to the indicated
            # geometry type
            match geometry_type_to_show:
                case GeometryType.TECHNOLOGICAL:
                    self.__build_regions()
                case GeometryType.SECTORIZED:
                    self.__build_sector_regions()
                case _:
                    raise ValueError(f"{geometry_type_to_show}: unhandled "
                                     "type of geometry to show.")
            # Assign the same color to all the regions having the same property
            # value, if any has been specified to show
            self.__associate_colors_to_regions(property_type_to_show)
            # Add the cell regions to the study, each with an assigned color
            for region in self.regions:
                # Set the region color in the viewer
                set_color_face(region.face, region.color)
                # Add the cell region to the study
                region_id = add_to_study_in_father(
                    self.face, region.face, region.name)
                # Display the region in the current view, if needed
                if update_view:
                    display_shape(region_id)
        except:
            # Show only the whole cell face before raising the caught exception
            display_shape(self.face_entry_id)
            raise
        finally:
            # Show everything on the SALOME application
            update_salome_study()
            # Update the geometry type
            self.displayed_geom = geometry_type_to_show

    def __associate_colors_to_regions(
            self, property_type: Union[PropertyType, None]) -> None:
        """
        Method that assign the same color to all the regions having the same
        value for the given property type.

        Parameters
        ----------
        property_type : Union[PropertyType, None]
            The type of property for which colors must be assigned to regions
            having the same property type value. If None, the regions color
            is reset to its default value
        """
        # If no colorset to display, reset the region colors
        if not property_type:
            for region in self.regions:
                # Set the region color to its default value
                region.reset_color_to_default()
            return

        # Extract the unique values of the given property type associated to
        # each cell region
        values = set()
        missing_regions_points = []
        for region, prop_vals in self.tech_geom_props.items():
            if property_type not in prop_vals:
                missing_regions_points.append(
                    get_point_coordinates(make_vertex_inside_face(region)))
                continue
            values.add(prop_vals[property_type])
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
        # Generate a specific amount of colors as the number of different
        # values for the same given property type
        colors = generate_unique_random_colors(len(values))
        # Build a dictionary of values for the given property type VS color
        property_vs_color = dict(zip(list(values), colors))
        # Loop through all the regions and assign a color corresponding to
        # the value of the given property type
        for region in self.regions:
            # Get the value of the given property type associated to the
            # region, if any
            if property_type not in region.properties:
                raise RuntimeError(
                    f"No {property_type.name} property type has been "
                    f"associated to the {region.name}.")
            value = region.properties[property_type]
            # Set the region color
            region.set_property_color(property_vs_color[value])

    def get_regions_info(self) -> None:
        """
        Method for retrieving descriptive information about any of the cell
        regions whose corresponding object has been selected in the SALOME
        study.
        """
        # Extract the geometrical object currently selected in the current
        # SALOME study
        shape = get_selected_object()
        if not shape:
            raise RuntimeError("Please, select a single region whose data to "
                               "show.")
        # Get the region that corresponds to the given shape
        for region in self.regions:
            if get_min_distance(region.inner_point, shape) < 1e-5:
                # Print info about the region name and its properties
                print(f"{region.name}:")
                if not region.properties:
                    print("   No associated properties.")
                    return
                for prop_type, value in region.properties.items():
                    print(f"   {prop_type.name}: {value}")

    def set_region_property(
            self, property_type: PropertyType, value: str) -> None:
        """
        Method that allows to set a value of a given type of property for
        the cell region which is currently selected in the SALOME study.

        Parameters
        ----------
        property_type : PropertyType
            The value of the 'PropertyType' enumeration indicating which
            type of property to assign to the selected region
        value         : str
            The value of the property type to assign to the selected region
        """
        # Check which cell geometry type is currently shown; if different
        # from the TECHNOLOGICAL one, raise an exception
        if self.displayed_geom != GeometryType.TECHNOLOGICAL:
            raise RuntimeError(
                f"Currently showing the {self.displayed_geom.name} type of "
                "geometry. To set cell regions properties, show the '"
                "TECHNOLOGICAL' geometry first.")
        # Extract the geometrical object currently selected in the current
        # SALOME study
        shape = get_selected_object()
        if not shape:
            raise RuntimeError("Please, select a single region to assign "
                               "a property to.")
        # Get the region that corresponds to the given shape
        for region in self.tech_geom_props:
            if get_min_distance(make_vertex_inside_face(region), shape) < 1e-5:
                # Either update or create an entry to the region dictionary
                # of properties
                self.tech_geom_props[region] = {property_type: value}

    def update_geometry(self, geo_type: GeometryType) -> None:
        """
        Method for updating the cell geometry representation according to
        the provided type.
        The cell corresponding face object is updated with the one currently
        selected in the SALOME study.

        Parameters
        ----------
        geo_type  : GeometryType
            The cell type of geometry to update, as value of the
            'GeometryType' enumeration
        """
        # Check if the displayed geometry coincides with the one to modify
        if geo_type != self.displayed_geom:
            raise RuntimeError(
                f"Currently displaying the cell {self.displayed_geom.name} "
                "geometry: this does not match with the indicated "
                f"{geo_type.name} geometry type to modify.")
        # Extract the geometrical object currently selected in the current
        # SALOME study
        shape = get_selected_object()
        if not shape:
            raise RuntimeError("Please, select a shape to update the cell "
                               "geometry with.")

        # Handle the different geometry type differently
        match geo_type:
            case GeometryType.SECTORIZED:
                # Extract the additional edges by comparing with the ones
                # already present in the sectorized face
                self.added_edges = []
                for e1 in extract_sub_shapes(shape, ShapeType.EDGE):
                    found = False
                    for e2 in extract_sub_shapes(self.sectorized_face,
                                                   ShapeType.EDGE):
                        # Ignore the edge if it is already present in the
                        # cell sectorized face
                        if tuple(get_kind_of_shape(e1)[1:]) == tuple(
                            get_kind_of_shape(e2)[1:]):
                            found = True
                            break
                    if not found:
                        self.added_edges.append(e1)
                # Update the face object result of the cell sectorization
                self.sectorized_face = shape
            case GeometryType.TECHNOLOGICAL:
                # Reset the cell face object
                self.face = self.figure.face
                # Clear the previously stored 'Circle' objects
                self.inner_circles.clear()
                # Loop through all the edge objects of the cell to identify
                # the presence of internal circles
                for edge in extract_sub_shapes(shape, ShapeType.EDGE):
                    # Get the information of the given edge object
                    data = get_kind_of_shape(edge)
                    if str(data[0]) == 'CIRCLE':
                        # Add each circle to the cell according to its position
                        self.add_circle(
                            data[7], get_point_coordinates(make_cdg(edge)))
            case _:
                raise ValueError(
                    f"{geo_type}: unhandled type of geometry to update.")

        # Update the dictionary storing the regions properties, given the
        # 'Region' objects; since the cell geometry is shown in the 3D
        # viewer, this list has been built and no check on its presence is
        # necessary.
        for region in self.regions:
            for zone in self.tech_geom_props:
                if get_min_distance(make_vertex_inside_face(zone),
                                    region.face) < 1e-7:
                    self.tech_geom_props[zone] = region.properties

class RectCell(GenericCell):
    """
    Class providing the means for building a cartesian (i.e. rectangular)
    cell.
    This cell can be characterised by a specified number of concentric
    circles, each one associated with a different property, to represent
    the different zones of a fuel pin cell.
    This technological geometry can be sectorized, i.e. subdivision points
    on the inner circles can be defined in order to partition the face
    into subfaces defining the cell sectors.
    Properties of the technological geometry zones are associated to each
    of the corresponding sectors.

    Parameters
    ----------
    center          : Union[Tuple[float, float, float], None] = None
                      The X-Y-Z coordinates of the cell center
    height_x_width  : Tuple[float, float] = (1, 1)
                      A tuple providing the cell height and width
    rounded_corners : Union[List[Tuple[int, float]], None]
                      A list of tuples providing the corner index and the
                      curvature radius of the same corner
    name            : str = "Cell"
                      The cell name in the current SALOME study

    Attributes
    ----------
    cell_type       : Union[CellType, None]
        The type of cell expressed as a value of the 'CellType' enumeration
    face            : Any
        The face object providing the cell technological geometry
        representation
    face_entry_id   : Union[str, None]
        The ID of the cell face object used in the current study
    figure          : GenericSurface
        The figure representing the cell shape
    height              : float
        The cell height
    inner_circles   : List[Circle]
        The list of 'Circle' objects representing the cell concentric circles
    is_windmill_applied : bool
        Boolean flag stating if the windmill sectorization is applied
    name            : str
        The name of the cell to be used in the current study
    regions         : List[Region]
        The list of 'Region' objects storing information about either the
        regions of the cell technological geometry or the ones resulting from
        its sectorization
    rotation        : float
        The rotation angle of the cell expressed in radians
    sectorized_face : Any
        The face object providing the cell sectorized geometry representation
        (used for calculations)
    tech_geom_props : Dict[Any, Dict[PropertyType, str]]
        A dictionary between the regions of the cell technological geometry
        VS the dictionary of property types VS value associated to each zone
    tech_geom_sect_opts : Dict[Any, Tuple[int, float]]
        A dictionary between the regions of the cell technological geometry
        VS the tuple identifying the associated sectorization options
    width               : float
        The cell width

    VALID_SECTOR_NO_VS_ANGLE  : Dict[int, List[float]]
        Dictionary providing the number of valid sectors each cell zone
        could be subdivided VS the accepted angles the sectorization can
        start from.
    """
    # Number of valid sectors each cell zone could be subdivided VS the
    # accepted angles the sectorization can start from
    VALID_SECTOR_NO_VS_ANGLE: Dict[int, List[float]] = {1 : [0],
                                                        4 : [0, 45],
                                                        8 : [0, 22.5],
                                                        16: [0]}

    def __init__(self,
                 center: Union[Tuple[float, float, float], None] = None,
                 height_x_width: Tuple[float, float] = (1, 1),
                 rounded_corners: Union[List[Tuple[int, float]], None] = None,
                 name: str = "Cell") -> None:
        # Initialize the superclass passing the rectangular shape
        shape = Rectangle(center=center,
                          height=height_x_width[0],
                          width=height_x_width[1],
                          rounded_corners=rounded_corners)
        super().__init__(shape, name)
        # Store the cell dimensions
        self.height: float = height_x_width[0]
        self.width: float = height_x_width[1]
        # Set the cell type
        self.cell_type: Union[CellType, None] = CellType.RECT

    def _initialize_specific_cell(self) -> None:
        """
        Method for initializing instance attributes that are specific for
        the cartesian cell.
        """
        self.windmills: List[Any] = []

    def _check_radius_vs_cell_dim(self, radius: float) -> None:
        """
        Method for assessing if the circle, whose radius is given as input,
        can be added to the cell. The radius is checked against the width
        and height of the cell.
        In case the radius is greater than the cell characteristic dimensions,
        an exception is raised.

        Parameters
        ----------
        radius  : float
                  The radius of the circle to check against the cell
                  dimensions
        """
        if 2*radius > self.width or 2*radius > self.height:
            raise RuntimeError("The circle cannot be added to the cell as its "
                               "dimensions exceedes the ones of the cell.")

    def sectorize(
        self, sectors_no: List[int], angles: List[float], **kwargs) -> None:
        """
        Method that subdivides the cell into sectors. Given the number of
        sectors for each cell zone and the values of the angles to start
        the sectorization from, points are built on the geometric surfaces
        of the cell.
        Lines are drawn between those points to define the sectors and the
        corresponding subfaces are extracted.

        Parameters
        ----------
        sectors_no  : List[int]
                      List of integers representing the number of subdivision
                      for each cell zone coming from the technological geometry
        angles      : List[float]
                      List of angles (in degree) the sectorization should start
                      from for each cell zone coming from the technological
                      geometry
        kwargs      : Any
                      Additional parameters specific to the cartesian cell type
                      to be sectorized (e.g., windmill)
        """
        # Check if the 'windmill' flag, stating that a windmill sectorization
        # of the cell needs to be considered, has been passed. If not, 'False'
        # is passed to the method that performs the cell sectorization
        if 'windmill' not in kwargs:
            self._sectorize_cell(sectors_no, angles, False)
        else:
            self._sectorize_cell(sectors_no, angles, kwargs['windmill'])
        # Clear the externally added edges as calling this method resets any
        # previously applied sectorization
        self.added_edges.clear()


class HexCell(GenericCell):
    """
    Class providing the means for building an hexagonal cell.
    This cell can be characterised by a specified number of concentric
    circles, each one associated with a different property, to represent
    the different zones of a fuel pin cell.
    This technological geometry can be sectorized, i.e. subdivision points
    on the inner circles can be defined in order to partition the face
    into subfaces defining the cell sectors.
    Properties of the technological geometry zones are associated to each
    of the corresponding sectors.

    Parameters
    ----------
    center      : Union[Tuple[float, float, float], None] = None
                  The X-Y-Z coordinates of the cell center
    edge_length : float = 1.0
                  The length of the hexagonal cell edge
    name        : str = "Hexagon"
                  The cell name in the current SALOME study

    Attributes
    ----------
    apothem         : float
        The length of the hexagonal cell apothem
    cell_type       : Union[CellType, None]
        The type of cell expressed as a value of the 'CellType' enumeration
    edge_length     : float
        The length of the hexagonal cell edge
    face            : Any
        The face object providing the cell technological geometry
        representation
    face_entry_id   : Union[str, None]
        The ID of the cell face object used in the current study
    figure          : GenericSurface
        The figure representing the cell shape
    height              : float
        The cell height
    inner_circles   : List[Circle]
        The list of 'Circle' objects representing the cell concentric circles
    is_windmill_applied : bool
        Boolean flag stating if the windmill sectorization is applied
    name            : str
        The name of the cell to be used in the current study
    regions         : List[Region]
        The list of 'Region' objects storing information about either the
        regions of the cell technological geometry or the ones resulting from
        its sectorization
    rotation        : float
        The rotation angle of the cell expressed in radians
    sectorized_face : Any
        The face object providing the cell sectorized geometry representation
        (used for calculations)
    tech_geom_props : Dict[Any, Dict[PropertyType, str]]
        A dictionary between the regions of the cell technological geometry
        VS the dictionary of property types VS value associated to each zone
    tech_geom_sect_opts : Dict[Any, Tuple[int, float]]
        A dictionary between the regions of the cell technological geometry
        VS the tuple identifying the associated sectorization options
    width               : float
        The cell width

    VALID_SECTOR_NO_VS_ANGLE  : Dict[int, List[float]]
        Dictionary providing the number of valid sectors each cell zone
        could be subdivided VS the accepted angles the sectorization can
        start from.
    """
    # Number of valid sectors each cell zone could be subdivided VS the
    # accepted angles the sectorization can start from
    VALID_SECTOR_NO_VS_ANGLE: Dict[int, List[float]] = {1 : [0],
                                                        6 : [0, 30]}

    def __init__(self,
                 center: Union[Tuple[float, float, float], None] = None,
                 edge_length: float = 1.0,
                 name: str = "Hexagon") -> None:
        # Initialize the superclass passing the hexagon shape
        shape = Hexagon(center=center, edge_length=edge_length)
        super().__init__(shape, name=name)
        # Store the cell dimensions
        self.edge_length: float = edge_length
        self.apothem: float = shape.apothem
        self.cell_type: Union[CellType, None] = CellType.HEX

    def _initialize_specific_cell(self) -> None:
        """
        Method for initializing instance attributes that are specific for
        the hexagonal cell.
        Up to now, nothing needs to be performed here.
        """

    def _check_radius_vs_cell_dim(self, radius: float) -> None:
        """
        Method for assessing if the circle, whose radius is given as input,
        can be added to the cell. The radius is checked against the apothem
        of the hexagonal cell.
        In case the radius is greater than the cell characteristic dimension,
        an exception is raised.

        Parameters
        ----------
        radius  : float
                  The radius of the circle to check against the cell
                  dimensions
        """
        if radius > self.apothem:
            raise ValueError("The circle cannot be added to the cell as its "
                             "dimensions exceedes the one of the cell")

    def sectorize(
        self, sectors_no: List[int], angles: List[float], **kwargs) -> None:
        """
        Method that subdivides the cell into sectors. Given the number of
        sectors for each cell zone and the values of the angles to start
        the sectorization from, points are built on the geometric surfaces
        of the cell.
        Lines are drawn between those points to define the sectors and the
        corresponding subfaces are extracted.

        Parameters
        ----------
        sectors_no  : List[int]
                      List of integers representing the number of subdivision
                      for each cell zone coming from the technological geometry
        angles      : List[float]
                      List of angles (in degree) the sectorization should start
                      from for each cell zone coming from the technological
                      geometry
        kwargs      : Any
                      Additional parameters specific to the hexagona cell type
                      to be sectorized. Up to now, nothing is passed
        """
        # Since the hexagonal cell type does not come with a 'windmill'
        # sectorization, 'False' is passed to the method that performs
        # the cell sectorization
        self._sectorize_cell(sectors_no, angles, False)
        # Clear the externally added edges as calling this method resets any
        # previously applied sectorization
        self.added_edges.clear()


class SalomeCell():
    """
    Class providing the means for interacting with a generic cell built from
    within SALOME using its GEOM module.
    Given the cell face ID in the current SALOME study, the face is retrieved
    and all the needed informations are extracted.
    This cell can be characterised by a specified number of concentric
    circles, each one associated with a different property, to represent
    the different zones of a fuel pin cell.
    This technological geometry can be sectorized, i.e. subdivision points
    on the inner circles can be defined in order to partition the face
    into subfaces defining the cell sectors.
    Properties of the technological geometry zones are associated to each
    of the corresponding sectors.

    Parameters
    ----------
    face_id   : str
                The cell face ID in the current SALOME study
    cell_type : CellType
                Value of the 'CellType' enumeration identifying the type
                of cell
    Attributes
    ----------
    ...
    """
    def __init__(self, face_id: str, cell_type: CellType) -> None:
        # Extract the face the given ID corresponds to in the SALOME study
        self.face : Any = get_object_from_id(face_id)
        print(f"The face has name '{get_shape_name(self.face)}'")
        # TODO continue implementation
        # super().__init__(None, self.face.GetName())


if __name__ == "__main__":
    # -----------------------------------
    # Testing the module functionalities.
    # -----------------------------------
    # Build a cartesian cell
    cell = RectCell(name="Cartesian cell")
    # Add three inner circles to the first cell
    radii = [0.65/2, 0.3, 0.75/2]
    for r in radii:
        cell.add_circle(r)
    # Add both cell base types to the SALOME study
    cell.show()

    # Assign the materials to each zone in the cell
    cell.set_properties({
        PropertyType.MATERIAL: ["TSTR.'TMIL_MOC'.'COMB0201'.1",
                                "TSTR.'TMIL_MOC'.'COMB0201'.2",
                                "TSTR.'TMIL_MOC'.'COMB0201'.3",
                                "TSTR.'TMIL_MOC'.'MODEXTCRAYON'"]
    })

    # Build the cell sectorization and show it in the 3D viewer
    cell.sectorize([4, 1, 8, 16], [0, 0, 22.5, 0], windmill=True)
    cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
    # Print the properties assigned to each cell sector
    print("No. regions after sectorization:", len(cell.regions))
    for region in cell.regions:
        print(f"Sector face {region.face_entry_id}, {region.name}, "
              f"{region.properties}")

    # Show everything on the SALOME application
    update_salome_study()
