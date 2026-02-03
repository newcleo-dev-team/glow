"""
Module containing classes providing the means for creating cells from basic
geometries in SALOME.
"""
import math

from typing import Any, Dict, List, Tuple

from glow.geometry_layouts.fillable_layouts import Fillable
from glow.geometry_layouts.geometries import Circle, Hexagon, Surface, \
    Rectangle, build_parallelogram, build_regular_triangle, \
    build_right_triangle, build_right_triangle_from_catheti
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_entities import Edge, wrap_shape
from glow.interface.geom_interface import get_bounding_box, get_min_distance, \
    get_point_coordinates, make_cdg, make_compound, make_edge, \
    make_intersection
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.support.utility import build_subdvision_vertices_on_edge, \
    get_vertices_on_edges, sort_shapes_from_vertex


class Cell(Fillable):
    """
    Class for representing any cell characterised in terms of its geometry
    layout made of regions and the associated properties.

    The cell geometry layouts can be represented in terms of its technological
    geometry. The surfaces that constitute this layout represent the regions
    associated with a different material property.
    In fuel assemblies, a unit cell represents a fuel pin characterised by
    concentric circular regions. However, this class offers no limitation on
    the type of surfaces and their positioning, provided they are compatible
    with the dimensions of the cell.
    In addition, since ``Cell`` is a subclass of ``Fillable``, it can contain
    not only ``Region`` objects, but also other ``Fillable`` instances (e.g.,
    a ``Lattice`` so to model an assembly).

    In any case, the associated ``Region`` objects describe the technological
    geometry layout of the cell. Additionally, this geometry can be refined by
    subdividing it into sectors or by adding circles within a circular region.
    The GEOM edge objects resulting from the refinement are stored in the map
    associating the compound to the type of geometry.

    This class is intended to represent a generic cell as it is not based on a
    pre-defined characteristic shape, which, in turns, is provided at its
    initialisation.
    Specialisations of this class provide the shape of the cell (i.e. a
    rectangle for a Cartesian cell or a hexagon for a hexagonal one).

    Parameters
    ----------
    shape : Surface
        The ``Surface`` object representing the characteristic shape of the
        cell.
    base_props : Dict[PropertyType, str] | None = None
        The mapping from ``PropertyType`` items to values associated to the
        characteristic shape of the cell.
    name : str = "Cell"
        The name of the cell when added to the current SALOME study.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the shape of the cell.
    entry_id : str | None
        The ID attributed by SALOME when the GEOM object is added to the
        study.
    geom_obj : Any | None
        The internal `GEOM_Object` representative of the layout this instance
        refers to.
    geometry_maps : Dict[GeometryType, Compound]
        A mapping from ``GeometryType`` values to ``Compound`` objects.
        Each entry provides a different representation for the geometry
        layout this instance refers to. It is used to switch between
        different visualisation types (e.g., technological, refined).
    layers : List[List[Region | Self]]
        A list of layers, each layer itself being a list of ``Region`` objects
        or ``Fillable`` instances. Layers represent the hierarchical structure
        of the geometry layout.
    name : str | None = None
        The name of the GEOM object identifying the layout this instance
        refers to. It is used when the cell is added to the current SALOME
        study.
    o : Vertex
        The ``Vertex`` object being the centre of the GEOM object which
        represents the geometry layout of the cell.
    regions : List[Region]
        A flat list of ``Region`` objects obtained by collapsing all the
        layers. Maintained in addition to the `layers` structure to allow
        the visualization of the layout with a property colour map.
    rot_angle : float
        The rotation angle (in degrees) of the cell's GEOM object wrt the
        X-axis.
    shape : Surface
        The ``Surface`` object representing the characteristic shape of the
        cell.
    state : LayoutState
        Providing the state of the layout in the SALOME study.
    symmetry_map : Dict[SymmetryType, Face]
        A mapping from ``SymmetryType`` values to ``Face`` objects. Each
        entry provides the characteristic shape of the corresponding symmetry
        type.
    """
    def __init__(
            self,
            shape: Surface,
            base_props: Dict[PropertyType, str] | None = None,
            name: str = "Cell"
        ) -> None:
        super().__init__()
        # Set the name of this instance and of the corresponding GEOM object
        # by appending the ID of this instance to the provided name
        self.name = f"{name}_{id(self)}"
        # Store the shape of cell and build the corresponding region
        self.shape = shape
        region = Region(shape.geom_obj, properties=base_props)
        # Set the region's name to include the cell's one
        region.name = f"{self.name}_{region.name}"
        # Initialize the list of layers by putting the region at the bottom
        # layer
        self.layers = [[region]]
        self.regions = [region.clone()]
        # Update the GEOM compound this instance refers to
        self.update(wrap_shape(make_compound(self.regions)))
        self.dimensions = self.shape.dimensions

    def get_centred_regions(self) -> List[Region]:
        """
        Method that returns a list of the cell-centred ``Region`` objects of
        the cell sorted in terms of their distance from the centre of the
        cell.

        Returns
        -------
        List[Region]:
            A list of ``Region`` objects whose centres coincide with the
            cell's one, sorted by the distance from the centre.
        """
        # Re-build the regions first, if an update is required
        if self.state.is_update_needed:
            self.update_hierarchical_structure()
        return sort_shapes_from_vertex(
            [
                region for region in self.get_regions()
                if get_min_distance(region.o, self.o) < 1e-5
            ],
            self.o
        )

    def sectorize(
            self, sectors_no: List[int], angles: List[float], **kwargs: Any
        ) -> None:
        """
        Method that subdivides the technological geometry of the cell into
        sectors.
        Given the number of sectors for each region and the values of the
        angles to start the sectorization from, edges are built so that they
        propagate radially from the centre of the cell.
        The result of the intersection between each region and the subdivision
        edges is collected and stored as mapping from
        ``GeometryType.SECTORIZED`` to ``Compound`` object.

        Parameters
        ----------
        sectors_no : List[int]
            List of integers representing the number of subdivision for each
            cell region coming from the technological geometry.
        angles : List[float]
            List of angles (in degrees) the sectorization should start from
            for each cell's region coming from the technological geometry.
        kwargs : Any
            Additional parameters specific to the type of cell to be
            sectorized. No parameters must be passed.

        Raises
        ------
        TypeError
            If additional keyword arguments are provided.
        """
        if kwargs:
            raise TypeError(
                "No additional keyword arguments are expected for a "
                f"'{self.__class__.__name__}' instance."
            )
        # Apply the sectorization
        self._sectorize_cell(sectors_no, angles, False)

    def _build_sectorization_edges(
            self,
            sectors_no: List[int],
            angles: List[float],
            windmill: bool = False
        ) -> List[Edge]:
        """
        Method that builds the ``Edge`` objects that subdivide the regions
        of the technological geometry of the cell into sectors.
        The returned list stores only the edges subdividing consecutive
        cell-centred regions into sectors.

        Parameters
        ----------
        sectors_no : List[int]
            List of integers representing the number of subdivision for
            each cell-centred region.
        angles : List[float]
            List of angles (in degrees) the sectorization should start
            from for each cell-centred region.
        windmill : bool = False
            Flag indicating if a `windmill` sectorization of the cell needs
            to be considered. It defaults to ``False`` as this parameter
            does not have any meaning for this class.

        Returns
        -------
        List[Edge]
            A list of ``Edge`` objects resulting from the sectorization
            operation.
        """
        # Declare the list storing all the cell sectorization edges
        sect_edges = []
        # Associate each cell-centred region of the technological geometry
        # to a circle that fully contains it
        regions_to_circles = self._map_region_to_circle()
        # Determine the parameter identifying where a point on the circles
        # will be added. It based on the rotation of the cell and ranges
        # from 0 to 1. In case the rotation is negative, the opposite value
        # is taken into account.
        u_rotation = self.rot_angle / (2 * math.pi)
        rotation = u_rotation if self.rot_angle >= 0 else 1 + u_rotation
        # Loop through all the cell-centred regions to build the edges of
        # the sectors
        for i, (region, circle) in enumerate(regions_to_circles.items()):
            # Extract the number of subdivision sectors and the angle to
            # start the subdivision from
            sector = sectors_no[i]
            angle = angles[i]
            if sector == 1:
                # Continue as no subdivision has to be performed on the
                # current region
                continue
            # Build subdivision points on the circle containing the region.
            # The intersections between the region and each edge from the cell
            # centre to the subdivision points are sectorization edges to
            # collect
            for pnt in build_subdvision_vertices_on_edge(
                sector, circle.borders[0], angle/360.0 + rotation
            ):
                sect_edges.append(
                    wrap_shape(
                        make_intersection(region, make_edge(self.o, pnt))
                    )
                )
        # Return the built list of edges
        return sect_edges

    def _check_sectorization_elements_len(
            self,
            no_regions: int,
            sectors_no_len: int,
            angles_len: int
        ) -> None:
        """
        Method for checking the correctness of the length of the lists for
        the sectorization operation. It verifies that:
        - the length of the lists is the same;
        - their size coincides with the number of cell-centred regions of
          the technological geometry.

        Parameters
        ----------
        no_regions : int
            The number of cell-centred regions of the technological geometry.
        sectors_no_len  : int
            The length of the list providing the number of sectors each of
            the cell-centred regions has to be subdivided into.
        angles_len : int
            The length of the list providing the angles the sectorization
            starts from for each cell-centred region.
        """
        if not sectors_no_len == angles_len:
            raise ValueError(
                f"The numbers of sectors ({sectors_no_len}) and angles to "
                f"start the sectorization from ({angles_len}) do not "
                "coincides."
            )
        if not sectors_no_len == no_regions:
            raise ValueError(
                "The number of elements in the lists for the sectorization "
                f"({sectors_no_len}) does not coincide with the number of "
                f"cell regions ({no_regions})."
            )

    def _map_region_to_circle(self) -> Dict[Region, Circle]:
        """
        Method that associates a ``Circle`` object to each cell-centred
        region.
        For each region, this method computes the diagonal of its bounding
        box and creates a circle centred in the cell's origin.
        The resulting mapping pairs each region with its corresponding
        circle.

        Returns
        -------
        Dict[Region, Circle]
            A dictionary mapping each region to a circle that fully contains
            it.
        """
        # Associate a 'Circle' object to each cell-centred region
        regions_to_circles: Dict[Region, Circle] = {}
        for r in self.get_centred_regions():
            # Get the diagonal of the bounding box for the current region
            x_min, x_max, y_min, y_max = get_bounding_box(r)
            width = x_max - x_min
            heigth = y_max - y_min
            diag = math.sqrt(width*width + heigth*heigth)
            # Build a circle containing the current region and associate it
            # to the region
            regions_to_circles[r] = Circle(
                get_point_coordinates(self.o), diag/2
            )
        return regions_to_circles

    def _sectorize_cell(
            self,
            sectors_no: List[int],
            angles: List[float],
            windmill: bool = False
        ) -> None:
        """
        Method that performs the sectorization of the cell-centred regions
        belonging to the technological geometry of the cell.
        Given the number of sectors for each region and the values of the
        angles to start the sectorization from, edges are built so that they
        propagate radially from the centre of the cell.
        The result of the intersection between each region and the subdivision
        edges is collected and stored as mapping from
        ``GeometryType.SECTORIZED`` to ``Compound`` object.

        Parameters
        ----------
        sectors_no : List[int]
            List of integers representing the number of sectors into which
            each cell-centred region has to be subdivided into.
        angles : List[float]
            List of angles (in degrees) the sectorization should start from
            for each cell-centred region.
        windmill : bool = False
            Flag indicating if a `windmill` sectorization of the cell needs
            to be considered.
        """
        # Apply the sectorization only if the cell contains 'Region' objects,
        # not 'Fillable' objects
        if any(
            isinstance(layout, Fillable)
            for layer in self.layers
            for layout in layer
        ):
            raise RuntimeError(
                "No sectorization available for a cell containing " \
                "other 'Fillable' objects."
            )
        # Check the correctness of the lists storing the information for
        # performing the sectorization
        self._check_sectorization_elements_len(
            len(self.get_centred_regions()),
            len(sectors_no),
            len(angles)
        )
        # Handle the sectorization operation: the resulting list holds all
        # the edges identifying:
        # . the sectors between consecutive cell-centred regions;
        # . the lines connecting points as the result of the intersection of
        #   the sector edges for the outmost region and the cell's borders (if
        #   the 'windmill' option is active)
        edges = self._build_sectorization_edges(sectors_no, angles, windmill)
        # Store the GEOM compound of the edges in the mapping
        self.geometry_maps[GeometryType.SECTORIZED] = make_compound(edges)

    def _translate_layout_specific_elems(
            self, new_cntr: Tuple[float, float, float]
        ) -> None:
        """
        Method for translating the shape of the cell, considering the centre
        of the translated layout is positioned at the given XYZ coordinates.

        Parameters
        ----------
        new_cntr : Tuple[float, float, float]
            The XYZ coordinates of the centre of the translated layout.
        """
        self.shape.translate(new_cntr)

    def __repr__(self):
        """
        Return a descriptive string of this instance containing a
        descriptive string with the name and the coordinates of the
        centre of the cell

        Returns
        -------
        str
            A descriptive string of the current instance.
        """
        return f"{self.name}, {wrap_shape(make_cdg(self.geom_obj))}"


class CartesianCell(Cell):
    """
    Class for representing a Cartesian cell characterised in terms of its
    geometry layouts and the associated properties.
    It is based on a rectangular shape whose corresponding ``Surface`` object
    is built from the geometric settings provided at its initialisation.

    The cell geometry layout can be represented in terms of its technological
    geometry. The surfaces that constitute this layout represent the regions
    associated with a different material property.
    In fuel assemblies, a unit cell represents a fuel pin characterised by
    concentric circular regions. However, this class offers no limitation on
    the type of surfaces and their positioning, provided they are compatible
    with the dimensions of the cell.
    In addition, since it is a subclass of ``Cell``, and consequently of
    ``Fillable``, it can contain not only ``Region`` objects, but also other
    ``Fillable`` instances (e.g., a ``Lattice`` so to model an assembly).

    In any case, the associated ``Region`` objects describe the technological
    geometry layout of the cell. Additionally, this geometry can be refined by
    subdividing it into sectors or by adding circles within a circular region.
    The GEOM edge objects resulting from the refinement are stored in the map
    associating the compound to the type of geometry.

    Parameters
    ----------
    center : Tuple[float, float, float] | None = None
        The XYZ coordinates of the centre of the cell.
    width_height : Tuple[float, float] = (1.0, 1.0)
        A tuple providing the width and the height of the rectangle
        representing the shape of the cell.
    rounded_corners : List[Tuple[int, float]] | None = None
        A list of tuples providing for each corner of the rectangle its index
        and the curvature radius, if any.
    base_props : Dict[PropertyType, str] | None = None
        The mapping from ``PropertyType`` items to values associated to the
        characteristic rectangular shape of the cell.
    name : str = "Cartesian_Cell"
        The name of the cell when added to the current SALOME study.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the shape of the cell.
    entry_id : str | None
        The ID attributed by SALOME when the GEOM object is added to the
        study.
    geom_obj : Any | None
        The internal `GEOM_Object` representative of the layout this instance
        refers to.
    geometry_maps : Dict[GeometryType, Compound]
        A mapping from ``GeometryType`` values to ``Compound`` objects.
        Each entry provides a different representation for the geometry
        layout this instance refers to. It is used to switch between
        different visualisation types (e.g., technological, refined).
    layers : List[List[Region | Self]]
        A list of layers, each layer itself being a list of ``Region`` objects
        or ``Fillable`` instances. Layers represent the hierarchical structure
        of the geometry layout.
    name : str | None = None
        The name of the GEOM object identifying the layout this instance
        refers to. It is used when the cell is added to the current SALOME
        study.
    o : Vertex
        The ``Vertex`` object being the centre of the GEOM object which
        represents the geometry layout of the cell.
    regions : List[Region]
        A flat list of ``Region`` objects obtained by collapsing all the
        layers. Maintained in addition to the `layers` structure to allow
        the visualization of the layout with a property colour map.
    rot_angle : float
        The rotation angle (in degrees) of the cell's GEOM object wrt the
        X-axis.
    shape : Surface
        The ``Surface`` object representing the characteristic shape of the
        cell.
    state : LayoutState
        Providing the state of the layout in the SALOME study.
    symmetry_map : Dict[SymmetryType, Face]
        A mapping from ``SymmetryType`` values to ``Face`` objects. Each
        entry provides the characteristic shape of the corresponding symmetry
        type.
    """
    def __init__(
            self,
            center: Tuple[float, float, float] | None = None,
            width_height: Tuple[float, float] = (1.0, 1.0),
            rounded_corners: List[Tuple[int, float]] | None = None,
            base_props: Dict[PropertyType, str] | None = None,
            name: str = "Cartesian_Cell"
        ) -> None:
        # Initialize the superclass passing the rectangular shape
        super().__init__(
            Rectangle(
                center=center,
                width=width_height[0],
                height=width_height[1],
                rounded_corners=rounded_corners
            ),
            base_props,
            name
        )

    def sectorize(
            self, sectors_no: List[int], angles: List[float], **kwargs: Any
        ) -> None:
        """
        Method that subdivides the technological geometry of the cell into
        sectors.
        Given the number of sectors for each region and the values of the
        angles to start the sectorization from, edges are built so that they
        propagate radially from the centre of the cell.
        The result of the intersection between each region and the subdivision
        edges is collected and stored as mapping from
        ``GeometryType.SECTORIZED`` to ``Compound`` object.

        Parameters
        ----------
        sectors_no : List[int]
            List of integers representing the number of subdivision for each
            cell region coming from the technological geometry.
        angles : List[float]
            List of angles (in degrees) the sectorization should start from
            for each cell's region coming from the technological geometry.
        kwargs : Any
            Additional parameters specific to the cartesian cell type to be
            sectorized (e.g., windmill).
        """
        # Check if the 'windmill' flag, stating that a windmill sectorization
        # of the cell needs to be considered, has been passed. If not, 'False'
        # is passed to the method that performs the cell sectorization
        if 'windmill' not in kwargs:
            self._sectorize_cell(sectors_no, angles, False)
        else:
            self._sectorize_cell(sectors_no, angles, kwargs['windmill'])

    def _apply_windmill_sectorization(
            self,
            last_sect_no: int,
            last_sect_angle: float,
            is_windmill: bool,
            sect_edges: List[Edge]
        ) -> List[Edge]:
        """
        Method that builds the wings of the windmill sectorization, if needed.
        The method connects two successive points of the sectorization edges
        for the outmost region which lay on the cell's borders.
        The edges are built only if required (``is_windmill`` parameter set
        to `True`) and the number of sectors for the outmost region of the
        technological geometry is either 8 or 16.
        Otherwise, nothing is performed and an empty list is returned.

        Parameters
        ----------
        last_sect_no : int
            The number of sectors for the outmost region of the cell's
            technological geometry.
        last_sect_angle : float
            The starting angle of the sectorization for the outmost region.
        is_windmill : bool
            Flag stating whether the windmill sectorization is needed.
        sect_edges : List[Edge]
            The list of edges resulting from applying the sectorization.

        Returns
        -------
        List[Edge]
            The list of ``Edge`` objects representing the wings of the
            windmill sectorization, if required, otherwise an empty list.
        """
        # Return immediately if the windmill flag is False
        if not is_windmill:
            return []
        # Declare the settings for supported windmill layouts and get the ones
        # corresponding to the number of sectors for the outmost region
        w_settings = {
            8: {"angles": [22.5 + j * 45 for j in range(8)], "offset": 1},
            16: {"angles": [0.0 + j * 45 for j in range(16)], "offset": 2},
        }
        current_config = w_settings.get(last_sect_no)
        # Return an empty list if no configuration is found or the starting
        # angle of the outmost region is not among the valid ones for the
        # windmill sectorization
        if current_config is None:
            return []
        if not any(
            math.isclose(last_sect_angle, a, abs_tol=1e-6)
            for a in current_config["angles"]
        ):
            return []

        # Get the vertices of the sectorization edges on the cell's
        # borders, sorted in counterclockwise order
        last_region_edges = sect_edges[-last_sect_no:]
        vertices = get_vertices_on_edges(
            last_region_edges, self.shape.borders, self.o
        )
        # Return the wing edges built between the indicated vertices
        offset = current_config["offset"]
        return [
            wrap_shape(make_edge(vertices[i], vertices[i + offset]))
            for i in range(0, last_sect_no, 2 * offset)
        ]

    def _build_sectorization_edges(
            self,
            sectors_no: List[int],
            angles: List[float],
            windmill: bool = False
        ) -> List[Edge]:
        """
        Method that builds the ``Edge`` objects that subdivide the regions
        into sectors.
        The returned list stores:
        - the sectors between consecutive cell-centred regions;
        - the lines connecting points as the result of the intersection of
          the sector edges for the outmost region and the cell's borders (if
          the ``windmill`` flag is ``True``).

        Parameters
        ----------
        sectors_no : List[int]
            List of integers representing the number of subdivision for
            each cell's region coming from the technological geometry.
        angles : List[float]
            List of angles (in degrees) the sectorization should start
            from for each cell's region coming from the technological
            geometry.
        windmill : bool = False
            Flag indicating if a `windmill` sectorization of the cell needs
            to be considered.

        Returns
        -------
        List[Edge]
            A list of ``Edge`` objects resulting from the sectorization
            operation.
        """
        # Call the method of the superclass to collect the sectorization edges
        sect_edges = super()._build_sectorization_edges(sectors_no, angles)
        # Build and store the wings of the sectorization, if required
        sect_edges += self._apply_windmill_sectorization(
            sectors_no[-1], angles[-1], windmill, sect_edges
        )
        # Return the list collecting all the edges of the sectorization
        return sect_edges

    def _build_symmetry_shape(
            self,
            symmetry: SymmetryType,
            x_min_max: Tuple[float, float],
            y_min_max: Tuple[float, float],
            o_xyz: Tuple[float, float, float]
        ) -> Surface:
        """
        Method that builds the geometric shape corresponding to the given
        symmetry type for a Cartesian cell.

        Supported symmetries are:
        - FULL: returns the characteristic shape of the cell.
        - HALF: builds a rectangle representing half of the cell.
        - QUARTER: builds a rectangle representing one quarter of the cell.
        - EIGHTH: builds a right triangle representing one eighth of the cell.
        - DIAG: builds a right triangle representing the diagonal symmetry of
          the cell.

        Parameters
        ----------
        symmetry : SymmetryType
            The symmetry type for which the characteristic shape is built.
        x_min_max : tuple[float, float]
            The minimum and maximum extension of the geometry bounding box
            along X-axis.
        y_min_max : tuple[float, float]
            The minimum and maximum extension of the geometry bounding box
            along Y-axis.
        o_xyz : tuple[float, float, float]
            The XYZ coordinates of the cell centre.

        Returns
        -------
        Surface
            A geometric shape representing the requested symmetry.

        Raises
        ------
        RuntimeError
            If the indicated symmetry type is not supported for a Cartesian
            cell.
        """
        # Get the XY dimensions of the bounding box for the geometry layout
        x_min, x_max = x_min_max
        y_min, y_max = y_min_max
        # Get the coordinates of the cell centre
        o_x, o_y, o_z = o_xyz
        # Match the shape with the type of symmetry
        match symmetry:
            case SymmetryType.DIAG:
                return build_right_triangle_from_catheti(
                    x_max - x_min, y_max - y_min, (x_min, y_min, o_z)
                )
            case SymmetryType.EIGHTH:
                return build_right_triangle_from_catheti(
                    x_max - o_x, y_max - o_y, (o_x, o_y, o_z)
                )
            case _:
                return super()._build_symmetry_shape(
                    symmetry, x_min_max, y_min_max, o_xyz
                )

class HexCell(Cell):
    """
    Class for representing a hexagonal cell characterised in terms of its
    geometry layouts and the associated properties.
    It is based on a hexagonal shape whose corresponding ``Surface`` object
    is built from the geometric settings provided at its initialisation.

    The cell geometry layout can be represented in terms of its technological
    geometry. The surfaces that constitute this layout represent the regions
    associated with a different material property.
    In fuel assemblies, a unit cell represents a fuel pin characterised by
    concentric circular regions. However, this class offers no limitation on
    the type of surfaces and their positioning, provided they are compatible
    with the dimensions of the cell.
    In addition, since it is a subclass of ``Cell``, and consequently of
    ``Fillable``, it can contain not only ``Region`` objects, but also other
    ``Fillable`` instances (e.g., a ``Lattice`` so to model an assembly).

    In any case, the associated ``Region`` objects describe the technological
    geometry layout of the cell. Additionally, this geometry can be refined by
    subdividing it into sectors or by adding circles within a circular region.
    The GEOM edge objects resulting from the refinement are stored in the map
    associating the compound to the type of geometry.

    Parameters
    ----------
    center : Tuple[float, float, float] | None = None
        The XYZ coordinates of the centre of the cell.
    side : float = 1.0
        The length of the side of the regular hexagon representing the shape
        of the cell.
    name : str = "Hexagonal_Cell"
        The name of the cell when added to the current SALOME study.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the shape of the cell.
    entry_id : str | None
        The ID attributed by SALOME when the GEOM object is added to the
        study.
    geom_obj : Any | None
        The internal `GEOM_Object` representative of the layout this instance
        refers to.
    geometry_maps : Dict[GeometryType, Compound]
        A mapping from ``GeometryType`` values to ``Compound`` objects.
        Each entry provides a different representation for the geometry
        layout this instance refers to. It is used to switch between
        different visualisation types (e.g., technological, refined).
    layers : List[List[Region | Self]]
        A list of layers, each layer itself being a list of ``Region`` objects
        or ``Fillable`` instances. Layers represent the hierarchical structure
        of the geometry layout.
    name : str | None = None
        The name of the GEOM object identifying the layout this instance
        refers to. It is used when the cell is added to the current SALOME
        study.
    o : Vertex
        The ``Vertex`` object being the centre of the GEOM object which
        represents the geometry layout of the cell.
    regions : List[Region]
        A flat list of ``Region`` objects obtained by collapsing all the
        layers. Maintained in addition to the `layers` structure to allow
        the visualization of the layout with a property colour map.
    rot_angle : float
        The rotation angle (in degrees) of the cell's GEOM object wrt the
        X-axis.
    shape : Surface
        The ``Surface`` object representing the characteristic shape of the
        cell.
    state : LayoutState
        Providing the state of the layout in the SALOME study.
    symmetry_map : Dict[SymmetryType, Face]
        A mapping from ``SymmetryType`` values to ``Face`` objects. Each
        entry provides the characteristic shape of the corresponding symmetry
        type.
    """

    def __init__(
            self,
            center: Tuple[float, float, float] | None = None,
            side: float = 1.0,
            base_props: Dict[PropertyType, str] | None = None,
            name: str = "Hexagonal_Cell"
        ) -> None:
        # Initialize the superclass passing the rectangular shape
        super().__init__(
            Hexagon(
                center=center,
                edge_length=side
            ),
            base_props,
            name
        )

    def sectorize(
            self, sectors_no: List[int], angles: List[float], **kwargs: Any
        ) -> None:
        """
        Method that subdivides the technological geometry of the cell into
        sectors.
        Given the number of sectors for each region and the values of the
        angles to start the sectorization from, edges are built so that they
        propagate radially from the centre of the cell.
        The result of the intersection between each region and the subdivision
        edges is collected and stored as mapping from
        ``GeometryType.SECTORIZED`` to ``Compound`` object.

        Parameters
        ----------
        sectors_no : List[int]
            List of integers representing the number of subdivision for each
            cell region coming from the technological geometry.
        angles : List[float]
            List of angles (in degrees) the sectorization should start from
            for each cell's region coming from the technological geometry.
        kwargs : Any
            Additional parameters specific to the type of cell to be
            sectorized. No parameters must be passed.

        Raises
        ------
        TypeError
            If additional keyword arguments are provided.
        """
        if kwargs:
            raise TypeError(
                "No additional keyword arguments are expected for a "
                f"'{self.__class__.__name__}' instance."
            )
        # Since the hexagonal cell type does not come with a 'windmill'
        # sectorization, 'False' is passed to the method that performs
        # the cell sectorization
        self._sectorize_cell(sectors_no, angles, False)

    def _build_symmetry_shape(
            self,
            symmetry: SymmetryType,
            x_min_max: Tuple[float, float],
            y_min_max: Tuple[float, float],
            o_xyz: Tuple[float, float, float]
        ) -> Surface:
        """
        Method that builds the geometric shape corresponding to the given
        symmetry type for a hexagonal cell.

        Supported symmetries are:
        - FULL: returns the characteristic shape of the cell.
        - HALF: builds a rectangle representing half of the cell.
        - QUARTER: builds a rectangle representing one quarter of the cell.
        - THIRD: builds a parallelogram representing one third of the cell.
        - SIXTH: builds a regular triangle representing one sixth of the cell.
        - TWELFTH: builds a right triangle representing one twelfth of the
          cell.

        Parameters
        ----------
        symmetry : SymmetryType
            The symmetry type for which the characteristic shape is built.
        x_min_max : tuple[float, float]
            The minimum and maximum extension of the geometry bounding box
            along X-axis.
        y_min_max : tuple[float, float]
            The minimum and maximum extension of the geometry bounding box
            along Y-axis.
        o_xyz : tuple[float, float, float]
            The XYZ coordinates of the cell centre.

        Returns
        -------
        Surface
            A geometric shape representing the requested symmetry.

        Raises
        ------
        RuntimeError
            If the indicated symmetry type is not supported for a hexagonal
            cell.
        """
        # Get the XY dimensions of the bounding box for the geometry layout
        x_min, x_max = x_min_max
        y_min, y_max = y_min_max
        # Get the coordinates of the cell centre
        o_x, o_y, o_z = o_xyz
        # Get the characteristic dimension of the hexagon
        side_len = (x_max - x_min) / 2
        apothem = (y_max - y_min) / 2
        # Match the shape with the type of symmetry
        match symmetry:
            case SymmetryType.THIRD:
                return build_parallelogram(
                    side_len,
                    side_len,
                    60,
                    (-(o_x + side_len/2), -(o_y + apothem), o_z)
                )
            case SymmetryType.SIXTH:
                return build_regular_triangle(
                    side_len,
                    (-(o_x + side_len/2), -(o_y + apothem), o_z)
                )
            case SymmetryType.TWELFTH:
                return build_right_triangle(
                    side_len, side_len*math.cos(math.pi/6), o_xyz
                )
            case _:
                return super()._build_symmetry_shape(
                    symmetry, x_min_max, y_min_max, o_xyz
                )
