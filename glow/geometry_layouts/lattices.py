"""
Module containing classes providing the means for creating lattices from base
cells types in SALOME.
"""
from copy import deepcopy
from typing import List, Tuple

from glow.geometry_layouts.cells import Cell
from glow.geometry_layouts.fillable_layouts import Fillable
from glow.geometry_layouts.geometries import GenericSurface, Hexagon, \
    Rectangle, Surface
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_entities import Edge, wrap_shape
from glow.interface.geom_interface import ShapeType, get_closed_free_boundary, \
    get_point_coordinates, make_compound, make_face, make_partition, \
    make_vertex
from glow.support.utility import are_same_shapes, build_compound_borders, \
    build_subdvision_vertices_on_edge, build_z_axis_from_vertex


class Lattice(Fillable):
    """
    Class for representing any lattice characterised in terms of its geometry
    layout made of cells that do not follow a specific pattern.
    This class can be used to model a portion of a generic lattice assembled
    by simply positioning the cells one by one. In this sense, it could be
    useful to represent a portion of a colorset without the need to build all
    the assemblies around the central one.
    This class does not support the placement of one or several rings of the
    same cell as it serves to describe a generic pattern.
    Subclasses of ``Lattice`` provide specialised methods to address this need
    according to the type of pattern, i.e. either Cartesian or hexagonal.

    Parameters
    ----------
    cells : List[Cell] = []
        The list of cells that constitute the lattice, as objects of the class
        ``Cell`` or of its subclasses.
    centre : Tuple[float, float, float] | None = None
        The coordinates of the lattice centre, if any.
    shape : Surface | None = None
        The ``Surface`` object representing the characteristic shape of the
        lattice.
    base_props : Dict[PropertyType, str] | None = None
        The mapping from ``PropertyType`` items to values associated to the
        characteristic shape of the lattice.
    name : str = "Lattice"
        The name of the lattice when added to the current SALOME study.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the shape of the lattice.
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
        refers to. It is used when the lattice is added to the current SALOME
        study.
    o : Vertex
        The ``Vertex`` object being the centre of the GEOM object which
        represents the geometry layout of the lattice.
    regions : List[Region]
        A flat list of ``Region`` objects obtained by collapsing all the
        layers. Maintained in addition to the `layers` structure to allow
        the visualization of the layout with a property colour map.
    rot_angle : float
        The rotation angle (in degrees) of the lattice's GEOM object wrt the
        X-axis.
    shape : Surface
        The ``Surface`` object representing the characteristic shape of the
        lattice.
    state : LayoutState
        Providing the state of the layout in the SALOME study.
    symmetry_map : Dict[SymmetryType, Face]
        A mapping from ``SymmetryType`` values to ``Face`` objects. Each
        entry provides the characteristic shape of the corresponding symmetry
        type.
    """
    def __init__(
            self,
            cells: List[Cell] = [],
            centre: Tuple[float, float, float] | None = None,
            name: str = "Lattice"
        ) -> None:
        super().__init__()
        # Set the name of this instance and of the corresponding GEOM object
        # by appending the ID of this instance to the provided name
        self.name = f"{name}_{id(self)}"
        # Initialize the list of layers by putting the given cells at the
        # bottom layer, if any are provided. Regions are extracted from each
        # cell
        if cells:
            self.layers = [[c.clone() for c in cells]]
            self.regions = self.get_regions()
            # Set the shape and the GEOM compound of the layout from the
            # compound of the regions
            # TODO the shape should update whenever the dimensions changes
            # i.e. a new ring is added
            self.shape = GenericSurface(
                make_face(build_compound_borders(make_compound(self.regions)))
            )
            # Update the GEOM compound this instance refers to
            self.update(self.shape)
            self.dimensions = self.shape.dimensions
        # Set the vertex of the layout centre, if any is provided, otherwise
        # the XYZ origin is used
        if centre is not None:
            self.o = make_vertex(centre if centre else (0.0, 0.0, 0.0))

    def add(self,
            layout: Region | Fillable,
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
        To collapse the layers and update all the contained layouts of this
        instance without displaying the geometry in the 3D viewer of SALOME,
        call the method ``update_hierarchical_structure``.
        To update the and display the geometry layout with all the regions of
        the technological geometry, call the method ``show``.
        """
        # Call the 'add()' method of the 'Fillable' superclass
        super().add(layout, position, layer_index)
        # Update the lattice characteristic shape and dimensions
        self.shape = GenericSurface(
            make_face(
                build_compound_borders(make_compound(self.get_regions()))
            )
        )
        self.dimensions = tuple(
            c_1 + c_2
            for c_1, c_2 in zip(self.dimensions, self.shape.dimensions)
        )

    def _compute_layer_index(self, layer_index: int | None) -> int:
        """
        Method that returns the layer index to which cells have to be added.
        Depending on the value of ``layer_index``, different actions are
        taken:
        - If ``None`` or equal to ``len(self.layers)``, a new empty layer is
          created and appended to the ``self.layers`` attribute.
        - If in ``[0, len(self.layers) - 1]``, the given index is returned.

        Parameters
        ----------
        layer_index : int | None
            The index of the layer to check and return.

        Returns
        -------
        int
            The validated layer index.

        Raises
        ------
        ValueError
            If ``layer_index`` is not ``None`` and does not fall within the
            valid range of existing layers.
        """
        if layer_index is None or layer_index == len(self.layers):
            # Initialise an empty sublist and return the previous length of
            # the layers list
            self.layers.append([])
            return len(self.layers) - 1
        elif 0 <= layer_index < len(self.layers):
            # The indicated layer index already exists
            return layer_index
        else:
            raise ValueError(f"Invalid layer index {layer_index}.")

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
        if self.shape is not None:
            self.shape.translate(new_cntr)


class CartesianLattice(Lattice):
    """
    Class for representing a Cartesian lattice characterised in terms of its
    geometry layout made of cells arranged according to a rectangular grid.
    This class can be used to model a full Cartesian lattice assembled either
    by providing its cells directly or successively calling the methods for
    adding one or more cells to the layout.

    Parameters
    ----------
    cells : List[Cell] = []
        The list of cells that constitute the lattice, as objects of the class
        ``Cell`` or of its subclasses.
    centre : Tuple[float, float, float] | None = None
        The coordinates of the lattice centre, if any.
    name : str = "CartesianLattice"
        The name of the lattice when added to the current SALOME study.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the shape of the lattice.
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
        refers to. It is used when the lattice is added to the current SALOME
        study.
    o : Vertex
        The ``Vertex`` object being the centre of the GEOM object which
        represents the geometry layout of the lattice.
    regions : List[Region]
        A flat list of ``Region`` objects obtained by collapsing all the
        layers. Maintained in addition to the `layers` structure to allow
        the visualization of the layout with a property colour map.
    rot_angle : float
        The rotation angle (in degrees) of the lattice's GEOM object wrt the
        X-axis.
    shape : Surface
        The ``Surface`` object representing the characteristic shape of the
        lattice.
    state : LayoutState
        Providing the state of the layout in the SALOME study.
    symmetry_map : Dict[SymmetryType, Face]
        A mapping from ``SymmetryType`` values to ``Face`` objects. Each
        entry provides the characteristic shape of the corresponding symmetry
        type.
    """
    def __init__(
            self,
            cells: List[Cell] = [],
            centre: Tuple[float, float, float] | None = None,
            name: str = "CartesianLattice"
        ) -> None:
        # Initialise the superclass
        super().__init__(cells=cells, centre=centre, name=name)
        # Update the state of the lattice
        self.state.is_update_needed = True

    def add_ring_of_cells(
            self,
            cell: Cell,
            ring_index: int,
            layer_index: int | None = None
        ) -> None:
        """
        Method that adds a ring of ``Cell`` objects to the lattice at the
        positions that correspond to the given ring index and layer.
        This method iteratively adds the provided ``Cell`` object at specific
        construction points determined by the `ring_index` parameter, which
        indicates the ring where the cells should be placed.
        Optionally, a ``layer_index`` can be provided to specify the layer to
        which the ring of cells is added; if not provided, the cells are
        added to a new layer.

        Parameters
        ----------
        cell : Cell
            The ``Cell`` instance to be added repeatedly to form a ring.
        ring_index : int
            The index indicating which ring to add the cells to.
        layer_index : int | None = None
            The index of the layer to which the cells are added. If ``None``,
            a new layer is created.

        Raises
        ------
        ValueError
            If indicating ``0`` as the ring index where cells should be added.
        ValueError
            If indicating a layer index not falling within the valid range
            of existing layers.

        Notes
        -----
        Index 0 refers to the lattice centre (for lattices with an even
        number of cells on one if the sides) or the central cell (valid for
        lattices with an odd number of cells on both sides); indices greater
        than 0 refer to rings around the centre.
        """
        # Check the validity of the indicated ring index
        ensure_not_zero(
            ring_index,
            "It is not possible to add a ring of cells at the given 0 index"
        )
        # Include the given layout at the end of the layer sublist specified
        # by the indicated index, if any; otherwise, either create a new
        # sublist or raise an exception
        layer_index = self._compute_layer_index(layer_index)
        # Get the cell's dimensions
        cell_width = cell.shape.dimensions[0]
        cell_heigth = cell.shape.dimensions[1]

        # Evaluate the multiplication factor for determining the construction
        # figure where the centers of the cells will be placed
        n = self._evaluate_ring_factor(ring_index)
        # Build the construction figure
        construction_fig = self._build_ring_construction_figure(
            cell_width, cell_heigth, n
        )
        # Compute the positions of the centres of the cells along the
        # indicated ring
        centres = compute_subdivision_points_on_borders(n, construction_fig)

        # If needed, update the cell's tree before adding the cell
        cell.update_hierarchical_structure()
        # Add the cell at the given centre positions
        self.layers[layer_index].extend(
            get_cell_at_centres(
                cell,
                centres,
                self.rot_angle,
                wrap_shape(build_z_axis_from_vertex(self.o))
            )
        )
        # Update the characteristic dimensions of the lattice, its shape and
        # its state
        self.dimensions = (
            self.dimensions[0] + (ring_index + 1) * cell_width,
            self.dimensions[1] + (ring_index + 1) * cell_heigth
        )
        self.state.is_update_needed = True

    def add_rings_of_cells(
            self,
            cell: Cell,
            no_rings: int,
            ring_index: int,
            layer_index: int | None = None
        ) -> None:
        """
        Method that adds several rings of ``Cell`` objects to the lattice at
        the positions calculated for each ring indicated by the total number
        of rings, starting from the `ring_index` parameter.
        Optionally, a ``layer_index`` can be provided to specify the layer to
        which the ring of cells is added; if not provided, the cells are
        added to a new layer.

        Parameters
        ----------
        cell : Cell
            The ``Cell`` instance to be iteratively added in order to build
            the rings of cells.
        no_rings : int
            The number of rings to add starting from the indicated ring index.
        ring_index : int
            The index indicating the starting ring index.
        layer_index : int | None = None
            The index of the layer to which the cells are added. If ``None``,
            a new layer is created.

        Raises
        ------
        ValueError
            If indicating a value less than ``1`` for the total number of
            rings to add.
        ValueError
            If indicating ``0`` as the ring index from which the rings of
            cells should be added.
        ValueError
            If indicating a layer index not falling within the valid range
            of existing layers.

        Notes
        -----
        Index 0 refers to the lattice centre (for lattices with an even
        number of cells on one if the sides) or the central cell (valid for
        lattices with an odd number of cells on both sides); indices greater
        than 0 refer to rings around the centre.
        """
        # Raise an exception if the number of rings to add is less than 1
        if no_rings < 1:
            raise ValueError(
                f"Wrong number ({no_rings}) of rings of cells to add has "
                "been indicated."
            )
        # Check the validity of the indicated ring index
        ensure_not_zero(
            ring_index,
            "It is not possible to add a ring of cells at the given 0 index"
        )
        # Include the layouts at the end of the layer sublist specified
        # by the indicated index, if any; otherwise, either create a new
        # sublist or raise an exception
        layer_index = self._compute_layer_index(layer_index)
        # Loop through the ring indices starting from the indicated one
        for i_ring in range(ring_index, ring_index+no_rings):
            # Add a ring of cells at the current index
            self.add_ring_of_cells(cell, i_ring, layer_index)
        # Update the characteristic dimensions of the lattice, its shape and
        # its state
        self.dimensions = (
            self.dimensions[0]
                + no_rings*(ring_index + 1) * cell.dimensions[0],
            self.dimensions[1]
                + no_rings*(ring_index + 1) * cell.dimensions[1]
        )
        self.state.is_update_needed = True

    def _build_ring_construction_figure(
            self, width: float, heigth: float, n: int
        ) -> Rectangle:
        """
        Method that builds and return a ``Rectangle`` object representing
        the construction figure on whose borders the cells belonging to the
        same ring will be placed. The dimensions of the rectangle depends on
        the given dimensions and the multiplication factor.
        The resulting figure is built with centre in the lattice centre, and,
        if needed, is rotated accordingly with the rotation angle of the
        lattice.

        Parameters
        ----------
        width : float
            The width of the ``Rectangle`` object.
        heigth : float
            The heigth of the ``Rectangle`` object.
        n : int
            The scaling factor for the dimensions of the ``Rectangle`` object.

        Returns
        -------
        Rectangle
            The ``Rectangle`` object built accordingly with the given
            dimensions.
        """
        construction_fig = Rectangle(
            get_point_coordinates(self.o), n*heigth, n*width
        )
        # Rotate the construction figure, if needed
        construction_fig.rotate(self.rot_angle)
        return construction_fig

    def _evaluate_ring_factor(self, ring_indx: int) -> int:
        """
        Method that evaluates the multiplication factor used to derive
        the construction figure along whose borders the cells belonging
        to the same ring index will be placed.
        The returned value depends on whether a central cell is present.

        Parameters
        ----------
        ring_indx : int
            The index of the ring for which the factor has to be derived.

        Returns
        -------
        int
            The multiplication factor dependent on the presence of a central
            cell in the Cartesian lattice.
        """
        for layer in self.layers:
            for layout in layer:
                if are_same_shapes(layout.o, self.o, ShapeType.VERTEX):
                    return 2*ring_indx
            else:
                continue
        else:
            return 2*ring_indx - 1


class HexLattice(Lattice):
    """
    Class for representing a hexagonal lattice characterised in terms of its
    geometry layout made of cells arranged according to a hexagonal grid.
    This class can be used to model a full hexagonal lattice assembled either
    by providing its cells directly or successively calling the methods for
    adding one or more cells to the layout.
    By default, the hexagonal lattice expects cells are provided so that the
    characteristic shape of the lattice follows a X-oriented hexagon. This
    means that cells should be rotated by 90° before being added.

    Parameters
    ----------
    cells : List[Cell] = []
        The list of cells that constitute the lattice, as objects of the class
        ``Cell`` or of its subclasses. Considered if no ``setup`` is provided.
    centre : Tuple[float, float, float] | None = None
        The coordinates of the lattice centre, if any.
    name : str = "Lattice"
        The name of the lattice when added to the current SALOME study.

    Attributes
    ----------
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the shape of the lattice.
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
        refers to. It is used when the lattice is added to the current SALOME
        study.
    o : Vertex
        The ``Vertex`` object being the centre of the GEOM object which
        represents the geometry layout of the lattice.
    regions : List[Region]
        A flat list of ``Region`` objects obtained by collapsing all the
        layers. Maintained in addition to the `layers` structure to allow
        the visualization of the layout with a property colour map.
    rot_angle : float
        The rotation angle (in degrees) of the lattice's GEOM object wrt the
        X-axis.
    shape : Surface
        The ``Surface`` object representing the characteristic shape of the
        lattice.
    state : LayoutState
        Providing the state of the layout in the SALOME study.
    symmetry_map : Dict[SymmetryType, Face]
        A mapping from ``SymmetryType`` values to ``Face`` objects. Each
        entry provides the characteristic shape of the corresponding symmetry
        type.
    """
    def __init__(
            self,
            cells: List[Cell] = [],
            centre : Tuple[float, float, float] | None = None,
            name: str = "Lattice"
        ) -> None:
        # Initialise the superclass
        super().__init__(cells=cells, centre=centre, name=name)
        # Update the state of the lattice
        self.state.is_update_needed = True

    def add_ring_of_cells(
            self,
            cell: Cell,
            ring_index: int,
            layer_index: int | None = None
        ) -> None:
        """
        Method that adds a ring of ``Cell`` objects to the lattice at the
        positions that correspond to the given ring index and layer.
        This method iteratively adds the provided ``Cell`` object at specific
        construction points determined by the `ring_index` parameter, which
        indicates the ring where the cells should be placed.
        Optionally, a ``layer_index`` can be provided to specify the layer to
        which the ring of cells is added; if not provided, the cells are
        added to a new layer.

        Parameters
        ----------
        cell : Cell
            The ``Cell`` instance to be added repeatedly to form a ring.
        ring_index : int
            The index indicating which ring to add the cells to.
        layer_index : int | None = None
            The index of the layer to which the cells are added. If ``None``,
            a new layer is created.

        Raises
        ------
        ValueError
            If indicating ``0`` as the ring index where cells should be added.
        ValueError
            If indicating a layer index not falling within the valid range
            of existing layers.

        Notes
        -----
        Index 0 refers to the central cell of the lattice; indices greater
        than 0 refer to rings around the centre.
        """
        # Check the validity of the indicated ring index
        ensure_not_zero(
            ring_index,
            "It is not possible to add a ring of cells at the given 0 index"
        )
        # Include the given layout at the end of the layer sublist specified
        # by the indicated index, if any; otherwise, either create a new
        # sublist or raise an exception
        layer_index = self._compute_layer_index(layer_index)
        # Get the cell's dimensions
        cell_side = cell.shape.dimensions[0]
        cell_apothem = cell.shape.dimensions[1]

        # Evaluate the multiplication factor for determining the construction
        # figure where the centers of the cells will be placed
        n = 2*ring_index
        # Build the construction figure
        construction_fig = self._build_ring_construction_figure(
            cell_apothem, n
        )
        # Compute the positions of the centres of the cells along the
        # indicated ring
        centres = compute_subdivision_points_on_borders(
            ring_index, construction_fig
        )

        # If needed, update the cell's tree before adding the cell
        cell.update_hierarchical_structure()
        # Add the cell at the given centre positions
        self.layers[layer_index].extend(
            get_cell_at_centres(
                cell,
                centres,
                self.rot_angle,
                wrap_shape(build_z_axis_from_vertex(self.o))
            )
        )
        # Update the characteristic dimensions of the lattice and its state
        self.dimensions = (
            self.dimensions[0] + (ring_index + 1) * cell_side,
            self.dimensions[1] + (ring_index + 1) * cell_apothem
        )
        self.state.is_update_needed = True

    def add_rings_of_cells(
            self,
            cell: Cell,
            no_rings: int,
            ring_index: int,
            layer_index: int | None = None
        ) -> None:
        """
        Method that adds several rings of ``Cell`` objects to the lattice at
        the positions calculated for each ring indicated by the total number
        of rings, starting from the `ring_index` parameter.
        Optionally, a ``layer_index`` can be provided to specify the layer to
        which the ring of cells is added; if not provided, the cells are
        added to a new layer.

        Parameters
        ----------
        cell : Cell
            The ``Cell`` instance to be iteratively added in order to build
            the rings of cells.
        no_rings : int
            The number of rings to add starting from the indicated ring index.
        ring_index : int
            The index indicating the starting ring index.
        layer_index : int | None = None
            The index of the layer to which the cells are added. If ``None``,
            a new layer is created.

        Raises
        ------
        ValueError
            If indicating a value less than ``1`` for the total number of
            rings to add.
        ValueError
            If indicating ``0`` as the ring index from which the rings of
            cells should be added.
        ValueError
            If indicating a layer index not falling within the valid range
            of existing layers.

        Notes
        -----
        Index 0 refers to the central cell of the lattice; indices greater
        than 0 refer to rings around the centre.
        """
        # Raise an exception if the number of rings to add is less than 1
        if no_rings < 1:
            raise ValueError(
                f"Wrong number ({no_rings}) of rings of cells to add has "
                "been indicated."
            )
        # Check the validity of the indicated ring index
        ensure_not_zero(
            ring_index,
            "It is not possible to add a ring of cells at the given 0 index"
        )
        # Include the layouts at the end of the layer sublist specified
        # by the indicated index, if any; otherwise, either create a new
        # sublist or raise an exception
        layer_index = self._compute_layer_index(layer_index)
        # Loop through the ring indices starting from the indicated one
        n0 = ring_index
        for i_ring in range(n0, n0+no_rings):
            # Add a ring of cells at the current index
            self.add_ring_of_cells(cell, i_ring, layer_index)
        # Set the need to update the lattice geometry
        self.state.is_update_needed = True

    def _build_ring_construction_figure(
            self, apothem: float, n: int
        ) -> Rectangle:
        """
        Method that builds and return a ``Hexagon`` object representing
        the construction figure on whose borders the cells belonging to the
        same ring will be placed. The dimensions of the hexagon depends on
        the given apothem and the scaling factor.
        The resulting figure is built with centre in the lattice centre, and,
        if needed, is rotated accordingly with the rotation angle of the
        lattice.

        Parameters
        ----------
        apothem : float
            The apothem of the ``Hexagon`` object.
        n : int
            The scaling factor for the dimensions of the ``Hexagon`` object.

        Returns
        -------
        Rectangle
            The ``Rectangle`` object built accordingly with the given
            dimensions.
        """
        construction_fig = Hexagon(
            get_point_coordinates(self.o), n*apothem
        )
        # Rotate the construction figure, if needed
        construction_fig.rotate(self.rot_angle)
        return construction_fig


# -------------------------------------------------------------------------- #
#                                FUNCTIONS                                   #
# -------------------------------------------------------------------------- #

def compute_subdivision_points_on_borders(
        no_vrtcs: int, surface: Surface
    ) -> List[Tuple[float, float, float]]:
    """
    Function that loops through the borders of the given ``Surface`` object
    and builds a list of evenly spaced vertices for each.

    Parameters
    ----------
    no_vrtcs : int
        The number of evenly spaced vertices to subdivide the edge into.
    surface : Surface
        The geometric surface on which borders vertices are built.

    Returns
    -------
    List[Tuple[float, float, float]]
        The list of XYZ coordinates for the vertices built on the borders of
        the given ``Surface`` object.
    """
    centres = []
    for border in surface.borders:
        centres.extend(
            [
                get_point_coordinates(p)
                for p in build_subdvision_vertices_on_edge(no_vrtcs, border)
            ]
        )
    return centres


def ensure_not_zero(value: int, message="Value must not be zero.") -> None:
    """
    Function that raises a ``ValueError`` exception if the given value is
    zero.

    Parameters
    ----------
    value : int
        The value to validate.
    message : str
        Error message to raise if ``value`` is zero.

    Raises
    ------
    ValueError
        If ``value`` is equal to zero.
    """
    if value == 0:
        raise ValueError(message)


def get_cell_at_centres(
        cell: Cell,
        centres: List[Tuple[float, float, float]],
        rot_angle: float,
        rot_axis: Edge
    ) -> List[Cell]:
    """
    Function that builds a list of the same given ``Cell`` instance where
    every element is placed according to the XYZ coordinates indicated by
    the elements of the input list.
    The ``Cell`` object is first cloned and rotated, then for each point
    in the list, the copied instance is cloned, translated and stored in the
    returned list.

    Parameters
    ----------
    cell : Cell
        The ``Cell`` instance to position according to the provided centres.
    centres : List[Tuple[float, float, float]]
        The list of the XYZ coordinates of the centres each cell should be
        placed at.
    rot_angle : float
        The rotation angle by which each cell should be rotated.
    rot_axis : Edge
        The rotation axis.

    Returns
    -------
    List[Cell]
        A list of ``Cell`` instances each positioned according to the
        provided centres.
    """
    # Initialise the list of cells to return
    cells = []
    # Clone and rotate the cell so that all the copies shares the same
    # rotation
    cell = cell.clone()
    cell.rotate(rot_angle, rot_axis)
    for centre in centres:
        # Clone, translate the cell in the given centre and append it
        # to the list to return
        cell = cell.clone()
        cell.translate(centre)
        cells.append(cell)
    return cells
