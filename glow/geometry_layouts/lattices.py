"""
Module containing classes providing the means for creating lattices from base
cells types in SALOME.
"""
import math

from copy import deepcopy
from typing import Any, Dict, List, Set, Tuple, Union

from glow.geometry_layouts.cells import Cell, HexCell, RectCell, Region, \
    check_cell_circle_are_cut, get_region_info
from glow.geometry_layouts.geometries import Hexagon, Rectangle, Surface, \
    build_hexagon
from glow.support.utility import are_same_shapes, build_compound_borders, \
    compute_point_by_reference, generate_unique_random_colors, \
    retrieve_selected_object, translate_wrt_reference
from glow.interface.geom_interface import ShapeType, add_to_study, \
    add_to_study_in_father, clear_view, display_shape, \
    extract_sorted_sub_shapes, extract_sub_shapes, get_basic_properties, \
    get_bounding_box, get_min_distance, get_object_from_id, \
    get_point_coordinates, is_point_inside_shape, make_cdg, make_common, \
    make_compound, make_cut, make_edge, make_face, make_partition, \
    make_rotation, make_translation, make_vector_from_points, \
    make_vertex, make_vertex_inside_face, make_vertex_on_curve, \
    remove_from_study, set_color_face, update_salome_study
from glow.support.types import CELL_VS_SYMM_VS_TYP_GEO, GeometryType, \
    LatticeGeometryType, PropertyType, SymmetryType, CellType


class Lattice():
    """
    Class that represents a lattice made by a group of cells.

    Parameters
    ----------
    cells : List[Cell]
            The list of cells that constitute the lattice, as objects
            of the 'Cell' subclasses
    name  : str = "Lattice"
            The lattice name in the current SALOME study

    Attributes
    ----------
    lattice_cells     : List[Cell]
                        The list of cells that constitute the lattice, as
                        objects of the 'Cell' subclasses
    lattice_cmpd      : Any
                        A GEOM compound object grouping all the faces of the
                        cells in the lattice
    lattice_entry_id  : Union[str, None]
                        An ID associated to the lattice surface in the current
                        SALOME study
    name              : str
                        The lattice name in the current SALOME study
    rings_no          : int
                        Index providing the total number of cells rings around
                        the lattice center
    distance          : float
                        Distance from the lattice origin to the outmost cell
                        CDG
    """
    # Admitted values for the rotation angle (in degrees) of the cells in
    # the lattice
    VALID_CELLS_ANGLES: List[float] = [0.0, 90.0, 180.0, 270.0, 360.0]

    def __init__(self,
                 cells: List[Cell] = [],
                 name: str = "Lattice",
                 center: Union[Tuple[float, float, float], None] = None,
                 boxes_thick: List[float] = []) -> None:
        # -------------------------
        # Attributes initialization
        # -------------------------
        self.lattice_center: Any = self.__evaluate_lattice_center(cells,
                                                                  center)
        self.cells_type: CellType | None = None
        self.__cells_rot: float | None = None
        self.lattice_cells: List[Cell] = deepcopy(cells)
        self.layers: List[List[Cell]] = [deepcopy(self.lattice_cells)]
        self.name: str = name
        self.lattice_entry_id: Union[str, None] = None
        self.rings_no: int = 0
        self.distance: float = 0.0
        self.__type_geo: LatticeGeometryType = LatticeGeometryType.ISOTROPIC
        self.symmetry_type: SymmetryType = SymmetryType.FULL
        self.lx: float = 0.0
        self.ly: float = 0.0
        self.box_layers: List[float] = boxes_thick
        self.__lattice_box: Union[Cell, None] = None
        self.lattice_symm: Union[Any, None] = None
        self.regions: List[Region] = []
        self.displayed_geom: GeometryType = GeometryType.TECHNOLOGICAL
        self.is_update_needed: bool = False

        # Update the instance attributes if any cells are provided
        self.__update_attributes(cells)

    @property
    def cells_rot(self) -> float | None:
        """
        Get or set the common rotation angle, in degrees, of the cells
        belonging to the main pattern of cells in the lattice.

        Parameters
        ----------
        cells_rot : float
            The rotation angle, in degrees, of the cells belonging to the
            main pattern of cells in the lattice.

        Returns
        -------
        float | None
            The rotation angle, in degrees, of the cells belonging to the
            main pattern of cells.
        """
        return self.__cells_rot

    def __get_main_pattern_rotation(self, new_val: float | None) -> float:
        """
        Method that determines the rotation angle (in degrees) of the cells
        representing the main pattern in the lattice, i.e. the most common
        value among the cells.

        Parameters
        ----------
        new_val : float
            An additional value of rotation angle to update the collection
            with.

        Returns
        -------
        float
            The rotation angle (in degrees) of the cells representing the
            main pattern in the lattice.
        """
        # Associate to each different rotation angle its frequency
        val_no: Dict[float, int] = {}
        for cell in self.lattice_cells:
            rot = math.degrees(cell.rotation)
            if rot in val_no:
                val_no[rot] += 1
            else:
                val_no[rot] = 1
        # Eventually update the dictionary with the given value
        if new_val != None and new_val in val_no:
            val_no[new_val] += 1
        elif new_val != None and new_val not in val_no:
            val_no[new_val] = 1
        # Return the rotation angle value of the main pattern of cells
        return max(val_no, key=val_no.get)

    @property
    def type_geo(self) -> LatticeGeometryType | None:
        """
        Get or set the lattice type of geometry as item of the enumeration
        `LatticeGeometryType`.

        Parameters
        ----------
        type_geo : LatticeGeometryType
            Item of the enumeration `LatticeGeometryType` to set the lattice
            type of geometry to

        Returns
        -------
        LatticeGeometryType | None
            Item of the enumeration `LatticeGeometryType` if get,
            `None` if set.

        Raises
        ------
        RuntimeError
            Invalid type of geometry for the current lattice geometry layout.
        """
        return self.__type_geo

    @type_geo.setter
    def type_geo(self, type_geo: LatticeGeometryType) -> None:
        self.set_type_geo(type_geo)

    @property
    def lattice_box(self) -> Cell | None:
        """
        Get or set the `Cell` subclass object providing the lattice's box.

        Parameters
        ----------
        cell : Cell | None
            An instance of the `Cell` subclass to set the lattice's box. If
            `None`, the previous box, if present, is removed.

        Returns
        -------
        Cell | None
            Instance the `Cell` subclass for the lattice's box, or `None` if
            get, `None` if set.

        Raises
        ------
        RuntimeError
            If the given cell's center does not coincide with the lattice's
            one.
        """
        return self.__lattice_box

    @lattice_box.setter
    def lattice_box(self, cell: Cell | None) -> None:
        # Check whether the cell's center coincide with the lattice's one
        if (cell is not None and
            get_min_distance(cell.figure.o, self.lattice_center) > 1e-5):
            raise RuntimeError(
                "The cell's center having coordinates "
                f"'{get_point_coordinates(cell.figure.o)}', "
                "does not coincide with the lattice's one being "
                f"'{get_point_coordinates(self.lattice_center)}'.")
        # Shape to update the lattice's characteristic dimensions with
        if cell is None:
            shape = make_compound([cell.face for cell in self.lattice_cells])
            # Set the lattice's box to None
            self.__lattice_box = None
        else:
            # Set the lattice's box
            self.__lattice_box = cell
            shape = cell.face
            # Update the lattice compounds with the box
            self.__update_lattice_compounds_with_box()

        # Update the characteristic dimensions of the lattice
        self.__update_lattice_dimensions(shape)
        # Set the need to update the lattice geometry
        self.is_update_needed = True

    def __update_lattice_dimensions(self, shape: Any) -> None:
        """
        Method that updates the lattice's XY characteristic dimensions
        from the given shape for which its min-max extension is calculated.

        Parameters
        ----------
        shape : Any
            The reference shape for calculating the lattice's characteristic
            dimensions.
        """
        # Get the min-max extension for the given shape
        x_min, x_max, y_min, y_max = get_bounding_box(shape)
        # Factor to evaluate the lattice's characteristic dimensions depending
        # on the cells' type
        n_type = 1
        if self.cells_type == CellType.HEX:
            n_type = 2
        self.lx = (x_max - x_min) / n_type
        self.ly = (y_max - y_min) / n_type

    def __evaluate_cells_rotation(self, cells: List[Cell]) -> float:
        """
        Method that checks if the rotation angle of the given cells is the
        same for all the given ones; if not, an exception is raised.
        A check on the validity of the cell rotation according to the admitted
        values is also present; an exception is raised if the value is not
        valid.
        The method returs the rotation angle of the cells, if all checks pass.

        Parameters
        ----------
        cells : List[Cell]
                The list of 'Cell' subclasses, representing the cells
                whose rotation angle must be checked

        Returns
        -------
        The value of the rotation angle (in degrees) shared by all the cells.
        """
        # Get the rotation angle of the first cell in the list
        cell_rot = math.degrees(cells[0].rotation)
        # If any of the cells has a different angle, an exception is raised
        if not all(
            [math.isclose(math.degrees(cell.rotation), cell_rot)
                for cell in cells]):
            raise RuntimeError("The given cells do not share the same "
                               "rotation angle.")
        # Check if the rotation angle of the cells is valid
        self.__check_cell_rotation_validity(cell_rot)
        # Return the rotation angle of all the cells
        return cell_rot

    def __check_cell_rotation_validity(self, cell_rot: float) -> None:
        """
        Method that checks if the provided cell rotation angle is valid.

        Parameters
        ----------
        cell_rot : float
            The rotation angle of the cell to validate.

        Raises
        ------
        RuntimeError
            If the provided cell rotation angle is not among the valid ones.
        """
        if not cell_rot in self.VALID_CELLS_ANGLES:
            raise RuntimeError("Cells can only have any of the "
                               f"{self.VALID_CELLS_ANGLES}° rotation angles")

    def __evaluate_lattice_center(
            self,
            cells: Union[List[Cell], None],
            center: Union[Tuple[float, float, float], None]) -> Any:
        """
        Method that evaluates the lattice center during the initialization:
        if no center is provided, the XYZ origin is returned. Otherwise, a
        check is performed to identify if there is a cell whose center
        coincides with the specified one.
        In case none is found or the center lays within any of the cells'
        faces, an exception is raised.

        Parameters
        ----------
        cells : List[Cell] | None
            A list of `Cell` objects constituting the lattice
        center : Tuple[float, float, float]
            The indicated X-Y-Z coordinates of the lattice center

        Returns
        -------
        The vertex object identifying the lattice center.

        Raises
        ------
        RuntimeError
            In case the specified center lays inside any of the cells.
        """
        if center:
            # Vertex object for the center's coordinates
            center_vrtx = make_vertex(center)
            if not cells:
                return center_vrtx
            # Check if the given lattice center corresponds to any of the
            # given cells centers
            for cell in cells:
                if get_min_distance(cell.figure.o, center_vrtx) < 1e-7:
                    return center_vrtx
                # Raise an exception if the center is within the cell
                if is_point_inside_shape(center_vrtx, cell.face) and not (
                   is_point_inside_shape(center_vrtx,
                                         make_compound(cell.figure.borders))):
                    raise RuntimeError(
                        f"The specified center {center} lays inside one "
                        "of the provided cells without coinciding with its "
                        "center or its borders.")
            else:
                return center_vrtx
        else:
            # No center has been specified: return the XYZ origin
            return make_vertex((0, 0, 0))

    def __build_lattice(
            self, geo_type: GeometryType = GeometryType.TECHNOLOGICAL) -> None:
        """
        Method that builds the lattice out of its cells according to the given
        type of geometry, i.e. the technological or the refined ones for the
        calculations. According to this, we have:
        - the GEOM compound object storing either the cells faces or those
          resulting from the face sectorization;
        - the GEOM compound object storing only the unique edges of all the
          cells comprising the lattice;

        In addition, if the lattice has a box, this container is rebuilt by
        means of the thicknesses of each layer stored as an instance
        attribute.
        """
        # Build the compound object out of the cells faces according to the
        # given type (either from the technological geometry or from the
        # sectorization)
        self.lattice_cmpd = get_compound_from_geometry(geo_type,
                                                       self.lattice_cells)
        # Handle the construction of the lattice container, if any is required
        self.__assemble_box(geo_type)

    def __assemble_box_with_lattice(
            self,
            inner_box: Any,
            lattice_cmpd: Any,
            geo_type: GeometryType) -> Any:
        """
        Method that assembles the lattice, in terms of its compound, with the
        container.
        Given the face built on the box inner borders, a common operation
        with the lattice compound is performed to take into account also for
        boxes being slighlty lesser than the lattice itself.
        The result is a lattice cut in correspondence of the intersection
        points.
        A partition operation between the lattice compound and the whole box
        face follows to assemble both objects into a whole new compound.

        Parameters
        ----------
        inner_box     : Any
            The face object representing the inner part of the lattice box
        lattice_cmpd  : Any
            The compound object representing the lattice the box should be
            assembled with
        geo_type : GeometryType
            The type of geometry indicating the layout to use for the box cell

        Returns
        -------
        The lattice compound assembled with its box.
        """
        # Extract the common part between the lattice compound and the inner
        # box.
        lattice_cmpd = make_common(inner_box, lattice_cmpd)
        # Update the lattice compound by partitioning it with the box face
        # and return it
        return make_partition(
            [lattice_cmpd] + extract_sub_shapes(
                get_compound_from_geometry(geo_type, [self.__lattice_box]),
                ShapeType.FACE),
            [],
            ShapeType.FACE)

    def build_lattice_box(self, box_thick: List[float] = []) -> None:
        """
        Method for building the box the lattice is inserted into, given the
        thickness of each layer.
        On the basis of the cells type (i.e. cartesian or hexagonal), a proper
        box is built and used to update the compound of the lattice.
        Since the method builds a new box for the lattice, any property
        associated to the regions of a previuosly set box are lost.

        N.B. The method accepts that the first value of layers thicknesses
        could be negative, meaning that the box is lesser than the lattice.
        The other values must be positive.

        Parameters
        ----------
        box_thick : List[float]
                    List storing the thickness of each layer of the lattice
                    box
        """
        # Return immediately if the method is called without specifying the
        # thicknesses of the box layers
        if not box_thick:
            print("WARNING! Method 'build_lattice_box' called without "
                  "specifying the thicknesses of the box layers")
            return
        if box_thick[0] < 0 and len(box_thick) > 1:
            for t in box_thick[1:-1]:
                if t < 0.0:
                    raise AssertionError("Only the first value of the box "
                                         "layers thicknesses can be negative.")

        # Store the thicknesses of the box layers: if the first is positive,
        # a zero value is put at posizion 0 in the list
        self.box_layers = list(box_thick)
        if box_thick[0] > 0:
            self.box_layers.insert(0, 0.0)

        # Build the lattice container according to the cells geometry
        self.__build_lattice_box_type()
        # Update the lattice compounds with the box
        self.__update_lattice_compounds_with_box()
        # Set the need to update the lattice geometry
        self.is_update_needed = True

    def __update_lattice_compounds_with_box(self) -> None:
        """
        Method that updates the lattice's compound object by assembling the
        cells with the inner area of the box so to take into account for
        boxes that cut the last ring of cells.
        In addition, if a symmetry has already been applied, the operation
        is performed again so to update the compound specific for the
        symmetry.
        """
        # Update the lattice compound from the cells' faces and the inner box
        self.lattice_cmpd = self.__assemble_box_with_lattice(
            self.__extract_inner_box(),
            make_compound([cell.face for cell in self.lattice_cells]),
            GeometryType.TECHNOLOGICAL)
        # Re-apply the symmetry operation if any has already been performed
        if self.symmetry_type != SymmetryType.FULL:
            self.apply_symmetry(self.symmetry_type)

    def __extract_inner_box(self) -> Any:
        """
        Method that extracts the inmost face object of the cell representing
        the lattice container.
        The box faces are extracted and then the one closest to the lattice
        center is returned.

        Returns
        -------
        The face object for the area of the lattice box closest to the lattice
        center.
        """
        # Extract the lattice box faces
        box_faces = extract_sub_shapes(self.__lattice_box.face,
                                       ShapeType.FACE)
        if not box_faces:
            box_faces = [self.__lattice_box.face]
        # Get the inner box face by sorting the faces according to the
        # distance from the lattice center
        return sorted(
            box_faces,
            key=lambda item: get_min_distance(self.lattice_center, item))[0]

    def __build_lattice_box_type(self) -> None:
        """
        Method that builds the lattice box as an instance of the `Cell` class.
        The container geometry is built accordingly with the type of geometry
        of the cells in the lattice.
        Either a rectangle or a hexagon is built for each layer of the box and
        a partition with all the figures that make up the box is performed.
        The lattice box is instantiated either as a `RectCell` or a `HexCell`
        and its face updated with the figure built herein.
        """
        # Declare a list storing the geometrical figures that constitute
        # the box
        box_surfaces: List[Surface] = []
        center = get_point_coordinates(self.lattice_center)
        # Get the bounding box of the lattice
        x_min, x_max, y_min, y_max = get_bounding_box(
            make_compound([cell.face for cell in self.lattice_cells]))
        # Build the container for the lattice according to the cells geometry
        # type
        match self.cells_type:
            case CellType.RECT:
                # Declare the starting dimensions of the box
                height = y_max - y_min
                width = x_max - x_min
                # Build a rectangle for each layer
                for thick in self.box_layers:
                    height += 2*thick
                    width += 2*thick
                    box_surfaces.append(Rectangle(center, height, width))
                # Perform a 'partition' operation to assemble the lattice box
                box_face = make_partition(
                    [rect.face for rect in box_surfaces], [], ShapeType.FACE)
                # Declare a 'Rectangle' instace from the built face
                self.__lattice_box = RectCell(center, (height, width))
            case CellType.HEX:
                # Calculate the apothem of the hexagon enclosing the lattice
                # according to the valid rotation angles of the cells (0° or
                # 90°)
                if math.isclose(self.__cells_rot, 0.0):
                    box_apothem = 0.5 * (x_max - x_min)
                else:
                    box_apothem = 0.5 * (y_max - y_min)
                # Build a hexagon for each layer
                for thick in self.box_layers:
                    box_apothem += thick
                    box_surfaces.append(
                        build_hexagon(box_apothem, center))
                # Perform a 'partition' operation to assemble the lattice box
                box_face = make_partition(
                    [hex.face for hex in box_surfaces], [], ShapeType.FACE)
                # Declare a 'Hexagon' instace from the built face
                self.__lattice_box = HexCell(center,
                                             box_apothem/math.sin(math.pi/3))
            case _:
                raise RuntimeError("Unhandled cell geometry type.")
        self.__lattice_box.update_geometry_from_face(
            GeometryType.TECHNOLOGICAL, box_face)
        # Rotate the box so to enclose the lattice (hexagonal lattice only)
        if self.__cells_rot == 0.0 and self.cells_type == CellType.HEX:
            self.__lattice_box.rotate(90)
        # Update the characteristic dimensions of the lattice
        self.__update_lattice_dimensions(self.__lattice_box.face)

    def __configure_lattice_type(self,
                                 cells_type: CellType,
                                 no_cells: int = 1) -> None:
        """
        Method that allows to configure the type of geometry of the lattice,
        as value of the 'LatticeGeometryType' enumeration, according to the
        type of the cells in the lattice and the number of cells that are
        present.
        The following convention is adopted, according to what required by
        the SALT module of DRAGON5:
        - lattice made of rectangular cells: the type of geometry depends on
          the number of cells:
            - only one cell - the type is RECTANGLE_TRAN
            - more than one cell - the type is RECTANGLE_SYM
        - lattice made of hexagonal cells: the type of geometry depends on the
          number of cells:
            - only one cell - the type is HEXAGON_TRAN
            - more than one cell - the type is ISOTROPIC

        The type of geometry is expressed as a value of the 'GeometryType'
        enumeration.
        Given the type of geometry, the corresponding type of boundary is
        set as well.

        Parameters
        ----------
        cells_type  : CellType
                      The value of the 'CellType' enumeration identifying the
                      type of cells in the lattice
        no_cells    : int
                      The number of cells in the lattice
        """
        print("The type of cells is", cells_type)
        match cells_type:
            case CellType.RECT:
                if no_cells > 1:
                    # Case of a lattice made of cartesian cells: this case is
                    # characterised by an isotropic geometry type with VOID
                    # BCs.
                    self.__type_geo = LatticeGeometryType.ISOTROPIC
                else:
                    # Only one cartesian cell is present: this case is
                    # characterised by an cartesian geometry type with
                    # translation on all sides
                    self.__type_geo = LatticeGeometryType.RECTANGLE_TRAN
            case CellType.HEX:
                if no_cells > 1:
                    # Case of a lattice made of hexagonal cells: this case is
                    # characterised by an isotropic geometry type with VOID
                    # BCs.
                    self.__type_geo = LatticeGeometryType.ISOTROPIC
                else:
                    # Only one hexagonal cell is present: this case is
                    # characterised by an hexagonal geometry type with
                    # translation on all sides
                    self.__type_geo = LatticeGeometryType.HEXAGON_TRAN

    def add_cell(
            self,
            cell: Cell,
            position: Tuple[float, float, float]) -> None:
        """
        Method that allows to add a new cell to the lattice at the specified
        position. The cell is added to a new layer.

        Parameters
        ----------
        cell : Cell
            A cell to add to the current lattice, as object of the 'Cell'
            subclasses.
        position : Tuple[float, float, float]
            The X-Y-Z coordinates of the position where the cell should be
            added, i.e. the cell center should be placed at those coordinates.
            If the tuple is empty, the cell is placed at the position
            indicated by its center.

        Raises
        ------
        RuntimeError
          - If the cell to add has a `CellType` which differs from the one
            of the cells in the lattice.
          - If the cell to add has an invalid rotation angle.
        """
        try:
            # Check that the cell to add has the same type of the ones already
            # present
            self.__check_cell_type(cell.cell_type)
            # Check the validity of the cell's rotation angle and update the
            # value representative of the main pattern of cells
            self.__set_main_rotation_angle(math.degrees(cell.rotation))
        except RuntimeError as e:
            raise RuntimeError(f"Error when adding a cells: {e}")
        # Initialize a new sub list identifying a new layer of cells
        self.layers.append([])
        # Add the cell to the newly created layer
        self.__add_cell_to_layer(cell, position, len(self.layers)-1)
        # Set the need to update the lattice geometry
        self.is_update_needed = True

    def __set_main_rotation_angle(self, cell_rot: float) -> None:
        """
        Method that sets the rotation angle of the main pattern of cells
        in the lattice based on the provided cell rotation.
        The rotation angle of the given cell is assessed first and the
        rotation angle of the main pattern of cells re-evaluated by taking
        into account also for the given value.

        Parameters
        ----------
        cell_rot : float
            The rotation angle, in degrees, to be evaluated.
        """
        # Check the validity of the cell's rotation angle
        self.__check_cell_rotation_validity(cell_rot)
        # Update the rotation angle of the main pattern of cells
        self.__cells_rot = self.__get_main_pattern_rotation(cell_rot)

    def __add_cell_to_layer(
            self, cell: Cell, position: Tuple[float], layer_indx: int) -> None:
        """
        Method that adds a given cell, as subclass of the 'Cell' class, to
        the lattice. The cell is added to the layer whose index is provided
        as input.

        Parameters
        ----------
        cell : Cell
            A cell to add to the current lattice, as object of the 'Cell'
            subclasses
        position : Tuple[float, float, float]
            The X-Y-Z coordinates of the position where the cell should be
            added, i.e. the cell center should be placed at those coordinates.
            If the tuple is empty, the cell is placed at the position
            indicated by its center
        """
        # Move the cell in the given position, if necessary and re-evaluate
        # the number of cells rings in the lattice
        print("Adding the cell to position", position)
        if position and not all(math.isclose(x, 0.0) for x in position):
            cell = cell.translate(position)
            self.__evaluate_no_rings(position)

        # Set the cell's name by appending a global index
        cell.name = f"Cell_{len(self.lattice_cells)+1}"
        # Add the cell to the given layer
        self.layers[layer_indx].append(deepcopy(cell))
        # Add the cell to the list of lattice cells
        self.lattice_cells.append(deepcopy(cell))
        # Re-evaluate the types of the geometry and BCs
        self.__configure_lattice_type(self.cells_type,
                                      len(self.lattice_cells))
        # Log
        print("Number of lattice cells is", len(self.lattice_cells))
        print("Number of lattice rings is", self.rings_no)

    def __evaluate_no_rings(
            self, position: Tuple[float, float, float]) -> None:
        """
        Method that evaluates the number of cells rings in the lattice, given
        the position of a cell, as the X-Y-Z coordinates of its center.

        If the distance of the cell from the lattice center is greater that
        the one stored in the corresponding attribute, the value is updated.
        The number of cells rings is updated as well only if this distance is
        greater than the one for the current number of rings: this value
        depends on the cells rotation.

        Parameters
        ----------
        position  : Tuple[float, float, float]
                    The X-Y-Z coordinates of the center of the cell used for
                    evaluating if the number of rings should be updated
        """
        origin_distance = math.dist(get_point_coordinates(self.lattice_center),
                                    position)
        if origin_distance - self.distance > 1e-7: # TODO assign a constant EPSILON
            # Update the distance to the outmost cell
            self.distance = origin_distance
            # Update the number of rings in the lattice only if the new
            # distance is greater than the maximum one for the current ring
            # number; this distance depends on the cells rotation.
            if math.isclose(self.__cells_rot, 0.0):
                # The distance is calculated using the cells characteristic
                # dimension along the X-axis
                ring_dist = 2*self.lx * self.rings_no
            elif math.isclose(self.__cells_rot, 90.0):
                # The distance is calculated using the cells characteristic
                # dimension along the Y-axis
                ring_dist = 2*self.ly * self.rings_no
            # Check if the rings number should be updated
            if self.distance > ring_dist + 1e-7: # TODO assign a constant EPSILON
                self.rings_no += 1

    def apply_symmetry(self, symmetry: SymmetryType) -> None:
        """
        Method that modifies the lattice in order to apply the given symmetry
        type. The symmetry types provided by the 'SymmetryType' enumeration
        are handled.
        According to the type, geometry transformations (i.e. cuts) by means
        of the SALOME GEOM module operations are performed on the lattice face.

        Parameters
        ----------
        symmetry  : SymmetryType
                    The type of symmetry to handle
        """
        # Restore the whole lattice if the 'FULL' symmetry type is provided
        if symmetry == SymmetryType.FULL:
            self.symmetry_type = SymmetryType.FULL
            # Re-evaluate the BC type according to the lattice
            self.__configure_lattice_type(self.cells_type,
                                          len(self.lattice_cells))
            # Update the lattice in the current SALOME study
            self.show()
            self.lattice_symm = self.lattice_cmpd
            return

        # Assemble all the lattice layers to build the lattice compounds
        self.__update_lattice_compounds(self.__assemble_layers())
        # Handle the symmetry operation according to the type of cells of
        # the lattice
        match self.cells_type:
            case CellType.RECT:
                self.__apply_rect_symmetry(symmetry)
            case CellType.HEX:
                self.__apply_hex_symmetry(symmetry)
            case _:
                raise AssertionError(
                    f"Unrecognized type '{self.cells_type}' for the lattice "
                    "cells.")
        # Update the lattice symmetry type
        self.symmetry_type = symmetry

        # Set the need to update the lattice geometry
        self.is_update_needed = True
        # Update the lattice in the current SALOME study
        self.show()

    def __apply_rect_symmetry(self, symmetry: SymmetryType):
        """
        Method that modifies the lattice, made of cartesian cells, to apply
        the given symmetry type. Only the symmetry types, provided by the
        'SymmetryType' enumeration, valid for rectangular shapes are handled.

        Parameters
        ----------
        symmetry  : SymmetryType
                    The type of symmetry to handle
        """
        # Get the lattice bounding box
        b_box = get_bounding_box(self.lattice_cmpd)
        # Handle the construction of the face object from which to derive the
        # lattice symmetry
        match symmetry:
            case SymmetryType.HALF:
                # ------------------
                # HALF symmetry case
                # ------------------
                # Build the face object for extracting the lattice symmetry
                face = self.__handle_rect_symmetry(b_box, 0)
                 # Update the lattice type of geometry
                self.__type_geo = LatticeGeometryType.RECTANGLE_SYM
            case SymmetryType.QUARTER:
                # ---------------------
                # QUARTER symmetry case
                # ---------------------
                # Build the face object for extracting the lattice symmetry
                face = self.__handle_rect_symmetry(b_box, 1)
                # Update the lattice type of geometry
                self.__type_geo = LatticeGeometryType.RECTANGLE_SYM
            case SymmetryType.EIGHTH:
                # --------------------
                # EIGHTH symmetry case
                # --------------------
                # Build the face object of the triangle representing the
                # symmetry case
                face = self.__build_face_from_vertices([
                    self.lattice_center,
                    make_vertex((
                        b_box[1],
                        get_point_coordinates(self.lattice_center)[1],
                        0.0)),
                    make_vertex((b_box[1], b_box[3], 0.0))])
                # Update the lattice type of geometry
                self.__type_geo = LatticeGeometryType.RECTANGLE_EIGHT
            case _:
                raise AssertionError(
                    f"The provided '{symmetry}' symmetry type is not "
                    "admitted for lattices with cartesian cells.")

        # Perform the 'common' operation on the lattice compound with the
        # face that identifies the symmetry to exctract the lattice portion
        # corresponding to the symmetry
        self.lattice_symm = make_common(self.lattice_cmpd, face)

    def __apply_hex_symmetry(self, symmetry: SymmetryType) -> None:
        """
        Method that modifies the lattice, made of hexagonal cells, to apply
        the given symmetry type. Only the symmetry types, provided by the
        'SymmetryType' enumeration, valid for hexagonal shapes are handled.

        Parameters
        ----------
        symmetry  : SymmetryType
                    The type of symmetry to handle
        """
        # Raise an exception if the lattice is not included within a box
        if not self.__lattice_box:
            raise AssertionError(
                "The hexagonal lattice is not included within a box: the "
                f"requested '{symmetry}' symmetry operation cannot be "
                "applied.")
        match symmetry:
            case SymmetryType.THIRD:
                # ------------------
                # R120 symmetry case
                # ------------------
                # Build the shape identifying the symmetry
                face = self.__handle_third_symmetry()
                # Update the lattice type of geometry
                self.__type_geo = LatticeGeometryType.R120
            case SymmetryType.SIXTH:
                # ------------------
                # SA60 symmetry case
                # ------------------
                # Build the shape identifying the symmetry
                face = self.__handle_sixth_symmetry()
                # Update the lattice type of geometry
                self.__type_geo = LatticeGeometryType.SA60
            case SymmetryType.TWELFTH:
                # -----------------
                # S30 symmetry case
                # -----------------
                # Build the shape identifying the symmetry
                face = self.__build_triangle_on_lattice(0.0, 1/12)
                # Update the lattice type of geometry
                self.__type_geo = LatticeGeometryType.S30
            case _:
                raise AssertionError(
                    f"The provided '{symmetry}' symmetry type is not "
                    "admitted for lattices with hexagonal cells.")
        # Perform the 'common' operation to extract the symmetry type
        # on the lattice and update the stored compound object
        self.lattice_symm = make_common(self.lattice_cmpd, face)

    def __build_vertices(self, u_params: List[float]) -> List[Any]:
        """
        Method that builds a list of vertices on the construction circle of
        the lattice box, with the first element being the lattice center.
        The position of the points along the curve is expressed in terms of
        a parameter in the [0-1] range, whose values are given as input.

        Parameters
        ----------
        u_params  : List[float]
            List of parameters in the [0-1] range identifying the positions
            of the points

        Returns
        -------
        A list of vertex objects built on the lattice box construction circle.
        """
        # Initialize the list with the lattice center
        points = [self.lattice_center]
        points += [
            make_vertex_on_curve(self.__lattice_box.figure.out_circle, u)
            for u in u_params]
        # Return the built list of vertex objects
        return points

    def __handle_sixth_symmetry(self) -> Any:
        """
        Method that handles the construction of the shape surface that
        identifies a sixth symmetry for lattices with hexagonal cells.
        The shape vertices are built depending on the rotation angle of
        the cells, allowing to handle both 0° and 90° cases.
        A face object is then built from the identified points.

        Returns
        -------
        A face object representing the sixth symmetry shape.
        """
        if self.__cells_rot == 0.0:
            # Build the triangle identifying the symmetry type
            return self.__build_triangle_on_lattice(1/12, 1 - 1/12)
        elif self.__cells_rot == 90.0:
            # Build the triangle identifying the symmetry type
            return self.__build_triangle_on_lattice(2/3, 5/6)
        else:
            raise AssertionError("The cells rotation of "
                                  f"{self.__cells_rot}° is not admitted.")

    def __handle_rect_symmetry(
            self, b_box: List[float], param: float) -> Any:
        """
        Method that handles the symmetry application for cartesian cells
        lattices. Given the X-Y min/max values of the bounding box of the
        lattice and a parameter that identifies either the 'HALF' or 'QUARTER'
        symmetry case, a rectangular face is built and returned.

        Parameters
        ----------
        b_box : List[float]
            List providing the X-Y min/max values of the bounding box of the
            lattice
        param : float
            Parameter for identifying either the 'HALF' or 'QUARTER' symmetry
            case

        Returns
        -------
        A face object representing the rectangle that identifies either the
        'HALF' or 'QUARTER' symmetry case.
        """
        # Get the X-Y coordinates of the lattice center
        o_x, o_y, _ = get_point_coordinates(self.lattice_center)
        # Get the symmetry rectangle center, depending on the symmetry case
        symm_center = (
            (b_box[1] - o_x) / 2 + o_x,
            param * ((b_box[3] - o_y) / 2) + o_y,
            0)
        # Build a 'Rectangle' object identifying the symmetry case
        rect = Rectangle(
            symm_center,
            (b_box[3] - b_box[2]) / (2 - (1 - param)/(1 + param)),
            (b_box[1] - b_box[0]) / 2)
        # Return the face
        return rect.face

    def __handle_third_symmetry(self) -> Any:
        """
        Method that handles the construction of the shape surface that
        identifies a third symmetry for lattices with hexagonal cells.
        The shape vertices are built depending on the rotation angle of
        the cells, allowing to handle both 0° and 90° cases.
        A face object is then built from the identified points.

        Returns
        -------
        A face object representing the third symmetry shape.
        """
        if self.__cells_rot == 0.0:
            points = self.__build_vertices([1/4, 7/12, 2/3])
        elif self.__cells_rot == 90.0:
            points = self.__build_vertices([0, 1-1/6, 1-1/3])
        else:
            raise AssertionError(
                f"The cells rotation of {self.__cells_rot}° is not admitted.")
        # Return a face built from the vertices
        return self.__build_face_from_vertices(points)

    def __build_triangle_on_lattice(
            self, vert2_u: float, vert3_u: float) -> Any:
        """
        Method that builds a triangle with a vertex in the lattice center.
        The other two vertices are positioned on the circle the lattice is
        inscribed into; their positions are determined by their U-parameters
        provided as parameters.
        Given the vertices, the edges and the corresponding face are built.

        Parameters
        ----------
        vert2_u : float
                  U-parameter for identifying the second triangle vertex
        vert3_u : float
                  U-parameter for identifying the third triangle vertex

        Returns
        -------
        The GEOM face object identifying the triangle
        """
        # Build the three vertices of the triangle
        points = self.__build_vertices([vert2_u, vert3_u])
        # Return the face object built on the vertices
        return self.__build_face_from_vertices(points)

    def __build_face_from_vertices(self, vertices: List[Any]) -> Any:
        """
        FIXME it can be transformed into a function
        Method that, given the point objects, builds and returns a face
        object.

        Parameters
        ----------
        vertices: List[Any]
            List of the vertices of the edges being the face borders

        Returns
        -------
        A face object built on the edges described by the given vertices.
        """
        # Build the list of edges from the given vertices
        edges = [make_edge(vertices[i],
                           vertices[(i+1) % len(vertices)])
                           for i in range(len(vertices))]
        # Return the face object built on the edges
        return make_face(edges)

    def __assemble_layers(self) -> List[Cell]:
        """
        Method that assembles all the lattice layers made by list of cells
        in a single compound. Each layer will cut all the layers below itself,
        if any overlapping occurs. Lastly, the lattice box, if present, is
        assembled with the whole compound.
        The method populate and returns a list of 'Cell' objects representing
        the cells currently present in the lattice after removal and cut
        operations due to layers overlapping.

        Returns
        -------
        A list of 'Cell' objects for all the cells currently present in the
        lattice.
        """
        # Store a copy of the layers in reverse order
        layers = self.layers[::-1]
        # Loop through all the layers to cut those that are overlapped
        for i, layer in enumerate(layers):
            # Extract a slice of the list of layers starting from the next
            # one, i.e. the layer below the current one and all the subsequent
            # ones
            sub_layers = layers[i + 1:]
            # Loop through all the layers below to cut the cells if the
            # current layer overlaps any inferior layer. Modifications to
            # the layers (cells cut) are kept as operating on shallow copies
            # of the two lists lists.
            for j, sub_layer in enumerate(sub_layers):
                self.__overlap_layer_to(layer, sub_layer)
        # Return the flattened list of list of cells
        return [cell for layer in layers for cell in layer]

    def __update_lattice_compounds(
            self,
            cells: List[Cell],
            geo_type: GeometryType = GeometryType.TECHNOLOGICAL) -> None:
        """
        Method that rebuilds the compound objects of the lattice grouping
        the given cells. According to the indicated geometry type, either
        the technological or the sectorized geometry layout of the cells
        is used.
        If any box is declared, the cells are assembled with it.

        Parameters
        ----------
        cells : List[Cell]
            List of 'Cell' objects to build a compound object from
        geo_type : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry indicating the cells' layout to use
        """
        self.lattice_cmpd = get_compound_from_geometry(geo_type, cells)

        # Handle the construction of the lattice container, if any
        self.__assemble_box(geo_type)

        # Update the compound object storing the applied symmetry
        if self.symmetry_type != SymmetryType.FULL:
            self.lattice_symm = make_common(self.lattice_cmpd,
                                            self.lattice_symm)

    def __assemble_box(self, geo_type: GeometryType) -> None:
        """
        Method that builds the lattice box, if it was not built yet, from the
        stored layers thicknesses. The container geometry depends on the
        type of geometry of the cells (i.e. rectangular or hexagonal).
        The box is then assembled with the lattice, eventually cutting the
        cells that are overlapped by it.

        Parameters
        ----------
        geo_type : GeometryType
            The type of geometry indicating the layout to use for the box cell
        """
        if not self.box_layers and not self.__lattice_box:
            return
        # If no lattice box has still been build, do it now
        if not self.__lattice_box:
            self.build_lattice_box(self.box_layers)
        else:
            # Assemble the lattice with the box subface closest to the lattice
            # center
            self.lattice_cmpd = self.__assemble_box_with_lattice(
                self.__extract_inner_box(), self.lattice_cmpd, geo_type)

    def __overlap_layer_to(
            self, layer: List[Cell], sub_layer: List[Cell]) -> None:
        """
        Method that overlaps a layer to an inferior one, both defined as
        lists of 'Cell' objects.
        If the two layers do not share a common part or their distance is
        greater than zero, the method returns without cutting any cell of the
        inferior layer. On the contrary, the common part is removed from the
        inferior layer by cutting its cells appropriately.
        In particular, cells that are completely overlapped are removed from
        their layer, whereas, for those partially overlapped, their geometry
        layout is updated with the corresponding cut version.

        Parameters
        ----------
        layer : List[Cell]
            The superior layer in terms of a list of 'Cell' objects
        sub_layer : List[Cell]
            The inferior layer in terms of a list of 'Cell' objects
        """
        # Return if any of the two layers do not have cells
        if len(layer) < 1 or len(sub_layer) < 1:
            return
        # Build a compound for each layer
        layer_cmpd = make_compound([cell.face for cell in layer])
        sub_layer_cmpd = make_compound([cell.face for cell in sub_layer])
        # Return immediately if the two layers do not overlap
        if get_min_distance(layer_cmpd, sub_layer_cmpd) > 0.0 or \
            not extract_sub_shapes(make_common(layer_cmpd, sub_layer_cmpd),
                                  ShapeType.FACE):
            return
        print(f"The current layer overlaps the inferior one.")

        # List storing the indices of the cells to remove, as completely
        # overlapped
        indices: List[int] = []
        for i, cell in enumerate(sub_layer):
            # Continue with the next cell if the common operation between
            # the cell's face and the superior layer does not return any
            # face, meaning there is no overlapping
            if not extract_sub_shapes(make_common(cell.face, layer_cmpd),
                                      ShapeType.FACE):
                continue
            # Cut the cell face with the layer and check if the result has
            # any face; if not, it means the cell is completely overlapped,
            # hence its index is stored to remove the corresponding cell
            cut_cell = make_cut(cell.face, layer_cmpd)
            if not extract_sub_shapes(make_compound([cut_cell]),
                                      ShapeType.FACE):
                indices.append(i)
                continue
            # Update the lattice cell face with the result of the cut
            # operation between the inferior layer cell face and the
            # superior layer
            print(f"Updating cell #{i} face ...")
            sub_layer[i].update_geometry_from_face(GeometryType.TECHNOLOGICAL,
                                                   cut_cell)
            # Update also the sectorized face, if any
            if cell.sectorized_face:
                sub_layer[i].update_geometry_from_face(
                    GeometryType.SECTORIZED,
                    make_cut(cell.sectorized_face, layer_cmpd))
        # Remove the cells completely overlapped by the superior layer
        for index in sorted(indices, reverse=True):
            sub_layer.pop(index)

    def show(self,
             property_type_to_show: Union[PropertyType, None] = None,
             geometry_type_to_show: GeometryType = GeometryType.TECHNOLOGICAL
             ) -> None:
        """
        Method that allows to show the lattice and its cell regions into
        the current SALOME study.
        If the whole lattice compound were already present in the study,
        it is fist removed, to keep only one entry at a time.
        It is then added to the study together with the edges compound,
        the lattice center point and all the lattice regions, colored
        according to the value of the property type provided as input to
        this method.
        If no property is selected, the regions are displayed without
        colouring.

        Parameters
        ----------
        property_type_to_show : Union[PropertyType, None] = None
            Either a value of the 'PropertyType' enumeration, indicating
            the property by which the regions are coloured, or None.
        geometry_type_to_show : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry to show: regions are shown accordingly
            with the selected type
        """
        # Erase all objects from the current view
        clear_view()
        # If already present in the current SALOME study, remove the compound
        # object
        if self.lattice_entry_id and get_object_from_id(self.lattice_entry_id):
            remove_from_study(self.lattice_entry_id)

        # Build regions and color them, if needed. If an exception is raised,
        # is caught and re-raised
        try:
            # Build the cell regions, as 'Region' objects, according to the
            # indicated geometry type. It also updates the lattice compound.
            self.build_regions(geometry_type_to_show)
            # Add the lattice compound to the current SALOME study as
            # representing the parent item for regions in the Object Browser
            self.lattice_entry_id = add_to_study(self.lattice_cmpd, self.name)
            # Assign the same color to all the regions having the same property
            # value, if any has been specified to show
            self.__associate_colors_to_regions(property_type_to_show)
            # Add all the lattice regions to the study, each with an assigned
            # color
            for region in self.regions:
                # Set the region color in the viewer
                set_color_face(region.face, region.color)
                # Add the lattice region to the study
                region.face_entry_id = add_to_study_in_father(
                    self.lattice_cmpd, region.face, region.name)
                if not region.face_entry_id:
                    raise RuntimeError(
                        f"Problem arose when adding the region {region.name}"
                        "to the SALOME study")
                # Display the region in the current view
                display_shape(region.face_entry_id)
        except:
            # Show only the whole lattice before raising the caught exception
            display_shape(self.lattice_entry_id)
            raise
        finally:
          # Show everything on the SALOME application
          update_salome_study()
          # Update the geometry type
          self.displayed_geom = geometry_type_to_show

    def __associate_colors_to_regions(
            self, property_type: Union[PropertyType, None]) -> None:
        """
        Method that assigns the same color to all the regions having the same
        value for the given property type.

        Parameters
        ----------
        property_type : PropertyType
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
        # each region
        values = self.__get_unique_values_for_property(property_type)
        # Generate a specific amount of colors as the number of different
        # values for the same given property type
        colors = generate_unique_random_colors(len(values))
        # Build a dictionary of values for the given property type VS color
        property_vs_color = dict(zip(list(values), colors))
        # Loop through all the regions and assign a color corresponding to
        # the value of the given property type
        for region in self.regions:
            # Get the value of the given property type associated to the
            # region
            value = region.properties[property_type]
            # Set the region color
            region.set_property_color(property_vs_color[value])

    def __get_unique_values_for_property(
            self, property_type: PropertyType) -> Set[str]:
        """
        Method that gets the unique values of the given property type for the
        lattice regions. If any 'Region' object does not have any property or
        the given property type is missing, a reference point for the region
        is stored for logging purposes. An exception showing the coordinates
        of the points of the problematic regions is raised.

        Parameters
        ----------
        property_type : PropertyType
            The type of property whose unique values to collect

        Raises
        ------
        RuntimeError
            Showing the coordinates of the points of the problematic regions

        Returns
        -------
        A set of the unique names for the given property type that have been
        associated to the lattice regions.
        """
        values = set()
        missing_regions_points = []
        for region in self.regions:
            if not region.properties or property_type not in region.properties:
                missing_regions_points.append(
                    get_point_coordinates(
                        make_vertex_inside_face(region.face)))
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
        return values

    def translate(self, new_pos: Tuple[float, float, float]) -> None:
        """
        Method that allows to translate the lattice compound in the XYZ space.
        It changes this class instance with a copy where the lattice compounds
        (face and edges) have been moved to the new position.

        Parameters
        ----------
        new_pos : Tuple[float, float, float]
                  The X-Y-Z coordinates of the position where the lattice
                  should be moved
        """
        # Store previous lattice center
        pre_center = self.lattice_center
        # Build a vector from the current center to the new one
        transl_vect = make_vector_from_points(pre_center,
                                              make_vertex(new_pos))
        # Update the lattice center
        self.lattice_center = make_vertex(new_pos)
        # Translate each cell of each layer in the lattice
        translated_layers: List[List[Cell]] = []
        # Loop through all the layers to translate the cells
        for layer in self.layers:
            translated_layers.append(
                self.__translate_cells(layer, new_pos, pre_center))
        # Update the data structure holding the translated cells for each layer
        self.layers = [
            [cell for cell in layer] for layer in translated_layers]
        # Update the list of lattice cells
        self.lattice_cells = [cell for layer in self.layers for cell in layer]

        # Translate the lattice compound
        self.lattice_cmpd = make_translation(self.lattice_cmpd, transl_vect)
        # Translate the lattice box, if any
        if self.__lattice_box:
            self.__lattice_box = self.__lattice_box.translate(new_pos)
        # Translate the lattice compound on which a symmetry operation has
        # been applied, if any
        if self.symmetry_type != SymmetryType.FULL:
            # Update the compound position relative to the shifted lattice
            # center and apply translation
            self.lattice_symm = translate_wrt_reference(
                self.lattice_symm, pre_center, new_pos)

        # Set the need to update the lattice geometry
        self.is_update_needed = True
        # Show the new lattice compound in the current SALOME study
        self.show()

    def __translate_cells(self,
                          cells: List[Cell],
                          ref_coords: Tuple[float, float, float],
                          original_ref_point: Any) -> List[Cell]:
        """
        Method that applies the translation operation to all the cells of the
        given list of `Cell` objects.
        The new position of each cell relative to the given coordinates is
        calculated. The translation is then applied.

        Parameters
        ----------
        cells : List[Cell]
            List of `Cell` objects to translate so to keep the relative
            position wrt to the given lattice center
        ref_coords : Tuple[float, float, float]
            XYZ coordinates of the point for which the relative distance of
            each cell must be kept
        original_ref_point : Any
            The vertex object each cell was relative to

        Returns
        -------
        List[Cell]
            The list of translated cells so that the relative distance from
            the new reference point is kept.
        """
        translated_cells = []
        for cell in cells:
            # Update the cell position relative to the shifted lattice center
            # and apply translation
            cell_to_center = compute_point_by_reference(
                cell.figure.o, original_ref_point, ref_coords)
            translated_cells.append(cell.translate(cell_to_center))
        return translated_cells

    def build_regions(
            self, geo_type: GeometryType = GeometryType.TECHNOLOGICAL) -> None:
        """
        Method that extracts all the subfaces from the lattice compound
        object corresponding to the given type of geometry and builds a
        'Region' object for each one.
        In case any symmetry operation has been applied, this method acts
        on that result (i.e. the symmetry compound object).

        Each lattice face object position is compared to the regions of
        each cell: if any match is found, the cell region properties are
        assigned to the 'Region' object being built for the lattice.

        Parameters
        ----------
        geo_type  : GeometryType
            The type of geometry identifying which lattice compound to
            use for building the lattice regions.
        """
        # Get the cells by assembling all the lattice layers
        self.lattice_cells = self.__assemble_layers()
        self.__update_lattice_compounds(self.lattice_cells, geo_type)
        # Get the lattice compound object, given the geometry type and the
        # current applied symmetry
        cmpd = self.__get_compound_from_type(geo_type, self.lattice_cells)
        # Extract the lattice subfaces
        subfaces = extract_sub_shapes(cmpd, ShapeType.FACE)

        # Re-initialize the lattice regions
        self.regions.clear()
        # Extract the 'Region' objects for the lattice cells only
        self.regions.extend(
            self.__build_regions_for_subfaces(0, subfaces, self.lattice_cells))
        print("No. cell regions in lattice:", len(self.regions))
        # Handle the construction of 'Region' objects for the lattice box
        # subfaces, if any
        if self.__lattice_box:
            box_subfaces = self.__extract_box_subfaces(cmpd, geo_type)
            print("No. lattice box subfaces", len(box_subfaces))
            # Add the 'Region' objects for the lattice box subfaces
            self.regions.extend(
                self.__build_regions_for_subfaces(
                    len(self.regions), box_subfaces, [self.__lattice_box]))
            # Update the lattice found subfaces
            subfaces += box_subfaces

        print("Total no. regions in lattice:", len(self.regions))
        print("Total no. subfaces in lattice:", len(subfaces))
        # Check if the number of subfaces coincides with the one of the built
        # regions
        if len(subfaces) != len(self.regions):
            raise AssertionError("Mismatch between the number of lattice "
                                 f"subfaces ({len(subfaces)}) and that of "
                                 f"found regions ({len(self.regions)})")
        # Set that there is no longer a need to update the lattice geometry
        self.is_update_needed = False

    def __build_regions_for_subfaces(
            self,
            indx: int,
            faces: List[Any],
            lattice_cells: List[Cell]
            ) -> List[Region]:
        """
        Method that builds a list of 'Region' objects for each of the given
        face objects.
        The names of the regions are set with an increasing index with a
        given starting value.
        Properties to associate to each 'Region' object are taken from the
        properties dictionary for the cell corresponding to the region. Cells
        are provided as third input to the method.

        Parameters
        ----------
        indx : int
            The starting value for the index used when assigning the region
            name
        faces : List[Any]
            A list of face objects for which 'Region' objects are built
        lattice_cells : List[Cell]
            A list of 'Cell' objects representing the lattice cells

        Returns
        -------
        A list of 'Region' objects each corresponding to a specific face
        object.
        """
        # List of regions to build
        regions = []
        # Loop through the faces
        for face in faces:
            found = False
            # Get a vertex inside the subface
            subface_point = make_vertex_inside_face(face)
            # Loop through all the regions and properties of each lattice cell
            for cell in lattice_cells:
                if not is_point_inside_shape(subface_point, cell.face):
                    continue
                for region, properties in cell.tech_geom_props.items():
                    # Check if the face point is within the lattice region
                    if is_point_inside_shape(subface_point, region):
                        # Update the index
                        indx += 1
                        # Build a 'Region' object
                        regions.append(
                            Region(
                                face=face,
                                inner_point=subface_point,
                                name=f"Region {indx}",
                                properties=deepcopy(properties)
                            ))
                        found = True
                        break
                # If a region has been found for the current cell, a new
                # lattice subface is evaluated
                if found:
                    break
            if not found:
                raise RuntimeError(
                    "No cell found for face whose inner point has "
                    f"coordinates: {get_point_coordinates(subface_point)}")
        # Return the list
        return regions

    def __get_compound_from_type(
            self, geo_type: GeometryType, lattice_cells: List[Cell]) -> Any:
        """
        Method that returns the compound object of the lattice according to
        the given type of geometry and the applied symmetry type.
        The compound for the given list of cells is extracted correspondingly
        to the geometry type.
        In case of a 'FULL' symmetry, the full lattice compound, eventually
        cut by the box, is returned. If any symmetry is applied, the shape
        of the symmetry is built and used to cut the corresponding part of
        the lattice compound.

        Parameters
        ----------
        geo_type : GeometryType
            Member of the 'GeometryType' enumeration indicating which cell
            faces to use for extracting the lattice compound
        lattice_cells : List[Cell]
            List of 'Cell' objects to extract the full compound from,
            according to the type of geometry

        Returns
        -------
        A compound object containing only the lattice cells. This is built
        from the technological geometry or the sectorized one; the currently
        applied symmetry determines the layout of the returned compound.
        """
        # Get the lattice compound that matches the indicated geometry type
        lattice_cmpd = get_compound_from_geometry(geo_type, lattice_cells)
        # Handle the different type of symmetry to build the final compound
        match self.symmetry_type:
            case SymmetryType.FULL:
                # If any lattice box is present, cut it out from the lattice
                # compound object
                return self.__cut_lattice_with_box(lattice_cmpd)
            case _:
                # Extract the shape of the symmetry by extracting its
                # borders and building a face object
                shape = make_face(
                    build_compound_borders(self.lattice_symm))
                # If any lattice box is present, cut it out from the lattice
                # compound object. Then, perform a common operation of the
                # result with the symmetry shape.
                return make_common(
                    self.__cut_lattice_with_box(lattice_cmpd), shape)

    def __cut_lattice_with_box(self, cmpd: Any) -> Any:
        """
        Method that cuts the given compound object with the box area closest
        to the lattice center, if any container has been declared.
        If not, the same compound is returned untouched.

        Parameters
        ----------
        cmpd : Any
            The compound object to cut with the inmost box area

        Returns
        -------
        The given compound object cut by the inmost box area or the same
        compound, if no container has been declared.
        """
        # If any lattice box is present, cut it out from the
        # lattice compound object (the one with only the cells
        # and no box)
        if self.__lattice_box:
            return make_common(cmpd, self.__extract_inner_box())
        return cmpd

    def add_ring_of_cells(
            self,
            cell: Cell,
            ring_indx: int,
            layer_indx: int | None = None) -> None:
        """
        Method that adds a ring of cells of the same type to the lattice at
        a given ring index and layer.
        This method iteratively adds the provided `Cell` object at specific
        construction points determined by the `ring_indx` parameter, which
        indicates the ring (distance from the center) where the cells should
        be placed. Optionally, a `layer_indx` can be provided to specify the
        layer to which the ring of cells is added; if not provided, the cells
        are added to a new layer.

        Notes
        -----
        - Index 0 refers to the center cell of the lattice (valid for lattices
          with hexagonal cells or cartesian ones, the latter if an odd number
          of cells is provided on the side); indices greater than 0 refer to
          rings around the center.
        - The method ensures that the cell to add has the same `CellType` of
          the ones in the lattice.

        Parameters
        ----------
        cell : Cell
            The `Cell` instance to be added repeatedly to form a ring
        ring_indx : int
            The index indicating which ring to add the cells to
        layer_indx : int | None = None
            The index of the layer to which the cells are added. If None, a
            new layer is created.

        Raises
        ------
        RuntimeError
        - If the cell's `CellType` does not match the one of the lattice
          cells.
        - If indicating `0` as the ring index where cells should be added.
        - If the cell to add has an invalid rotation angle.
        - If the type of cells of the lattice is not supported.
        """
        try:
            # Check that the cell to add has the same type of the ones already
            # present
            self.__check_cell_type(cell.cell_type)
            # Check the validity of the indicated ring index
            if ring_indx == 0:
                raise RuntimeError(
                    "It is not possible to add a ring of cells at the "
                    "indicated 0 index.")
            # Check the validity of the cell's rotation angle and update the
            # value representative of the main pattern of cells
            self.__set_main_rotation_angle(math.degrees(cell.rotation))
        except RuntimeError as e:
            raise RuntimeError(f"Error when adding a ring of cells: {e}")
        # If no layer index is provided, the ring of cells is added to a new
        # layer
        if layer_indx is None:
            layer_indx = len(self.layers)
            self.layers.append([])
        # Evaluate the multiplication factor for determining the construction
        # figure where the centers of the cells will be placed
        if not cell.cell_type == CellType.HEX:
            for cell in self.lattice_cells:
                if get_min_distance(cell.figure.o,
                                    self.lattice_center) < 1e-5:
                    n = 2*ring_indx
                    break
            else:
                n = 2*ring_indx - 1
        else:
            n = 2*ring_indx
        # Build the construction figure accordingly with the type of lattice
        # cells (either a rectangle or a hexagon figure)
        match self.cells_type:
            case CellType.RECT:
                construction_fig = Rectangle(
                    get_point_coordinates(self.lattice_center),
                    n*cell.figure.ly,
                    n*cell.figure.lx)
                # Rotate by 90° the construction figure if the cells have
                # a rotation angle of 90°
                if math.isclose(90.0, self.__cells_rot):
                    construction_fig.rotate(90.0)
            case CellType.HEX:
                construction_fig = Hexagon(
                    get_point_coordinates(self.lattice_center),
                    n*cell.figure.ly)
                # Rotate by 90° the construction figure if the cells have
                # a rotation angle of 0°
                if math.isclose(0.0, self.__cells_rot):
                    construction_fig.rotate(90.0)
                # Parameter for identifying the subdivision points
                n = int(n/2)
            case _:
                raise RuntimeError(f"The {self.cells_type} is not handled.")
        # Loop through the borders of the construction figure and place a
        # cell on each subdivision point, given by the ring index and the
        # subdivision parameter
        for border in construction_fig.borders:
            for i in range(0, n):
                # Build a subdivision point onto the figure borders
                p = get_point_coordinates(make_vertex_on_curve(border, i/(n)))
                # Add the cell at the position of the subdivision point
                self.__add_cell_to_layer(cell, p, layer_indx)
        # Set the need to update the lattice geometry
        self.is_update_needed = True

    def __check_cell_type(self, cell_type: CellType):
        """
        Method that checks the consistency of the given cell type with the
        one associated to the lattice cells.

        If the lattice does not yet have a defined cell type, it is set to
        the provided one. If the cell type for the lattice is defined, an
        exception is raised if the provided one differs from it.

        Parameters
        ----------
        cell_type : CellType
            The geometric type of the cell to check or set.

        Raises
        ------
        RuntimeError
            If the provided `cell_type` does not match the existing lattice
            cell type.
        """
        # Store the value, if not already set
        if not self.cells_type:
            self.cells_type = cell_type
        else:
            if cell_type != self.cells_type:
                raise RuntimeError(
                    f"The geometric type '{cell_type}' differs from the "
                    f"'{self.cells_type}' one of the cells in the lattice.")

    def add_rings_of_cells(self, cell: Cell, no_rings: int) -> None:
        """
        Method that adds to the lattice several rings of cells of the same
        type. The cell is provided as an object of one of the 'Cell'
        subclasses, which is iteratively added, for each index of rings, at
        specific construction points.
        The rings of cells are added starting from the current maximum value
        of rings for the lattice and to a new layer.

        Parameters
        ----------
        cell      : Cell
                    The cell instance to be iteratively added in order to
                    build the rings of cells
        no_rings  : int
                    The number of rings to add starting from the current
                    maximum value of rings for the lattice
        """
        # Raise an exception if the number of rings to add is less than 1
        if no_rings < 1:
            raise AssertionError(f"Wrong number ({no_rings}) of rings of "
                                 "cells to add has been indicated.")
        # Initialize a new sub list identifying a new layer of cells
        self.layers.append([])
        # Loop through the ring indices starting from the current number of
        # rings of cells
        for i_ring in range(self.rings_no+1, no_rings+1):
            # Add a ring of cells at the current index
            self.add_ring_of_cells(cell, i_ring, len(self.layers) - 1)
        # Set the need to update the lattice geometry
        self.is_update_needed = True

    def set_lattice_box_properties(
            self, properties: Dict[PropertyType, List[str]]) -> None:
        """
        Method that allows to set the properties of the lattice box itself and
        to the area between the cells and the box contour, if any.
        The assumed convention for the elements in the list is outwards from
        the lattice center. If any area between the cells and the box is
        present, the properties occupies the first position in the list for
        each type, then the ones for the box layers follow always moving from
        the inner to the outer layer wrt the lattice center.

        Parameters
        ----------
        properties  : Dict[PropertyType, List[str]]
                      Dictionary collecting the properties for each region
                      of the lattice box; different types, provided by the
                      'PropertyType' enumeration, can be provided
        """
        if not properties:
            raise RuntimeError("No properties have been provided")
        if not self.box_layers and self.__lattice_box:
            raise RuntimeError(
                "The lattice box has been declared by directly assigning the "
                "cell instance instead of providing its layers. This method "
                "only works when building the box from the layers by calling "
                "the method 'build_lattice_box().'")
        # Check if the number of box regions coincides with the number of
        # property elements for all the given property types: for hexagonal
        # lattices it is needed to consider also the area between the cells
        # and the box, hence the number of regions are 1 more than the case
        # for cartesian cells.
        no_regions = len(self.box_layers)
        if self.cells_type == CellType.RECT:
            no_regions -= 1
        for type, values in properties.items():
            if not no_regions == len(values):
                message = "There is no correspondence between the number " +\
                          f"of box regions ({no_regions}) and that of the " +\
                          f"properties with type '{type.name}' (no. " +\
                          f"{len(values)})"
                raise RuntimeError(message)

        # Clear any previously set entry in the dictionary associating the
        # regions to the properties, as this method sets the properties for
        # all.
        self.__lattice_box.tech_geom_props.clear()
        # Extract the lattice box subfaces
        box_subfaces = extract_sorted_sub_shapes(self.__lattice_box.face,
                                                 ShapeType.FACE)
        # Sort the box regions according to their distance from the lattice
        # center and their perimeter in reverse order: this so that the
        # outmost layer of the box is the first element and so on with the
        # others.
        box_subfaces = sorted(
            box_subfaces,
            key=lambda subface: get_min_distance(self.lattice_center,
                                                 make_cdg(subface))
                                and get_basic_properties(subface)[0],
            reverse=True)
        indx = 0
        # Assign the value for each property type to all the box regions
        for region in box_subfaces:
            prop_types = {}
            for type, values in properties.items():
                # Reverse-sort the values
                values = values[::-1]
                if indx < len(values) - 1:
                    prop_types[type] = values[indx]
                else:
                    # Properties for regions between cells and layers
                    prop_types[type] = values[-1]
            # Associate a region with its properties
            self.__lattice_box.tech_geom_props[region] = prop_types
            # Update the index
            indx += 1

    def __extract_box_subfaces(
            self, lattice_cmpd: Any, geom_type: GeometryType) -> List[Any]:
        """
        Method that extracts the lattice box subfaces, i.e. the box layers
        and the areas between the cells and the container layers, if any
        (as depending on the cells geometry).
        If any symmetry type, other than `FULL`, is applied, the common part
        between the box and the shape of the symmetry is extracted first.
        By cutting the box with the compound made of lattice cells only, the
        remaining parts are the ones only strictly belonging to the box.

        Parameters
        ----------
        lattice_cmpd : Any
            The lattice compound object made by the cells only
        geom_type : GeometryType
            Indicating the geometry layout of the box cell from which face
            objects are extracted

        Raises
        ------
        RuntimeError
            - If no face objects are available after cutting the box face
              with the compound made from the lattice cells only.
            - In case the geometry type is different from `TECHNOLOGICAL`
              or `SECTORIZED` cases.

        Returns
        -------
        A list of face objects of the container layers and the areas between
        the cells and the box. Elements are sorted according to their distance
        from the lattice center and their perimeter in reverse order. The
        result is that the outmost layer of the box occupies the first
        position in the returned list and so on with the others.
        """
        # Get the box geometry layout according to the given geometry type
        box = get_compound_from_geometry(geom_type, [self.__lattice_box])
        # Handle symmetry types other than FULL
        if self.symmetry_type != SymmetryType.FULL:
            # Get the shape of the symmetry
            shape = make_face(build_compound_borders(self.lattice_symm))
            # Extract the common part between the box and the symmetry shape
            box = make_common(box, shape)
        # Cut the box with the lattice cells, so that only the box areas
        # remain, and extract its face objects, if any
        box_faces = extract_sorted_sub_shapes(
            make_cut(box, lattice_cmpd), ShapeType.FACE)
        if not box_faces:
            raise RuntimeError(
                "Error in extracting the box areas of the lattice.")
        # Return the sorted box regions
        return sorted(
            box_faces,
            key=lambda subface: get_min_distance(self.lattice_center,
                                                 make_cdg(subface))
                                and get_basic_properties(subface)[0],
            reverse=True)

    def get_regions_info(self) -> None:
        """
        Method for retrieving descriptive information about any of the lattice
        regions whose corresponding object has been selected in the SALOME
        study.
        """
        # Extract the geometrical objects the given ID corresponds to in the
        # current SALOME study
        shape = retrieve_selected_object(
            "Please, select a single region whose data to show.")
        # Get the region that corresponds to the given shape and print the
        # corresponding data
        print(get_region_info(shape, self.regions))

    def rotate(self, angle: float = 0.0) -> None:
        """
        Method for rotating the whole lattice by the given angle in degrees.

        Parameters
        ----------
        angle : float
                The rotation angle in degrees
        """
        # Return immediately if the rotation angle is 0.0
        if math.isclose(angle, 0.0):
            print("No rotation is performed as the given angle is 0.0°")
            return
        # Convert the rotation angle in radians
        rotation = math.radians(angle)

        # Build the Z-axis of rotation
        x, y, _ = get_point_coordinates(self.lattice_center)
        z_axis = make_vector_from_points(
            self.lattice_center, make_vertex((x, y, 1)))
        # Rotate the lattice compound
        self.lattice_cmpd = make_rotation(self.lattice_cmpd, z_axis, rotation)
        # Loop through all the layers to rotate the cells
        rotated_layers: List[List[Cell]] = []
        for layer in self.layers:
            rotated_layers.append(
                self.__rotate_cells(layer, angle, z_axis))
        # Update the data structure holding the rotated cells for each layer
        self.layers = [
            [cell for cell in layer] for layer in rotated_layers]
        # Update the list of lattice cells
        self.lattice_cells = [cell for layer in self.layers for cell in layer]
        # Rotate the lattice box, if any
        if self.__lattice_box:
            self.__lattice_box.rotate_from_axis(angle, z_axis)
        # Rotate the lattice compound on which a symmetry operation has been
        # applied, if any
        if self.symmetry_type != SymmetryType.FULL:
            self.lattice_symm = make_rotation(
                self.lattice_symm, z_axis, rotation)
        # Update the rotation angle of the main pattern of cells
        self.__cells_rot = self.__evaluate_cells_rotation(self.lattice_cells)

        # Set the need to update the lattice geometry
        self.is_update_needed = True
        # Show the new lattice compound in the current SALOME study
        self.show()

    def __rotate_cells(self,
                       cells: List[Cell],
                       angle: float,
                       axis: Any) -> List[Cell]:
        """
        Method that applies the rotation operation to all the cells of the
        given list of `Cell` objects. It returns a list of the rotated
        cells.

        Parameters
        ----------
        cells : List[Cell]
            List of `Cell` objects to rotate
        angle : float
            The rotation angle in degrees
        axis : Any
            The vector object being the rotation axis

        Returns
        -------
        List[Cell]
            The list of rotated cells.
        """
        rotated_cells = []
        for cell in cells:
            cell.rotate_from_axis(angle, axis)
            rotated_cells.append(deepcopy(cell))
        return rotated_cells

    def set_region_property(
            self,
            property_type: PropertyType,
            value: str,
            region: Any | None = None) -> None:
        """
        Method that allows to set the value of the given type of property for
        the lattice region which is currently selected in the SALOME study.
        After applying the modification, the lattice is shown again in the
        SALOME viewer with the colorset associated to the values of the type
        of property.

        Parameters
        ----------
        property_type : PropertyType
            The member of the 'PropertyType' enumeration indicating which
            type of property to modify for the selected region
        value : str
            The value of the property type to assign to the selected region
        region : Any | None = None
            The lattice's region of the technological geometry whose property
            to change. When not provided, the region currently selected is
            considered.
        """
        # Check which cell geometry type is currently shown; if different
        # from the TECHNOLOGICAL one, raise an exception
        if self.displayed_geom != GeometryType.TECHNOLOGICAL:
            raise RuntimeError(
                f"Currently showing the {self.displayed_geom.name} type of "
                "geometry. To set lattice regions properties, show the '"
                "TECHNOLOGICAL' geometry first.")
        # Extract the geometrical object currently selected in the current
        # SALOME study, if no one is provided as input
        if region is None:
            region = retrieve_selected_object(
                "Please, select a single region to assign a property to.")

        # Get the list of cells, eventually cut by the lattice box
        cells = self.__get_cells_and_box()
        # Get the region that corresponds to the given shape and set the
        # value for the indicated property
        self.__set_region_property(region, cells, property_type, value)
        # Set the need to update the lattice geometry
        self.is_update_needed = True

    def __get_cells_and_box(self) -> List[Cell]:
        """
        Method that provides a copy of the list of `Cell` objects comprising
        all the cells in the lattice including the box, if present.
        The `Cell` object for the lattice box is cut by removing all the area
        occupied by the cells.

        Returns
        -------
        List[Cell]
            A list of `Cell` objects made by the lattice cells and the box
            cell, if any.
        """
        # Get the compound corresponding to the technological geometry of the
        # cells and the currently applied symmetry
        cmpd = self.__get_compound_from_type(GeometryType.TECHNOLOGICAL,
                                             self.lattice_cells)
        # Store all the cells to look for the region to update
        cells: List[Cell] = []
        # If a box is present, add it to the list of cells
        if self.__lattice_box:
            # Copy the lattice box cell and update its face by cutting out the
            # cells; this to prevent to select the box region instead of a
            # cell as the box is not empty
            lattice_box = deepcopy(self.__lattice_box)
            lattice_box.face = make_compound(
                self.__extract_box_subfaces(cmpd, GeometryType.TECHNOLOGICAL))
            cells.append(lattice_box)
        # Add all the lattice cells traversing the layers in reverse order
        cells.extend([cell for layer in self.layers[::-1] for cell in layer])
        return cells

    def __set_region_property(
            self,
            shape: Any,
            cells: List[Cell],
            property_type: PropertyType,
            value: str) -> None:
        """
        Method that looks for the given region in the dictionary of cells'
        regions vs properties and sets the value for the indicated property
        type.

        Parameters
        ----------
        shape : Any
            The lattice's region of the technological geometry whose property
            to change.
        cells : List[Cell]
            The list of `Cell` objects made by all the lattice's cells and
            the box.
        property_type : PropertyType
            Indicating the type of property the value associated to the
            region needs to be changed.
        value : str
            The value of the property type to set.

        Raises
        ------
        RuntimeError
            If no region can be found among the ones stored in the properties
            dictionary for each of the analysed cells.
        """
        # Point for idendifying the shape in the geometry
        point = make_vertex_inside_face(shape)
        # Get the region that corresponds to the given shape
        for i, cell in enumerate(cells):
            # Continue with another cell if the shape is not contained within
            # the cell's area (reference point inside cell's face and common
            # area between cell and shape is the shape itself)
            if (not is_point_inside_shape(point, cell.face)):
                continue
            # Search for the region among the faces stored in the cell
            # dictionary of regions VS properties (box case, if any)
            if i == 0 and self.__lattice_box is not None:
                for zone in self.__lattice_box.tech_geom_props:
                    if (is_point_inside_shape(point, zone)):
                        self.__lattice_box.tech_geom_props[zone][
                            property_type] = value
                        return
                else:
                    continue
            # Search for the region among the cells (lattice's cells)
            for region in cell.tech_geom_props:
                if (is_point_inside_shape(point, region) and
                    are_same_shapes(shape, region, ShapeType.FACE)):
                    cell.tech_geom_props[region][property_type] = value
                    return
        # Raise an exception if no region has been found
        else:
            raise RuntimeError(
                "No lattice region could be found for the selected shape.")

    def restore_cells(self,
                      cells: List[Cell],
                      properties: Dict[PropertyType, str],
                      ignore_not_cut: bool = True) -> None:
        """
        Method that restores the geometry layout of the given cells by
        removing any circular region, while setting properties accordingly
        with the provided dictionary.
        If any cell has no centered circular regions, the restore operation
        is not performed.
        If the `ignore_not_cut` boolean flag is `True`, this operation is
        performed only if the compound made of the circular regions is cut
        by the cell's face. Otherwise, if `True`, all the given cells are
        restored.

        Parameters
        ----------
        cells: List[Cell]
            List of shallow copies of the ones in the lattice that should
            be restored.
        properties : Dict[PropertyType, str]
            Providing the values for each type of property to assign to the
            cells being restored.
        ignore_not_cut : bool = True
            Boolean flag indicating whether the cells whose circular regions
            are not cut should also be restored or not.

        Notes
        -----
        This method should be called after building the lattice
        regions by overlapping all the layers of cells, i.e. after calling
        either the `show` or the `build_regions` methods.
        This method has an effect only if the given list of cells contains
        shallow copies of the `Cell` objects included in the present `Lattice`
        instance.
        """
        # Filter the cells to be restored
        if ignore_not_cut:
            # Loop through all the given cells to get the ones to restore
            cells_to_restore = [
                cell for cell in cells if check_cell_circle_are_cut(cell)
            ]
        else:
            cells_to_restore = cells
        # Restore the geometry and properties of the cells
        for cell in cells_to_restore:
            cell.restore()
            cell.set_properties({k: [v] for k, v in properties.items()})

        # Set the need to update the lattice geometry
        self.is_update_needed = True

    def set_type_geo(self, type_geo: LatticeGeometryType) -> None:
        """
        Method that sets the geometry type for the lattice, as item of the
        `LatticeGeometryType` enumeration.
        It checks if the provided geometry type is valid for the type of
        cells of the lattice and the currently applied type of symmetry.
        In case of types involving TRAN-type of BCs, the assignement is done
        only if the lattice is either made by a single cell or if enclosed in
        a box.
        An exception is raised if the geometry type is not compatible.

        Parameters
        ----------
        type_geo : LatticeGeometryType
            The geometry type to set for the lattice.

        Raises
        ------
        RuntimeError
            If the given geometry type is not compatible with the type of
            cells and applied symmetry type.
        """
        # Check if the given type of geometry is valid for the type of cells
        # in the lattice and the currently applied type of symmetry.
        try:
            # Get the list of types of geometry available for the lattice
            types_geo = CELL_VS_SYMM_VS_TYP_GEO[self.cells_type][
                self.symmetry_type]
            if type_geo not in types_geo:
                raise KeyError
        except KeyError:
            raise RuntimeError(
                f"The given type of geometry '{type_geo}' is not compatible "
                "with the type of cells of the lattice (i.e. "
                f"'{self.cells_type}') and the applied symmetry type '"
                f"{self.symmetry_type}'. Expected values are {types_geo}.")

        # Check if the given type of geometry is compatible with the lattice
        # layout; types involving TRAN-type of BCs can be applied only if the
        # lattice has either a single cell or is enclosed in a box
        if type_geo in [LatticeGeometryType.HEXAGON_TRAN,
                        LatticeGeometryType.RECTANGLE_TRAN]:
            if len(self.lattice_cells) > 1 and self.__lattice_box is None:
                raise RuntimeError(
                    f"The given type of geometry '{type_geo}' is not "
                    "compatible with the current lattice geometry layout as "
                    "made by more than one cell without being enclosed in a "
                    "box."
                )
        # Assign the lattice type of geometry
        self.__type_geo = type_geo

    def __update_attributes(self, cells: List[Cell]):
        """
        Method that updates the instance attributes according to the provided
        list of `Cell` objects, unless it is empty.

        This method performs the following checks:
          - verifies that all provided cells are of the same type;
          - evaluates and verifies that all the cells have the same rotation
            angle.
        It updates:
          - the `CellType` of the cells in the lattice;
          - the common rotation angle of the given cells;
          - the `LatticeGeometryType` indicating the lattice type of geometry;
          - the characteristic dimensions of the lattice using the ones of the
            first cell in the given list;
          - the compound object for the whole lattice;
          - the total number of rings of cells based on the provided cells.

        Parameters
        ----------
        cells : List[Cell]
            List of `Cell` objects to use for updating the lattice attributes.

        Raises
        ------
        RuntimeError
          - If the provided cells have different geometry types.
          - If there is no common rotation angle or any cell does not have
            one of the admitted values.
        """
        # If no cells are provided, return immediately without assembling the
        # lattice compound
        if not cells:
            return
        # Check that the input cells are all of the same type
        self.cells_type = cells[0].cell_type
        if not all(cell.cell_type == self.cells_type for cell in cells):
            raise RuntimeError(
                "The lattice presents cells with different geometry types.")
        # Check if the rotation angle of the given cells is the same for all.
        # If so, store the value.
        self.__cells_rot = self.__evaluate_cells_rotation(cells)
        # Set the lattice type of geometry according to the cells' type and
        # the number of cells
        self.__configure_lattice_type(self.cells_type,
                                      len(self.lattice_cells))
        # Update the characteristic dimensions of the lattice with the ones
        # of the first of the given cells
        self.lx = cells[0].figure.lx
        self.ly = cells[0].figure.ly

        # Build the GEOM compound object assembling all the cells' faces
        # of the lattice, eventually with the box
        self.__build_lattice()

        # Update the total number of cells' rings
        for cell in cells:
            self.__evaluate_no_rings(get_point_coordinates(cell.figure.o))

        # Set the need to update the lattice geometry
        self.is_update_needed = True


def get_compound_from_geometry(
        geo_type: GeometryType, lattice_cells: List[Cell]) -> Any:
    """
    Function that gets the lattice compound that matches the indicated
    geometry type.

    Parameters
    ----------
    geo_type : GeometryType
        The type of geometry of the cells to build a compound from
    lattice_cells : List[Cell]
        The list of cells whose face objects to use to build a compound

    Returns
    -------
    A compound object collecting all the faces, either from the technological
    geometry or the sectorized one, of the given cells.
    """
    cell_faces = []
    match geo_type:
        case GeometryType.TECHNOLOGICAL:
            cell_faces = [cell.face for cell in lattice_cells]
        case GeometryType.SECTORIZED:
            for cell in lattice_cells:
                if not cell.sectorized_face:
                    print(f"The {cell.name} cell has not been sectorized: "
                          "the non sectorized face will be used.")
                    cell_faces.append(cell.face)
                    continue
                cell_faces.append(cell.sectorized_face)
        case _:
            raise ValueError(f"{geo_type}: unhandled type.")
    return make_compound(cell_faces)

def get_changed_cells(lattice: Lattice) -> List[Cell]:
    """
    Function that returns a list of `Cell` objects belonging to the given
    `Lattice` instance. These cells have their geometry layout changed
    compared to their original one.
    For each cell in each layer, the current shape of the cell is built and
    its area compared with the area of its specific characteristic figure (as
    a `Surface` instance).
    Those showing a different value for the area indicates a change in their
    geometry layout has occurred and are collected into the returned list.

    This function can be used to retrieve those cells that have been modified
    (e.g., by overlap with a superior layer of cells) within the lattice.

    Parameters
    ----------
    lattice : Lattice
        The lattice instance to check for cells whose geometry layout has
        changed.

    Returns
    -------
    List[Cell]
        A list of `Cell` objects whose geometry layout differs from their
        original one.
    """
    cells = []
    for layer in lattice.layers:
        for cell in layer:
            # Build a face object over the borders of the cell
            cell_shape = make_face(build_compound_borders(cell.face))
            # Compare the area of the current cell's face with the one of its
            # original shape
            if get_basic_properties(cell_shape)[1] != \
                get_basic_properties(cell.figure.face)[1]:
                cells.append(cell)
    # Return the list of changed cells
    return cells
