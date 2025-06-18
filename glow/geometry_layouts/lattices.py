"""
Module containing classes providing the means for creating lattices from base
cells types in SALOME.
"""
import math

from copy import deepcopy
from typing import Any, Dict, List, Set, Tuple, Union

from glow.geometry_layouts.cells import Cell, HexCell, RectCell, Region
from glow.geometry_layouts.geometries import GenericSurface, Hexagon, \
    Rectangle, build_hexagon
from glow.geometry_layouts.utility import build_compound_borders, \
    update_relative_pos
from glow.interface.geom_interface import ShapeType, add_to_study, \
    add_to_study_in_father, clear_view, display_shape, \
    extract_sorted_sub_shapes, extract_sub_shapes, get_basic_properties, \
    get_bounding_box, get_min_distance, get_object_from_id, \
    get_point_coordinates, get_selected_object, is_point_inside_shape, \
    make_cdg, make_common, make_compound, make_cut, make_edge, make_face, \
    make_partition, make_rotation, make_translation, make_vector_from_points, \
    make_vertex, make_vertex_inside_face, make_vertex_on_curve, \
    remove_from_study, set_color_face, update_salome_study
from glow.generator.support import GeometryType, LatticeGeometryType, \
    PropertyType, SymmetryType, BoundaryType, CellType, TYPEGEO_VS_BC, \
    generate_unique_random_colors


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
    lattice_edges     : Any
                        A GEOM compound object grouping all the unique edges
                        of the cells in the lattice
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
    VALID_CELLS_ANGLES: List[float] = [0.0, 90.0]

    def __init__(self,
                 cells: List[Cell],
                 name: str = "Lattice",
                 center: Union[Tuple[float, float, float], None] = None,
                 boxes_thick: List[float] = []) -> None:
        # --------------
        # Initial checks
        # --------------
        # Check the input list is not empty
        if not cells:
            raise AssertionError("The lattice cannot be built without "
                                 "specifying any cell.")
        # Check if the lattice center is coherent with the given cells, if any
        # is provided. If not, the XYZ origin is selected.
        self.lattice_center: Any = self.__evaluate_lattice_center(cells,
                                                                  center)
        # Check that the input cells have all the same type
        self.cells_type: CellType = cells[0].cell_type
        if not all(cell.cell_type == self.cells_type for cell in cells):
            raise AssertionError("The lattice presents cells with different "
                                 "geometry types.")
        # Check if the rotation angle of the given cells is the same for all.
        # If so, store the value.
        self.cells_rot: float = self.__evaluate_cells_rotation(cells)

        # Initialize the instance attributes
        self.lattice_cells: List[Cell] = deepcopy(cells)
        self.layers: List[List[Cell]] = [deepcopy(self.lattice_cells)]
        self.name: str = name
        self.lattice_entry_id: Union[str, None] = None
        self.rings_no: int = 0
        self.distance: float = 0.0
        # Set the variables related to the geometry and BC types according to
        # the cells type
        self.type_geo: LatticeGeometryType = LatticeGeometryType.HEXAGON_TRAN
        self.boundary_type: BoundaryType = TYPEGEO_VS_BC[self.type_geo]
        self.__configure_lattice_types(self.cells_type,
                                       len(self.lattice_cells))

        self.symmetry_type: SymmetryType = SymmetryType.FULL # TODO Add method for setting
        # self.boundary_type: int = TYPEGEO_VS_BC[self.type_geo][0]  # TODO Add method for setting
        # TODO check if each boundary can be set to a different BC type
        # FIXME the 'lx' and 'ly' characteristic dimensions should be set in
        # accordance to the number of cells in the lattice. Here they are
        # initialized to the values for the first cell
        self.lx: float = cells[0].figure.lx
        self.ly: float = cells[0].figure.ly
        self.box_layers: List[float] = boxes_thick
        self.lattice_box: Union[Cell, None] = None
        self.lattice_symm: Union[Any, None] = None
        self.lattice_tech: Union[Any, None] = None # FIXME not necessary?
        self.regions: List[Region] = []
        self.type_vs_property_vs_color: Dict[ # FIXME not necessary?
            PropertyType, Dict[str, Tuple[int, int, int]]] = {}
        self.displayed_geom: GeometryType = GeometryType.TECHNOLOGICAL

        # Build the GEOM compound objects for the lattice storing the cells
        # faces and the unique edges respectively, as well as the lattice
        # box, if any should be present
        self.__build_lattice()

        # Update the total number of cells rings
        for cell in cells:
            # self.__add_cell_regions(cell)
            self.__evaluate_no_rings(get_point_coordinates(cell.figure.o))

        # Show the lattice in the current SALOME study
        self.show()

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
        if not cell_rot in self.VALID_CELLS_ANGLES:
            raise RuntimeError("Cells can only have any of the "
                               f"{self.VALID_CELLS_ANGLES}° rotation angles")
        # Return the rotation angle of all the cells
        return cell_rot

    def __evaluate_lattice_center(
            self,
            cells: List[Cell],
            center: Union[Tuple[float, float, float], None]) -> Any:
        """
        Method that evaluates the lattice center during the initialization:
        if no center is provided at class instantiation, the XYZ origin is
        returned. In the opposite case, a check is performed to identify if
        there is a cell whose center coincides with the specified one.
        In case none is found, an exception is raised.

        Parameters
        ----------
        cells   : List[Cell]
                  A list of cell objects constituting the lattice
        center  : Tuple[float, float, float]
                  The X-Y-Z coordinates of the lattice center

        Returns
        -------
        The vertex object identifying the lattice center.
        """
        if center:
            # Check if the given lattice center corresponds to any of the given
            # cells centers
            for cell in cells:
                if get_min_distance(cell.figure.o, make_vertex(center)) < 1e-7:
                    # Return the center of the cell being the lattice center
                    return cell.figure.o
            # Handle the case where the specified center does not corresponds
            # to any of the cells ones
            raise Exception(f"The specified center {center} does not "
                            "correspond to any of the given cells centers.")
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
        # Build a GEOM compound storing the lattice unique edges of its cells
        self.lattice_edges: Any = make_partition( # FIXME not necessary?
            [self.lattice_cmpd], [], ShapeType.EDGE)
        # Store the original lattice compound (either technological or
        # sectorized geometry) without any box
        self.lattice_org = deepcopy(self.lattice_cmpd) # FIXME not necessary?

        # Handle the construction of the lattice container, if any is required
        self.__assemble_box()

    def __assemble_box_with_lattice(
            self, inner_box: Any, lattice_cmpd: Any) -> Any:
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

        Returns
        -------
        The lattice compound assembled with its box.
        """
        # Extract the common part between the lattice compound and the inner
        # box.
        lattice_cmpd = make_common(inner_box, lattice_cmpd)
        # Update the lattice compound by partitioning it with the box face
        lattice_cmpd = make_partition(
            [lattice_cmpd, self.lattice_box.face], [], ShapeType.FACE)
        # Update the lattice edges by partitioning the lattice compound
        self.lattice_edges = make_partition(
            [lattice_cmpd], [], ShapeType.EDGE)

        # Return the lattice assembled with the box
        return lattice_cmpd

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
        # Get the box subface closest to the lattice center
        inner_box = self.__extract_inner_box()
        # Update both the current lattice compound and the one coming from
        # the cells technological geometry
        self.lattice_cmpd = self.__assemble_box_with_lattice(
            inner_box, self.lattice_org)
        self.lattice_tech = self.__assemble_box_with_lattice(
            inner_box, self.lattice_org)
        # Re-apply the symmetry operation if any has already been performed
        if self.symmetry_type != SymmetryType.FULL:
            self.apply_symmetry(self.symmetry_type)
        # Update the 'lx' characteristic dimension of the lattice
        self.lx = self.lattice_box.figure.lx

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
        box_faces = extract_sub_shapes(self.lattice_box.face, ShapeType.FACE)
        # Get the inner box face by sorting the faces according to the
        # distance from the lattice center
        return sorted(
            box_faces,
            key=lambda item: get_min_distance(self.lattice_center, item))[0]

    def __build_lattice_box_type(self) -> None:
        """
        Method that builds the lattice box as an instance of the 'Cell'
        class, which is stored as an instance attribute.
        The container geometry is built accordingly with the type of geometry
        of the cells in the lattice.
        Either a rectangle or a hexagon is built for each layer of the box and
        a partition with all the figures that make up the box is performed.
        Either a 'RectCell' or a 'HexCell' is built for the lattice box and
        with its face updated with the figure one built herein.
        """
        # Declare a list storing the geometrical figures that constitute
        # the box
        box_surfaces: List[GenericSurface] = []
        center = get_point_coordinates(self.lattice_center)
        # Build the container for the lattice according to the cells geometry
        # type
        match self.cells_type:
            case CellType.RECT:
                # Declare the starting dimensions of the box
                height = self.lattice_cells[0].figure.ly*(2*self.rings_no + 1)
                width = self.lattice_cells[0].figure.lx*(2*self.rings_no + 1)

                # Build a rectangle for each layer
                for thick in self.box_layers:
                    height += 2*thick
                    width += 2*thick
                    box_surfaces.append(Rectangle(center, height, width))
                # Perform a 'partition' operation to assemble the lattice box
                box_face = make_partition(
                    [rect.face for rect in box_surfaces], [], ShapeType.FACE)
                # Declare a 'Rectangle' instace from the built face
                self.lattice_box = RectCell(center, (height, width))
            case CellType.HEX:
                # Calculate the apothem of the hexagon enclosing the lattice
                box_apothem = 0.5 * self.lattice_cells[0].figure.lx * (
                    math.cos(math.pi/3) + (1 + math.cos(math.pi/3)) *
                    (1 + self.rings_no*2))

                # Build a hexagon for each layer
                for thick in self.box_layers:
                    box_apothem += thick
                    box_surfaces.append(
                        build_hexagon(box_apothem, center))
                # Perform a 'partition' operation to assemble the lattice box
                box_face = make_partition(
                    [hex.face for hex in box_surfaces], [], ShapeType.FACE)
                # Declare a 'Hexagon' instace from the built face
                self.lattice_box = HexCell(center,
                                           box_apothem/math.cos(math.pi/3))
            case _:
                raise RuntimeError("Unhandled cell geometry type.")
        self.lattice_box.update_geometry_from_face(GeometryType.TECHNOLOGICAL,
                                                   box_face)
        # Rotate the box so to enclose the lattice (hexagonal lattice only)
        if self.cells_rot == 0.0 and self.cells_type == CellType.HEX:
            self.lattice_box.rotate(90)

    def __configure_lattice_types(self,
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
                    # characterised by specular reflection on each side
                    self.type_geo = LatticeGeometryType.RECTANGLE_SYM
                else:
                    # Only one cartesian cell is present: this case is
                    # characterised by an cartesian geometry type with
                    # translation on all sides
                    self.type_geo = LatticeGeometryType.RECTANGLE_TRAN
            case CellType.HEX:
                if no_cells > 1:
                    # Case of a lattice made of hexagonal cells: this case is
                    # characterised by an isotropic geometry type with VOID
                    # BCs.
                    self.type_geo = LatticeGeometryType.ISOTROPIC
                else:
                    # Only one hexagonal cell is present: this case is
                    # characterised by an hexagonal geometry type with
                    # translation on all sides
                    self.type_geo = LatticeGeometryType.HEXAGON_TRAN

        # Set the BC type according to the geometry type
        self.boundary_type: BoundaryType = TYPEGEO_VS_BC[self.type_geo][0]

    def add_cell(self,
                 cell_to_add: Cell,
                 position: Tuple[float, float, float]) -> None:
        """
        Method that allows to add a new cell to the lattice at the specified
        position. The cell is added to a new layer.

        Parameters
        ----------
        cell_to_add   : Cell
                        A cell to add to the current lattice, as object of
                        the 'Cell' subclasses
        position      : Tuple[float, float, float]
                        The X-Y-Z coordinates of the position where the cell
                        should be added, i.e. the cell center should be placed
                        at those coordinates. If the tuple is empty, the cell
                        is placed at the position indicated by its center
        """
        # Check that the added cell has the same type of the others
        if cell_to_add.cell_type != self.cells_type:
            raise AssertionError(
                "The cell to add has a geometry type "
                f"'{cell_to_add.cell_type}' which differs from the '"
                f"{self.cells_type}' one of the lattice.")
        # Initialize a new sub list identifying a new layer of cells
        self.layers.append([])
        # Add the cell to the newly created layer
        self.__add_cell_to_layer(cell_to_add, position, len(self.layers)-1)

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
        if position or not all(math.isclose(x, 0) for x in position):
            cell = cell.translate(position)
            self.__evaluate_no_rings(position)

        # Add the cell to the given layer
        self.layers[layer_indx].append(cell)
        # Add the cell to the list of lattice cells
        self.lattice_cells.append(cell)
        # Re-evaluate the types of the geometry and BCs
        self.__configure_lattice_types(self.cells_type,
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
            if math.isclose(self.lattice_cells[0].rotation, 0.0):
                # The distance is calculated using the cells characteristic
                # dimension along the X-axis
                ring_dist = 2*self.lattice_cells[0].figure.lx * self.rings_no
            elif math.isclose(self.lattice_cells[0].rotation, math.pi/2):
                # The distance is calculated using the cells characteristic
                # dimension along the Y-axis
                ring_dist = 2*self.lattice_cells[0].figure.ly * self.rings_no
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
            self.__configure_lattice_types(self.cells_type,
                                           len(self.lattice_cells))
            # Update the lattice in the current SALOME study
            self.show()
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

        # Update the type of BCs
        # FIXME how to handle the case where multiple BC types are associated
        # to the same geometry type?? A method for setting it should be added
        self.boundary_type = TYPEGEO_VS_BC[self.type_geo][0]
        # Update the lattice symmetry type
        self.symmetry_type = symmetry

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
                self.type_geo = LatticeGeometryType.RECTANGLE_SYM
            case SymmetryType.QUARTER:
                # ---------------------
                # QUARTER symmetry case
                # ---------------------
                # Build the face object for extracting the lattice symmetry
                face = self.__handle_rect_symmetry(b_box, 1)
                # Update the lattice type of geometry
                self.type_geo = LatticeGeometryType.RECTANGLE_SYM
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
                self.type_geo = LatticeGeometryType.RECTANGLE_EIGHT
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
        if not self.lattice_box:
            raise AssertionError(
                "The hexagonal lattice is not included within a box: the "
                f"requested '{symmetry}'symmetry operation cannot be "
                "applied.")
        match symmetry:
            case SymmetryType.THIRD:
                # ------------------
                # R120 symmetry case
                # ------------------
                # Build the shape identifying the symmetry
                face = self.__handle_third_symmetry()
                # Update the lattice type of geometry
                self.type_geo = LatticeGeometryType.R120
            case SymmetryType.SIXTH:
                # ------------------
                # SA60 symmetry case
                # ------------------
                # Build the shape identifying the symmetry
                face = self.__handle_sixth_symmetry()
                # Update the lattice type of geometry
                self.type_geo = LatticeGeometryType.SA60
            case SymmetryType.TWELFTH:
                # -----------------
                # S30 symmetry case
                # -----------------
                # Build the shape identifying the symmetry
                face = self.__build_triangle_on_lattice(0.0, 1/12)
                # Update the lattice type of geometry
                self.type_geo = LatticeGeometryType.S30
            case _:
                raise AssertionError(
                    f"The provided '{symmetry}' symmetry type is not "
                    "admitted for lattices with hexagonal cells.")
        # Perform the 'common' operation to extract the symmetry type
        # on the lattice and update the stored compound object
        self.lattice_symm = make_common(self.lattice_cmpd, face)

    def __build_vertices(self, u_params: List[float]) -> List[Any]:
        """
        FIXME it can be transformed into a function
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
            make_vertex_on_curve(self.lattice_box.figure.out_circle, u)
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
        if self.cells_rot == 0.0:
            # Build the triangle identifying the symmetry type
            return self.__build_triangle_on_lattice(1/12, 1 - 1/12)
        elif self.cells_rot == 90.0:
            # Build the triangle identifying the symmetry type
            return self.__build_triangle_on_lattice(2/3, 5/6)
        else:
            raise AssertionError("The cells rotation of "
                                  f"{self.cells_rot}° is not admitted.")

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
        if self.cells_rot == 0.0:
            points = self.__build_vertices([1/4, 7/12, 2/3])
        elif self.cells_rot == 90.0:
            points = self.__build_vertices([0, 1-1/6, 1-1/3])
        else:
            raise AssertionError(
                f"The cells rotation of {self.cells_rot}° is not admitted.")
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

    def __update_lattice_compounds(self, cells: List[Cell]) -> None:
        """
        Method that rebuilds the compound objects of the lattice grouping the
        given cells.
        If any box is declared, the cells are assembled with it.

        Parameters
        ----------
        cells : List[Cell]
            List of 'Cell' objects to build a compound object from
        """
        self.lattice_cmpd = make_partition(
            [c.face for c in cells], [], ShapeType.FACE)
        self.lattice_tech = self.lattice_cmpd
        # Build a GEOM compound storing the lattice unique edges of its cells
        self.lattice_edges: Any = make_partition(
            [self.lattice_cmpd], [], ShapeType.EDGE)
        # Store the original lattice compound (either technological or
        # sectorized geometry) without any box
        self.lattice_org = deepcopy(self.lattice_cmpd)

        # Handle the construction of the lattice container, if any
        self.__assemble_box()

    def __assemble_box(self) -> None:
        """
        Method that builds the lattice box, if it was not built yet, from the
        stored layers thicknesses. The container geometry depends on the
        type of geometry of the cells (i.e. rectangular or hexagonal).
        The box is then assembled with the lattice, eventually cutting the
        cells that are overlapped by it.
        """
        if not self.box_layers:
            return
        # If no lattice box has still been build, do it now
        if not self.lattice_box:
            self.build_lattice_box(self.box_layers)
        else:
            # Assemble the lattice with the box subface closest to the lattice
            # center
            self.lattice_cmpd = self.__assemble_box_with_lattice(
                self.__extract_inner_box(), self.lattice_cmpd)

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
                id = add_to_study_in_father(
                    self.lattice_cmpd, region.face, region.name)
                if not id:
                    raise RuntimeError(
                        f"Problem arose when adding the region {region.name}"
                        "to the SALOME study")
                # Display the region in the current view
                display_shape(id)
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

        # Translate each cell in the lattice
        translated_cells = self.__translate_cells(
            self.lattice_cells, new_pos, pre_center)
        # Update the list of lattice cells based on the translated ones
        self.lattice_cells = [cell for cell in translated_cells]

        # Translate each cell of each layer in the lattice
        translated_layers: List[List[Cell]] = []
        # Loop through all the layers to translate the cells
        for i, layer in enumerate(self.layers):
            translated_layers.append(
                self.__translate_cells(layer, new_pos, pre_center))
        # Update the data structure holding the translated cells for each layer
        self.layers = [
            [cell for cell in layer] for layer in translated_layers]

        # Translate the lattice compounds
        self.lattice_cmpd = make_translation(self.lattice_cmpd, transl_vect)
        self.lattice_tech = make_translation(self.lattice_tech, transl_vect)
        self.lattice_org = make_translation(self.lattice_org, transl_vect)
        # Rebuild the compound storing the lattice unique edges of its cells
        self.lattice_edges = make_partition(
            [self.lattice_cmpd], [], ShapeType.EDGE)

        # Translate the lattice box, if any
        if self.lattice_box:
            self.lattice_box = self.lattice_box.translate(new_pos)

        # Translate the lattice compound on which a symmetry operation has
        # been applied, if any
        if self.symmetry_type != SymmetryType.FULL:
            # Update the compound position relative to the shifted lattice
            # center and apply translation
            face_to_center = update_relative_pos(
                make_cdg(self.lattice_symm), pre_center, new_pos)
            self.lattice_symm = make_translation(
                self.lattice_symm, make_vector_from_points(
                    make_cdg(self.lattice_symm), make_vertex(face_to_center)))

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
            cell_to_center = update_relative_pos(
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
        self.lattice_cells = deepcopy(self.__assemble_layers())
        self.__update_lattice_compounds(self.lattice_cells)
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
        if self.lattice_box:
            box_subfaces = self.__extract_box_subfaces(cmpd)
            print("No. lattice box subfaces", len(box_subfaces))
            # Add the 'Region' objects for the lattice box subfaces
            self.regions.extend(
                self.__build_regions_for_subfaces(
                    len(self.regions), box_subfaces, [self.lattice_box]))
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
        if self.lattice_box:
            return make_common(cmpd, self.__extract_inner_box())
        return cmpd

    def add_ring_of_cells(self,
                          cell: Cell,
                          ring_indx: int,
                          layer_indx: int | None = None) -> None:
        """
        Method that adds a ring of cells of the same type. The cell is
        provided as an object of one of the 'Cell' subclasses, which
        is iteratively added at specific construction points identified by
        the given index of the lattice ring where the cells have to be added.
        Additionaly, the layer index can be provided: it specifies to which
        layer the ring of cells is added; if none, the cells are added to a
        new layer.

        **N.B.** Index 0, means the center cell of the lattice, whereas
        index greater than 0 indicates one of the rings of cells around
        this center cell.

        Parameters
        ----------
        cell : Cell
            The cell instance to be iteratively added in order to build a
            ring of cells
        ring_indx : int
            The index indicating where the ring of cells should be added
        layer_indx : int | None = None
            Identifying the layer index where cells are added; if None, cells
            are added to a new index
        """
        # Check that the cell to add has the same type of the ones already
        # present
        if cell.cell_type != self.cells_type:
            raise AssertionError(
                "The cell to add has a geometric surface type "
                f"'{cell.cell_type}' which differs from the '"
                f"{self.cells_type}' one of the lattice.")
        # If no layer index is provided, the ring of cells is added to a new
        # layer
        if layer_indx is None:
            layer_indx = len(self.layers)
            self.layers.append([])
        # Handle the addition differently according to the type of lattice
        # cells
        match self.cells_type:
            case CellType.RECT:
                # Build a construction rectangle where cartesian cells have
                # to be placed along its edges
                construction_fig = Rectangle(
                    get_point_coordinates(self.lattice_center),
                    2*cell.figure.ly * ring_indx,
                    2*cell.figure.lx * ring_indx)
                # Parameter for identifying the subdivision points
                sbdv_param = 2
            case CellType.HEX:
                # Build a construction hexagon where hexagonal cells have to
                # be placed along its edges. Use the 'ly' attribute as it
                # identifies the hexagon apothem
                construction_fig = Hexagon(
                    get_point_coordinates(self.lattice_center),
                    2*cell.figure.ly * ring_indx)
                # Rotate the construction hexagon if the lattice is rotated
                # by 90°
                if math.isclose(0.0,
                                math.degrees(self.lattice_cells[0].rotation)):
                    construction_fig.rotate(90.0)
                # Parameter for identifying the subdivision points
                sbdv_param = 1
            case _:
                raise RuntimeError(f"The {self.cells_type} is not handled.")

        # Loop through the borders of the construction figure and place a
        # cell on each subdivision point, given by the ring index and the
        # subdivision parameter
        for border in construction_fig.borders:
            # Build a subdivision point onto the figure borders
            for i in range(0, sbdv_param*ring_indx):
                p = make_vertex_on_curve(border, i/(sbdv_param*ring_indx))
                # Add the cell at the position of the subdivision point
                self.__add_cell_to_layer(
                    cell, get_point_coordinates(p), layer_indx)

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
            raise Exception("No properties have been provided")
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
                raise AssertionError(message)

        # Clear any previously set entry in the dictionary associating the
        # regions to the properties, as this method sets the properties for
        # all.
        self.lattice_box.tech_geom_props.clear()
        # Extract the lattice box subfaces
        box_subfaces = extract_sorted_sub_shapes(self.lattice_box.face,
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
            self.lattice_box.tech_geom_props[region] = prop_types
            # Update the index
            indx += 1

    def __extract_box_subfaces(self, lattice_cmpd: Any) -> List[Any]:
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

        Raises
        ------
        RuntimeError
            If no face objects are available after cutting the box face
            with the compound made from the lattice cells only

        Returns
        -------
        A list of face objects of the container layers and the areas between
        the cells and the box. Elements are sorted according to their distance
        from the lattice center and their perimeter in reverse order. The
        result is that the outmost layer of the box occupies the first
        position in the returned list and so on with the others.
        """
        box = self.lattice_box.face
        if self.symmetry_type != SymmetryType.FULL:
            # Get the shape of the symmetry
            shape = make_face(build_compound_borders(self.lattice_symm))
            # Extract the common part between the box and the symmetry shape
            box = make_common(self.lattice_box.face, shape)
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
        shape = get_selected_object()
        if not shape:
            raise Exception("Please, select a single region whose data to "
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

    def update_geometry(self, geo_type: GeometryType) -> None:
        """
        Method for updating the lattice geometry with the given compound
        that has been externally manipulated due to some operation (e.g.
        cut, scaling, etc.).
        The SALOME ID of the compound is provided, given which the GEOM
        object is retrieved.
        A common operation between the stored lattice compound and the
        new one is performed. Different operations are performed depending
        on the result:

        - If None, it means that the lattice position has changed. However,
          if the geometrical properties are different, a warning is raised
          as the center of the new geometry cannot be univocally determined.
          In this case, the method `update_geometry_from_center` should be
          called. If no changes in the geometrical properties are present,
          the lattice (and all the related attributes) is translated in the
          new position. A rotation is performed as well, if needed.
        - If not None, it means that the two compounds shares all the lattice
          or a part of it. If the geometric properties of the common result
          and the new compound are different, the same check is performed on
          the previous and current lattice:

          - If different, nothing can be said about its center: a warning
            is raised notifying the need to call the method
            `update_geometry_from_center`.
          - If equal, the new compound is just placed in a different
            position: the current one is translated and a rotation is
            performed as well, if needed.

          If the common result and the new compound are the same, then a
          rotation is performed, if needed. A comparison on the previous and
          actual lattice area is performed to check if any symmetry could be
          applied. In any case, if the areas are different, the lattice
          compound is simply updated with the new one.

        Parameters
        ----------
        geo_type  : GeometryType
            The lattice type of geometry to update, as value of the
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
        # Update the lattice compounds with the given shape
        self.lattice_cmpd = shape
        self.lattice_tech = shape
        self.lattice_edges: Any = make_partition(
            [self.lattice_cmpd], [], ShapeType.EDGE)

        # Update the lattice in the current SALOME study
        self.show()

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
            # FIXME put into log
            print("No rotation is performed as the given angle is 0.0°")
            return
        # Convert the rotation angle in radians
        rotation = math.radians(angle)

        # Build the Z-axis of rotation
        x, y, _ = get_point_coordinates(self.lattice_center)
        z_axis = make_vector_from_points(
            self.lattice_center, make_vertex((x, y, 1)))
        # Rotate the lattice compounds
        self.lattice_cmpd = make_rotation(self.lattice_cmpd, z_axis, rotation)
        self.lattice_tech = make_rotation(self.lattice_tech, z_axis, rotation)
        self.lattice_org = make_rotation(self.lattice_org, z_axis, rotation)
        # Re-build the lattice edges from the compound
        self.lattice_edges = make_partition(
            [self.lattice_cmpd], [], ShapeType.EDGE)
        # Rotate each cell of the lattice
        for cell in self.lattice_cells:
            cell.rotate_from_axis(angle, z_axis)
        # Rotate each cell of each layer
        for layer in self.layers:
            for cell in layer:
                cell.rotate_from_axis(angle, z_axis)
        # Rotate the lattice box, if any
        if self.lattice_box:
            self.lattice_box.rotate_from_axis(angle, z_axis)
        # Rotate the lattice compound on which a symmetry operation has been
        # applied, if any
        if self.symmetry_type != SymmetryType.FULL:
            self.lattice_symm = make_rotation(
                self.lattice_symm, z_axis, rotation)

        # Show the new lattice compound in the current SALOME study
        self.show()

    def sectorize(self, sect_indexes: Tuple[int, int]) -> None:
        """
        Method for applying the sectorization to all the cells in the
        lattice at once.

        Parameters
        ----------
        sect_indexes  : List[Tuple[int, int]]
            Containing the sectorization index, 'isect', and 'jsect', the
            number of embedded tubes not sectorized (0 < jsect < no_tubes).
            These values are provided for each different type of cells in
            the lattice, i.e. having different number of tubes.
        """
        # Group similar cells according to the number of inner tubes
        grouped_objects = {}
        # FIXME to remove as receiving a list of tuples -> only for demo
        sect_indexes = [sect_indexes]

        # Group the objects by the number attribute
        for cell in self.lattice_cells:
            if len(cell.inner_circles) not in grouped_objects:
                grouped_objects[len(cell.inner_circles)] = []
            grouped_objects[len(cell.inner_circles)].append(cell)

        # Check for a mismatch between the number of cell types and the number
        # of sectorization indices
        if len(grouped_objects.keys()) != len(sect_indexes):
            raise AssertionError("Mismatch between the number of cell types "
                                 f"({len(grouped_objects.keys())}) and the "
                                 f"number of sectorization indices ("
                                 f"{len(sect_indexes)})")

        # FIXME this implementation must be applied to all cells grouped by
        # the number of inner circles

        # Handle the sectorization operation according to the type of cells
        # in the lattice
        match self.cells_type:
            case CellType.RECT:
                raise Exception("Sectorization not yet implemented for "
                                "lattices made by rectangular cells.")
            case CellType.HEX:
                # FIXME only one couple of indices is handled for now
                isect = sect_indexes[0][0]
                jsect = sect_indexes[0][1]
                # Check if isect = 0, -1
                if isect != 0 and isect != -1:
                    raise AssertionError(
                        "Wrong 'isect' index: it must be either 0 or -1 "
                        "for hexagonal cells.")
                # Check if 0 < jsect < no_tubes
                if jsect < 0 or jsect > len(
                    self.lattice_cells[0].inner_circles):
                    raise AssertionError(
                        "Wrong 'jsect' index: it must be between 0 and "
                        f"{len(self.lattice_cells[0].inner_circles)}")
                # Loop through all the cells and apply the sectorization
                for cell in self.lattice_cells:
                    sector_indx = [6] * (len(
                        self.lattice_cells[0].inner_circles) + 1)
                    for j in range(jsect):
                        sector_indx[j] = 1
                    angles = [0]*len(sector_indx)
                    cell.sectorize(sector_indx, angles)
            case _:
                raise Exception(f"The cell type {self.cells_type} is not "
                                "valid for a sectorization operation.")
        # Re-build the lattice compound and regions
        self.__build_lattice(GeometryType.SECTORIZED)
        # Update the lattice in the current SALOME study
        self.show(geometry_type_to_show=GeometryType.SECTORIZED)

    def set_region_property(
            self, property_type: PropertyType, value: str) -> None:
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
        """
        # Check which cell geometry type is currently shown; if different
        # from the TECHNOLOGICAL one, raise an exception
        if self.displayed_geom != GeometryType.TECHNOLOGICAL:
            raise RuntimeError(
                f"Currently showing the {self.displayed_geom.name} type of "
                "geometry. To set lattice regions properties, show the '"
                "TECHNOLOGICAL' geometry first.")
        # Extract the geometrical object currently selected in the current
        # SALOME study
        shape = get_selected_object()
        if not shape:
            raise RuntimeError("Please, select a single region to assign "
                               "a property to.")

        # Point for idendifying the shape in the geometry
        point = make_vertex_inside_face(shape)
        # Get the compound corresponding to the technological geometry of the
        # cells and the currently applied symmetry
        cmpd = self.__get_compound_from_type(GeometryType.TECHNOLOGICAL,
                                             self.lattice_cells)
        # Store all the cells to look for the region to update
        cells: List[Cell] = []
        # If a box is present, add it to the list of cells
        if self.lattice_box:
            # Copy the lattice box cell and update its face by cutting out the
            # cells; this to prevent to select the box region instead of a
            # cell as the box is not empty
            lattice_box = deepcopy(self.lattice_box)
            lattice_box.face = make_compound(
                self.__extract_box_subfaces(cmpd))
            cells.append(lattice_box)
        # Add all the lattice cells
        cells.extend(self.lattice_cells)

        # Get the region that corresponds to the given shape
        for i, cell in enumerate(cells):
            # Continue with another cell if the shape is not contained within
            # the cell's area (reference point inside cell's face and common
            # area between cell and shape is the shape itself)
            if (not is_point_inside_shape(point, cell.face)):
                continue
            # Search for the region among the faces stored in the cell
            # dictionary of regions VS properties
            for region in cell.tech_geom_props:
                if (is_point_inside_shape(point, region)):
                    cell.tech_geom_props[region][property_type] = value
                    break
            else:
                continue
            # Update the cell properties of the lattice box, if it has
            # been modified (first position in the list of cells)
            if i == 0:
                for zone in self.lattice_box.tech_geom_props:
                    if is_point_inside_shape(point, zone):
                        self.lattice_box.tech_geom_props[zone][
                            property_type] = value
                        break
            break
        # Raise an exception if no region has been found
        else:
            raise RuntimeError(
                "No cell region could be found for the selected shape.")
        # Update the lattice in the SALOME viewer
        self.show(property_type)


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
