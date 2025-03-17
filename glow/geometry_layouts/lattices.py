"""
Module containing classes providing the means for creating lattices from base
cells types in SALOME.
"""
import math

from copy import deepcopy
from typing import Any, Dict, List, Tuple, Union
from warnings import WarningMessage

from glow.geometry_layouts.cells import GenericCell, HexCell, RectCell, Region
from glow.geometry_layouts.geometries import GenericSurface, Hexagon, \
    Rectangle, build_hexagon, build_rectangle
from glow.interface.geom_interface import ShapeType, add_to_study, \
    add_to_study_in_father, clear_view, display_shape, \
    extract_sorted_sub_shapes, extract_sub_shapes, fuse_edges_in_wire, \
    get_basic_properties, get_bounding_box, get_closed_free_boundary, \
    get_inertia_matrix, get_min_distance, get_point_coordinates, \
    get_selected_object, is_point_inside_shape, make_cdg, make_common, \
    make_compound, make_cut, make_edge, make_face, make_fuse, \
    make_partition, make_rotation, make_translation, make_vector, \
    make_vector_from_points, make_vertex, make_vertex_inside_face, \
    make_vertex_on_curve, remove_from_study, set_color_face, \
    update_salome_study
from glow.generator.support import GeometryType, LatticeGeometryType, \
    PropertyType, SymmetryType, BoundaryType, CellType, TYPEGEO_VS_BC, \
    generate_unique_random_colors


class Lattice():
    """
    Class that represents a lattice made by a group of cells.

    Parameters
    ----------
    cells : List[GenericCell]
            The list of cells that constitute the lattice, as objects
            of the 'GenericCell' subclasses
    name  : str = "Lattice"
            The lattice name in the current SALOME study

    Attributes
    ----------
    lattice_cells     : List[GenericCell]
                        The list of cells that constitute the lattice, as
                        objects of the 'GenericCell' subclasses
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
                 cells: List[GenericCell],
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
        self.lattice_cells: List[GenericCell] = deepcopy(cells)
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
        self.lattice_box: Union[GenericCell, None] = None
        self.lattice_symm: Union[Any, None] = None
        self.lattice_tech: Union[Any, None] = None
        self.regions: List[Region] = []
        self.type_vs_property_vs_color: Dict[
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
        self.show(update_view=False)

    def __evaluate_cells_rotation(self, cells: List[GenericCell]) -> float:
        """
        Method that checks if the rotation angle of the given cells is the
        same for all the given ones; if not, an exception is raised.
        A check on the validity of the cell rotation according to the admitted
        values is also present; an exception is raised if the value is not
        valid.
        The method returs the rotation angle of the cells, if all checks pass.

        Parameters
        ----------
        cells : List[GenericCell]
                The list of 'GenericCell' subclasses, representing the cells
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
            raise AssertionError("The given cells do not share the same "
                                 "rotation angle.")
        # Check if the rotation angle of the cells is valid
        if not cell_rot in self.VALID_CELLS_ANGLES:
            raise AssertionError("Cells can only be rotated by any of the "
                                 f"{self.VALID_CELLS_ANGLES}° angles")
        # Return the rotation angle of all the cells
        return cell_rot

    def __evaluate_lattice_center(
            self,
            cells: List[GenericCell],
            center: Union[Tuple[float, float, float], None]) -> Any:
        """
        Method that evaluates the lattice center during the initialization:
        if no center is provided at class instantiation, the XYZ origin is
        returned. In the opposite case, a check is performed to identify if
        there is a cell whose center coincides with the specified one.
        In case none is found, an exception is raised.

        Parameters
        ----------
        cells   : List[GenericCell]
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
        cell_faces = []
        match geo_type:
            case GeometryType.TECHNOLOGICAL:
                cell_faces = [cell.face for cell in self.lattice_cells]
                # Store the lattice coming from the cells technological
                # geometries
                self.lattice_tech = make_compound(cell_faces)
            case GeometryType.SECTORIZED:
                for cell in self.lattice_cells:
                    if not cell.sectorized_face:
                        print(f"The {cell.name} cell has not been sectorized: "
                              "the non sectorized face will be used.")
                        cell_faces.append(cell.face)
                        continue
                    cell_faces.append(cell.sectorized_face)
            case _:
                raise ValueError(f"{geo_type}: unhandled type.")
        self.lattice_cmpd = make_compound(cell_faces)

        # Build a GEOM compound storing the lattice unique edges of its cells
        self.lattice_edges: Any = make_partition(
            [self.lattice_cmpd], [], ShapeType.EDGE)
        # Store the original lattice compound (either technological or
        # sectorized geometry) without any box
        self.lattice_org = deepcopy(self.lattice_cmpd)

        # Handle the construction of the lattice container, if any is required
        if self.box_layers:
            # If no lattice box has still been build, do it now
            if not self.lattice_box:
                # In case, the first layer is negative, remove the second
                # value as inserted it on purpose
                if self.box_layers[0] < 0:
                    self.box_layers.pop(1)
                self.build_lattice_box(self.box_layers)
            else:
                box_layers = []
                # Extract the box layers only (i.e. the ones centered in the
                # lattice center)
                for subface in self.lattice_box.extract_subfaces():
                    if get_min_distance(make_cdg(subface),
                                        self.lattice_center) <= 1e-6:
                        box_layers.append(subface)
                # Sort the layers surfaces and extract the inmost one
                inner_layer = sorted(
                    box_layers,
                    key=lambda item: get_min_distance(item,
                                                      self.lattice_center))[0]
                # Build a face for the inmost wire of the extracted layer
                inner_box = make_face(
                    extract_sorted_sub_shapes(inner_layer, ShapeType.WIRE)[0]
                )
                # Assemble the lattice with the box
                self.lattice_cmpd = self.__assemble_box_with_lattice(
                    inner_box, self.lattice_cmpd)

    def __add_cell_regions(self, cell: GenericCell):
        """
        """
        # Get the current number of regions
        indx = len(self.regions)
        # Extract the lattice subfaces, either from the whole compound or
        # the one after applying the symmetry
        cmpd = self.lattice_cmpd if self.symmetry_type == SymmetryType.FULL \
            else self.lattice_symm
        subfaces = extract_sub_shapes(cmpd, ShapeType.FACE)

        # Loop through all the cell regions
        for region in cell.regions:
            # Loop through all the subfaces
            for subface in subfaces:
                # Get a vertex inside the subface
                subface_point = make_vertex_inside_face(subface)
                # Check if the lattice subface overlaps the cell one
                if get_min_distance(region.face, subface_point) < 1e-7:
                    # Update the index
                    indx += 1
                    # Build a region for the lattice
                    self.regions.append(Region(
                        face=subface,
                        inner_point=subface_point,
                        name=f"Region {indx}",
                        properties=region.properties
                    ))
                    break

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

    def __cut_box_with_cells(self, inner_box: Any) -> None:
        """
        Method that extracts from the whole lattice box face only the part
        not belonging to the cells; this is performed by cut boolean operation
        that removes the cells shape from the lattice box face.

        Parameters
        ----------
        inner_box : Any
            Face object identifying the inner region of the lattice box
        """
        # Extract the common part between the inner box and the original
        # lattice
        cut_lattice = make_common(inner_box, self.lattice_org)
        # Update the lattice box face by cutting out the cells
        self.lattice_box.face = make_cut(self.lattice_box.face, cut_lattice)
        self.lattice_box.figure.update_from_face(self.lattice_box.face)

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

        # Store the thickness of the first box layer: if positive, a zero
        # value is stored at posizion 0 in the list, if negative (box
        # overlapping the cells), the opposite value is stored at position
        # 1 in the list
        self.box_layers = list(box_thick)
        if box_thick[0] > 0:
            self.box_layers.insert(0, 0.0)
        else:
            self.box_layers.insert(1, -box_thick[0])

        # Build the outer container for the lattice according to the cells
        # geometry
        self.__build_lattice_box_type()

        # Extract the lattice box faces
        box_faces = extract_sub_shapes(self.lattice_box.face, ShapeType.FACE)
        # Get the inner box face by sorting the faces according to the
        # distance from the lattice center
        inner_box = sorted(
            box_faces,
            key=lambda item: get_min_distance(self.lattice_center, item))[0]
        # Cut out the cells from the lattice box
        self.__cut_box_with_cells(inner_box)

        # Update the lattice box face and its underlying figure
        # FIXME to put all this part into a method for the 'GenericCell' class
        subfaces = self.lattice_box.extract_subfaces()
        if not subfaces:
            subfaces = [self.lattice_box.face]
        # Initialize the regions properties dictionary
        self.lattice_box.tech_geom_props = {
                region: {} for region in subfaces}

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

        # FIXME If the lattice has a symmetry, regions are built because of
        # the previous call to the 'apply_symmetry()' method, then they are
        # removed from the view. Do not build regions when calling the
        # 'apply_symmetry()' method but at the end of this method.
        # Update the lattice representation in the current SALOME study
        # without showing it the viewer, i.e. only the object browser is
        # updated.
        self.show(update_view=False)

    def __build_lattice_box_type(self) -> None:
        """
        Method that builds the lattice box as an instance of the 'GenericCell'
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
                    box_surfaces.append(
                        build_rectangle(height, width, center))
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
        self.lattice_box.face = box_face
        self.lattice_box.figure.update_from_face(box_face)
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
                 cell_to_add: GenericCell,
                 position: Tuple[float, float, float],
                 update_study: bool = True) -> None:
        """
        Method that allows to add a new cell to the current lattice
        representation.

        Parameters
        ----------
        cell_to_add   : GenericCell
                        A cell to add to the current lattice, as object of
                        the 'GenericCell' subclasses
        position      : Tuple[float, float, float]
                        The X-Y-Z coordinates of the position where the cell
                        should be added, i.e. the cell center should be placed
                        at those coordinates
        update_study  : bool
                        Flag stating if the study should be updated both in
                        the viewer and in the object browser; it defaults to
                        True
        """
        # Check that the added cell has the same type of the others
        # cell_type = type(cell_to_add)
        # cells_type = type(self.lattice_cells[0])
        if cell_to_add.cell_type != self.lattice_cells[0].cell_type:
            raise AssertionError(
                "The cell to add has a geometry type "
                f"'{cell_to_add.cell_type}' which differs from the '"
                f"{self.lattice_cells[0].cell_type}' one of the lattice.")
        print("Adding the cell to position", position)
        # Move the cell in the given position
        cell = cell_to_add.translate(position)
        # Check whether the position of the cell to add has already been
        # occupied by another cell. The 'common' operation between the whole
        # lattice and the face of the cell to add is used to determine this.
        common = make_common(self.lattice_cmpd, cell.face)
        # Extract a list of the faces resulting from the 'common' operation
        common_face = extract_sub_shapes(common, ShapeType.FACE)
        # If the list is not empty, it means the cell to add intersect the
        # lattice; hence, an exception is raised.
        if common_face:
            # Add the result of the 'common' operation for 'logging' purposes
            add_to_study(common, "Common face")
            update_salome_study()
            raise Exception(f"The cell, whose position is {position}, cannot "
                            "be added to the lattice since this cell overlaps "
                            "the lattice.")

        # Add the cell to the list of lattice cells
        self.lattice_cells.append(cell)
        # Re-evaluate the types of the geometry and BCs
        self.__configure_lattice_types(self.lattice_cells[0].cell_type,
                                       len(self.lattice_cells))
        # Re-evaluate the number of cells rings in the lattice
        self.__evaluate_no_rings(position)
        # Rebuild the lattice GEOM compound objects as a new cell is present
        self.__build_lattice()
        # print("2) DISTANZA O-outmost cell", self.distance)
        # FIXME this should be the method to call so to add only a region each
        # time a cell is added to the lattice
        # self.__add_cell_regions(cell)

        # FIXME add to log
        print("Number of lattice cells is", len(self.lattice_cells))
        print("Number of lattice rings is", self.rings_no)

        # Update the lattice representation in the current SALOME study, if
        # required
        if update_study:
            self.show()

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
        # print("1) DISTANCE O-outmost cell", self.distance, "O-cell", origin_distance)
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
                self.type_geo = LatticeGeometryType.SYMMETRIES_TWO
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
            case SymmetryType.SIXTH:
                # ------------------
                # SA60 symmetry case
                # ------------------
                # Extract the lattice portion corresponding to the symmetry
                self.__handle_hex_symmetry(6)
                 # Update the lattice type of geometry
                self.type_geo = LatticeGeometryType.SA60
            case SymmetryType.TWELFTH:
                # -----------------
                # S30 symmetry case
                # -----------------
                # Extract the lattice portion corresponding to the symmetry
                self.__handle_hex_symmetry(12)
                # Update the lattice type of geometry
                self.type_geo = LatticeGeometryType.SYMMETRIES_TWO
            case _:
                raise AssertionError(
                    f"The provided '{symmetry}' symmetry type is not "
                    "admitted for lattices with hexagonal cells.")

    def __handle_hex_symmetry(self, param: float) -> None:
        """
        Method that handles the symmetry application for hexagonal cells
        lattices. Given the position of two of the vertices of a triangle,
        with the third being the lattice center, a triangular shape is built.
        This shape identifies the type of symmetry to apply.
        A 'common' operation is performed between the shape and the lattice
        compound so to extract the part of interest for the symmetry.

        Parameters
        ----------
        param : float
            Parameter for identifying the triangle vertices (6 for 'SIXTH',
            12 for 'TWELFTH' symmetry type)
        """
        if self.cells_rot == 0.0:
            # Build the triangle identifying the symmetry type
            face = self.__build_triangle_on_lattice(0.75 - 1/(2*param),
                                                    0.75 + 1/(2*param))
        elif self.cells_rot == 90.0:
            # Build the triangle identifying the symmetry type
            face = self.__build_triangle_on_lattice(0.0, 1/param)
        else:
            raise AssertionError("The cells rotation of "
                                  f"{self.cells_rot}° is not admitted.")
        # Perform the 'common' operation to extract the symmetry type
        # on the lattice and update the stored compound object
        self.lattice_symm = make_common(self.lattice_cmpd, face)

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
        # Build the rectangle face object
        rect.build_face()
        # Return the face
        return rect.face

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
        points = [self.lattice_center,
                  make_vertex_on_curve(
                      self.lattice_box.figure.out_circle, vert2_u),
                  make_vertex_on_curve(
                      self.lattice_box.figure.out_circle, vert3_u)]
        # Return the face object built on the vertices
        return self.__build_face_from_vertices(points)

    def __build_face_from_vertices(self, vertices: List[Any]) -> Any:
        """
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

    def show(self,
             property_type_to_show: Union[PropertyType, None] = None,
             geometry_type_to_show: GeometryType = GeometryType.TECHNOLOGICAL,
             update_view: bool = True) -> None:
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
        property_type_to_show : Union[PropertyType, None]
            Either a value of the 'PropertyType' enumeration, indicating
            the property by which the regions are coloured, or None.
        geometry_type_to_show : GeometryType = GeometryType.TECHNOLOGICAL
            The type of geometry to show: regions are shown accordingly
            with the selected type
        update_view           : bool
            Flag stating if the viewer should be updated by displaying the
            lattice regions and edges
        """
        # Erase all objects from the current view
        clear_view()

        # If already present in the current SALOME study, remove the GEOM
        # compound object
        if self.lattice_entry_id:
            remove_from_study(self.lattice_entry_id)

        # Add the GEOM compound objects to the current SALOME study only if
        # the geometry type to display is not changed
        if self.displayed_geom == geometry_type_to_show:
            self.lattice_entry_id = add_to_study(self.lattice_cmpd, self.name)
            # self.lattice_edges_id = add_to_study_in_father(
            #     self.lattice_cmpd, self.lattice_edges, self.name + "_edges")
            # # Show the edges of the complete lattice, if needed
            # if update_view and self.symmetry_type != SymmetryType.FULL:
            #     display_shape(self.lattice_edges_id)
            # # Show the lattice center point
            # add_to_study_in_father(
            #     self.lattice_cmpd, self.lattice_center, self.name + "_center")

        # Build regions and color them, if needed. If an exception is raised,
        # is caught and re-raised
        try:
            # If there is no need to update the study, return without building
            # the lattice regions
            if not update_view:
                return
            # Build the cell regions, as 'Region' objects, according to the
            # indicated geometry type
            self.build_regions(geometry_type_to_show)
            # If the geometry to show is different from the previously
            # displayed one, the lattice compound is added to the current
            # SALOME study
            if self.displayed_geom != geometry_type_to_show:
                self.lattice_entry_id = add_to_study(self.lattice_cmpd,
                                                     self.name)
            # Assign the same color to all the regions having the same property
            # value, if any has been specified to show
            self.__associate_colors_to_regions(property_type_to_show)
            # Add all the lattice regions to the study, each with an assigned
            # color
            for region in self.regions:
                # Set the region color in the viewer
                set_color_face(region.face, region.color)
                # Add the lattice region to the study
                id = add_to_study_in_father(self.lattice_cmpd,
                                            region.face,
                                            region.name)
                if not id:
                    raise RuntimeError("Problem arose when adding the region "
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
        Method that assign the same color to all the regions having the same
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
        translated_cells = []
        for cell in self.lattice_cells:
            # Update the cell position relative to the shifted lattice center
            # and apply translation
            cell_to_center = update_shape_position_on_lattice(
                cell.figure.o, pre_center, new_pos)
            translated_cells.append(cell.translate(cell_to_center))
        # Update the list of lattice cells based on the translated ones
        self.lattice_cells = [cell for cell in translated_cells]

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
            face_to_center = update_shape_position_on_lattice(
                make_cdg(self.lattice_symm), pre_center, new_pos)
            self.lattice_symm = make_translation(
                self.lattice_symm, make_vector_from_points(
                    make_cdg(self.lattice_symm), make_vertex(face_to_center)))

        # Re-build the lattice 'Region' objects
        # self.build_regions()

        # Show the new lattice compound in the current SALOME study
        self.show()

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
        # Get the lattice compound object, given the geometry type and the
        # current applied symmetry
        cmpd = self.__get_compound_from_type(geo_type)
        # Extract the lattice subfaces
        subfaces = extract_sub_shapes(cmpd, ShapeType.FACE)

        # Re-initialize the lattice regions
        self.regions.clear()
        # Extract the 'Region' objects for the lattice cells only
        self.regions.extend(
                self.__build_regions_for_subfaces(
                    0, subfaces, self.lattice_cells))
        # Handle the construction of 'Region' objects for the lattice box
        # subfaces, if any
        if self.lattice_box:
            box_subfaces = self.__extract_box_subfaces()
            # Add the 'Region' objects for the lattice box subfaces
            self.regions.extend(
                self.__build_regions_for_subfaces(
                    len(self.regions), box_subfaces, [self.lattice_box]))
            # Update the lattice found subfaces
            subfaces += box_subfaces

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
            lattice_cells: List[GenericCell]
            ) -> List[Region]:
        """
        Method that builds a list of 'Region' objects for each of the given
        face objects.
        The names of the regions are set with an increasing index with a
        given starting value.
        Properties to associate to each 'Region' object are taken from each
        given cell dictionary.

        Parameters
        ----------
        indx : int
            The starting value for the index used when assigning the region
            name
        faces : List[Any]
            A list of face objects for which 'Region' objects are built
        lattice_cells : List[GenericCell]
            A list of 'GenericCell' objects representing the lattice cells

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
                                properties=properties
                            ))
                        found = True
                        break
                # If a region has been found for the current cell, a new
                # lattice subface is evaluated
                if found:
                    break
        # Return the list
        return regions

    def __get_compound_from_type(self, geo_type: GeometryType) -> Any:
        """
        Method that returns the compound object of the lattice according to
        the given type of geometry, as value of the 'GeometryType'
        enumeration:
        - TECHNOLOGICAL: in case of a 'FULL' symmetry, the corresponding
          instance attribute is returned, otherwise, all the subfaces of
          lattice to which the symmetry has been applied are fused to extract
          the characteristic shape of the symmetry. A 'common' operation with
          the whole lattice (technological geometry) provides the compound to
          return.
        - SECTORIZED: according to the symmetry, the corresponding instance
          attribute is returned.

        Parameters
        ----------
        geo_type  : GeometryType
            Value of the 'GeometryType' enumeration indicating the lattice
            compound type to return

        Returns
        -------
        A compound object representing the lattice (either full or not) built
        from the technological geometry or the one for the calculations (i.e.
        the sectorized one).
        """
        # Rebuild the lattice compound using the sectorized cells, where
        # defined, or the technological ones, if the given geometry type is
        # changed
        if self.displayed_geom != geo_type:
            self.__build_lattice(geo_type)
        # Handle the different type of geometry to build regions from
        match geo_type:
            case GeometryType.TECHNOLOGICAL:
                if self.symmetry_type == SymmetryType.FULL:
                    # If any lattice box is present, cut it out from the
                    # lattice compound object (the one with only the cells
                    # and no box)
                    if self.lattice_box:
                        return make_cut(self.lattice_tech,
                                        self.lattice_box.face)
                    # Return the lattice without box
                    return self.lattice_org
                else:
                    # Extract the shape of the symmetry by extracting its
                    # borders and building a face object
                    shape = make_face(self.__build_borders(self.lattice_symm))
                    # If any lattice box is present, cut it out from the
                    # lattice compound object (the one with only the cells
                    # and no box). Then, perform a common operation of the
                    # result with the symmetry shape.
                    if self.lattice_box:
                        cut = make_cut(self.lattice_org,
                                       self.lattice_box.face)
                        return make_common(cut, shape)
                    # Return the common part between the lattice and the shape
                    # of the symmetry
                    return make_common(self.lattice_tech, shape)
            case GeometryType.SECTORIZED:
                if self.symmetry_type == SymmetryType.FULL:
                    # If any lattice box is present, cut it out from the
                    # lattice compound object (the one with only the cells
                    # and no box)
                    if self.lattice_box:
                        return make_cut(self.lattice_cmpd,
                                        self.lattice_box.face)
                    return self.lattice_cmpd
                else:
                    # Extract the shape of the symmetry by extracting its
                    # borders and building a face object
                    shape = make_face(self.__build_borders(self.lattice_symm))
                    # If any lattice box is present, cut it out from the
                    # lattice compound object (the one with only the cells
                    # and no box). Then, perform a common operation of the
                    # result with the symmetry shape.
                    if self.lattice_box:
                        cut = make_cut(self.lattice_org,
                                       self.lattice_box.face)
                        return make_common(cut, shape)
                    # Return the common part between the lattice and the shape
                    # of the symmetry
                    return make_common(self.lattice_cmpd, shape)
            case _:
                raise ValueError(f"{geo_type}: unhandled type of geometry.")

    def __build_borders(self, lattice_cmpd: Any) -> List[Any]:
        # FIXME To transform into a function to put into a support module?
        """
        Method that extracts the borders of the lattice from the GEOM
        compound object storing the faces and edges of the lattice.

        Parameters
        ----------
        lattice_cmpd  : Any
                        The lattice GEOM compound object

        Returns
        -------
        The list of GEOM edge objects that represent the lattice borders.
        """
        # Extract a list of closed boundaries from the lattice compound
        closed_boundaries = get_closed_free_boundary(lattice_cmpd)
        # Handle the case where more than a closed wire is extracted from
        # the lattice compound
        if len(closed_boundaries) > 1:
            # Build a face for each of the extracted closed wires after
            # fusing adjacent edges
            shapes =  []
            for wire in closed_boundaries:
                wire_mod = fuse_edges_in_wire(wire)
                shapes.append(make_face(wire_mod))
            # Fuse all the faces into a single shape
            shape = make_fuse(shapes)
            # Return the edges of the fused shape
            return extract_sorted_sub_shapes(shape, ShapeType.EDGE)

        # Suppress vertices internal to the edges of the wire
        borders_wire = fuse_edges_in_wire(closed_boundaries[0])
        # Extract the edge objects of the border wire
        return extract_sorted_sub_shapes(borders_wire, ShapeType.EDGE)

    def add_ring_of_cells(self, cell: GenericCell, ring_indx: int) -> None:
        """
        Method that adds a ring of cells of the same type. The cell is
        provided as an object of one of the 'GenericCell' subclasses, which
        is iteratively added at specific construction points identified by
        the given index of the lattice ring where the cells have to be added.

        **N.B.** Index 0, means the center cell of the lattice, whereas
        index greater than 0 indicates one of the rings of cells around
        this center cell.

        Parameters
        ----------
        cell      : GenericCell
                    The cell instance to be iteratively added in order to
                    build a ring of cells
        ring_indx : int
                    The index indicating where the ring of cells should be
                    added
        """
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
                self.add_cell(cell, get_point_coordinates(p), False)

    def add_rings_of_cells(self, cell: GenericCell, no_rings: int) -> None:
        """
        Method that adds to the lattice several rings of cells of the same
        type. The cell is provided as an object of one of the 'GenericCell'
        subclasses, which is iteratively added, for each index of rings, at
        specific construction points.
        The rings of cells are added starting from the current maximum value
        of rings for the lattice.

        Parameters
        ----------
        cell      : GenericCell
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
        # Loop through the ring indices starting from the current number of
        # rings of cells
        for i_ring in range(self.rings_no+1, no_rings+1):
            # Add a ring of cells at the current index
            self.add_ring_of_cells(cell, i_ring)

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

        # Re-build all the lattice regions so to assign the box properties
        # to the corresponding ones
        # FIXME: if a box is already present, only the corresponding regions
        # should be updated by modifying the property attribute, so to avoid
        # re-building all the regions each time
        # self.build_regions()

    def __extract_box_subfaces(self) -> List[Any]:
        """
        Method that extracts the lattice box subfaces, i.e. the box layers
        and the areas between the cells and the container layers, if any
        (as depending on the cells geometry).

        Returns
        -------
        A list of face objects of the container layers and the areas between
        the cells and the box. Elements are sorted according to their distance
        from the lattice center and their perimeter in reverse order. The
        result is that the outmost layer of the box occupies the first
        position in the returned list and so on.
        """
        # Extract the lattice box subfaces: if any symmetry is applied, they
        # result from the common operation with the compound storing the cut
        # lattice
        box = self.lattice_box.face
        if self.symmetry_type != SymmetryType.FULL:
            box = make_common(self.lattice_symm, box)
            if not box:
                raise Exception(
                    "Error in extracting the box part of the lattice with "
                    f"applied {self.symmetry_type.name} symmetry type.")

        # Extract the lattice box subfaces
        box_subfaces = extract_sorted_sub_shapes(box, ShapeType.FACE)

        # Sort the box regions according to their distance from the lattice
        # center and their perimeter in reverse order: this so that the
        # outmost layer of the box is the first element and so on with the
        # others.
        return sorted(
            box_subfaces,
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

        # Get the new compound CDG
        cdg = make_cdg(shape)

        # Perform a 'common' operation between the previous lattice compound
        # and the modified one. If the operation returns the same element, it
        # means no changes in the lattice center is needed.
        common = make_common(self.lattice_cmpd, shape)
        if common:
            # Check if they have the same geometric properties
            if get_basic_properties(common) != get_basic_properties(shape):
                print("The compound is different from the common part!!")
                # Verify if the geometric properties have changed
                if (get_basic_properties(self.lattice_cmpd) !=
                    get_basic_properties(shape)):
                    print("The new compound is different from the previous "
                          "one!!")
                    # Raise a warning, since no center can be clearly
                    # identified
                    raise WarningMessage(
                        "Warning: the updated compound appears to be placed "
                        "in a different position. The area is different as "
                        "well. Please, call the "
                        "'update_geometry_from_center' method to specify "
                        "also the new center.")
                # Same area but different location: translate the lattice
                self.translate(get_point_coordinates(cdg))
                # Check if any rotation is present
                delta = compare_compounds_rotation(self.lattice_cmpd, shape)
                # Rotate the lattice, if required
                self.rotate(delta)

            # The new compound is part of the previous compound
            # Check if any rotation is present
            delta = compare_compounds_rotation(common, shape)
            # Rotate the lattice, if required
            self.rotate(delta)
            # --------------------------------------------------------------
            # Check if any symmetry applies by comparing the area of the new
            # compound wrt the previous one
            pre_area = get_basic_properties(self.lattice_cmpd)[1]
            new_area = get_basic_properties(shape)[1]
            # Apply the symmetry if the result of the division is:
            # 1/2, 1/4, 1/6, 1/8, 1/12
            res = new_area / pre_area
            if math.isclose(res, 1/2):
                self.apply_symmetry(SymmetryType.HALF)
            elif math.isclose(res, 1/4):
                self.apply_symmetry(SymmetryType.QUARTER)
            elif math.isclose(res, 1/6):
                self.apply_symmetry(SymmetryType.SIXTH)
            elif math.isclose(res, 1/8):
                self.apply_symmetry(SymmetryType.EIGHTH)
            elif math.isclose(res, 1/12):
                self.apply_symmetry(SymmetryType.TWELFTH)

            # If none of the previous operations have been performed, simply
            # update the compound
            if not math.isclose(pre_area, new_area):
                self.lattice_cmpd = shape
                # Rebuild the compound storing the lattice unique edges
                self.lattice_edges: Any = make_partition(
                    [self.lattice_cmpd], [], ShapeType.EDGE)
                print("Other elements need to be updated??")
        else:
            print("The compound position has changed")
            # Check if they have the same geometric properties
            if (get_basic_properties(self.lattice_cmpd) !=
                get_basic_properties(shape)):
                print("The new compound is different from the previous one!!")
                # Raise a warning, since no center can be clearly identified
                raise WarningMessage("Warning: the updated compound appears "
                                     "to be placed in a different position. "
                                     "The area is different as well. Please, "
                                     "call the 'update_geometry_from_center' "
                                     "method to specify also the new center.")
            # The new compound has the same properties: the lattice can be
            # safely translated to the new position given by the new CDG
            self.translate(get_point_coordinates(cdg))
            # Check if any rotation is present
            delta = compare_compounds_rotation(self.lattice_cmpd, shape)
            # Rotate the lattice, if required
            self.rotate(delta)

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
        z_axis = make_vector((0, 0, 1))
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
        # Rotate the lattice box, if any
        if self.lattice_box:
            self.lattice_box.rotate_from_axis(angle, z_axis)
        # Rotate the lattice compound on which a symmetry operation has been
        # applied, if any
        if self.symmetry_type != SymmetryType.FULL:
            self.lattice_symm = make_rotation(
                self.lattice_symm, z_axis, rotation)

        # FIXME re-evaluate lattice characteristic dimensions??

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
        Method that allows to set a value of a given type of property for
        the lattice region which is currently selected in the SALOME study.

        Parameters
        ----------
        property_type : PropertyType
            The value of the 'PropertyType' enumeration indicating which
            type of property to assign to the selected region
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
        found = False
        # Point for idendifying the shape in the geometry
        point = make_vertex_inside_face(shape)
        cells: List[GenericCell] = self.lattice_cells + [self.lattice_box]
        # Get the region that corresponds to the given shape
        for cell in cells:
            # Continue with another cell if the shape is not contained within
            # the cell's area
            if (not is_point_inside_shape(point, cell.face) and not
                get_basic_properties(make_common(cell.face, shape)) ==
                    get_basic_properties(shape)):
                continue
            # Search for the region among the faces stored in the cell
            # dictionary of regions VS properties
            for region in cell.tech_geom_props:
                if (is_point_inside_shape(point, region) and
                    get_basic_properties(make_common(region, shape)) ==
                        get_basic_properties(shape)):
                    found = True
                    cell.tech_geom_props[region] = {property_type: value}
                    break
            if found:
                break

        # Raise an exception if no region has been found
        if not found:
            raise RuntimeError(
                "No cell region could be found for the selected shape.")


def compare_compounds_rotation(cmpd1: Any, cmpd2: Any) -> float:
    """
    Function for comparing the rotation angle of two given compounds. If they
    have a different angle, their difference is calculated and checked against
    the allowed values for the lattice.

    Parameters
    ----------
    cmpd1 : Any
            One of the geometrical compounds to compare
    cmpd2 : Any
            One of the geometrical compounds to compare

    Returns
    -------
    The difference (in degrees) of the rotation angles of the two compounds,
    if they have different rotations; otherwise, the common rotation angle.
    """
    # Get the principal axis angle for both compounds
    angle1 = get_principal_axis_angle(cmpd1)
    angle2 = get_principal_axis_angle(cmpd2)
    # Check if the rotation angle is the same
    if angle1 != angle2:
        print("The rotation angle has changed")
        # Calculate the rotation angle difference
        d_angle = angle2 - angle1
        # Check if the rotation angle is one of the allowed ones
        if not abs(d_angle) in Lattice.VALID_CELLS_ANGLES:
            raise ValueError("The lattice compound has been rotated by "
                             f"an invalid angle of {d_angle}. The allowed "
                             f"ones are {Lattice.VALID_CELLS_ANGLES}°.")
        # Return the difference of the compounds rotation angles
        return d_angle
    # Return the common rotation angle of the compounds
    return angle1


def get_principal_axis_angle(shape: Any) -> float:
    """
    Function for getting the principal axis angle of the given geometrical
    shape.
    The rotation angle of the shape in the XY plane is given from the atan
    of the XY components of this axis.

    Parameters
    ----------
    shape : Any
            The geometrical shape in the XY plane whose rotation angle has
            to be determined

    Returns
    -------
    The rotation angle (in degrees) of the shape in the XY plane.
    """
    # Get the shape inertia matrix
    inertia_matrix = get_inertia_matrix(shape)
    # Extract the principal X-Y axes coordinates (2D shape)
    px = inertia_matrix[-3]
    py = inertia_matrix[-2]
    # Return the rotation angle in degrees
    return math.degrees(math.atan2(py, px))


def update_shape_position_on_lattice(
        center: Any,
        pre_center: Any,
        new_pos: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Function that updates a shape center which is relative to a given point.
    The new position of the point is given, so that the returned XYZ
    coordinates indicates the position of the shape relative to that.

    Parameters
    ----------
    center      : Any
                  The current shape center point
    pre_center  : Any
                  The previous point the shape is relative to
    new_pos     : Tuple[float, float, float]
                  The XYZ coordinates of the new point the shape must be
                  relative to

    Returns
    -------
    The XYZ coordinates of the new shape center so that is still relative to
    the new center position.
    """
    # Get the shape position wrt the previous point position
    shape_to_point = tuple(o - o_shape for o, o_shape in zip(
        get_point_coordinates(pre_center),
        get_point_coordinates(center)))
    # Return the updated shape position relative to the new point position
    return tuple(o - o_shape for o, o_shape in zip(
        new_pos,
        shape_to_point))
