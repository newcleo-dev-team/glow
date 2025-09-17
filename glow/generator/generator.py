"""
This module contains the classes that support the generation of the output
TDT file containing the geometry representation for further analysis in
DRAGON.
"""
import math
import os

from dataclasses import dataclass, field
from io import TextIOWrapper
from pathlib import Path
from typing import List, Tuple

from glow.generator.geom_extractor import Boundary, Edge, Face
from glow.support.types import EDGE_NAME_VS_TYPE, BoundaryType, EdgeType, \
    LatticeGeometryType, SymmetryType


# Precision in terms of number of digits after the decimal
PRECISION : int = 7

# Format string for floating-point values in scientific notation with
# defined precision
FORMAT: str = f"{{:.{PRECISION}E}}"


@dataclass
class TdtData():
    """
    Dataclass storing the geometric information, in terms of subfaces, edges
    and boundaries of the geometry layout, as well as other data providing
    the type of geometry, of symmetry, and the properties associated to the
    regions of the geometry.
    """
    filename: str = os.path.join(Path(__file__).resolve().parent.parent,
                                       "tdt_lattice.dat")
    """Name of the file in the TDT format to be generated (without the `.dat`
        extension)."""
    edges: List[Edge] = field(default_factory=list)
    """List of the geometry layout edges, as ``Edge`` objects."""
    faces: List[Face] = field(default_factory=list)
    """List of the geometry layout regions, as ``Face`` objects."""
    boundaries: List[Boundary] = field(default_factory=list)
    """List of the geometry layout borders, as ``Boundary`` objects."""
    type_geo: LatticeGeometryType = LatticeGeometryType.HEXAGON_TRAN
    """
    The type of geometry applied to the geometry layout, as element of the
    ``LatticeGeometryType`` enumeration.
    """
    type_sym: SymmetryType = SymmetryType.FULL
    """
    The type of the symmetry applied to the geometry layout, as element of
    the ``SymmetryType`` enumeration.
    """
    albedo: float | None = None
    """Identifying the value for the `albedo` applied to the lattice's BCs."""
    impressions: Tuple[int, int] = (0, 0)
    """Options for printing the geometric data."""
    precisions: Tuple[float, float] = (1e-5, 1e-5)
    """Options for the geometric precision of the data."""
    properties: List[str] = field(init=False)
    """
    List of the names of the properties the regions of the geometry layout
    are associated with.
    """
    property_ids : List[int] = field(init=False)
    """List of the IDs of the properties in the geometry layout."""
    nb_folds      : int = field(init=False)
    """
    The number of times the geometry layout has to be unfolded to replicate
    the full geometry, if any symmetry is applied.
    """

    def __post_init__(self) -> None:
        """
        Method that is automatically run after the dataclass initialization
        for setting all the attributes that depends on others.
        """
        # Set the number of folds for the lattice according to the type of
        # geometry and of symmetry. For geometries already counting a symmetry
        # (i.e. type_geo > 2) this value is set to 0.
        self.nb_folds = 0
        if self.type_geo.value <= LatticeGeometryType.ROTATION.value:
            self.nb_folds = self.type_sym.value
        # Set the albedo for the lattice's BCs according to the type of
        # geometry, i.e. by default is 1.0 if ISOTROPIC, 0.0 for the other
        # types
        if self.type_geo == LatticeGeometryType.ISOTROPIC:
            self.albedo = 1.0 if self.albedo is None else self.albedo
        else:
            if self.albedo is not None and self.albedo > 0.0:
                raise RuntimeError(
                    f"A value of {self.albedo} for the albedo is not "
                    f"compatible with the '{self.type_geo}' type of "
                    "geometry.")
            self.albedo = 0.0

        # Set the list of property names and IDs
        self.__build_properties_id()

    def __build_properties_id(self) -> None:
        """
        Method that builds two lists for the properties associated to the
        regions of the lattice: one containing the names of the properties,
        the other containing the corresponding IDs, as integer indices, so
        that they appear only once.
        """
        # Initialize the list of properties names
        self.properties = list()
        # Initialize the list of property IDs as a list of '-1' with
        # dimension being the size of the list of subfaces in the lattice
        self.property_ids = [-1]*len(self.faces)
        # Loop through the 'Face' objects
        for face in self.faces:
            # Add the property name to the corresponding list, if not
            # already present
            if face.property not in self.properties:
                self.properties.append(face.property)
                # Update the unique index for the properties
                prop_indx = len(self.properties)
            else:
                # Get the index that corresponds to the property name and
                # increment it by 1 (as indices start from 0)
                prop_indx = self.properties.index(face.property) + 1

            # Add the property ID in the list of properties: the index at
            # which it is set corresponds to the 'no' attribute of the
            # associated 'Face' object
            self.property_ids[face.no - 1] = prop_indx


def write_tdt_file(tdt_data: TdtData) -> None:
    """
    Function that writes the output file with the characteristics of
    the geometry layout in the TDT format.

    Parameters
    ----------
    tdt_data  : TdtData
        The instance of the ``TdtData`` class storing the information of
        the lattice geometry.
    """
    # Open the file for writing the geometry information in the TDT format
    print(tdt_data.filename)
    with open(tdt_data.filename, "w") as file:
        # Write the TDT header of the file
        _write_header(file, tdt_data)
        # Write the list of regions in the TDT file
        _write_regions(file, tdt_data)
        # Write the list of edges in the TDT file
        _write_edges(file, tdt_data)
        # Write the list of boundaries in the TDT file
        _write_boundary_conditions(file, tdt_data)
        # Write the material information, in terms of IDs (associated
        # to each face in the lattice) and names, to the TDT file
        _write_properties(file, tdt_data)
        # Write the ending lines in the TDT file
        file.write("-" * 60 + "\n")
        file.write("   this is the last line generated by conversion\n")


def _write_header(file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for writing to file the header for TDT-format file.

    Parameters
    ----------
    file : TextIOWrapper
        Handle for the opened file to write.
    tdt_data : TdtData
        The instance of the ``TdtData`` class storing the information of
        the lattice geometry.
    """
    # Declare the number of nodes and regions as equal to the faces number
    nbnodes    = len(tdt_data.faces)
    nbregions  = len(tdt_data.faces)
    # Declare the number of elements to be equal to the number of edges
    nbelements = len(tdt_data.edges)
    # Get the index for the geometry type of the lattice
    typegeom   = tdt_data.type_geo.value
    # Get the number of folds of the lattice
    nb_folds   = tdt_data.nb_folds
    file.writelines([
         "\n",
         "	dat input file for DRAGON5\n",
         "------------------------------------------------------------\n\n\n",
         "* typge, nbfo, node, elem, macr, nreg,    z, mac2\n",
        f"  {typegeom:5d},{nb_folds:5d}, {nbnodes:5d}, {nbelements:5d}, "
        f"{1:5d},{nbregions:5d}, {0:5d}, {1:5d}\n",
         "* index  kindex\n",
        f"  {tdt_data.impressions[0]:5d}  {tdt_data.impressions[1]:6d}  1\n",
         "*     eps    eps0\n",
        f"  {tdt_data.precisions[0]:7E}   {tdt_data.precisions[1]:7E}\n"])


def _write_regions(file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for writing to file the list of regions and macros.

    Parameters
    ----------
    file : TextIOWrapper
        Handle for the opened file to write.
    tdt_data : TdtData
        The instance of the ``TdtData`` class storing the information of
        the lattice geometry.
    """
    # Declare the number of regions to be equal to the number of faces
    nbregions = len(tdt_data.faces)
    # Write a line in the file
    file.writelines(["*   flux region number per geometry region (mesh)\n"])

    # Write the list of region numbers
    for i in range(1, nbregions):
        file.write(f"{i:4d},")
        # Write on a new line after 12 values
        if i % 12 == 0: file.write("\n")
    # Write the last region number without ','
    file.write(f"{nbregions:4d}")
    file.write("\n")

    # Write the information about the macros
    file.writelines([
        "*   names of macros\n",
        "mac\n",
        "*   macro order number per flux region\n",
        f"{nbregions}*1\n"])


def _write_edges(file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for writing to file the list of edges of the lattice. These
    elements are sorted by their number.

    Parameters
    ----------
    file : TextIOWrapper
        Handle for the opened file to write.
    tdt_data : TdtData
        The instance of the ``TdtData`` class storing the information of
        the lattice geometry.
    """
    # Write the header line for this section
    file.writelines([ "   elements\n" ])
    # Loop through the 'Edge' objects sorted by their number attribute
    for edge in sorted(tdt_data.edges):
        # Write the information about the current edge
        # Get the element type index and its descriptive string
        type_indx, type_descr = EDGE_NAME_VS_TYPE[edge.kind.name]

        # Get the index of the left/right faces if any is associated
        # to the edge
        left_indx = edge.left.no if edge.left else 0
        right_indx = edge.right.no if edge.right else 0

        # Write the information about the edge element, its number and
        # descriptive string
        file.write(f"* ELEM  {edge.no}  {type_descr}\n")
        # Write the information about the index of the edge type and those
        # identifying the associated left and right faces
        file.write(f" {type_indx.value}, {right_indx}, {left_indx}\n")
        file.write("*\n")

        # Write the geometric data of the edge according to its type
        if type_indx == EdgeType.SEGMENT:
            # Extract the edge data as 'x1, y1, z1, x2, y2, z2'
            x1, y1, _, x2, y2, _ = edge.data[1:]
            dx = x2-x1
            dy = y2-y1
            if abs(dx) < 1e-7:
                dx = 0.0
            if abs(dy) < 1e-7:
                dy = 0.0
            # Write the info about the X-Y coordinates of the first point
            # and the X-Y distances between the segment vertices
            file.write(f"  {FORMAT.format(x1)}, {FORMAT.format(y1)}, " + \
                       f"{FORMAT.format(dx)}, {FORMAT.format(dy)}\n")
            continue
        if type_indx == EdgeType.CIRCLE:
            # Extract the edge data as 'xc, yc, zc, dx, dy, dz, R'
            xc, yc, _, _, _, _, R = edge.data[1:]
            # Write the info about the X-Y coordinates of the circle center,
            # the radius and a 4th data (value '0.0')
            file.write(f"  {FORMAT.format(xc)}, {FORMAT.format(yc)}, " + \
                       f"{FORMAT.format(R)}, {FORMAT.format(0.0)}\n")
            continue
        if type_indx == EdgeType.ARC_CIRCLE:
            # Extract the edge data as 'xc, yc, zc, dx, dy, dz, R,
            #                           x1, y1, z1, x2, y2, z2'
            xc, yc, _, _, _, _, R, x1, y1, _, x2, y2, _ = edge.data[1:]
            # Calculate the angles (in degree) of the vertices wrt the
            # circle center
            angle_1 = math.atan2((y1-yc), (x1-xc)) * (180/math.pi) % 360
            angle_2 = math.atan2((y2-yc), (x2-xc)) * (180/math.pi) % 360
            # Since it is necessary to have positive value for the angles
            # difference, the absolute value is considered
            delta_angle = (angle_2 - angle_1) % 360
            if abs(delta_angle) < 1e-7:
                delta_angle = 0.0

            # Write the info about the X-Y coordinates of the arc circle
            # center, its radius, the angle of the first point of the arc
            # and the angle difference between the two arc points
            file.write(f"  {FORMAT.format(xc)}, {FORMAT.format(yc)}, " + \
                       f"{FORMAT.format(R)}, {FORMAT.format(angle_1)}, " + \
                       f"{FORMAT.format(delta_angle)}\n")
            continue


def _write_boundary_conditions(
        file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for determining and writing to file the list of boundary
    conditions.

    Parameters
    ----------
    file : TextIOWrapper
        Handle for the opened file to write.
    tdt_data : TdtData
        The instance of the ``TdtData`` class storing the information of
        the lattice geometry.
    """
    # Write the header line for this section
    file.write("* boundaries conditions: defaul nbbcda allsur\n")
    # Declare the default number of BCs (0) and the non-default one, equal to
    # the size of the list of 'Boundary' objects identifying the lattice
    # borders
    default_bc = 0
    non_default_bc_no = 0
    # Update the number of BCs only if a specific geometry type is set
    if tdt_data.type_geo != LatticeGeometryType.ISOTROPIC:
        non_default_bc_no = len(tdt_data.boundaries)

    # Write the information about the number of BCs (both default and non) and
    # the albedo
    file.write(f"  {default_bc}, {non_default_bc_no}, 0\n")
    file.write("* albedo\n")
    file.write(f"  {tdt_data.albedo:.1f}\n")

    # Nothing more to write if the geometry type is ISOTROPIC
    if tdt_data.type_geo == LatticeGeometryType.ISOTROPIC:
        return
    # Loop through all the 'Boundary' objects
    for bc in tdt_data.boundaries:
        # Write the BCs type (as an index) and number of lattice border edges
        file.write("* type  number of elements\n")
        file.write(f"  {bc.get_bc_type_number()}, {len(bc.edge_indxs)}\n")
        file.write("*   elements\n")
        # Loop through all the indexes of the lattice border edges
        for edge_no in bc.edge_indxs:
            # Write the index of the lattice border edge
            file.write(f"{edge_no}\n")
        # Check if the BC type is allowed
        if bc.type not in [BoundaryType.AXIAL_SYMMETRY,
                           BoundaryType.ROTATION,
                           BoundaryType.TRANSLATION]:
            raise RuntimeError(
                f"The '{bc.type}' BC type cannot be treated by the SALT "
                "module of DRAGON5.")
        # Write the X-Y coordinates of the border axes
        file.write("* tx, ty, angle\n")
        file.write(f"  {FORMAT.format(bc.tx)} {FORMAT.format(bc.ty)} " + \
                   f"{FORMAT.format(bc.angle)}\n")


def _write_properties(file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for writing the indices of the properties associated with each
    region of the lattice to the TDT-format file.

    Parameters
    ----------
    file : TextIOWrapper
        Handle for the opened file to write.
    tdt_data : TdtData
        The instance of the ``TdtData`` class storing the information of
        the lattice geometry.
    """
    # Write the names of the properties that are present in the lattice prior
    # to the header line. Each line starts with a '#' so to be ignored.
    for id, name in enumerate(tdt_data.properties):
        file.write(f"# {(id+1):2d} - {name}\n")
    # Write the header line for this section
    file.write("* medium number per region\n")
    # Loop through the IDs of the materials associated to a face
    for material_id in tdt_data.property_ids:
        # Write the ID of the material associated to a face to the TDT file
        file.write(f"  {material_id}\n")
