"""
This module contains the classes that support the generation of the output
TDT file containing the geometry representation for further analysis in
DRAGON.
"""


from dataclasses import dataclass, field
from enum import Enum
from io import TextIOWrapper
import math
import os
from pathlib import Path
from typing import List, Tuple
from glow.generator.support import BoundaryType, LatticeGeometryType, \
    SymmetryType
from glow.generator.geom_extractor import Boundary, Edge, Face


class Element (Enum):
    """
    Enumeration assigning each geometric element type an index.
    """
    SEGMENT = 1
    CIRCLE  = 2
    ARC     = 3


# Dictionary of element type VS a tuple containing the corresponding
# attribute of the 'Element' class and a descriptive string
TYPE_ELEM = {"SEGMENT"     : (Element.SEGMENT, "line segment"),
             "CIRCLE"      : (Element.CIRCLE, "circle"),
             "ARC_CIRCLE"  : (Element.ARC, "circular arc")}


@dataclass
class TdtData():
    """
    Class that uses the lattice definition in terms of its subfaces and edges
    to generate a data file containing the lattice geometric description in
    the TDT format.

    Attributes
    ----------
    filename      : str
                    Name of the TDT format data file to be generated (.dat)
    edges         : List[Edge]
                    List of edges as 'Edge' objects
    faces         : List[Face]
                    List of faces as 'Face' objects
    boundaries    : List[Boundary]
                    List of the 'Boundary' objects for the lattice borders
    type_sym      : GeometryType
                    The type of the lattice symmetry
    impressions   : Tuple[int, int]
                    Options for printing data
    precisions    : Tuple[float, float]
                    Options for geometric precision
    properties    : List[str]
                    List of the names of the properties the lattice zones are
                    associated with
    properties_id : List[int]
                    List of the IDs of the properties in the lattice
    """
    filename      : str = os.path.join(Path(__file__).resolve().parent.parent,
                                       "tdt_lattice.dat")
    edges         : List[Edge] = field(default_factory=list)
    faces         : List[Face] = field(default_factory=list)
    boundaries    : List[Boundary] = field(default_factory=list)
    type_geo      : LatticeGeometryType = LatticeGeometryType.HEXAGON_TRAN
    type_sym      : SymmetryType = SymmetryType.FULL
    impressions   : Tuple[int, int] = (0, 0)
    precisions    : Tuple[float, float] = (1e-5, 1e-5)
    properties    : List[List[str]] = field(init=False)
    properties_id : List[int] = field(init=False)
    nb_folds      : int = field(init=False)

    def __post_init__(self):
        """
        Method that is automatically run after the dataclass initialization
        for setting all the attributes that depends on others.
        """
        # Set the number of folds for the lattice according to the type of
        # geometry and of symmetry. For geometries already counting a symmetry
        # (i.e. type_geo > 1) this value is set to 0.
        self.nb_folds = 0
        if self.type_geo.value <= LatticeGeometryType.SYMMETRIES_TWO.value:
            self.nb_folds = self.type_sym.value
        # Set the list of property names and IDs
        self.__build_properties_id()

    def __build_properties_id(self):
        """
        Method that builds two lists, one containing the names of the
        properties, the other containing the ID of one of the properties
        that is associated to a subface of the lattice.
        """
        # Initialize the index for the number of properties
        prop_no = 0
        # Initialize the list of properties names
        self.properties = list()
        # Initialize the list of property IDs as a list of '-1' with
        # dimension being the size of the list of subfaces in the lattice
        self.properties_id = [-1]*len(self.faces)
        # Loop through the 'Face' objects
        for face in self.faces:
            # Add the property name to the list, if not already present
            for property in face.properties:
                if not property in self.properties:
                    self.properties.append(property)
                    # Update the index for the number of properties
                    prop_no = prop_no + 1
                    prop = prop_no
                else:
                    # Get the index corresponding to the property name in
                    # the list and increment it by 1
                    prop = self.properties.index(property) + 1

                ### Association face number - property
                # Add the property ID in the list of properties associated
                # to a face at position given by the 'Face' object number
                self.properties_id[face.no - 1] = prop


def write_tdt_file(tdt_data: TdtData) -> None:
    """
    Function that allows to write the TDT output file of the geometry
    conversion.

    Parameters
    ----------
    tdt_data  : TdtData
                The instance of the `TdtData` class storing the information
                of the lattice geometry
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
    file      : TextIOWrapper
                Handle for the opened file to write
    tdt_data  : TdtData
                The instance of the `TdtData` class storing the information
                of the lattice geometry
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
        f"  {typegeom:5d},{nb_folds:5d}, {nbnodes:5d}, {nbelements:5d}, {1:5d},"
        f"{nbregions:5d}, {0:5d}, {1:5d}\n",
         "* index  kindex\n",
        f"  {tdt_data.impressions[0]:5d}  {tdt_data.impressions[1]:6d}  1\n",
         "*     eps    eps0\n",
        f"  {tdt_data.precisions[0]:7E}   {tdt_data.precisions[1]:7E}\n"])

def _write_regions(file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for writing to file the list of regions and macros.

    Parameters
    ----------
    file      : TextIOWrapper
                Handle for the opened file to write
    tdt_data  : TdtData
                The instance of the `TdtData` class storing the information
                of the lattice geometry
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
    file      : TextIOWrapper
                Handle for the opened file to write
    tdt_data  : TdtData
                The instance of the `TdtData` class storing the information
                of the lattice geometry
    """
    # Write the header line for this section
    file.writelines([ "   elements\n" ])
    # Loop through the 'Edge' objects sorted by their number attribute
    for edge in sorted(tdt_data.edges):
        # Write the information about the current edge
        # Get the element type index and its descriptive string
        type_indx, type_descr = TYPE_ELEM[edge.kind]

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
        if type_indx == Element.SEGMENT:
            # Extract the edge data as 'x1, y1, z1, x2, y2, z2'
            x1, y1, _, x2, y2, _ = edge.data[1:]
            # Write the info about the X-Y coordinates of the first point
            # and the X-Y distances between the segment vertices
            file.write(f"  {x1:6f}, {y1:6f}, {(x2-x1):6f}, {(y2-y1):6f}\n")
        elif type_indx == Element.CIRCLE:
            # Extract the edge data as 'xc, yc, zc, dx, dy, dz, R'
            xc, yc, _, _, _, _, R = edge.data[1:]
            # Write the info about the X-Y coordinates of the circle center,
            # the radius and a 4th data (value '0.0') which is mandatory
            # despite not being indicated in the APOLLO2 doc of 15/06/09
            file.write(f"  {xc:6f}, {yc:6f}, {R:6f}, 0.0\n")
        elif type_indx == Element.ARC:
            # Extract the edge data as 'xc, yc, zc, dx, dy, dz, R,
            #                           x1, y1, z1, x2, y2, z2'
            xc, yc, _, _, _, _, R, x1, y1, _, x2, y2, _ = edge.data[1:]
            # Calculate the angles (in degree) of the vertices wrt the
            # circle center
            angle_1 = math.atan2((y1-yc), (x1-xc)) * (180/math.pi)
            angle_2 = math.atan2((y2-yc), (x2-xc)) * (180/math.pi)
            # Since it is necessary to have positive value for the angles
            # difference, the absolute value is considered
            delta_angle = abs(abs(angle_2) - abs(angle_1))

            # Write the info about the X-Y coordinates of the arc circle
            # center, its radius, the angle of the first point of the arc
            # and the angle difference between the two arc points
            file.write(f"  {xc:6f}, {yc:6f}, {R:6f}, {angle_1:6f}, "
                       f"{delta_angle:6f}\n")

def _write_boundary_conditions(file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for determining and writing to file the list of boundary
    conditions.

    Parameters
    ----------
    file      : TextIOWrapper
                Handle for the opened file to write
    tdt_data  : TdtData
                The instance of the `TdtData` class storing the information
                of the lattice geometry
    """
    # Write the header line for this section
    file.write("* boundaries conditions: defaul nbbcda allsur\n")
    # Declare the default number of BCs (0) and the non-default one, equal to
    # the size of the list of 'Boundary' objects identifying the lattice
    # borders
    default_bc = 0 # original line + angle??
    non_default_bc_no = 0
    # Update the number of BCs only if a specific geometry type is set
    # TODO check if also typgeo 1 and 2 should be included
    if not tdt_data.type_geo == LatticeGeometryType.ISOTROPIC:
        non_default_bc_no = len(tdt_data.boundaries)

    # Write the information about the number of BCs (both default and non) and
    # the albedo
    file.write(f"  {default_bc}, {non_default_bc_no}, 0\n")
    file.write("* albedo\n")
    file.write(f"  {1.0:.1f}\n")

    # Nothing more to write if no specific geometry type is set
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
        # Handle the different BC types
        if bc.type == BoundaryType.AXIAL_SYMMETRY:
            # Write the X-Y coordinates of the border edge origin and the
            # angle between its vertices
            file.write("* ax, ay, angle\n")
            file.write(f"  {bc.tx:6f} {bc.ty:6f} {bc.angle:6f}\n")
        elif bc.type == BoundaryType.TRANSLATION:
            # Write the X-Y coordinates of the borders axes
            file.write("* tx, ty, angle\n")
            file.write(f"  {bc.tx:6f} {bc.ty:6f} {0:6f}\n")
        else:
            raise AssertionError(
                f"Treatment of type {bc.type} not implemented")

def _write_properties(file: TextIOWrapper, tdt_data: TdtData) -> None:
    """
    Function for writing the indices of the properties associated to each
    subface of the lattice in the TDT-format file.

    Parameters
    ----------
    file      : TextIOWrapper
                Handle for the opened file to write
    tdt_data  : TdtData
                The instance of the `TdtData` class storing the information
                of the lattice geometry
    """
    # Write the names of the properties that are present in the lattice prior
    # to the header line. Each line starts with a '#' so to be ignored.
    for id, name in enumerate(tdt_data.properties):
        file.write(f"# {(id+1):2d} - {name}\n")
    # Write the header line for this section
    file.write("* medium number per region\n")
    # Loop through the IDs of the materials associated to a face
    for material_id in tdt_data.properties_id:
        # Write the ID of the material associated to a face to the TDT file
        file.write(f"  {material_id}\n")

