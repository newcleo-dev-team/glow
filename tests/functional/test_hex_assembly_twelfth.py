"""
Modelling of a Hexagonal Assembly made of a central pin fuel cell sorrounded by
6 rings of pin fuel cells. Symmetry is exploited, one-twelfth of the complete
hexagon is considered.
"""
import os
import sys

from glow.geometry_layouts.cells import *
from glow.geometry_layouts.lattices import *
from glow.main import analyse_and_generate_tdt
from glow.generator.generator import *

# Declare the values of the hexagonal cells geometrical characteristics
edge_length = 0.7852193995
radii_f=[0.1, 0.45, 0.465, 0.525] #fuel cell
#Build two types of hexagonal cells
fuel_cell= HexCell(edge_length=edge_length,
                    name="Fuel Cell")
# Rotate the cells
fuel_cell.rotate(90)
# Add the circles representing the different zones for pellet and cladding
for r in radii_f:
    fuel_cell.add_circle(r)

# Assign the properties to each cell region
fuel_cell.set_properties(
    {PropertyType.MATERIAL: [
        "HELIUM", "MOX", "HELIUM", "CLADDING", "COOLANT"]}
)

# Update the viewer showing a color for the MATERIAL property type
fuel_cell.show(PropertyType.MATERIAL)
# Build the assembly made of a central dummy cell and 6 rings of fuel cells around
# the central one
lattice = Lattice([fuel_cell],"Fuel assembly")
# Add a specific number of rings of cells to the lattice
lattice.add_rings_of_cells(fuel_cell, 6)

# Add the lattice box
lattice.build_lattice_box([0.05, 0.05])
# Assign the properties to the lattice box areas
lattice.set_lattice_box_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CLADDING", "COOLANT"]}
)

# Assign a boundary condition
lattice.type_geo = LatticeGeometryType.HEXAGON_TRAN

# Assign a symmetry
lattice.apply_symmetry(SymmetryType.TWELFTH)

# Perform the lattice faces and edges analysis and generate the output
# TDT file
analyse_and_generate_tdt(lattice, os.path.join(os.path.dirname(sys.argv[0]),'test_hex_assembly_twelfth'))




