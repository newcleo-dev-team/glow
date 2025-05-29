"""
Testing an assembly representing the ALFRED case, where the lattice is
made of several rings of hexagonal cells. This lattice is enclosed within an
hexagon-shaped box with two layers of same thickness.

The cell geometrical characteristics are taken from the ALFRED test case:
- edge: 0.7852193995 cm
- pellet inside (hollow) radius: 0.1 cm
- pellet outside radius: 0.45 cm
- cladding inside radius: 0.465 cm
- cladding outside radius: 0.525 cm

The cell properties are assigned starting from the inner circle moving
outwards; they are the followings: 'HOLLOW', 'MOX', 'GAP', 'CLADDING',
'COOLANT'. The last one is also assigned to the area comprised between the
box and the cells.
"""
# Import the classes and functions for performing the geometry conversion
import os
import sys
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice
from glow.generator.support import PropertyType
from glow.main import analyse_and_generate_tdt

# Import the module for measuring the execution time
import time

# Get the initial time
start_time = time.time()

# Declare the values of the hexagonal cell geometrical characteristics
edge_length = 0.7852193995
radii = [0.1, 0.45, 0.465, 0.525]
# Build an hexagonal cell
hex_cell = HexCell(edge_length=edge_length,
                    name="ALFRED Hexagonal Cell")
# Rotate the cell
hex_cell.rotate(90)
# Add the circles representing the different zones for pellet and cladding
for r in radii:
    hex_cell.add_circle(r)

# Assign the properties to each cell region
hex_cell.set_properties(
    {PropertyType.MATERIAL: [
        "HOLLOW", "MOX", "HOLLOW", "CLADDING", "COOLANT"]}
)
# Update the viewer showing a color for the MATERIAL property type
hex_cell.show(PropertyType.MATERIAL)

# Build a lattice made of a central cell and several rings of cells around
# the first one
lattice = Lattice([hex_cell], 'ALFRED Hexagonal Lattice')
# Add a specific number of rings of cells to the lattice
lattice.add_rings_of_cells(hex_cell, 6)

# Add the lattice box
lattice.build_lattice_box([0.05, 0.05])
# Assign the properties to the lattice box areas
lattice.set_lattice_box_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CLADDING", "COOLANT"]}
)

# hex_lattice.apply_symmetry(SymmetryType.SIXTH)
lattice.show(PropertyType.MATERIAL)

print(f"--- Generation executed in {time.time() - start_time} seconds ---")

# Perform the lattice faces and edges analysis and generate the output
# TDT file
analyse_and_generate_tdt(lattice, os.path.join(os.path.dirname(sys.argv[0]),
                                               'alfred_stress_test'))

# Print the execution time
print(f"--- Code executed in {time.time() - start_time} seconds ---")
