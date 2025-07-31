"""
Testing an assembly where the lattice is made of a central hexagonal
cell surronded by other six. This lattice is enclosed within an
hexagon-shaped box with two layers of same thickness. In addition, the
first layer is defined so that it overlaps the cells regions (negative
thickness value).

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

A symmetry is applied to the lattice to test the TDT file generation under
different cases.
"""
# Import the classes and functions for performing the geometry conversion
import math
import os
import sys
import time
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice
from glow.support.types import PropertyType, SymmetryType
from glow.main import analyse_and_generate_tdt

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

# Add both cell face and its subfaces from the technological geometry
# to the SALOME study
hex_cell.show()

# Assign the properties to each zone in the cell
hex_cell.set_properties(
    {PropertyType.MATERIAL: ["HOLLOW",
                             "MOX",
                             "GAP",
                             "CLADDING",
                             "COOLANT"]})

# Build a lattice made of a central cell and six cells around the first
lattice = Lattice([hex_cell], 'ALFRED Hexagonal Lattice')
dx = 2*hex_cell.apothem*math.cos(math.radians(60))
dy = 2*hex_cell.apothem*math.sin(math.radians(60))

lattice.add_cell(hex_cell, (dx, dy, 0))
lattice.add_cell(hex_cell, (-dx, dy, 0))
lattice.add_cell(hex_cell, (-dx, -dy, 0))
lattice.add_cell(hex_cell, (dx, -dy, 0))
lattice.add_cell(hex_cell, (2*dx, 0, 0))
lattice.add_cell(hex_cell, (-2*dx, 0, 0))

# Add the lattice box
lattice.build_lattice_box([-0.05, 0.05, 0.05])
# Show the partitioned box face
lattice.lattice_box.show()
# Assign the properties to the lattice box areas
lattice.set_lattice_box_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CLADDING", "COOLANT"]}
)

# Apply the symmetry
lattice.apply_symmetry(SymmetryType.TWELFTH)

# Show the lattice with regions colored according to the 'MATERIAL'
# property
lattice.show(PropertyType.MATERIAL)

# --------------------------------------------------------------------
# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt(lattice, os.path.join(os.path.dirname(sys.argv[0]),
                                               'lattice_symmetry_test'))

# Print the execution time
print(f"--- Code executed in {time.time() - start_time} seconds ---")
