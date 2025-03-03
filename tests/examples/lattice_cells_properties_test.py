"""
Testing an assembly where the lattice is made of a central hexagonal cell
surronded by other six.

The central cell has the following characteristics:
- side: 1 cm
- radii: 0.1, 0.2, 0.3, 0.5 cm
- properties: 'MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5'

The six outer cells has the following characteristics:
- side: 1 cm
- radii: 0.1, 0.45, 0.465, 0.525 cm
- properties: 'HOLLOW', 'MOX', 'GAP', 'CLADDING', 'COOLANT'

The lattice is enclosed within an hexagon-shaped box with different
thicknesses. The 'MATERIAL' property is assigned to each layer; the
direction of assignment is outwards from the lattice center.

The built lattice is shown in the 3D viewer of SALOME and its surface geometry
representation exported to an output TDT file.
"""
# Import the classes and functions for performing the geometry conversion
import math
import os
import sys
import time
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice
from glow.generator.support import PropertyType
from glow.main import analyse_and_generate_tdt


# Get the initial time
start_time = time.time()

# Declare the values of the hexagonal cells geometrical characteristics
radii_1 = [0.1, 0.2, 0.3, 0.5]
radii_2 = [0.1, 0.45, 0.465, 0.525]
# Build the two types of hexagonal cells
hex_cell_1 = HexCell(edge_length=1.0, name="Hexagonal Cell 1")
hex_cell_2 = HexCell(edge_length=1.0, name="Hexagonal Cell 2")

# Add the circles representing the different zones
for r in radii_1:
    hex_cell_1.add_circle(r)
for r in radii_2:
    hex_cell_2.add_circle(r)

# Add both cell face and its subfaces from the technological geometry
# to the SALOME study
hex_cell_1.show()
hex_cell_2.show()

# Assign the properties to each zone in the two cell types
hex_cell_1.set_properties(
    {PropertyType.MATERIAL: ["MAT1",
                             "MAT2",
                             "MAT3",
                             "MAT4",
                             "MAT5",]})
hex_cell_2.set_properties(
    {PropertyType.MATERIAL: ["HOLLOW",
                             "MOX",
                             "GAP",
                             "CLADDING",
                             "COOLANT"]})

# Build a lattice made of a central cell and six cells around the first
lattice = Lattice([hex_cell_1], 'Mixed Cells Hexagonal Lattice')
dx = 2*hex_cell_1.apothem*math.sin(math.radians(60))
dy = 2*hex_cell_1.apothem*math.cos(math.radians(60))

lattice.add_cell(hex_cell_2, (dx, dy, 0))
lattice.add_cell(hex_cell_2, (-dx, dy, 0))
lattice.add_cell(hex_cell_2, (-dx, -dy, 0))
lattice.add_cell(hex_cell_2, (dx, -dy, 0))
lattice.add_cell(hex_cell_2, (0, 2*dy, 0))
lattice.add_cell(hex_cell_1, (0, -2*dy, 0))

# Add the lattice box
lattice.build_lattice_box([0.05, 0.075])
# Assign the properties to the lattice box areas
lattice.set_lattice_box_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CLADDING", "COOLANT"]}
)

# Show the lattice with regions colored according to the 'MATERIAL'
# property
lattice.show(PropertyType.MATERIAL)

# --------------------------------------------------------------------
# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt(lattice, os.path.join(os.path.dirname(sys.argv[0]),
                                               'Mixed_lattice_7cells'))

# Print the execution time
print(f"--- Code executed in {time.time() - start_time} seconds ---")
