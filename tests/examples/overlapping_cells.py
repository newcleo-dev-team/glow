"""
Testing an assembly where the lattice is made of hexagonal cells of different
dimensions: several rings of hexagonal cells characterised by a smaller size
are overlapped by others of a greater size.

The smaller hexagonal cells have the following characteristics (ALFRED test
case):
- edge: 0.7852193995 cm
- pellet inside (hollow) radius: 0.1 cm
- pellet outside radius: 0.45 cm
- cladding inside radius: 0.465 cm
- cladding outside radius: 0.525 cm
- properties: 'HOLLOW', 'MOX', 'HOLLOW', 'CLADDING', 'COOLANT'

The greater hexagonal cells have the following characteristics:
- edge: 2.0 cm
- cladding inside radius: 1.0 cm
- cladding outside radius: 1.25 cm
- properties: 'COOLANT', 'CLADDING', 'COOLANT'

The lattice is enclosed within a hexagon-shaped box with different thicknesses.
It has the following characteristics:
- thicknesses: -0.05, 0.05, 0.05 cm
- properties: 'COOLANT', 'CLADDING', 'COOLANT'

The built lattice is shown in the 3D viewer of SALOME and its surface geometry
representation exported to an output TDT file.
"""
# Import the module for measuring the execution time
import math
import sys
import time

from pathlib import Path

# Import the GLOW geometry classes
from glow.main import analyse_and_generate_tdt
from glow.interface.geom_interface import *
from glow.geometry_layouts.geometries import Hexagon
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice
from glow.generator.support import *

# Get the initial time
start_time = time.time()

# ----------------------------------------------------------------------
# Build the hexagonal cell with greater size
cell1 = HexCell(edge_length=2.0, name="Cell 1")
radii = [1.0, 1.25]
# Add the circles delimiting the different regions
for r in radii:
    cell1.add_circle(r)

# Assign the properties to each cell region
cell1.set_properties(
    {PropertyType.MATERIAL: ['COOLANT', 'CLADDING', 'COOLANT']}
)
# Show the result in the SALOME 3D view
cell1.show(PropertyType.MATERIAL)

# ----------------------------------------------------------------------
# Build the hexagonal cell with smaller size
cell2 = HexCell(edge_length=0.7852193995, name="Cell2")
cell2.rotate(90)
radii = [0.1, 0.45, 0.465, 0.525]
# Add the circles delimiting the different regions
for r in radii:
    cell2.add_circle(r)

# Assign the properties to each cell region
cell2.set_properties(
    {PropertyType.MATERIAL: [
        "HOLLOW", "MOX", "HOLLOW", "CLADDING", "COOLANT"]}
)
# Show the result in the SALOME 3D view
cell2.show(PropertyType.MATERIAL)

# ----------------------------------------------------------------------
# Build the lattice
lattice = Lattice(cells=[cell2])
lattice.add_rings_of_cells(cell2, 6)
# Add a container for the assembly and assign properties
lattice.build_lattice_box([-0.05, 0.05, 0.05])
lattice.set_lattice_box_properties(
    {PropertyType.MATERIAL: ['COOLANT', 'CLADDING', 'COOLANT']})

# ----------------------------------------------------------------------
# Add the cells with greater size
x = 3.4001
y = 3.533487
lattice.add_cell(cell1, ())
lattice.add_cell(cell1, (x,y,0))
lattice.add_cell(cell1, (-x,y,0))
lattice.add_cell(cell1, (x,-y,0))
lattice.add_cell(cell1, (-x,-y,0))

# ----------------------------------------------------------------------
# Display the assembly in the SALOME 3D viewer
lattice.show(PropertyType.MATERIAL)
print(f"--- Geometry generated in {time.time() - start_time} seconds ---")

# ----------------------------------------------------------------------
# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt(
    lattice,
    str(Path(sys.argv[0]).parent.parent / 'outputs/overlapping_cells'))

print(f"--- Code executed in {time.time() - start_time} seconds ---")
