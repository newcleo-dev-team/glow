"""
Demo showing some of the functionalities of GLOW. A simple geometry layout
made by hexagonal cells with characteristics taken from the ALFRED case is
built.
"""
# Import the module for measuring the execution time
import time
# Get the initial time
start_time = time.time()

from pathlib import Path
import sys

# Import the GLOW geometry classes
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice
from glow.interface.geom_interface import *
from glow.support.types import *
from glow.main import analyse_and_generate_tdt

# ----------------------------------------------------------------------
# Declare the values of the hexagonal cell geometrical characteristics
edge_length = 0.7852193995
radii = [0.1, 0.45, 0.465, 0.525]
# Build an hexagonal cell
cell = HexCell(edge_length=edge_length, name="ALFRED Hexagonal Cell")
# Rotate the cell
cell.rotate(90)
# Add the circles representing the different zones for pellet and cladding
for r in radii:
    cell.add_circle(r)

# Assign the properties to each cell region
cell.set_properties(
    {PropertyType.MATERIAL: [
        "HOLLOW", "MOX", "HOLLOW", "CLADDING", "COOLANT"]}
)

# Build a lattice made of a central cell and several rings of cells around
# the first one
lattice = Lattice([cell], 'ALFRED Hexagonal Lattice')
# Add a specific number of rings of cells to the lattice
lattice.add_rings_of_cells(cell, 2)

# Add the lattice box
lattice.build_lattice_box([0.05, 0.05])
# Assign the property labels to the lattice box regions
lattice.set_lattice_box_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CLADDING", "COOLANT"]}
)

# Perform a SIXTH symmetry operation (SA60 case)
lattice.apply_symmetry(SymmetryType.SIXTH)

# Show the lattice with colors associated to the MATERIAL property
lattice.show(PropertyType.MATERIAL)

# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt(
    lattice, str(Path(sys.argv[0]).parent.parent / 'outputs/demo'))

print(f"--- Code executed in {time.time() - start_time} seconds ---")
