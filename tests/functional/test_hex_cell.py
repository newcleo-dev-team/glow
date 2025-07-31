"""
Testing the construction of a hexagonal cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry.
"""
import os
import sys

file_path = os.path.abspath(__file__)
glow_path = os.path.abspath(os.path.join(file_path, "..", "..", ".."))  

if glow_path not in sys.path:
    sys.path.insert(0, glow_path)

from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice
from glow.support.types import PropertyType
from glow.main import analyse_and_generate_tdt

edge_length = 0.7852193995
# Build a hexagonal cell
hex_cell = HexCell(edge_length=edge_length,
                   name="Hexagonal cell")

radii = [0.1, 0.45, 0.465, 0.525]
# Add three inner circles to the cell
for r in radii:
    hex_cell.add_circle(r)

# Assign the materials to each zone in the cell
hex_cell.set_properties(
    {PropertyType.MATERIAL: ["HOLLOW",
                             "FUEL",
                             "HOLLOW",
                             "CLADDING",
                             "LEAD"]}
)
# Update the viewer showing a color for the MATERIAL property type
hex_cell.show(PropertyType.MATERIAL)

# Build a lattice made of a central cell 
lattice = Lattice([hex_cell], 'Hexagonal Lattice')

# Perform the lattice faces and edges analysis and generate the output
# TDT file
analyse_and_generate_tdt(lattice, os.path.join(os.path.dirname(sys.argv[0]),'test_hex_cell'))

