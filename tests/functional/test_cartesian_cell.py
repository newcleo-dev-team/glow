"""
Testing the construction of a cartesian cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry.
"""

import os
import sys

file_path = os.path.abspath(__file__)
glow_path = os.path.abspath(os.path.join(file_path, "..", "..", "..")) 

if glow_path not in sys.path:
    sys.path.insert(0, glow_path)

from glow.geometry_layouts.cells import RectCell
from glow.geometry_layouts.lattices import Lattice
from glow.support.types import PropertyType
from glow.main import analyse_and_generate_tdt

# Build a cartesian cell
rect_cell = RectCell(name="Cartesian cell")
# Add three inner circles to the cell
radii = [0.65/2, 0.3, 0.75/2]
for r in radii:
    rect_cell.add_circle(r)

# Assign the materials to each zone in the cell
rect_cell.set_properties(
    {PropertyType.MATERIAL: ["FUEL",
                             "HOLLOW",
                             "CLADDING",
                             "COOLANT"]}
)

# Show the sectorized cell with regions colored according to the 'MATERIAL'
# property
rect_cell.show(PropertyType.MATERIAL)
# Build a lattice made of a central cell 
lattice = Lattice([rect_cell], 'Cartesian Lattice')

# Perform the lattice faces and edges analysis and generate the output
# TDT file
analyse_and_generate_tdt(lattice, os.path.join(os.path.dirname(sys.argv[0]),'test_cartesian_cell'))

