"""
Testing the construction of a lattice made by cartesian cells. The 'MATERIAL'
property is assigned to each region of the cell technological geometry.
A lattice is built by including 36 of the built cells. Symmetry is exploited,
one-eighth of the complete cartesian lattice is considered.
"""
import os
import sys

from glow.geometry_layouts.cells import RectCell
from glow.support.types import *
from glow.geometry_layouts.lattices import Lattice
from glow.main import analyse_and_generate_tdt
from glow.interface.geom_interface import *

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

# Show the cell with regions colored according to the 'MATERIAL'
# property
rect_cell.show(PropertyType.MATERIAL)

# Build a lattice made of cartesian cells
lattice = Lattice([], 'Cartesian Lattice')

lattice.add_rings_of_cells(rect_cell, 3)

lattice.apply_symmetry(SymmetryType.EIGHTH)

# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt(
    [lattice],
    os.path.join(
        os.path.dirname(sys.argv[0]), 'test_cartesian_assembly_eighth')
)
