"""
Use case showing the construction of an assembly where the lattice is made of
hexagonal cells of different dimensions: several rings of hexagonal cells
characterised by a smaller size are overlapped by others of a greater size.
The result is that the geometry layout of the small cells being overlapped is
cut, which is a scenario that cannot occur in real-case situations.
Hence, those cells are restored by removing all their circular regions and by
setting their 'MATERIAL' property to the same value.
The lattice is enclosed within a hexagon-shaped box with different thicknesses.
Its type of geometry is set so that the resulting surface geometry representation
applies BCs of type TRAN on the borders of the assembly; this implies the use
of a cycling tracking in SALT.
The built lattice is shown in the 3D viewer of SALOME and its surface geometry
representation exported to an output TDT file.
"""
from glow.main import analyse_and_generate_tdt
from glow.interface.geom_interface import *
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice, get_changed_cells
from glow.support.types import *


# ----------------------------------------------------------------------
# Build the hexagonal cell that constitutes the lattice. It is rotated
# by 90° as needed for tracking
cell_1 = HexCell(name="Cell 1")
cell_1.rotate(90)
radii = [0.1, 0.6, 0.625, 0.70]
for radius in radii:
    cell_1.add_circle(radius)
cell_1.set_properties(
    {PropertyType.MATERIAL: [
        "GAP", "FUEL", "GAP", "CLADDING", "COOLANT"]}
)

# ----------------------------------------------------------------------
# Build the second hexagonal cell
cell_2 = HexCell(edge_length=2.0, name="Cell 2")
radii = [1.0, 1.25]
for radius in radii:
    cell_2.add_circle(radius)
cell_2.set_properties(
    {PropertyType.MATERIAL: ['COOLANT', 'CLADDING', 'COOLANT']}
)

# ----------------------------------------------------------------------
# Build the lattice and add both types of cells
lattice = Lattice(cells=[cell_1])
lattice.add_rings_of_cells(cell_1, 6)
# XY coordinates of the centres of the cells with greater size
x = 4.330127
y = 4.5
lattice.add_cell(cell_2, ())
lattice.add_cell(cell_2, (x, y, 0.0))
lattice.add_cell(cell_2, (-x, y, 0.0))
lattice.add_cell(cell_2, (x, -y, 0.0))
lattice.add_cell(cell_2, (-x, -y, 0.0))
# Show the lattice's technological geometry with the 'MATERIAL' colorset
lattice.show(PropertyType.MATERIAL)

# Get the cells whose geometry layout has been cut and restore them by
# assigning a specific property type
lattice.restore_cells(
    get_changed_cells(lattice),
    {PropertyType.MATERIAL: 'COOLANT'}
)

# Add a container for the assembly and assign properties
lattice.build_lattice_box([0.15, 0.15])
lattice.set_lattice_box_properties(
    {PropertyType.MATERIAL: ['COOLANT', 'CLADDING', 'COOLANT']})
# Show the lattice's technological geometry
lattice.show(PropertyType.MATERIAL)

# ----------------------------------------------------------------------
# Change the lattice type of geometry to use 'TRANSLATION' BCs and cycling
# tracking type
lattice.type_geo = LatticeGeometryType.HEXAGON_TRAN

# ----------------------------------------------------------------------
# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt([lattice], 'hexagonal_assembly')
