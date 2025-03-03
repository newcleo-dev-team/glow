"""
Testing the construction of a lattice made by cartesian cells. The 'MATERIAL'
property is assigned to each region of the cell technological geometry. A
sectorization is applied to the cell and the result graphically shown in the
SALOME 3D viewer.
A lattice is built by including 9 of the built cells; the result is shown in
the SALOME 3D viewer and its surface geometry representation exported to an
output TDT file.
"""
from glow.geometry_layouts.cells import RectCell
from glow.generator.support import GeometryType, PropertyType, SymmetryType
from glow.geometry_layouts.lattices import Lattice
from glow.main import analyse_and_generate_tdt
from glow.interface.geom_interface import *

# -----------------------
# TEST THE CARTESIAN CELL
# -----------------------
# Build a cartesian cell
rect_cell = RectCell(name="Cartesian cell")
# Add three inner circles to the cell
radii = [0.65/2, 0.3, 0.75/2]
for r in radii:
    rect_cell.add_circle(r)

# Assign the materials to each zone in the cell
rect_cell.set_properties(
    {PropertyType.MATERIAL: ["TSTR.'TMIL_MOC'.'COMB0201'.1",
                             "TSTR.'TMIL_MOC'.'COMB0201'.2",
                             "TSTR.'TMIL_MOC'.'COMB0201'.3",
                             "TSTR.'TMIL_MOC'.'MODEXTCRAYON'"]}
)

# Build the cell sectorization
rect_cell.sectorize([1, 1, 4, 8], [0, 0, 0, 22.5], windmill=True)
# Show the sectorized cell with regions colored according to the 'MATERIAL'
# property
rect_cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
# Print the properties assigned to each cell sector
for sec in rect_cell.regions:
    print(f"Sector region {sec.name}, {sec.properties}")

# -----------------------------
# TEST THE LATTICE CONSTRUCTION
# -----------------------------
# Build a lattice made of sectorized cartesian cells
lattice = Lattice([rect_cell], 'Cartesian Lattice')
dx = rect_cell.width
dy = rect_cell.height
# lattice.add_cell(rect_cell, (dx, 0, 0))
# lattice.add_cell(rect_cell, (dx, dy, 0))
# lattice.add_cell(rect_cell, (0, dy, 0))
# lattice.add_cell(rect_cell, (-dx, dy, 0))
# lattice.add_cell(rect_cell, (-dx, 0, 0))
# lattice.add_cell(rect_cell, (-dx, -dy, 0))
# lattice.add_cell(rect_cell, (0, -dy, 0))
# lattice.add_cell(rect_cell, (dx, -dy, 0))
lattice.add_rings_of_cells(rect_cell, 2)

# Apply the symmetry to the cartesian-type lattice
lattice.apply_symmetry(SymmetryType.EIGHTH)

# Build a box for the lattice
lattice.build_lattice_box([-0.075, 0.075])

# Add MATERIAL property for the box layers
lattice.set_lattice_box_properties({
    PropertyType.MATERIAL: ["MAT1", "MAT2"]
})

# Show the result in the SALOME 3D viewer
lattice.show(PropertyType.MATERIAL)

# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt(lattice, "cartesian_lattice")
