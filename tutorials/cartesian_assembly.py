"""
Use case showing the construction of a lattice made by the same cartesian cell.
The 'MATERIAL' property is assigned to each region of the cell's technological
geometry. A sectorization is applied to the cell with the 'windmill' option
enabled.
The lattice is built by adding several rings of the same cell around a central
one and a box, which encloses the lattice, is declared as the union of different
rectangular shapes.
An eighth symmetry is applied to the whole assembly and the resulting geometry
layout is shown in the SALOME 3D viewer.
In the end, the surface geometry representation of one eighth of the assembly
is exported to an output TDT file.
"""
from glow.geometry_layouts.cells import RectCell
from glow.geometry_layouts.geometries import Rectangle
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.geometry_layouts.lattices import Lattice
from glow.main import TdtSetup, analyse_and_generate_tdt
from glow.interface.geom_interface import *

# Build the cell's geometry layout by adding three circular regions
cell = RectCell(name="Cartesian cell")
radii = [0.2, 0.3, 0.4]
for radius in radii:
    cell.add_circle(radius)
# Assign the materials to each zone in the cell
cell.set_properties(
      {PropertyType.MATERIAL: ["MAT_1", "MAT_2", "MAT_3", "MAT_4"]}
)
# Apply the cell's sectorization
cell.sectorize([1, 1, 4, 8], [0, 0, 0, 22.5], windmill=True)

# --------------------
# LATTICE CONSTRUCTION
# --------------------
# Build the lattice with several rings of the same cartesian cell
lattice = Lattice([cell], 'Cartesian Lattice')
lattice.add_rings_of_cells(cell, 4)
lattice.show(PropertyType.MATERIAL)

# Build the cell representing the lattice's box so that it sligthly cuts
# the outmost ring of cells; the box is subdivided by means of squares at
# its corners. The dimensions of the lattice are extracted to get the box
# dimensions.
x_min, x_max, y_min, y_max = get_bounding_box(lattice.lattice_cmpd)
thickness = 0.1
box = RectCell(
    height_x_width=((y_max-y_min) + thickness, (x_max-x_min) + thickness)
)
box.set_properties({PropertyType.MATERIAL: ["MAT_2"]})
layer_1 = Rectangle(
    height=(y_max-y_min) - thickness,
    width=(x_max-x_min) - thickness
)
corners = [
    Rectangle((x_max, y_max, 0.0), thickness, thickness),
    Rectangle((x_max, y_min, 0.0), thickness, thickness),
    Rectangle((x_min, y_min, 0.0), thickness, thickness),
    Rectangle((x_min, y_max, 0.0), thickness, thickness),
]
# Assemble all the geometric shapes together
box_face = make_partition(
    [box.face],
    [layer.face for layer in [layer_1] + corners],
    ShapeType.COMPOUND
)
# Update the box cell's technological geometry with the assembled one
box.update_geometry_from_face(GeometryType.TECHNOLOGICAL, box_face)

# Assemble the box's cell with the whole lattice and show the result in the
# SALOME 3D viewer
lattice.lattice_box = box
lattice.show(PropertyType.MATERIAL)

# Apply the eighth symmetry type to the cartesian lattice
lattice.apply_symmetry(SymmetryType.EIGHTH)
# Show the resulting layout with the 'MATERIAL' colorset
lattice.show(PropertyType.MATERIAL)

# Perform the geometry analysis and export the TDT file of the surface
# geometry
analyse_and_generate_tdt(
    lattice, "cartesian_lattice", TdtSetup(GeometryType.SECTORIZED))
