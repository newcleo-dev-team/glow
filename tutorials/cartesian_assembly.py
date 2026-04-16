"""
Use case showing the construction of an assembly made by a lattice of
Cartesian cells framed within a box cutting the cells of the outer cells.
The 'MATERIAL' property is assigned to each region of the cell's technological
geometry. A sectorisation is applied to the cell with the 'windmill' option
enabled.
The lattice is built by adding several rings of the same cell around a central
one.
The assembly is modelled from a 'RectCell' instance built by cutting out the
central region, by adding the lattice and by overlapping the box contour.
A one eighth symmetry is applied to the whole assembly, while a Cartesian mesh
is built. The geometry layout resulting from applying the refinement of its
regions is shown in the SALOME 3D viewer according to the 'MATERIAL' color map.
Lastly, the surface geometry representation of the layout is exported to file
according to the TDT format.
"""
from math import pi

from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.geometries import Circle, Rectangle
from glow.geometry_layouts.lattices import CartesianLattice
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_entities import wrap_shape
from glow.interface.geom_interface import *
from glow.main import TdtSetup, export_layout_to_tdt
from glow.support.types import GeometryType, LayoutGeometryType, PropertyType, \
    SymmetryType

# -----------------
# CELL CONSTRUCTION
# -----------------
# Intialise three lists one storing the circular regions radii, the other the
# names of the materials, sorted from the inner to the outer region
radii = [0.2, 0.3, 0.4]
materials = ["MAT_1", "MAT_2", "MAT_3"]
# Build the geometry layout of the cell by adding three circular regions, from
# external to internal
cell = CartesianCell(
    name="Cartesian cell", base_props={PropertyType.MATERIAL: "MAT_4"}
)
for radius, mat in zip(radii[::-1], materials[::-1]):
    cell.add(
        Region(Circle(radius=radius), properties={PropertyType.MATERIAL: mat})
    )
# Apply the cell's sectorisation
cell.sectorize([1, 1, 4, 8], [0, 0, 0, 22.5], windmill=True)

# ---------------------
# ASSEMBLY CONSTRUCTION
# ---------------------
# Build the lattice with several rings of the same Cartesian cell
lattice = CartesianLattice([cell], name='Cartesian Lattice')
lattice.add_rings_of_cells(cell, 4, 1)

# Build the cell representing the entire assembly made by the lattice and a
# box that sligthly cuts the outmost ring of cells.
# The dimensions of the lattice are extracted to get the box dimensions.
box_w, box_h = (lattice.dimensions[0], lattice.dimensions[1])
thickness = 0.1

assembly = CartesianCell(
    width_height=(box_w + 2*thickness, box_h + 2*thickness),
    base_props={PropertyType.MATERIAL: "MAT_2"}
)
# Build the inner region of the assembly that is cut out from the box shape
inner_area = Rectangle(
    height=(box_h - thickness),
    width=(box_w - thickness)
)
# Build the box contour as a 'Region' obtained by cutting the entire assembly
# with the inner region and assigning a material
box_contour = Region(
    assembly - inner_area,
    properties={PropertyType.MATERIAL: "MAT_2"}
)

# Add the lattice to the assembly cell, then apply the assembly layer region
assembly.add(lattice)
assembly.add(box_contour)

# Update the hierarchical tree of the layout and collapse all the layers into
# one, while cutting the refined geometry due to the box layer overlapping
# the outer cells
assembly.update_hierarchical_structure(True)
# Apply the eighth symmetry type to the assembly
assembly.apply_symmetry(SymmetryType.EIGHTH)

# -----------------
# MESH CONSTRUCTION
# -----------------
# Get the X dimension of the pitch
pitch_x = cell.dimensions[0]
# Build XYZ axis vectors
o_x = make_vector((1, 0, 0))
o_y = make_vector((0, 1, 0))
o_z = make_vector((0, 0, 1))
# Build the outer edges of the mesh along Y by relying on the min/max lattice
# points
x_min, x_max, y_min, y_max = get_bounding_box(lattice.shape)
mesh_edge_outer_1 = make_edge(
    make_vertex((x_min + 0.5*thickness, y_min - thickness, 0.0)),
    make_vertex((x_min + 0.5*thickness, y_max + thickness, 0.0)),
)
mesh_edge_outer_2 = make_edge(
    make_vertex((x_max - 0.5*thickness, y_min - thickness, 0.0)),
    make_vertex((x_max - 0.5*thickness, y_max + thickness, 0.0)),
)
# Build the starting edges from which the Y-mesh is built with a 1D
# multi-translation
mesh_edge_1 = make_edge(
    make_vertex((x_min + pitch_x, y_min - thickness, 0.0)),
    make_vertex((x_min + pitch_x, y_max + thickness, 0.0)),
)
mesh_y = make_multi_translation_1d(mesh_edge_1, o_x, 1, 8)
# Join the Y-mesh edges into a single compound
mesh_y = make_compound([mesh_y, mesh_edge_outer_1, mesh_edge_outer_2])
# Rotate Y-mesh compound to obtain the X-mesh and join both in a compound
mesh_x = make_rotation(mesh_y, o_z, pi/2)
mesh = make_compound([mesh_x, mesh_y])

# Collect the refined geometry of the cells and join it with the built mesh
# by performing a partition
assembly.geometry_maps[GeometryType.SECTORIZED] = \
    assembly.get_geometry_map(GeometryType.SECTORIZED) // wrap_shape(mesh)

# Show the resulting refined layout with the 'MATERIAL' colour map
assembly.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)

# Export the surface representation of the layout according to the TDT format
export_layout_to_tdt(
    assembly,
    "cartesian_assembly",
    TdtSetup(
        geom_type=GeometryType.SECTORIZED,
        type_geo=LayoutGeometryType.RECTANGLE_EIGHT,
        symmetry_type=SymmetryType.EIGHTH
    )
)

