"""
Use case showing the construction of a colorset that reproduces the case of a
control rod assembly surrounded by fuel assemblies.
Given the available symmetry of one twelfth of the colorset, only the two
assemblies that concur in identifying the portion to study are built herein.

The fuel assembly is built first as a lattice with hexagonal cells framed in
a box which partially cuts the outmost ring of cells.
The control rod assembly is a lattice with a hexagonal cell, framed in a box,
that replicates the layout. The control rods pins are built as 'GenericCell'
instances from a compound made by three overlapping circles.
These cells are added to the lattice by placing them around the circumference
of two circles centred in the control rod assembly.

To replicate the colorset, the fuel assembly is translated to the upper-right
side of the control rod assembly.
The resulting compound made by the two assemblies is built and shown in the 3D
viewer of SALOME.
The shape that replicates the S30 type of symmetry is built and a common
operation extracts the portion of the colorset compound to study.

To enable the visualization of the regions of the colorset portion that belong
to the two assemblies at once with the same 'MATERIAL' color map, the regions,
and the names of the materials, are collected.
Given the number of materials, colors are generated and associated to the
regions according to their specific material name.
The regions are then added to the SALOME study so that they can be displayed
in the 3D viewer of SALOME.

Lastly, the surface geometry representation of the colorset portion is
exported to an output TDT file by setting the 'TdTSetup' attributes so that
the TDT file describes an S30 symmetry for an isotropic tracking (TISO) in
SALT.
"""
import math
import time

from glow import *
from glow.geometry_layouts.layouts import associate_colors_to_regions, build_compound_regions
from glow.support.types import *

# Get start time
t0 = time.time()

# -------------------------------------------------------------------------- #
#                                 FUNCTIONS                                  #
# -------------------------------------------------------------------------- #
def add_circular_regions(
        cell: Cell, radii: List[float], materials: List[float]
    ) -> None:
    """
    Function that adds cell-centred circular ``Region`` objects to the given
    ``Cell`` instance. Regions are characterised in terms of the radius and
    the material property.

    Parameters
    ----------
    cell : Cell
        The ``Cell`` instance the circular regions are added to.
    radii : List[float]
        The list of radii of the circular regions in ascending order.
    materials : List[str]
        The list of material names of the circular regions ordered from the
        inner to the outer region.
    """
    for radius, mat in zip(radii[::-1], materials[::-1]):
        cell.add(
            Region(
                Circle(radius=radius),
                properties={PropertyType.MATERIAL: mat}
            )
        )


def create_vertices_list(circle_radius: float, n_vertices: int) -> List[Any]:
    """
    Function that creates a list of vertex objects laying on the same
    circumference with the given radius. The number of vertices is provided
    as input.

    Parameters
    ----------
    circle_radius : float
        The radius of the circle on whose circumference the vertices are
        placed.
    n_vertices : int
        The number of vertices to build.

    Returns
    -------
    List[Any]
        The list of vertex objects laying on the same circumference.
    """
    # Build the circle
    circle = Circle(
        radius=circle_radius, name=f"Circle with radius {circle_radius}")
    vertices = []
    # Build the vertices on the circle's circumference
    for n in range(n_vertices):
        vertex = make_vertex_on_curve(circle.borders[0], n/n_vertices)
        vertices.append(vertex)
    return vertices


# -------------------------------------------------------------------------- #
# FUEL ASSEMBLY CONSTRUCTION                                                 #
# -------------------------------------------------------------------------- #
# Build the hexagonal cells of the fuel assembly
fuel_cell = HexCell(
    name="Fuel cell", base_props={PropertyType.MATERIAL: "COOLANT"}
)
fuel_cell.rotate(90)
# Add the circular regions with their materials
add_circular_regions(
    fuel_cell, [0.2, 0.6, 0.62, 0.68], ["GAP", "FUEL", "GAP", "CLADDING"]
)

central_cell = HexCell(
    name="Central cell", base_props={PropertyType.MATERIAL: "COOLANT"}
)
central_cell.rotate(90)
add_circular_regions(central_cell, [0.6, 0.65], ["GAP", "CLADDING"])

# Update the viewer showing the two cells with the MATERIAL color map
fuel_cell.show(PropertyType.MATERIAL)
central_cell.show(PropertyType.MATERIAL)

# Build the fuel assembly lattice of the colorset
fuel_lattice = HexLattice([central_cell], name="Fuel Assembly Lattice")
fuel_lattice.add_rings_of_cells(fuel_cell, 5)
# Build the cell framing the fuel lattice into a fuel assembly and add regions
# for layers and the lattice
layers_t = [-0.1, 0.3, 0.3]
fuel_assembly = HexCell(
    side=fuel_lattice.dimensions[0] + sum(layers_t),
    base_props={PropertyType.MATERIAL: "COOLANT"}
)
cladding_layer = Region(
    Hexagon(edge_length=fuel_lattice.dimensions[0] + layers_t[1])
    - Hexagon(edge_length=fuel_lattice.dimensions[0] + layers_t[0]),
    properties={PropertyType.MATERIAL: "CLADDING"}
)
fuel_assembly.add(fuel_lattice)
fuel_assembly.add(cladding_layer)

# Display the fuel assembly with the MATERIAL colour map
fuel_assembly.show(PropertyType.MATERIAL)

# -------------------------------------------------------------------------- #
# CONTROL ROD ASSEMBLY CONSTRUCTION                                          #
# -------------------------------------------------------------------------- #
# Data
pitch = fuel_assembly.dimensions[1] * 2
edge_bypass = (pitch) / math.sqrt(3)
edge_ext_wrap_o = (pitch - 0.4) / math.sqrt(3)
edge_ext_wrap_i = (pitch - 0.6) / math.sqrt(3)
r_cr_circles_i = 3.2
r_cr_circles_o = 5.2
cr_pin_radii = [0.68, 0.7, 0.78]
cr_wrapper_radii = [7, 7.25]
int_shaft_ir = 1.4
int_shaft_or = 1.7

# Build the control rod assembly as a hexagonal cell
cr_assembly = HexCell(
    side=edge_bypass,
    name= "Control Rod Assembly",
    base_props={PropertyType.MATERIAL: "COOLANT"}
)
# Add the box layers
cr_assembly.add(
    Region(
        Hexagon(edge_length=edge_ext_wrap_o),
        properties={PropertyType.MATERIAL: "CR_CLADDING"}
    )
)
cr_assembly.add(
    Region(
        Hexagon(edge_length=edge_ext_wrap_i),
        properties={PropertyType.MATERIAL: "CR_MIX"}
    )
)

# Add the circles representing the different zones of the wrapper
add_circular_regions(
    cr_assembly, cr_wrapper_radii, ["COOLANT", "CR_CLADDING"]
)
# Build the vertices at which the control rod regions are placed
cr_vertices_i = create_vertices_list(r_cr_circles_i, 6)
cr_vertices_o = create_vertices_list(r_cr_circles_o, 12)
# Build the circular control rod regions placed along two circumferences by
# specifying the same layer index (the last one) to reduce tree complexity
for v in cr_vertices_i + cr_vertices_o:
    cr_assembly.add(
        Region(
            Circle(radius=cr_pin_radii[2]),
            properties={PropertyType.MATERIAL: "CR_CLADDING2"}
        ),
        get_point_coordinates(v),
        len(cr_assembly.layers)
    )
    cr_assembly.add(
        Region(
            Circle(radius=cr_pin_radii[1]),
            properties={PropertyType.MATERIAL: "GAP"}
        ),
        get_point_coordinates(v),
        len(cr_assembly.layers)
    )
    cr_assembly.add(
        Region(
            Circle(radius=cr_pin_radii[0]),
            properties={PropertyType.MATERIAL: "ABSORBER"}
        ),
        get_point_coordinates(v),
        len(cr_assembly.layers)
    )
# Build the central shaft as made by two overlapping circular regions
add_circular_regions(
    cr_assembly, [int_shaft_ir, int_shaft_or], ["COOLANT", "CR_CLADDING"]
)

# Display the control rod assembly
cr_assembly.show(PropertyType.MATERIAL)

# -------------------------------------------------------------------------- #
# COLORSET CONSTRUCTION                                                      #
# -------------------------------------------------------------------------- #
# Translate the fuel assembly to the right of the control rod assembly
fuel_assembly.translate(
    (3/2*cr_assembly.dimensions[0], cr_assembly.dimensions[1], 0)
)

# Build the colorset as a 'Lattice' with the two assembly cells
colorset = Lattice([cr_assembly, fuel_assembly], (0.0, 0.0, 0.0), "Colorset")
# Display the colorset
colorset.show(PropertyType.MATERIAL)

# -------------------------------------------------------------------------- #
# COLORSET S30 SYMMETRY CONSTRUCTION                                         #
# -------------------------------------------------------------------------- #
# Extract the S30 symmetry portion out of the entire colorset
edges = build_contiguous_edges(
    [
        make_vertex((0.0, 0.0, 0.0)),
        make_vertex((
            3/2*fuel_assembly.dimensions[0],
            fuel_assembly.dimensions[1],
            0.0
        )),
        make_vertex((2*fuel_assembly.dimensions[0], 0.0, 0.0))
    ]
)
cutting_face = make_face(edges)
colorset_portion = make_common(colorset, cutting_face)
add_to_study(colorset_portion, "Colorset - S30 Symmetry")

t1 = time.time()
print(f"--- Geometry generated in {t1 - t0} s. ---")

# -------------------------------------------------------------------------- #
# COLORSET REGIONS VISUALIZATION                                             #
# -------------------------------------------------------------------------- #
# For each face in the colorset compound, recover the corresponding 'Region'
# in the colorset 'Lattice' and build the corresponding region to display
colorset_regions = build_compound_regions(colorset_portion, colorset.regions)
# Associate a unique colour to each region according to the material name
associate_colors_to_regions(PropertyType.MATERIAL, colorset_regions)
# Display the regions of the colorset portion with the material colour map
for region in colorset_regions:
    set_color_face(region, region.color)
    add_to_study_in_father(colorset_portion, region, region.name)

update_salome_study()

t2 = time.time()
print(f"--- Geometry displayed in {t2 - t1} s. ---")

# -------------------------------------------------------------------------- #
# COLORSET S30 PORTION TDT EXPORT                                            #
# -------------------------------------------------------------------------- #
# Generate the TDT file from the colorset portion using a specific typgeo and
# symmetry type
export_layout_to_tdt(
    colorset,
    "colorset_s30",
    tdt_setup=TdtSetup(
        type_geo=LayoutGeometryType.SYMMETRIES_TWO,
        symmetry_type=SymmetryType.TWELFTH
    ),
    compound_to_export=colorset_portion
)

print(f"--- Script executed in {time.time() - t0} s. ---")
