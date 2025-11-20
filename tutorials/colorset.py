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
from glow.support.types import *

# Get start time
t0 = time.time()

# -------------------------------------------------------------------------- #
# FUNCTIONS DECLARATION                                                      #
# -------------------------------------------------------------------------- #
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


def make_circular_cells_list(
        vertices: List[Any], radii: List[float]) -> List[GenericCell]:
    """
    Function that creates a list of circular cells, as ``GenericCell``
    instances. The cells centres are given by the input list of vertices.

    Parameters
    ----------
    vertices : List[Any]
        The list of vertex objects where each circular cell is placed.
    radii : List[float]
        The list of radii for the circles contained in the circular cell.

    Returns
    -------
    List[GenericCell]
        The list of circular cells.
    """
    cells = []
    # Loop through the given vertices and build the layout of each cell as
    # made by three concentric circles.
    for i, vertex in enumerate(vertices):
        cr_circle_inn = Circle(
            get_point_coordinates(vertex),
            radius=radii[0],
            name=f"Inner CR Circle vertex {get_point_coordinates(vertex)}"
        )
        cr_circle_mid = Circle(
            get_point_coordinates(vertex),
            radius=radii[1],
            name=f"Middle CR Circle vertex {get_point_coordinates(vertex)}"
        )
        cr_circle_out = Circle(
            get_point_coordinates(vertex),
            radius=radii[2],
            name=f"Outer CR Circle vertex {get_point_coordinates(vertex)}"
        )
        cr_pin_compound = make_partition(
            [cr_circle_out.face, cr_circle_mid.face, cr_circle_inn.face],
            [],
            ShapeType.FACE
        )
        # Set the cell's name and instantiate the corresponding 'GenericCell'
        set_shape_name(cr_pin_compound, f"CR_pin_compound {i}")
        cr_pin_cell = GenericCell(cr_pin_compound)
        cr_pin_cell.set_properties(
            {PropertyType.MATERIAL: ["ABSORBER", "GAP", "CR_CLADDING2"]}
        )
        # Dummy assignement of the cell's type so that it can be added to a
        # hexagonal lattice
        cr_pin_cell.cell_type = CellType.HEX
        # Add the cell to the returned list
        cells.append(cr_pin_cell)
    return cells


# -------------------------------------------------------------------------- #
# FUEL ASSEMBLY CONSTRUCTION                                                 #
# -------------------------------------------------------------------------- #
# Build the hexagonal cells of the fuel assembly
fuel_cell = HexCell(name="Cartesian cell")
fuel_cell.rotate(90)
radii = [0.2, 0.6, 0.62, 0.68]
for radius in radii:
    fuel_cell.add_circle(radius)
# Assign the materials to each zone in the cell
fuel_cell.set_properties(
      {PropertyType.MATERIAL: ["GAP", "FUEL", "GAP", "CLADDING", "COOLANT"]}
)
central_cell = HexCell(name="Central cell")
central_cell.rotate(90)
for radius in [0.6, 0.65]:
    central_cell.add_circle(radius)
# Assign the materials to each zone in the cell
central_cell.set_properties(
      {PropertyType.MATERIAL: ["GAP", "CLADDING", "COOLANT"]}
)
# Update the viewer showing the two cells with the MATERIAL color map
fuel_cell.show(PropertyType.MATERIAL)
central_cell.show(PropertyType.MATERIAL)

# Build the fuel assembly lattice of the colorset
fuel_assembly = Lattice([central_cell], "Fuel Assembly")
fuel_assembly.add_rings_of_cells(fuel_cell, 5)
# Build the fuel assembly box
fuel_assembly.build_lattice_box([-0.1, 0.3, 0.3])
fuel_assembly.set_lattice_box_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CLADDING", "COOLANT"]}
)
# Display the fuel assembly with the MATERIAL color map
fuel_assembly.show(PropertyType.MATERIAL)

# -------------------------------------------------------------------------- #
# CONTROL ROD ASSEMBLY CONSTRUCTION                                          #
# -------------------------------------------------------------------------- #
# Data
pitch = fuel_assembly.lattice_box.figure.ly * 2
edge_bypass = (pitch) / math.sqrt(3)
edge_ext_wrap_o = (pitch - 0.4) / math.sqrt(3)
edge_ext_wrap_i = (pitch - 0.6) / math.sqrt(3)
r_cr_circles_i = 3.2
r_cr_circles_o = 5.2
cr_pin_radii = [0.68, 0.7, 0.78]
cr_wrapper_radii = [7, 7.25]
int_shaft_ir = 1.4
int_shaft_or = 1.7

# Build the control rod assembly cells
cr_cell = HexCell(edge_length=edge_ext_wrap_i, name= "Control Rod cell")
# Add the circles representing the different zones
for r in cr_wrapper_radii:
    cr_cell.add_circle(r)
cr_cell.set_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CR_CLADDING", "CR_MIX"]}
)

# Build the control rod box cell
box_cell = HexCell(edge_length=edge_bypass, name='Box cell')
wrapper_i = Hexagon(
    edge_length=edge_ext_wrap_i, name="External Wrapper Inner")
wrapper_o = Hexagon(
    edge_length=edge_ext_wrap_o, name="External Wrapper Outer")
box_face = make_partition(
    [box_cell.face], [wrapper_i.face, wrapper_o.face], ShapeType.FACE)
box_cell.update_geometry_from_face(GeometryType.TECHNOLOGICAL, box_face)
box_cell.set_properties(
    {PropertyType.MATERIAL: ["CR_MIX", "CR_CLADDING", "COOLANT"]})
box_cell.show(PropertyType.MATERIAL)

# Build the control rod assembly
cr_assembly = Lattice([cr_cell], "Control Rod Assembly")
cr_assembly.lattice_box = box_cell
cr_assembly.show(PropertyType.MATERIAL)

# Build the vertices at which the control rod circles are placed
cr_vertices_i = create_vertices_list(r_cr_circles_i, 6)
cr_vertices_o = create_vertices_list(r_cr_circles_o, 12)
# Build the circular cells and add them to the control rod assembly
cr_cells_i = make_circular_cells_list(cr_vertices_i, cr_pin_radii)
cr_cells_o = make_circular_cells_list(cr_vertices_o, cr_pin_radii)
for cell in cr_cells_i + cr_cells_o:
    cr_assembly.add_cell(cell, ())

# Build the central shaft cell as made by three concentric circles
circle_shaft_i = Circle(
    radius=int_shaft_ir, name="Inner Shaft Circle"
)
circle_shaft_o = Circle(
    radius=int_shaft_or, name="Outer Shaft Circle"
)
shaft_compound = make_partition(
    [circle_shaft_o.face, circle_shaft_i.face], [], ShapeType.FACE)
set_shape_name(shaft_compound, "Shaft Cell")
shaft_cell = GenericCell(shaft_compound)
shaft_cell.set_properties(
    {PropertyType.MATERIAL: ["COOLANT", "CR_CLADDING"]}
)
# Dummy assignement of the cell's type so that it can be added to a hexagonal
# lattice
shaft_cell.cell_type = CellType.HEX
cr_assembly.add_cell(shaft_cell, ())
# Display the control rod assembly
cr_assembly.show(PropertyType.MATERIAL)

# -------------------------------------------------------------------------- #
# COLORSET CONSTRUCTION                                                      #
# -------------------------------------------------------------------------- #
# Translate the fuel assembly to the right of the control rod assembly
fuel_assembly.translate((3/2*cr_assembly.lx, cr_assembly.ly, 0))

# Build the colorset as a list of the two assemblies
colorset = [cr_assembly, fuel_assembly]
# Build the colorset compound and display it in the SALOME 3D viewer
colorset_cmpd = make_compound([lattice.lattice_cmpd for lattice in colorset])
add_to_study(colorset_cmpd, "Colorset")

# -------------------------------------------------------------------------- #
# COLORSET S30 SYMMETRY CONSTRUCTION                                         #
# -------------------------------------------------------------------------- #
# Extract the S30 symmetry portion out of the entire colorset
edges = build_contiguous_edges(
    [
        make_vertex((0.0, 0.0, 0.0)),
        make_vertex((3/2*fuel_assembly.lx, fuel_assembly.ly, 0.0)),
        make_vertex((2*fuel_assembly.lx, 0.0, 0.0))
    ]
)
cutting_face = make_face(edges)
colorset_portion = make_common(colorset_cmpd, cutting_face)
add_to_study(colorset_portion, "Colorset - S30 Symmetry")

t1 = time.time()
print(f"--- Geometry generated in {t1 - t0} s. ---")

# -------------------------------------------------------------------------- #
# COLORSET REGIONS VISUALIZATION                                             #
# -------------------------------------------------------------------------- #
# Display all the regions of the two assemblies by using the same colormap
colorset_regions: List[Region] = []
material_names: List[str] = []
for lattice in colorset:
    for region in lattice.regions:
        if get_min_distance(colorset_portion, region.face) > 1e-6:
             continue
        new_region = make_common(region.face, colorset_portion)
        # Continue with the next region if the common shape does not hold
        # any face, meaning that the region and the portion do not overlap
        if not extract_sub_shapes(make_compound([new_region]),
                                  ShapeType.FACE):
           continue
        # Store the material name, if not present
        mat_name = region.properties.get(PropertyType.MATERIAL)
        if mat_name is None:
            raise RuntimeError(f"No material for region {region}.")
        if mat_name not in material_names:
            material_names.append(mat_name)
        # If the result is a compound or a shell, extract the contained faces
        if get_shape_type(new_region) in [ShapeType.COMPOUND,
                                          ShapeType.SHELL]:
            for new_region in extract_sub_shapes(new_region, ShapeType.FACE):
                colorset_regions.append(
                    Region(
                        new_region,
                        name=region.name,
                        properties=deepcopy(region.properties)
                    )
                )
            continue
        colorset_regions.append(
            Region(
                new_region,
                name=region.name,
                properties=deepcopy(region.properties)
            )
        )

# Generate a specific amount of colors as the number of different
# values for the same given property type
colors = generate_unique_random_colors(len(material_names))
# Join the material names with the colors
materials_vs_color = dict(zip(material_names, colors))

# Display the regions of the colorset portion with the material color map
for region in colorset_regions:
    # Get the color according to the material name of the region
    region.color = materials_vs_color[
        region.properties[PropertyType.MATERIAL]
    ]
    set_color_face(region.face, region.color)
    add_to_study_in_father(colorset_portion, region.face, region.name)

update_salome_study()

t2 = time.time()
print(f"--- Geometry displayed in {t2 - t1} s. ---")

# Generate the TDT file from the colorset portion using a specific typgeo and
# symmetry type
analyse_and_generate_tdt(
    colorset,
    "colorset_s30",
    tdt_config=TdtSetup(
        type_geo=LatticeGeometryType.SYMMETRIES_TWO,
        symmetry_type=SymmetryType.TWELFTH),
    compound_to_export=colorset_portion)

print(f"--- Script executed in {time.time() - t0} s. ---")
