"""
Use case showing the construction of a cartesian cell and how to customise its
geometry layout. The 'MATERIAL' property is assigned to each region of the
cell technological geometry. A sectorisation is applied to the cell and the
result graphically shown in the SALOME 3D viewer.
In addition, the sectorised geometry is modified by means of the functions
that wrap the ones of the *GEOM* module of *SALOME*. The result is a sectorised
geometry costituted by more circles between the regions of the technological
geometry.
"""
from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_entities import wrap_shape
from glow.support.types import GeometryType, PropertyType
from glow.geometry_layouts.geometries import Circle
from glow.interface.geom_interface import *


# Intialise three lists one storing the circular regions radii, the other the
# names of the materials, sorted from the inner to the outer region
radii = [0.2, 0.3, 0.4]
materials = ["MAT_1", "MAT_2", "MAT_3"]
# Build the cell's geometry layout by adding three circular regions from the
# outer to the inner
cell = CartesianCell(
    name="Cartesian cell", base_props={PropertyType.MATERIAL: "MAT_4"}
)
for radius, mat in zip(radii[::-1], materials[::-1]):
    cell.add(
        Region(Circle(radius=radius), properties={PropertyType.MATERIAL: mat})
    )

# Build the cell's sectorised geometry with 'windimill' option enabled
cell.sectorize([1, 1, 4, 8], [0, 0, 0, 22.5], windmill=True)
# Show the sectorised cell with regions colored according to the 'MATERIAL'
# property
cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)

# ---------------------------------------------------------------------
# Update the cell's sectorised geometry with a face built with SALOME's
# functions
# ---------------------------------------------------------------------
# Setup the XYZ coordinates of the centres of the circles
centres = [(0, 0.1, 0), (0, -0.1, 0), (0.1, 0, 0), (-0.1, 0, 0)]
# Build the corresponding 'Circle' objects, all with the same radius
circles = [Circle(centre, radius=0.05) for centre in centres]
# Build circles positioned in the cell centre
center_circles = [Circle(radius=r) for r in [0.32, 0.34, 0.36, 0.38]]
# Update the list of 'Circle' objects
circles += center_circles

# Build a compound from the edges of the circles
circles_cmpd = make_compound([c.borders[0] for c in circles])
add_to_study(circles_cmpd, "")
# Update the cell's sectorised geometry with the compound of edges
cell.geometry_maps[GeometryType.SECTORIZED] = wrap_shape(
    make_compound(
        [cell.geometry_maps[GeometryType.SECTORIZED], circles_cmpd]
    )
)
# Show the result in the 3D viewer
cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
