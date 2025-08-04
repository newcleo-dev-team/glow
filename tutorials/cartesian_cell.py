"""
Use case showing the construction of a cartesian cell and how to customize its
geometry layout. The 'MATERIAL' property is assigned to each region of the
cell technological geometry. A sectorization is applied to the cell and the
result graphically shown in the SALOME 3D viewer.
In addition, the sectorized geometry is modified by means of the functions
that wrap the ones of the *GEOM* module of *SALOME*. The result is a sectorized
geometry costituted by more circles between the regions of the technological
geometry.
"""
from glow.geometry_layouts.cells import RectCell
from glow.support.types import GeometryType, PropertyType
from glow.geometry_layouts.geometries import Circle
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

# Build the cell's sectorized geometry with 'windimill' option enabled
cell.sectorize([1, 1, 4, 8], [0, 0, 0, 22.5], windmill=True)
# Show the sectorized cell with regions colored according to the 'MATERIAL'
# property
cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)

# ---------------------------------------------------------------------
# Update the cell's sectorized geometry with a face built with SALOME's
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

# Partition the original cell's technological geometry with all the circles
updated_face = make_partition(
    [cell.face], [c.face for c in circles], ShapeType.FACE)
# Update the cell's sectorized geometry with the just built shape
cell.update_geometry_from_face(GeometryType.SECTORIZED, updated_face)
# Show the result in the 3D viewer
cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
