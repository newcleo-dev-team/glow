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
import math
import time

from glow.geometry_layouts.geometries import Circle, Hexagon
from glow.geometry_layouts.layouts import Region
from glow.main import TdtSetup, export_layout_to_tdt
from glow.interface.geom_interface import *
from glow.geometry_layouts.cells import Cell, HexCell
from glow.geometry_layouts.lattices import HexLattice, Lattice
from glow.support.types import *
from glow.support.utility import build_compound_borders

# -------------------------------------------------------------------------- #
#                                 FUNCTIONS                                  #
# -------------------------------------------------------------------------- #
def add_circular_regions(
        cell: Cell, radii: List[float], materials: List[float]
    ) -> None:
    """
    Function that adds circular ``Region`` objects to the given ``Cell``
    instance. Regions are characterised in terms of the radius and the
    material property.

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


def get_modified_cells(lattice: Lattice) -> List[Cell]:
    """
    Function that returns a list of ``Cell`` objects belonging to the given
    ``Lattice`` instance. These cells have their geometry layout changed
    compared to their original one.
    For each cell in each layer, the current shape of the cell is built and
    its area compared with the area of its specific characteristic figure (as
    a ``Surface`` instance).
    Those showing a different value for the area indicates a change in their
    geometry layout has occurred and are collected into the returned list.

    This function retrieves cells that have been modified within the lattice,
    such as by overlap with a superior layer of cells.

    Parameters
    ----------
    lattice : Lattice
        The lattice instance to check for cells whose geometry layout has
        changed.

    Returns
    -------
    List[Cell]
        A list of ``Cell`` objects whose geometry layout differs from their
        original one.
    """
    cells = []
    lattice.update_hierarchical_structure()
    for layer in lattice.layers:
        for layout in layer:
            if isinstance(layout, Region):
                continue
            # Build a face object over the borders of the cell
            cell_shape = make_face(build_compound_borders(layout))
            # Compare the area of the current cell's face with the one of its
            # original shape
            if not math.isclose(
                get_basic_properties(cell_shape)[1],
                get_basic_properties(layout.shape)[1],
                abs_tol=1e-6
            ):
                cells.append(layout)
    # Return the list of changed cells
    return cells


t0 = time.time()
# ----------------------------------------------------------------------
# Build the hexagonal cell that constitutes the lattice.
cell_1 = HexCell(
    name="Cell 1",
    base_props={PropertyType.MATERIAL: "COOLANT"}
)
cell_1.rotate(90)
add_circular_regions(
    cell_1, [0.1, 0.6, 0.625, 0.70], ["GAP", "FUEL", "GAP", "CLADDING"]
)

# ----------------------------------------------------------------------
# Build the second hexagonal cell
cell_2 = HexCell(
    side=2.0,
    name="Cell 2",
    base_props={PropertyType.MATERIAL: "COOLANT"}
)
add_circular_regions(cell_2, [1.0, 1.25], ['COOLANT', 'CLADDING'])

# ----------------------------------------------------------------------
# Build the lattice and add both types of cells
lattice = HexLattice(cells=[cell_1], name="HexLattice")
lattice.add_rings_of_cells(cell_1, 6, 1)
# XY coordinates of the centres of the cells with greater size
x = 4.330127
y = 4.5
lattice.add(cell_2, ())
lattice.add(cell_2, (x, y, 0.0))
lattice.add(cell_2, (-x, y, 0.0))
lattice.add(cell_2, (x, -y, 0.0))
lattice.add(cell_2, (-x, -y, 0.0))
# Show the lattice's technological geometry with the 'MATERIAL' colour map
lattice.show(PropertyType.MATERIAL)

# Get the cells whose geometry layout has been cut and restore them by
# assigning a specific property type
for cell in get_modified_cells(lattice):
    cell.restore()
    cell.regions[0].properties = {PropertyType.MATERIAL: 'COOLANT'}

# Build the assembly cell from the lattice's dimensions + the sum of the box's
# layers thickness
box_thicknesses = [0.15, 0.15]
assembly = HexCell(
    side=lattice.dimensions[0] + sum(box_thicknesses),
    name="Hexagonal Assembly",
    base_props={PropertyType.MATERIAL: 'COOLANT'}
)
# Add the box layers first, then the lattice
assembly.add(
    Region(
        Hexagon(edge_length=lattice.dimensions[0]+box_thicknesses[0]),
        properties={PropertyType.MATERIAL: 'CLADDING'}
    )
)
assembly.add(
    Region(
        Hexagon(edge_length=lattice.dimensions[0]),
        properties={PropertyType.MATERIAL: 'COOLANT'}
    )
)
assembly.add(lattice)

# Show the lattice's technological geometry
assembly.show(PropertyType.MATERIAL)

# ----------------------------------------------------------------------
# Perform the geometry analysis and export the TDT file of the surface
# geometry
export_layout_to_tdt(
    assembly,
    'hexagonal_assembly',
    TdtSetup(type_geo=LayoutGeometryType.HEXAGON_TRAN)
)

print(f"Script executed in {time.time() - t0} s.")
