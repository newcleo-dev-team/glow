from copy import deepcopy
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import Lattice
from glow.main import TdtSetup, analyse_and_generate_tdt
from glow.support.types import LatticeGeometryType, PropertyType, SymmetryType
from glow.interface.geom_interface import *
from glow.support.utility import build_contiguous_edges

# Build a hexagonal cell
cell = HexCell(name="Cartesian cell")
radii = [0.2, 0.5, 0.6, 0.7]
for radius in radii:
    cell.add_circle(radius)
# Assign the materials to each zone in the cell
cell.set_properties(
      {PropertyType.MATERIAL: ["GAP", "FUEL", "GAP", "CLADDING", "COOLANT"]}
)

# Build the lattices of the colorset as copies of the same one translated
# accordingly
lattice = Lattice([cell], 'Hexagonal Lattice')
lattice.add_rings_of_cells(cell, 2)
lattice.build_lattice_box([0.15, 0.05])
lattice.set_lattice_box_properties({PropertyType.MATERIAL: ["COOLANT", "CLADDING", "COOLANT"]})

lattice2 = deepcopy(lattice)
lattice2.translate((lattice.lx, 3/2*lattice.ly, 0))
lattice3 = deepcopy(lattice)
lattice3.translate((lattice.lx, -3/2*lattice.ly, 0))
lattice4 = deepcopy(lattice)
lattice4.translate((2*lattice.lx, 0, 0))
# Build a list of the lattices belonging to the colorset
blocks = [lattice, lattice2, lattice3, lattice4]
# Assemble the compound of the colorset so that it can be shown in the SALOME
# viewer
colorset = make_compound([lattice.lattice_cmpd for lattice in blocks])
add_to_study(colorset, "Colorset")
update_salome_study()

# Build the GEOM elements to extract the portion of the colorset to export
edges = build_contiguous_edges(
    [
        make_vertex((0.0, 0.0, 0.0)),
        make_vertex((2*lattice.lx, 0.0, 0.0)),
        make_vertex((2*lattice.lx, lattice.ly, 0.0, 0.0))
    ]
)
cutting_face = make_face(edges)
add_to_study(cutting_face, "Cutting face")
# Extract the colorset portion as the common part between the colorset and the
# cutting shape
colorset_portion = make_common(colorset, cutting_face)
add_to_study(colorset_portion, "Colorset portion")

# Generate the TDT file from the colorset portion using a specific typgeo and
# symmetry type
analyse_and_generate_tdt(
    blocks,
    "colorset_s30",
    tdt_config=TdtSetup(
        type_geo=LatticeGeometryType.SYMMETRIES_TWO,
        symmetry_type=SymmetryType.TWELFTH),
    compound_to_export=colorset_portion)
