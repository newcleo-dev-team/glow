import os
from pathlib import Path

from glow.generator.support import GeometryType, LatticeGeometryType, PropertyType
from glow.geometry_layouts.lattices import Lattice


def analyse_and_generate_tdt(
        lattice: Lattice,
        filename: str,
        geom_type: GeometryType = GeometryType.TECHNOLOGICAL,
        property_type: PropertyType = PropertyType.MATERIAL) -> None:
    """
    Function that analyses the given lattice, as instance of the `Lattice`
    class, to extract information about the characteristics of its geometry
    and the properties associated to its regions.
    A TDT file, whose name is provided as second parameter, is generated,
    collecting all this information.
    Users can also specify to analyse the lattice according to:

    - the geometry type of its cells;
    - the type of property associated to the lattice regions.

    Parameters
    ----------
    lattice : Lattice
        The object storing the information about the geometry and the
        properties of the lattice
    filename : str
        The name of the output TDT file
    geom_type : GeometryType = GeometryType.TECHNOLOGICAL
        The type of geometry of the lattice cells to use in the analysis
    property_type : PropertyType = PropertyType.MATERIAL
        The type of property associated to the lattice regions to use in
        the analysis
    """
    # Import the 'time' module for evaluating the analysis performance
    import time
    # Get the start time
    start_time = time.time()

    # Import the classes and functions for performing the geometry conversion
    from glow.generator.geom_extractor import analyse_lattice
    from glow.generator.generator import TdtData, write_tdt_file

    # Perform the lattice faces and edges analysis for the given geometry and
    # property types
    data_extractor = analyse_lattice(lattice, geom_type, property_type)

    print("--- Lattice analysis executed in " + \
          f"{time.time() - start_time} seconds ---")

    # Instantiate the dataclass storing the needed lattice data
    tdt = TdtData(filename=os.path.join(
                      Path(os.getcwd()).resolve().parent,
                      filename + ".dat"),
                  edges=data_extractor.edges,
                  faces=data_extractor.subfaces,
                  boundaries=data_extractor.boundaries,
                  type_geo=data_extractor.lattice.type_geo,
                  type_sym=data_extractor.lattice.symmetry_type)

    print("--- TdtData class instantiation executed in " + \
          f"{time.time() - start_time} seconds ---")

    # Write the data to the output TDT-format file
    write_tdt_file(tdt)

    print("--- TDT file generation executed in " + \
          f"{time.time() - start_time} seconds ---")
