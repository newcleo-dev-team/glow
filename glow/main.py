import os
from pathlib import Path

from glow.geometry_layouts.lattices import Lattice


def analyse_and_generate_tdt(lattice: Lattice, filename: str):
    """
    Function that performs the analysis onto the given instance of the
    'Lattice' class, representing the geometry and the properties of the
    lattice for which the TDT file must be generated.

    Parameters
    ----------
    lattice   : Lattice
                The object storing the geometry and the properties of the
                lattice
    filename  : str
                The name of the output TDT file
    """
    # TODO Add a parameter that indicates the type of lattice geometry to
    # analyse (TECHNOLOGICAL or SECTORIZED); according to its value, the
    # regions are extracted
    import time
    # Get the start time
    start_time = time.time()

    # Import the classes and functions for performing the geometry conversion
    from generator.geom_extractor import analyse_lattice
    from generator.generator import TdtData, write_tdt_file

    # Perform the lattice faces and edges analysis
    data_extractor = analyse_lattice(lattice)

    t1 = time.time()
    print(f"--- Lattice analysis executed in {t1 - start_time} seconds ---")

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
          f"{time.time() - t1} seconds ---")

    # Write the data to the output TDT-format file
    write_tdt_file(tdt)

    print("--- TDT file generation executed in " + \
          f"{time.time() - start_time} seconds ---")
