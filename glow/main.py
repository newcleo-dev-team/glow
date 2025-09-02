import math
import os

from dataclasses import dataclass
from pathlib import Path
from typing import Any, List

from glow.support.types import GeometryType, PropertyType
from glow.geometry_layouts.lattices import Lattice


@dataclass
class TdtSetup:
    """
    Dataclass holding the settings for configuring how to export the TDT
    file for the lattice's geometry layout.

    Notes
    -----
    The `albedo` attribute can have values between ``0.0`` and ``1.0``, with
    the latter case indicating the `ALBE 1.0` BC used in DRAGON5 with a
    uniform tracking (i.e. ``LatticeGeometryType.ISOTROPIC``). If ``None``,
    the value that corresponds to the type of geometry of the lattice is used,
    i.e. ``0.0`` for values of ``LatticeGeometryType`` greater than zero,
    ``1.0`` otherwise.
    """
    geom_type: GeometryType = GeometryType.TECHNOLOGICAL
    """Identifying the type of geometry of the lattice's cells."""
    property_type: PropertyType = PropertyType.MATERIAL
    """Identifying the type of property associated to lattice's regions."""
    albedo: float | None = None
    """Identifying the value for the albedo applied to the lattice's BCs."""

    def __post_init__(self) -> None:
        """
        Method run after the dataclass initialization. It checks whether the
        value of the ``albedo`` attribute is in the 0.0-1.0 range. If not, an
        exception is raised.

        Raises
        ------
        RuntimeError
            If the ``albedo`` attribute is not in the 0.0-1.0 range.
        """
        if self.albedo is not None:
            if not (
                0.0 <= self.albedo <= 1.0 or
                math.isclose(self.albedo, 0.0) or
                math.isclose(self.albedo, 1.0)
            ):
                raise RuntimeError(
                    f"The value {self.albedo} for the albedo is out of "
                    "bounds [0.0, 1.0]")


def analyse_and_generate_tdt(
        lattices: List[Lattice],
        filename: str,
        tdt_config: TdtSetup = TdtSetup(
            GeometryType.TECHNOLOGICAL,
            PropertyType.MATERIAL,
            None),
        compound_to_export: Any | None = None
    ) -> None:
    """
    Function that analyses the given lattices, as instance of the ``Lattice``
    class, to extract information about the characteristics of its geometry
    and the properties associated to its regions.
    A TDT file, whose name is provided as second parameter, is generated,
    collecting all this information.
    Users can also specify to analyse the lattice according to:

    - the geometry type of its cells;
    - the type of property associated to the lattice regions;
    - the value for the albedo applied to the lattice's BCs. If ``None``, a
      default value that corresponds to the lattice's geometry type is
      adopted.

    If the ``compound_to_export`` parameter is provided, it will be the one
    to be analysed and exported, according to the information stored in the
    provided lattices. In this way, a colorset of the lattices can be treated.

    Parameters
    ----------
    lattices : List[Lattice]
        The object storing the information about the geometry and the
        properties of the lattices.
    filename : str
        The name of the output TDT file.
    tdt_config : TdtSetup
        Dataclass providing the settings for exporting the TDT representation
        of the geometry layout of the lattice.
    compound_to_export: Any | None = None
        The compound object to analyse and export to TDT, if present. If
        ``None`` is given, the lattices are considered instead.
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
    data_extractor = analyse_lattice(lattices, tdt_config, compound_to_export)

    t1 = time.time()
    print(f"--- Lattice analysis executed in {t1 - start_time} seconds ---")

    # Instantiate the dataclass storing the needed lattice data
    tdt = TdtData(filename=os.path.join(
                      Path(os.getcwd()).resolve().parent,
                      filename + ".dat"),
                  edges=data_extractor.edges,
                  faces=data_extractor.subfaces,
                  boundaries=data_extractor.boundaries,
                  type_geo=data_extractor.lattices[0].type_geo,
                  type_sym=data_extractor.lattices[0].symmetry_type,
                  albedo=tdt_config.albedo)

    print("--- TdtData class instantiation executed in " + \
          f"{time.time() - t1} seconds ---")

    # Write the data to the output TDT-format file
    write_tdt_file(tdt)

    print("--- TDT file generation executed in " + \
          f"{time.time() - start_time} seconds ---")
