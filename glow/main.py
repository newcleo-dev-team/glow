import math
import os

from dataclasses import dataclass
from pathlib import Path
from typing import Any, List

from glow.support.types import GeometryType, LatticeGeometryType, \
    PropertyType, SymmetryType
from glow.geometry_layouts.lattices import Lattice
from glow.support.utility import check_type_geo_consistency


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
    property_types: PropertyType | List[PropertyType] = PropertyType.MATERIAL
    """Identifying the type(s) of property associated to lattice's regions."""
    albedo: float | None = None
    """Identifying the value for the albedo applied to the lattice's BCs."""
    type_geo: LatticeGeometryType = LatticeGeometryType.ISOTROPIC
    """Identifying the value for the typegeo related to the layout."""
    symmetry_type: SymmetryType = SymmetryType.FULL
    """Identifying the value for the symmetry type applied to the lattice."""

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
        # Transform the given property type to a list
        if not isinstance(self.property_types, List):
            self.property_types = [self.property_types]


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
      adopted;
    - the value for the `typegeo` parameter;
    - the type of symmetry applied to the single lattice or the colorset.

    When exporting one or more ``Lattice`` instances provided by the input
    list, the `typegeo` and the symmetry type values of the ``TdtSetup``
    instance are neglected. The values of the corresponding attributes for
    the first lattice in the list (taken as reference) are considered instead.

    If the ``compound_to_export`` parameter is provided, it will be the one
    to be analysed and exported, according to the property information stored
    in the provided lattices. The indicated compound object must be a portion
    of the lattices, otherwise the successive steps of the analysis will fail.
    The values for the `typegeo` parameter and the symmetry type provided in
    the ``TdtSetup`` instance are considered regardless of what set in the
    list of involved lattices, even if only one lattice is present.
    Users should note that these two values must match with the provided
    compound. If values that do not match with the shape of the compound are
    provided, the validity of the results in DRAGON cannot be assured.

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

    Raises
    ------
    RuntimeError
        When multiple lattices, and no compound, are provided and they do
        not have the same ``SymmetryType.FULL`` symmetry.
        When a compound is provided and the corresponding lattices do not
        have the same ``SymmetryType.FULL`` symmetry.
        In case of inconsistencies in the values of the ``TdTSetup``
        settings.
        If the analysis fails due to the compound not been part of the
        lattices, no property found for a region, the impossibility to get
        the borders of the layout to export.
    """
    # Import the 'time' module for evaluating the analysis performance
    import time
    # Get the start time
    start_time = time.time()

    # Import the classes and functions for performing the geometry conversion
    from glow.generator.geom_extractor import analyse_lattice
    from glow.generator.generator import TdtData, write_tdt_file

    # The reference lattice is the first instance in the given list
    ref_lattice = lattices[0]
    # Choose whether to use data from the 'TdtSetup' instance or the lattice;
    # in the latter case, update the 'TdtSetup' instance values accordingly
    if compound_to_export is None:
        # If multiple lattices are provided, they must not have any symmetry
        if len(lattices) > 1 and not all(
            lattice.symmetry_type == SymmetryType.FULL
                for lattice in lattices
        ):
            raise RuntimeError(
                "When considering a colorset, the type of symmetry of all "
                "the involved lattices must be 'SymmetryType.FULL'. If a "
                "portion of the colorset is meant to be considered, run the "
                "analysis with a compound from the lattices.")
        # Modify the typegeo and symmetry type values of the 'TdtSetup'
        # instance
        tdt_config.symmetry_type = ref_lattice.symmetry_type
        tdt_config.type_geo = ref_lattice.type_geo
    else:
        # If a compound is provided, the lattices must not have any symmetry
        if not all(
            lattice.symmetry_type == SymmetryType.FULL
                for lattice in lattices
        ):
            raise RuntimeError(
                "When a portion of the lattices is provided, the type of "
                "symmetry of all the involved lattices must be "
                "'SymmetryType.FULL'.")
        # Check the correctness of the 'TdtSetup' settings
        check_type_geo_consistency(
            tdt_config.type_geo,
            ref_lattice.cells_type,
            tdt_config.symmetry_type
        )

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
                  type_geo=tdt_config.type_geo,
                  type_sym=tdt_config.symmetry_type,
                  albedo=tdt_config.albedo)

    print("--- TdtData class instantiation executed in " + \
          f"{time.time() - t1} seconds ---")

    # Write the data to the output TDT-format file
    write_tdt_file(tdt)

    print("--- TDT file generation executed in " + \
          f"{time.time() - start_time} seconds ---")
