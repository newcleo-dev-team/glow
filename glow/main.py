import math
import os

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict

from glow.geometry_layouts.cells import CartesianCell, Cell, HexCell
from glow.geometry_layouts.fillable_layouts import Fillable
from glow.support.types import GeometryType, LayoutGeometryType, LayoutType, \
    PropertyType, SymmetryType
from glow.geometry_layouts.lattices import CartesianLattice, HexLattice, \
    Lattice
from glow.support.utility import check_type_geo_consistency


# Dictionary associating the type of layout to each of the available classes
# describing either a cell or a lattice.
CLASS_VS_LAYOUT_TYPE: Dict[Fillable, LayoutType] = {
    CartesianCell: LayoutType.RECT,
    CartesianLattice: LayoutType.RECT,
    Cell: LayoutType.GENERIC,
    HexCell: LayoutType.HEX,
    HexLattice: LayoutType.HEX,
    Lattice: LayoutType.GENERIC
}

@dataclass
class TdtSetup:
    """
    Dataclass holding the settings for configuring how to export the TDT
    file for the current geometry layout.

    Notes
    -----
    The `albedo` attribute can have values between ``0.0`` and ``1.0``, with
    the latter case indicating the `ALBE 1.0` BC used in DRAGON5 with a
    uniform tracking (i.e. ``LayoutGeometryType.ISOTROPIC``). If ``None``,
    the value that corresponds to the type of geometry of the layout is used,
    i.e. ``0.0`` for values of ``LayoutGeometryType`` greater than zero,
    ``1.0`` otherwise.
    """
    geom_type: GeometryType = GeometryType.TECHNOLOGICAL
    """Identifying the type of geometry of the layout."""
    property_type: PropertyType = PropertyType.MATERIAL
    """Identifying the type of property associated to layout's regions."""
    albedo: float | None = None
    """Identifying the value for the albedo applied to the layout's BCs."""
    type_geo: LayoutGeometryType = LayoutGeometryType.ISOTROPIC
    """Identifying the value for the typegeo related to the layout."""
    symmetry_type: SymmetryType = SymmetryType.FULL
    """Identifying the value for the symmetry type applied to the layout."""
    layout_type: LayoutType = field(init=False, repr=False)
    """Identifying the type of the layout."""

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


def export_layout_to_tdt(
        layout: Fillable,
        filename: str,
        tdt_setup: TdtSetup = TdtSetup(
            GeometryType.TECHNOLOGICAL,
            PropertyType.MATERIAL,
            None),
        compound_to_export: Any | None = None
    ) -> None:
    """
    Function that analyses the given layout, as instance of the ``Fillable``
    class, to extract information about the characteristics of its geometry
    and the properties associated to its regions.
    A TDT file, whose name is provided as second parameter, is generated,
    collecting all this information.

    By properly configuring the ``TdtSetup`` instance, provided as third
    parameter, users can indicate which information about the geometry needs
    to be extracted from the layout and the tracking setup. In particular,
    the available options are:

    - the geometry type of the layout (either the technological or the refined
      geometry);
    - the type(s) of property associated to the regions of the layout;
    - the value for the albedo applied to the BCs of the layout. If ``None``,
      a default value that corresponds to the geometry type of the layout is
      adopted;
    - the value for the `type_geo` attribute which drives the type of tracking
      (either TISO or TSPC) to adopt accordingly with what requested by the
      `SALT:` module of DRAGON5;
    - the type of symmetry applied to the layout.

    When specifying values for the type of symmetry which differ from the one
    currently applied to the layout, and specified in the the `state`
    attribute of the given ``Fillable`` object, the indicated symmetry is
    applied, if the corresponding shape of symmetry has already been built
    by calling the ``apply_symmetry`` method.

    If the ``compound_to_export`` parameter is provided, it will be the one
    to be analysed and exported, according to the property information stored
    in the layout. The indicated compound object must be a portion of the
    entire layout, otherwise the successive steps of the analysis will fail.

    Parameters
    ----------
    layout : Fillable
        The layout, as instance of the ``Fillable`` subclasses, storing the
        information about the geometry and the assigned properties.
    filename : str
        The name of the output TDT file.
    tdt_setup : TdtSetup
        Dataclass providing the settings for exporting the TDT representation
        of the geometry layout.
    compound_to_export: Any | None = None
        The compound object (as portion of the layout) to analyse and export
        to TDT, if provided.

    Raises
    ------
    RuntimeError
        In case of inconsistencies in the values of the ``TdTSetup`` settings.
        If the analysis fails due to the compound not being part of the
        layout, no properties being found for any region of the layout, or the
        layout borders being impossible to obtain.

    Notes
    -----
    Users should note that:

    - The configuration values provided in the ``TdtSetup`` instance are
      considered regardless of what set in the given ``Fillable`` object.
    - The configuration values provided in the ``TdtSetup`` instance must
      match with the indicated compound object. If values that do not match
      with the shape of the compound are provided, the validity of the results
      in DRAGON cannot be assured.
    """
    # Import the 'time' module for evaluating the analysis performance
    import time
    # Get the start time
    start_time = time.time()

    # Import the classes and functions for performing the geometry conversion
    from glow.generator.geom_extractor import analyse_layout
    from glow.generator.generator import TdtData, write_tdt_file

    # Get the type of the layout and set the corresponding attribute in the
    # given 'TdtSetup' instance
    tdt_setup.layout_type = CLASS_VS_LAYOUT_TYPE[layout.__class__]

    # Check the correctness of the 'TdtSetup' settings
    check_type_geo_consistency(
        tdt_setup.type_geo, tdt_setup.layout_type, tdt_setup.symmetry_type
    )

    # Perform the analysis on the layout according to the given settings
    data_extractor = analyse_layout(layout, tdt_setup, compound_to_export)

    t1 = time.time()
    print(f"--- Lattice analysis executed in {t1 - start_time} seconds ---")

    # Instantiate the dataclass storing the needed data of the layout
    tdt = TdtData(
        filename=Path(os.getcwd()).resolve().parent / f"{filename}.dat",
        edges=data_extractor.edges,
        faces=data_extractor.subfaces,
        boundaries=data_extractor.boundaries,
        type_geo=tdt_setup.type_geo,
        type_sym=tdt_setup.symmetry_type,
        albedo=tdt_setup.albedo
    )

    print("--- TdtData class instantiation executed in " + \
          f"{time.time() - t1} seconds ---")

    # Write the data to the output TDT-format file
    write_tdt_file(tdt)

    print("--- TDT file generation executed in " + \
          f"{time.time() - start_time} seconds ---")
