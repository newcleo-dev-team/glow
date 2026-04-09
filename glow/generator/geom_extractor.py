"""
Module containing the class that deals with the extraction of the geometric
information from the geometry layout. The functionalities in this module serve
for preparing all the data for the output TDT file generation.
"""
import logging
import math

from typing import Any, Dict, List, Tuple

from glow.generator.export_data import BoundaryData, EdgeData, FaceData, \
    build_edge_id, classify_layout_edges
from glow.geometry_layouts.cells import Region
from glow.geometry_layouts.fillable_layouts import Fillable
from glow.geometry_layouts.layouts import build_compound_regions
from glow.interface.geom_entities import Compound, wrap_shape
from glow.interface.geom_interface import ShapeType, add_to_study, \
    extract_sub_shapes, get_bounding_box, get_in_place, \
    get_in_place_by_hystory, get_point_coordinates, get_shape_name, \
    get_shape_type, make_compound, make_face, make_partition, \
    make_partition_non_self_intersecting, make_vertex, update_salome_study
from glow.main import TdtSetup
from glow.support.types import GeometryType, LayoutGeometryType, LayoutType, \
    PropertyType, SymmetryType
from glow.support.utility import are_same_shapes, build_compound_borders, \
    translate_wrt_reference


class LayoutDataExtractor():
    """
    Class that extracts the geometric and properties data from a layout,
    either provided as an instance of the ``Fillable`` class only, or by
    additionally indicating a compound object representing a generic portion
    of the ``Fillable`` object.

    The extraction process relies on determining:

    - the association of GEOM faces with properties;
    - the association between GEOM edges and the GEOM faces they belong to;
    - the information about the borders of the layout in terms of their
      geometric data and applied BC types.

    Parameters
    ----------
    layout : Fillable
        A ``Fillable`` instance storing the geometric and properties data
        to extract.
    tdt_setup : TdtSetup
        Dataclass providing the settings for extracting the geometric and
        properties information from the given layout.
    compound_to_analyse: Any | None = None
        The compound object to analyse. If ``None`` is provided, the layout
        given as first parameter is considered instead.

    Attributes
    ----------
    geometry_layout : Fillable
        A ``Fillable`` instance storing the geometric and properties data
        to extract.
    borders : List[Any]
        The list of edge objects representing the borders of the layout.
    boundaries : List[BoundaryData]
        The list of ``BoundaryData`` objects storing the geometric and BCs
        information about the layout's borders.
    subfaces : List[FaceData]
        The list of ``FaceData`` objects storing the geometric and properties
        information about each region of the layout.
    edges : List[EdgeData]
        The list of ``EdgeData`` objects storing the geometric information
        about each edge of the layout and the faces they belong to.
    id_vs_edge : Dict[str, Any]
        A dictionary of the IDs for the layout's edges VS the corresponding
        edge objects.
    layout_edges : List[Any]
        The list of the edge objects contained in the layout.
    layout_centre : Tuple[float, float, float]
        The XYZ coordinates of the whole layout centre.
    dimensions : Tuple[float, float]
        The XY characteristic dimensions of the layout.
    """
    CASES_FOR_TRANSLATION = {
        LayoutType.GENERIC: [],
        LayoutType.RECT: [
            SymmetryType.FULL, SymmetryType.HALF, SymmetryType.DIAG],
        LayoutType.HEX: [SymmetryType.THIRD, SymmetryType.SIXTH]
    }
    """
    Identifying the symmetry types for each type of layout for which the
    layout translation should be evaluated.
    """

    def __init__(
            self,
            layout: Fillable,
            tdt_setup: TdtSetup,
            compound_to_analyse: Any | None = None
        ) -> None:
        # Raise an exception if the layout is empty
        if all(not layer for layer in layout.layers):
            raise RuntimeError(
                "No data extraction can be performed from an empty layout."
            )
        # Initialize the instance attributes
        self.geometry_layout: Fillable = layout.clone()
        self.borders: List[Any] = []
        self.layout_edges: List[Any] = []
        self.boundaries: List[BoundaryData] = []
        self.subfaces: List[FaceData] = []
        self.edges: List[EdgeData] = []
        self.regions: List[Region] = []
        self.layout_centre: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        self.dimensions: Tuple[float, float] = (0.0, 0.0)
        # Extract the information to be stored about borders and edges of the
        # layout according to the applied symmetry and geometry
        self._preprocess(tdt_setup, compound_to_analyse)
        # Associate each edge with an index and build a dictionary
        self.id_vs_edge: Dict[str, Any] = classify_layout_edges(
            self.layout_edges
        )

    def build_boundaries(self, type_geo: LayoutGeometryType) -> None:
        """
        Method that constructs a list of ``BoundaryData`` objects representing
        the layout boundary edges, i.e. those connected to a single face.
        All the GEOM edges being part of each boundary are associated to the
        same ``BoundaryData`` object.

        Parameters
        ----------
        type_geo : LatticeGeometryType
            Identifying the characteristic value associated to the type of
            symmetry and tracking.
        """
        # No boundaries to extract if an 'ISOTROPIC' type of geometry, meaning
        # 'VOID' or 'ALBE 1.0' BCs in DRAGON5.
        if type_geo == LayoutGeometryType.ISOTROPIC:
            return
        # Initialize the list of 'Boundary' objects
        boundary_edges = []

        # Loop through all the stored 'Edge' objects
        print("LEN EDGES", len(self.edges))
        for edge in self.edges:
            # Check if the current edge has no faces on the left or on the
            # right, i.e. it represents a boundary edge
            if not edge.left or not edge.right:
                # Append the corresponding GEOM edge object to the list
                boundary_edges.append(edge.edge)

        # Build a compound from all the edges placed on the layout borders
        boundary_edgs_cmpd = make_compound(boundary_edges)
        # Loop through all the edge objects representing the layout borders
        print("LEN BORDERS:", len(self.borders))
        for border in self.borders:
            # Build an object of the 'BoundaryData' class
            boundary = BoundaryData(
                border=border,
                type_geo=type_geo,
                layout_o=make_vertex(self.layout_centre),
                dimensions=self.dimensions
            )
            # Store all the indices of the edges belonging to the border
            boundary.find_edges_on_border(boundary_edgs_cmpd, self.id_vs_edge)
            # Append the built 'BoundaryData' object to the corresponding list
            self.boundaries.append(boundary)

    def build_edges(
            self, edge_names_vs_faces: Dict[str, List[Any | FaceData]]
        ) -> None:
        """
        Method that builds a list of ``EdgeData`` objects from the given
        dictionary.
        It associates for each edge name a list containing the corresponding
        GEOM edge and the ``FaceData`` objects; the latter represent the faces
        sharing the same edge.

        Parameters
        ----------
        edge_names_vs_faces : Dict[str, List[Any | FaceData]]
            Dictionary of edge names VS the list with the corresponding
            GEOM edge and the connected ``FaceData`` objects.
        """
        # Loop through all the lists of objects associated to each edge
        for shapes in edge_names_vs_faces.values():
            # Instantiate an object of the 'EdgeData' class and append to the
            # corresponding list
            self.edges.append(EdgeData(*shapes))

    def build_edges_and_faces_association(self) -> Dict[str, List[FaceData]]:
        """
        Method that associates the faces sharing the same edge with the edge
        name. These names contain a global index to identify the edge in the
        layout.

        Returns
        -------
        Dict[str, List[FaceData]]
            A dictionary whose entries associate a list of adjacent faces (as
            ``FaceData`` objects) to the name of the corresponding shared
            edge.
        """
        # Initialize the dictionary storing the edges names VS the list of
        # connected faces
        edges_name_vs_faces : Dict[str, List[Any | FaceData]] = {}
        # Loop through all the layout subfaces ('FaceData' objects)
        print("LEN SUBFACES", len(self.subfaces))
        for subface in self.subfaces:
            # Log the 'FaceData' characteristics
            print(subface)
            # Loop through all the edges of the current subface
            for subface_edge, edge_id in subface.edge_vs_id.items():
                # Extract the corresponding GEOM edge object(s)
                unique_edges = self._get_unique_edges(subface_edge, edge_id)
                # Update the dictionary of edge names VS connected faces
                self._update_edge_face_association(
                    edges_name_vs_faces, subface, unique_edges)
        # Return the dictionary of edge names VS the list of connected faces
        return edges_name_vs_faces

    def build_faces(self, property_types: List[PropertyType]) -> None:
        """
        Method that builds a list of ``FaceData`` objects from the GEOM face
        objects, extracted from the layout regions, and the property
        values. Each parsed ``Region`` object corresponds to a region in the
        analysed layoutaccording to the considered type of geometry.
        Each region must be associated with a value for the given property
        type. If not the case, an exception is raised.

        Parameters
        ----------
        property_type : List[PropertyType]
            The list of property types associated to the layout regions.

        Raises
        ------
        RuntimeError
            If no properties or the indicated one are associated to a layout
            region.
        RuntimeError
            If no value for the given property type is associated to a
            layout region.
        """
        # Loop through all the layout regions and build the corresponding
        # data structure storing the face object and the value of the given
        # type of property
        print("LEN REGIONS:", len(self.regions))
        for indx, region in enumerate(self.regions):
            # Build a 'FaceData' object and append to the corresponding list
            try:
                self.subfaces.append(
                    FaceData(region, indx+1, property_types)
                )
            except RuntimeError as e:
                raise RuntimeError("The layout analysis failed.") from e

    def print_log_analysis(
            self, edge_name_vs_faces: Dict[str, List[Any]]) -> None:
        """
        Method that prints on the stdout the log of the data extraction from
        the layout.

        Parameters
        ----------
        edge_name_vs_faces : Dict[str, List[Any]]
            A dictionary of the name of the layout edges VS the GEOM faces
            connected to it.
        """
        # Displays the number of faces, the number of edges and the edges
        # associated with one face and two faces.
        n0 = 0
        n1 = 0
        n2 = 0
        for f in edge_name_vs_faces.values():
            if len(f) == 2:
                n1 += 1
            elif len(f) == 3:
                n2 += 1
            else:
                n0 += 1
        print("\t# subparts (faces) : ", len(self.subfaces))
        print("\t# edges found : ", len(self.edges))
        print("\t# edges w/one face :", n1)
        print("\t# edges w/two faces :", n2)
        print("\t# edges w errors :", n0)

    def _apply_layout_elements_translation(
            self,
            layout_cmpd: Any,
            cur_centre: Any,
            new_centre: Tuple[float, float, float]
        ) -> Any:
        """
        Method that translates the layout regions and the given compound so
        that they are positioned according to the provided new centre.

        If the current layout centre differs from the provided one, all
        the regions of the layout are translated so to keep their relative
        distance from the new position of the layout centre.
        The same translation is applied also to the given layout compound.
        In addition, a vertex object with the given XYZ coordinates is
        assigned to the centre of the stored ``Lattice`` instance.

        Parameters
        ----------
        layout_cmpd : Any
            The layout compound object to be translated.
        cur_centre : Any
            The vertex object representing the current layout centre.
        new_centre : Tuple[float, float, float]
            The coordinates for the new layout centre.

        Returns
        -------
        Any
            The translated layout compound, or the same compound if the
            centre has not changed.
        """
        # Procede only if the layout centre has changed
        if all(math.isclose(c, nc) for c, nc in zip(
            get_point_coordinates(cur_centre), new_centre)):
            return layout_cmpd
        # Translate the regions of the layout
        for region in self.regions:
            region.update(
                wrap_shape(
                    translate_wrt_reference(region, cur_centre, new_centre)
                )
            )
        # Translate the given layout compound and return it
        return translate_wrt_reference(
            layout_cmpd, cur_centre, new_centre)

    def _evaluate_layout_centre(self) -> Tuple[float, float, float]:
        """
        Method that evaluates and returns the coordinates of the centre of the
        layout so that its lower-left corner is positioned in the XYZ space
        origin.

        The method constructs a face from the layout's borders, extracts
        its vertices, and determines the lower-left corner vertex based on
        the minimum X, Y, and Z coordinates.
        If the lower-left corner does not coincide with the XYZ space origin,
        it computes and returns the coordinates the centre should have so
        that the lower-left corner is placed in the origin.
        Otherwise, it returns the current layout centre coordinates.

        Returns
        -------
        Tuple[float, float, float]
            A tuple providing the coordinates of the layout centre so that
            the lower-left corner coincides with the XYZ space origin.
        """
        # Get the layout face vertices
        layout_vertices = extract_sub_shapes(
            make_face(self.borders), ShapeType.VERTEX
        )
        # Get the lower-left corner vertex as the one having the minimum value
        # for the XYZ coordinates
        coords = [get_point_coordinates(p) for p in layout_vertices]
        lower_left = min(coords, key=lambda c: (c[0], c[1], c[2]))
        # Check if the lower-left corner coincides with the XYZ origin; if
        # not, evaluate the new layout centre to fulfill the condition
        if any(not math.isclose(c, 0.0, abs_tol=1e-6) for c in lower_left):
            # Return the coordinates of the new centre
            return (0.0 - lower_left[0]), (0.0 - lower_left[1]), 0.0
        # Return the current layout centre
        return get_point_coordinates(self.geometry_layout.o)

    def _build_compound_regions(
            self, compound: Any, layout_regions: List[Region]
        ) -> None:
        """
        Method that builds a ``Region`` object for each face of the given
        compound. Properties are assigned by identifying the corresponding
        ``Region`` object from the given list.

        Parameters
        ----------
        compound : Any
            The compound object for whose faces ``Region`` objects are built.
        layout_regions : List[Region]
            The list of ``Region`` objects to use as reference.

        Raises
        ------
        RuntimeError
            If any of the face objects of the compound does not have a
            corresponding ``Region`` object among the given ones.
        """
        # Clear out any previously built regions
        self.regions.clear()
        logging.info(
            f"Got no. {len(layout_regions)} regions from the technological "
            "geometry."
        )
        # Build the 'Region' objects corresponding to the faces of the given
        # compound
        try:
            self.regions = build_compound_regions(compound, layout_regions)
        except RuntimeError as e:
            raise RuntimeError(
                f"No region could be found along the hierarchical structure "
                f"of the layout named '{self.geometry_layout.name}' that "
                "matches a face object extracted from the given compound "
                "object. Please ensure the compound is a portion of the "
                "indicated layout."
            ) from e

    def _build_refined_regions(self, refined_cmpd: Any) -> None:
        """
        Method that builds the ``Region`` objects to store from the partition
        of the technological geometry of the layout with the compound object
        collecting the edges that contribute to defining the refined geometry
        layout.
        Given the ``Region`` objects of the technological geometry, the face
        objects in which each region is subdivided are retrieved and
        corresponding ``Region`` objects sharing the same properties are
        built and stored in the instance attribute.

        Parameters
        ----------
        refined_cmpd : Any
            The compound object collecting the edges for the refined geometry
            layout.

        Notes
        -----
        The face objects in which each region is partitioned are retrieved by
        relying on the "in place" concept. This provides:
        - a FACE, if the region was not partitioned;
        - a COMPOUND/SHELL of faces if it has been split.
        """
        # Perform a non-intersecting partition between the regions and the
        # edges of the refined geometry layout; in case of an exception,
        # fallback to the standard partition operation
        try:
            part = make_partition_non_self_intersecting(
                [r.geom_obj for r in self.regions],
                [refined_cmpd],
                ShapeType.FACE
            )
        except RuntimeError:
            part = make_partition(
                [r.geom_obj for r in self.regions],
                [refined_cmpd],
                ShapeType.FACE
            )
        # For each region, recover its image in the partition, i.e. the faces
        # in which it has been partitioned
        refined_regions = []
        idx = 0
        for r in self.regions:
            in_place = get_in_place_by_hystory(part, r.geom_obj)
            # Extract faces from that in-place shape
            sub_faces = extract_sub_shapes(in_place, ShapeType.FACE)
            # If no subfaces are present, keep the in_place as it is a face
            if not sub_faces:
                sub_faces = [in_place]
            # Build a new Region object with the same properties for each face
            for f in sub_faces:
                idx += 1
                refined_regions.append(
                    Region(f, f"Region {idx}", r.properties)
                )
        # Update the stored regions with the refined ones
        self.regions = refined_regions

    def _get_layout_compound(
            self, tdt_setup: TdtSetup, compound_to_analyse: Any | None
        ) -> Any:
        """
        Method that returns the compound object based on the provided setup
        and compound to analyse, if any is indicated.
        This method either uses an existing compound and builds the
        corresponding ``Region`` objects, or extracts regions from the
        technological geometry of the layout considering the symmetry type
        indicated in the ``TdtSetup`` instance.
        In this case, the compound built from the extracted regions is
        returned.

        Parameters
        ----------
        tdt_setup : TdtSetup
            The ``TdtSetup`` instance containing symmetry type information.
        compound_to_analyse : Any | None
            A compound object representing a portion of the entire layout.
            If provided, the corresponding ``Region`` objects are built. If
            ``None``, regions will be extracted from the geometry layout.

        Returns
        -------
        Any
            A compound object representing either the input compound (if
            provided) or a newly created compound built from the regions
            extracted from the layout.

        Raises
        ------
        RuntimeError
            In case a face contained in the given compound does not correspond
            to any of the ``Region`` objects of the layout.
        RuntimeError
            In case the indicated symmetry type does not correspond to any
            symmetry shape in the layout.
        """
        if compound_to_analyse is not None:
            # Build the 'Region' objects corresponding to the given compound
            # by considering the entire layout
            self._build_compound_regions(
                compound_to_analyse, self.geometry_layout.get_regions()
            )
            layout_cmpd = compound_to_analyse
        else:
            # Get the regions of the technological geometry of the layout
            # according to indicated symmetry type, limiting the tolerances
            # of the regions to avoid issues on the result of the common
            # operations with the shape of the symmetry
            self.regions = self.geometry_layout.get_regions_with_symmetry(
                tdt_setup.symmetry_type, True
            )
            # Get the GEOM compound identifying either the full layout of a
            # part of it
            layout_cmpd = make_compound(self.regions)
            logging.info("Built a compound from the extracted regions.")
        return layout_cmpd

    def _get_refinement_edges(self, tdt_setup: TdtSetup) -> Compound:
        """
        Method that returns the ``Compound`` object of the edges of the
        refined geometry layout that matches the geometry type indicated
        in the given ``TdtSetup`` instance.
        An empty compound is returned if the ``GeometryType.TECHNOLOGICAL``
        value is present. The compound object is taken from the corresponding
        mapping in the stored ``Fillable`` instance.

        Parameters
        ----------
        tdt_setup : TdtSetup
            Dataclass providing the settings for extracting the geometric and
            properties information from the given layout.

        Returns
        -------
        Compound
            The ``Compound`` object containing the edges of the refined
            geometry layout. An empty ``Compound`` object, if the
            ``GeometryType.TECHNOLOGICAL`` is indicated.

        Raises
        ------
        RuntimeError
            If no ``Compound`` object is associated to the geometry type
            indicated in the given ``TdtSetup`` instance.
        """
        if tdt_setup.geom_type == GeometryType.TECHNOLOGICAL:
            return wrap_shape(make_compound([]))
        if tdt_setup.geom_type not in self.geometry_layout.geometry_maps:
            raise RuntimeError(
                f"Missing '{tdt_setup.geom_type.name}' type of geometry "
                f"for the layout named '{self.geometry_layout.name}'. "
                "Call the method 'show()' first with the desired type of "
                "geometry to enable the construction of the corresponding "
                "edges."
            )
        return self.geometry_layout.geometry_maps[tdt_setup.geom_type]

    def _get_unique_edges(
            self, subface_edge: Any, edge_id: str
        ) -> List[Any]:
        """
        Method that retrieves the GEOM edge objects associated to the given
        ID (second argument) from the attribute dictionary of IDs VS GEOM
        edges.
        In case of an exception, due to a missing entry in the ``id_vs_edge``
        dictionary, a further analysis is performed. This could happen for
        edges shared by two adjacent faces: the same edge could have a
        different orientation in the two faces, i.e. starting and ending
        points are inverted, or an edge having a single face from one side
        can be associated to different faces from the other.
        In both cases, it is important to retrieve all the corresponding sub
        edges by exploiting the GEOM function ``GetInPlace``; this is used to
        extract the sub-shape(s) of the layout unique edges, which are
        coincident with, or could be a part of, the GEOM edge provided as
        first argument.
        If any edge is retrieved, the corresponding ID is built and used to
        get the corresponding GEOM edges stored in the ``id_vs_edge``
        attribute.

        Parameters
        ----------
        subface_edge : Any
            The GEOM edge object to retrieve the corresponding sub-edges from.
        edge_id : str
            The ID of the edge to look for in the dictionary of IDs VS edges.

        Raises
        ------
        RuntimeError
            If no corresponding edge is found in the layout.

        Returns
        -------
        List[Any]
            A list of GEOM edge objects that is either directly associated to
            the given ID, or representing the sub-edges the edge can be
            subdivided into.
        """
        try:
            return [self.id_vs_edge[edge_id]]
        except KeyError as exc:
            # Get the sub-edges the given edge can be subdivided into: these
            # are associated to a different face of the layout
            edge = get_in_place(make_compound(self.layout_edges),
                                subface_edge)
            # Only edge and compound of edges are treated
            edge_type = get_shape_type(edge)
            error_message = "No corresponding edge in the layout could " +\
                f"be retrieved for the subface edge whose data is {edge_id}"
            if not edge or edge_type not in [ShapeType.COMPOUND,
                                             ShapeType.EDGE]:
                raise RuntimeError(error_message) from exc
            # EDGE-type case
            if not edge_type == ShapeType.COMPOUND:
                return [self.id_vs_edge[build_edge_id(edge)]]
            # COMPOUND-type case
            edges_in_place = extract_sub_shapes(edge, ShapeType.EDGE)
            if edges_in_place:
                return [
                    self.id_vs_edge[build_edge_id(e)] for e in edges_in_place]
            # Raise an exception if the compound does not have edges
            raise RuntimeError(error_message) from exc

    def _handle_refinement(
            self, layout_cmpd: Any, tdt_setup: TdtSetup
        ) -> Any:
        """
        Method that refines the given compound object with the edges stored
        in the mapping that corresponds to the type of geometry indicated by
        the given ``TdtSetup`` instance. The refinement is based on a
        partition operation between the layout's compound object and the
        ``Compound`` instance collecting the edges of the considered geometry
        type. In case of a ``GeometryType.TECHNOLOGICAL`` value, the partition
        is performed on the compound itself, an operation that removes any
        duplicate edges.
        In case of a refined geometry type, the stored ``Region`` objects are
        rebuilt from the refined compound by keeping the same properties of
        the current regions.

        Parameters
        ---------
        layout_cmpd : Any
            The layout compound to be refined, if needed.
        tdt_setup : TdtSetup
            Dataclass providing the settings for extracting the geometric and
            properties information from the given layout.

        Returns
        -------
        Any
            The refined input compound or the same one, if no refinement is
            needed.

        Raises
        ------
        RuntimeError
            If no ``Compound`` object is associated to the geometry type
            indicated in the given ``TdtSetup`` instance.
        """
        # Get the compound of edges related to the indicated geometry type,
        # if any
        edges = self._get_refinement_edges(tdt_setup)
        logging.info(
            f"Got the edges of the {tdt_setup.geom_type} mapping of the "
            f"layout {self.geometry_layout.name}."
        )
        # Update the compound by partitioning the content of the compound,
        # possibly with the refinement edges, if any
        layout_cmpd = make_partition([layout_cmpd], [edges], ShapeType.FACE)
        logging.info(
            f"Updated the compound of the layout {self.geometry_layout.name} "
            f"with the edges of the {tdt_setup.geom_type} mapping."
        )
        # Update the regions, if a refined geometry type is adopted
        if tdt_setup.geom_type != GeometryType.TECHNOLOGICAL:
            self._build_refined_regions(edges)
            logging.info(
                f"Re-extracted no. {len(self.regions)} regions from the "
                f"compound of the layout {self.geometry_layout.name} that "
                f"matches the {tdt_setup.geom_type} mapping."
            )
        return layout_cmpd

    def _handle_translation(
            self, layout_cmpd: Any, tdt_setup: TdtSetup
        ) -> Any:
        """
        Method that handles the translation of the layout compound and its
        regions so that the lower-left corner is in the XYZ space origin.
        This translation is performed only for the combination of specific
        symmetries and layout types or if the layout centre does not coincide
        with the XYZ origin.
        If a geometry type other than ``GeometryType.TECHNOLOGICAL`` is
        indicated, the corresponding compound object is translated as well.
        Lastly, the edges constituting the borders of the translated layout
        compound are built.

        Parameters
        ---------
        layout_cmpd : Any
            The layout compound to be translated, if needed.
        tdt_setup : TdtSetup
            Dataclass providing the settings for extracting the geometric and
            properties information from the given layout.

        Returns
        -------
        Any
            The translated input compound or the same one, if no translation
            is needed.
        """
        # Build and evaluate the two conditions
        symm_condition = (
            tdt_setup.symmetry_type
            in self.CASES_FOR_TRANSLATION[tdt_setup.layout_type]
        )
        centre_condition = (
            not are_same_shapes(
                self.geometry_layout.o,
                make_vertex((0.0, 0.0, 0.0)),
                ShapeType.VERTEX
            )
        )
        # Return if both conditions are not met
        if not symm_condition and not centre_condition:
            return layout_cmpd
        # Evaluate the new centre of the layout, if it has not been
        # translated yet, and apply the translation to the regions
        # and the layout compound
        self.layout_centre = (
            self._evaluate_layout_centre() if symm_condition
            else (0.0, 0.0, 0.0)
        )
        logging.info(f"Evaluated new centre '{self.layout_centre}'.")
        layout_cmpd = self._apply_layout_elements_translation(
            layout_cmpd, self.geometry_layout.o, self.layout_centre
        )
        # Translate the compound of edges corresponding to the indicated
        # geometry type
        if tdt_setup.geom_type != GeometryType.TECHNOLOGICAL:
            self.geometry_layout.geometry_maps[tdt_setup.geom_type] = \
                translate_wrt_reference(
                    self._get_refinement_edges(tdt_setup),
                    self.geometry_layout.o,
                    self.layout_centre
                )
        logging.info(
            "Translated layout elements so to have the lower-left corner "
            "of the layout in the XYZ origin."
        )
        # Re-evaluate the layout borders
        self.borders = build_compound_borders(layout_cmpd)
        logging.info(
            f"Re-extracted no. {len(self.borders)} borders from the "
            f"compound of the layout {self.geometry_layout.name}"
        )
        # Return the translated compound
        return layout_cmpd

    def _preprocess(
            self,
            tdt_setup: TdtSetup,
            compound_to_analyse: Any | None
        ) -> None:
        """
        Method that initialises the geometric information to be stored by
        considering either the ``Fillable`` or the compound objects provided
        when the class ``LayoutDataExtractor`` is instantiated.
        In case no compound is provided, this information is determined
        according to the symmetry and the geometry types of the layout found
        in the provided ``TdtSetup`` instance.

        This method updates the hierarchical structure of the layout, if it
        is needed, and builds the ``Region`` objects of the layout that
        correspond to the faces the layout compound can be subdivided into.
        The layout compound is either provided as input or determined from the
        ``Fillable`` object according to the indicated symmetry type.
        If the layout falls in one of the cases identified by the attribute
        ``CASES_FOR_TRANSLATION``, the determined compound object, and the
        associated regions, is translated so that the lower-left corner
        coincides with the origin of the XYZ space.

        In addition, if the geometry type differs from the
        ``GeometryType.TECHNOLOGICAL`` one, the layout compound is partitioned
        with the corresponding refinement edges; the same is done for the
        ``Region`` objects.
        Lastly, the edges of the considered layout are extracted and stored
        as an instance attribute.

        Parameters
        ----------
        tdt_setup : TdtSetup
            Dataclass providing the settings for extracting the geometric and
            properties information from the given layout.
        compound_to_analyse: Any | None = None
            The compound object to process, if present. If ``None`` is given,
            the stored ``Fillable`` object is considered instead.

        Raises
        ------
        RuntimeError
            If no ``Compound`` object is associated to the geometry type
            indicated in the given ``TdtSetup`` instance.
        RuntimeError
            During the extraction of the ``Region`` objects of the layout's
            technological geometry, if there is no match between a face of the
            given compound and the regions of the ``Fillable`` object or the
            indicated symmetry type is not available.
        RuntimeError
            If no closed boundaries can be extracted from the compound object
            to analyse.
        """
        # Update the hierarchical tree of the stored layout, if needed
        self.geometry_layout.update_hierarchical_structure()
        # Get the compound to export based on either the regions or the given
        # one
        layout_cmpd = self._get_layout_compound(
            tdt_setup, compound_to_analyse
        )
        logging.info(
            f"Extracted no. {len(self.regions)} regions of the technological "
            f"geometry from the layout {self.geometry_layout.name}."
        )
        # Extract the borders of the layout compound
        self.borders = build_compound_borders(layout_cmpd)
        logging.info(
            f"Extracted no. {len(self.borders)} borders from the compound "
            f"of the layout {self.geometry_layout.name}"
        )
        # Extract the layout characteristic dimensions depending on the layout
        # type
        x_min, x_max, y_min, y_max = get_bounding_box(layout_cmpd)
        self.dimensions = (x_max - x_min, y_max - y_min)
        if tdt_setup.layout_type == LayoutType.HEX:
            self.dimensions = tuple([d / 2 for d in self.dimensions])

        # Handle the layout translation so that the lower-left corner is in
        # the XYZ space origin; this is valid for specific symmetries and
        # layout types or if the layout centre does not coincide with the XYZ
        # origin
        layout_cmpd = self._handle_translation(layout_cmpd, tdt_setup)

        # Rebuild the regions and the layout compound in case the indicated
        # geometry type is not the technological one
        layout_cmpd = self._handle_refinement(layout_cmpd, tdt_setup)

        # Extract the layout edges from the layout compound to analyse
        self.layout_edges = extract_sub_shapes(layout_cmpd, ShapeType.EDGE)
        logging.info(
            f"Extracted no. {len(self.layout_edges)} edges from the "
            f"compound of the layout {self.geometry_layout.name}"
        )

        add_to_study(layout_cmpd, "CMPD TO ANALYSE")
        add_to_study(make_compound(self.layout_edges), "EDGES TO ANALYSE")
        add_to_study(make_compound(self.borders), "BOUNDARIES TO ANALYSE")
        update_salome_study()

    def _update_edge_face_association(
            self,
            edges_name_vs_faces: Dict[str, List[Any | FaceData]],
            subface: FaceData,
            edges: List[Any]
        ) -> None:
        """
        Method that updates the given dictionary of edge names VS connected
        faces with the provided ``FaceData`` object and the list of edges.
        For each edge, its ID is checked for its presence among the keys of
        the given dictionary: if so, the corresponding list is updated with
        the ``FaceData`` object, otherwise a new entry is created. The entry
        has as key the edge name, as value a list with the GEOM edge object
        as first element, followed by the given ``FaceData`` object.

        Parameters
        ----------
        edges_name_vs_faces : Dict[str, List[Any | FaceData]]
            A dictionary of edge names VS connected faces.
        subface : FaceData
            A ``FaceData`` object to associate the given edges to.
        edges : List[Any]
            A list of GEOM edge objects each associated to the given face.
        """
        for edge in edges:
            # Get the edge name
            edge_name = get_shape_name(edge)
            # Check if the edge ID is already stored, if so, append
            # the current face, otherwise add a new entry
            if edge_name in edges_name_vs_faces:
                edges_name_vs_faces[edge_name].append(subface)
            else:
                edges_name_vs_faces[edge_name] = [edge, subface]


def analyse_layout(
        layout: Fillable,
        tdt_setup: TdtSetup,
        compound_to_analyse: Any | None = None
    ) -> LayoutDataExtractor:
    """
    Function that analyses the given layout to extract the necessary data
    about the regions and their associated properties, the edges and their
    association with the faces they belong to. It also extracts information
    about the edges representing the boundaries of the layout.

    If the ``compound_to_analyse`` parameter is provided, it will be the one
    to be analysed, according to the information stored in the provided
    layout. In this way, a generic portion of the given ``Fillable`` object
    can be treated.

    Parameters
    ----------
    layout : Fillable
        The ``Fillable`` object storing the information about the geometry
        and the properties of the entire layout.
    tdt_setup : TdtSetup
        Dataclass providing the settings for extracting the geometric and
        properties information from the given layout.
    compound_to_analyse: Any | None = None
        The compound object to analyse, if present. If ``None`` is given, the
        given ``Fillable`` object is considered instead.

    Returns
    -------
    LayoutDataExtractor
        Object collecting all the information about the geometry and the
        properties extracted from the layout.
    """
    # Instantiate the class for extracting the geometric data from the layout
    # according to the given type of geometry
    data_extractor = LayoutDataExtractor(
        layout,
        tdt_setup,
        compound_to_analyse
    )
    # Call its method for performing the analysis
    data_extractor.build_faces(tdt_setup.property_types)
    edge_name_vs_faces = data_extractor.build_edges_and_faces_association()
    data_extractor.build_edges(edge_name_vs_faces)
    data_extractor.build_boundaries(tdt_setup.type_geo)
    data_extractor.print_log_analysis(edge_name_vs_faces)

    # Return the instance
    return data_extractor
