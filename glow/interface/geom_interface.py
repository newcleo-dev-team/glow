"""
Module containing functions providing an interface towards the GEOM functions
of SALOME.
"""
from enum import Enum
from typing import Any, Dict, List, Tuple, Union

# Import SALOME element
from glow.interface import geompy, gst, salome, SALOMEDS


class ShapeType(Enum):
    """
    Enumeration for defining the GEOM topological types of shapes.
    """
    COMPOUND : int = geompy.ShapeType['COMPOUND']
    """Indicating a compound shape type."""
    COMPSOLID : int = geompy.ShapeType['COMPSOLID']
    """Indicating a compound solid shape type."""
    SOLID : int = geompy.ShapeType['SOLID']
    """Indicating a solid shape type."""
    SHELL : int = geompy.ShapeType['SHELL']
    """Indicating a shell shape type."""
    FACE : int = geompy.ShapeType['FACE']
    """Indicating a face shape type."""
    WIRE : int = geompy.ShapeType['WIRE']
    """Indicating a wire shape type"""
    EDGE : int = geompy.ShapeType['EDGE']
    """Indicating an edge shape type."""
    VERTEX : int = geompy.ShapeType['VERTEX']
    """Indicating a vertex shape type."""
    SHAPE : int = geompy.ShapeType['SHAPE']
    """Indicating a generic shape type."""
    FLAT : int = geompy.ShapeType['FLAT']
    """Indicating a flat shape type."""


# Dictionary associating the name of the GEOM type of shape VS the
# corresponding enumeration element of 'ShapeType'
NAME_VS_SHAPE_TYPE: Dict[str, ShapeType] = {
    'COMPOUND': ShapeType.COMPOUND,
    'COMPSOLID': ShapeType.COMPSOLID,
    'SOLID': ShapeType.SOLID,
    'SHELL': ShapeType.SHELL,
    'FACE': ShapeType.FACE,
    'WIRE': ShapeType.WIRE,
    'EDGE': ShapeType.EDGE,
    'VERTEX': ShapeType.VERTEX,
    'SHAPE': ShapeType.SHAPE,
    'FLAT': ShapeType.FLAT,
    'DISK': ShapeType.FACE,
    'PLANAR': ShapeType.FACE,
    'POLYGON': ShapeType.FACE,
    'ARC_CIRCLE': ShapeType.EDGE,
    'CIRCLE': ShapeType.EDGE,
    'SEGMENT': ShapeType.EDGE,
    'DISK_CIRCLE': ShapeType.FACE
}


def add_to_study(shape: Any, name: str) -> str:
    """
    Function that adds the given shape in the current SALOME study and
    returns its ID entry assigned by SALOME.

    Parameters
    ----------
    shape : Any
        The shape to be added to the study.
    name : str
        The name of the shape to be visualized in the SALOME object browser.

    Returns
    -------
    str
        A string representing the ID entry assigned by SALOME to the shape
        when added to the current study.
    """
    return geompy.addToStudy(shape, name)


def add_to_study_in_father(father_shape: Any, shape: Any, name: str) -> str:
    """
    Function that adds the given shape in the current SALOME study and
    returns its ID entry assigned by SALOME.

    Parameters
    ----------
    father_shape : Any
        The father shape in the object browser of the study.
    shape : Any
        The shape to be added to the study.
    name : str
        The name of the shape to be visualized in the SALOME object browser.

    Returns
    -------
    str
        A string representing the ID entry assigned by SALOME to the shape
        when added to the current study under the given father shape.
    """
    return geompy.addToStudyInFather(father_shape, shape, name)


def clear_view() -> None:
    """
    Function that clears out any geometrical shape currently displayed in
    the SALOME 3D viewer.
    """
    if is_gui_available():
        salome.sg.EraseAll()


def display_shape(entry_id: str) -> None:
    """
    Function that displays the geometrical shape, whose entry ID is provided
    as input, in the SALOME 3D viewer.

    Parameters
    ----------
    id : str
        The entry ID of the shape to display in the SALOME 3D viewer.
    """
    if is_gui_available():
        salome.sg.Display(entry_id)


def extract_sorted_sub_shapes(shape: Any,
                              sub_shapes_type: ShapeType) -> List[Any]:
    """
    Function that explodes a shape on its sub-shapes of the given type.
    Sub-shapes are sorted by taking into account their gravity centers.

    Parameters
    ----------
    shape : Any
        The shape object to explode.
    sub_shapes_type : ShapeType
        The type of the sub-shapes to extract as value of the ``ShapeType``
        enumeration.

    Returns
    -------
    List[Any]
        A list of the sub-shapes of the given type that are contained in the
        provided shape.
    """
    return geompy.SubShapeAllSortedCentres(shape, sub_shapes_type.value)


def extract_sub_shapes(shape: Any, sub_shapes_type: ShapeType) -> List[Any]:
    """
    Function that explodes a shape on its sub-shapes of the given type
    without applying any sorting algorithm.

    Parameters
    ----------
    shape : Any
        The shape object to explode.
    sub_shapes_type : ShapeType
        The type of the sub-shapes to extract as item of the ``ShapeType``
        enumeration.

    Returns
    -------
    List[Any]
        A list of the sub-shapes of the given type that are contained in the
        provided shape.
    """
    return geompy.ExtractShapes(shape, sub_shapes_type.value)


def fuse_edges_in_wire(wire: Any) -> Any:
    """
    Function that modifies the given wire object by suppressing all the
    vertices in the C1 continuous adjacent edges of the wire.

    Parameters
    ----------
    wire : Any
        The wire object whose vertices to suppress.

    Returns
    -------
    Any
        A modified wire obtained by suppressing vertices in C1 continuous
        adjacent edges in the wire.
    """
    return geompy.FuseCollinearEdgesWithinWire(wire, [])


def get_angle_between_shapes(shape1: Any, shape2: Any) -> float:
    """
    Function that computes the angle in degrees between two EDGE-type shapes,
    which must both be linear edges, i.e. their type name in SALOME is
    ``SEGMENT``.

    Parameters
    ----------
    shape1 : Any
        The first shape object, expected to be a linear edge.
    shape2 : Any
        The second shape object, expected to be a linear edge.

    Returns
    -------
    float
        The angle in degrees between the two segment-type shapes.

    Raises
    ------
    RuntimeError
        If either ``shape1`` or ``shape2`` is not of type ``SEGMENT``.
    """
    # Get the specific type name of the two shapes
    shape1_type = str(get_kind_of_shape(shape1)[0])
    shape2_type = str(get_kind_of_shape(shape2)[0])
    # Check if the given shapes are of the expected 'SEGMENT' type
    if shape1_type != "SEGMENT" or shape2_type != "SEGMENT":
        raise RuntimeError(
            "Both shapes must be objects of type 'SEGMENT': however, shape1 "
            f"is of type '{shape1_type}' and shape2 is of type "
            f"'{shape2_type}'")
    # Return the angle in degrees between the two EDGE-type shapes
    return geompy.GetAngle(shape1, shape2)


def get_basic_properties(shape: Any) -> List[float]:
    """
    Function that returns the sum of the lengths of all the wires, the
    area and volume of the given shape.

    Parameters
    ----------
    shape : Any
        The shape whose geometrical properties to extract.

    Returns
    -------
    List[float]
        A list providing the sum of the lengths of all the wires, the area
        and volume of the given shape.
    """
    return geompy.BasicProperties(shape)


def get_bounding_box(shape: Any) -> List[float]:
    """
    Function that returns the bounding box extension of the given shape.

    Parameters
    ----------
    shape : Any
        The shape whose bounding box to extract.

    Returns
    -------
    List[float]
        The `[Xmin, Xmax, Ymin, Ymax]` values representing the shape bounding
        box extension.
    """
    return geompy.BoundingBox(shape)[:4]


def get_closed_free_boundary(compound: Any) -> List[Any]:
    """
    Function that evaluates the closed free boundaries of the given compound
    object. An exception is raised if any error occurred during the operation.

    Parameters
    ----------
    compound : Any
        The compound object to extract the closed boundary from.

    Returns
    -------
    List[Any]
        A list of the wire objects each one forming a closed boundary within
        the given compound.

    Raises
    ------
    RuntimeError
        If no closed boundaries could be extracted from the given compound
        object.
    """
    (isDone, closed, _) = geompy.GetFreeBoundary(compound)
    if not isDone or not closed:
        raise RuntimeError("No closed free boundaries could be extracted "
                           "from the given compound object")
    return closed


def get_id_from_object(shape: Any) -> str:
    """
    Function that retrieves the unique SALOME ID associated with a given
    shape object, if shown in the current SALOME study.

    Parameters
    ----------
    shape : Any
        The shape object for which to retrieve the SALOME ID.

    Returns
    -------
    str
        The unique SALOME ID corresponding to the provided shape.

    Raises
    ------
    RuntimeError
        If the shape is not present in the current SALOME study.
    """
    id = salome.ObjectToID(shape)
    if id is None:
        raise RuntimeError(
            f"The shape {get_shape_name(shape)} is not present in the "
            "current SALOME study.")
    return id


def get_in_place(shape1: Any, shape2: Any) -> Any:
    """
    Function that extracts the sub-shape(s) of first shape, which are
    coincident with, or could be a part of, the second shape.

    Parameters
    ----------
    shape1  : Any
        The shape to find sub-shapes of.
    shape2  : Any
        The shape specifying what to find in the first one.

    Returns
    -------
    Any
        A compound object including all the found sub-shapes in the first
        given shape.
    """
    return geompy.GetInPlace(shape1, shape2)


def get_inertia_matrix(shape: Any) -> List[float]:
    """
    Function that returns the inertia matrix of the given shape.

    Parameters
    ----------
    shape : Any
        The shape whose inertia matrix to extract.

    Returns
    -------
    List[float]
        A list providing, for the given shape, the components of its inertia
        matrix and the moments of inertia in the XYZ directions.
    """
    return geompy.Inertia(shape)


def get_kind_of_shape(shape: Any) -> List[Any]:
    """
    Function that returns the list of geometrical information about the
    given shape which depends on the kind of shape.
    Shapes of interest returns the following information:

      - `CIRCLE xc yc zc dx dy dz R` (X-Y-Z center coordinates, X-Y-Z normal
        vector elements, circle radius).
      - `ARC_CIRCLE xc yc zc dx dy dz R x1 y1 z1 x2 y2 z2`
        (X-Y-Z center coordinates, X-Y-Z normal vector elements, arc
        radius, X-Y-Z coordinates of arc starting and ending points).
      - `SEGMENT x1 y1 z1 x2 y2 z2`
        (X-Y-Z coordinates of segment starting and ending points).

    Notes
    -----
    Values lesser than the tolerance of 1e-10 are substituted with ``0`` in
    the returned list.

    Parameters
    ----------
    shape : Any
        The shape whose geometrical information to extract.

    Returns
    -------
    List[Any]
        A list providing, for the given shape, its geometrical information.
    """
    # Extract the list of information about the given shape
    kind_of_shape = geompy.KindOfShape(shape)
    # Loop over all the retrieved values and substitute with '0' values
    # lesser than the 1e-15 tolerance
    for i, info in enumerate(kind_of_shape):
        if isinstance(info, float):
            if abs(info) < 1e-10:
                kind_of_shape[i] = 0
    return kind_of_shape


def get_min_distance(shape1: Any, shape2: Any) -> float:
    """
    Function that returns the minimum distance between the given shapes.

    Parameters
    ----------
    shape1 : Any
        The first shape to get the minimum distance.
    shape2 : Any
        The second shape to get the minimum distance.

    Returns
    -------
    float
        The minimum distance between the two given shapes.
    """
    return geompy.MinDistance(shape1, shape2)


def get_object_from_id(entry_id: str) -> Any | None:
    """
    Function that returns the geometrical object associated to the entry ID
    declared in the current SALOME study.

    Parameters
    ----------
    entry_id : str
        The value of the entry ID associated to the geometrical object in the
        current SALOME study.

    Returns
    -------
    Any | None
        The geometrical object associated with the given entry ID in the
        current SALOME study or ``None``, if no object could be found.
    """
    try:
        return gst.getGeomObjectFromEntry(entry_id)
    except:
        return None


def get_point_coordinates(point: Any) -> Tuple[float, float, float]:
    """
    Function that returns the XYZ coordinates of the given vertex object.

    Parameters
    ----------
    point : Any
        The vertex object whose coordinates are returned.

    Returns
    -------
    Tuple[float, float, float]
        The XYZ coordinates of the given vertex object.
    """
    return geompy.PointCoordinates(point)


def get_selected_object() -> Any:
    """
    Function that returns the geometrical object being currently selected
    in the study.

    Returns
    -------
    Any
        The selected geometrical object in the study.
    """
    return gst.getGeomObjectSelected()


def get_shape_name(shape: Any) -> str:
    """
    Function that returns the ``name`` attribute assigned to the given shape.

    Parameters
    ----------
    shape : Any
        The geometrical shape whose name to retrieve.

    Returns
    -------
    str
        The value of the ``name`` attribute assigned to the given shape.
    """
    return shape.GetName()


def get_shape_type(shape: Any) -> ShapeType:
    """
    Function that returns the type of the given shape as value of the
    ``ShapeType`` enumeration.

    Parameters
    ----------
    shape : Any
        The shape whose type to determine.

    Returns
    -------
    ShapeType
        The shape type as value of the ``ShapeType`` enumeration.
    """
    # Get the name of the shape type
    type = str(get_kind_of_shape(shape)[0])
    # Return the corresponding value of the 'ShapeType' enumeration
    try:
        return NAME_VS_SHAPE_TYPE[type]
    except KeyError:
        print(f"WARNING, type {type} not recognized")


def get_subshape_id(shape: Any, subshape: Any) -> str:
    """
    Function that returns the entry ID of a given subshape contained in the
    parent shape.

    Parameters
    ----------
    shape : Any
        The parent shape.
    subshape : Any
        The subshape contained in the parent shape.

    Returns
    -------
    str
        A string representing the entry ID associated to the subshape in the
        current study.
    """
    return geompy.GetSubShapeID(shape, subshape)


def is_point_inside_shape(point: Any, shape: Any) -> bool:
    """
    Function that checks if the given point object is within the boundaries
    of the given geometrical shape.

    Parameters
    ----------
    point : Any
        The point object to check if inside the shape.
    shape : Any
        The shape object the point position has to be evaluated.

    Returns
    -------
    bool
        ``True``, if the point is inside the shape, ``False`` otherwise.
    """
    return geompy.AreCoordsInside(shape,
                                  list(get_point_coordinates(point)))[0]


def is_gui_available() -> None:
    """
    Function that returns a boolean flag indicating whether the SALOME GUI
    is available when running a script.

    Returns
    -------
    bool
        ``True`` if the SALOME GUI is available, ``False`` otherwise.
    """
    return salome.sg.hasDesktop()


def make_arc_edge(point1: Any, point2: Any, point3: Any) -> Any:
    """
    Function that returns the arc edge object built from the given three
    vertex objects.

    Parameters
    ----------
    point1 : Any
        The vertex object being the arc's center.
    point2 : Any
        The vertex object being the arc's start point.
    point3 : Any
        The vertex object being the arc's end point.

    Returns
    -------
    Any
        The arc edge built from the given three construction points.
    """
    return geompy.MakeArcCenter(point1, point2, point3, False)


def make_cdg(shape: Any) -> Any:
    """
    Function that returns the vertex object of the given shape CDG.

    Parameters
    ----------
    shape : Any
        The geometric shape whose CDG is returned.

    Returns
    -------
    Any
        The shape CDG as a vertex object.
    """
    return geompy.MakeCDG(shape)


def make_circle(center: Any, axis: Union[Any, None], radius: float) -> Any:
    """
    Function that returns a circle object, given its center, axis and radius.

    Parameters
    ----------
    center : Any
        The vertex object being the center of the circle.
    axis : Any
        The vector object being the normal axis of the circle.
    radius : float
        The value of the circle radius.

    Returns
    -------
    Any
        The circle object from the given center, axis and radius.
    """
    return geompy.MakeCircle(center, axis, radius)


def make_common(shape1: Any, shape2: Any) -> Any:
    """
    Function that performs the `common` boolean operation between the two
    given shapes.

    Parameters
    ----------
    shape1 : Any
        The first shape of the `common` operation.
    shape2 : Any
        The second shape of the `common` operation.

    Returns
    -------
    Any
        A shape object resulting from the `common` operation.
    """
    return geompy.MakeCommon(shape1, shape2)


def make_compound(shapes: List[Any]) -> Any:
    """
    Function that creates a compound object from the given list of shapes.

    Parameters
    ----------
    shapes : List[Any]
        The list of shapes to be put into the returned compound object.

    Returns
    -------
    Any
        A compound object made from the given list of shapes.
    """
    return geompy.MakeCompound(shapes)


def make_cut(shape1: Any, shape2: Any) -> Any:
    """
    Function that performs the `cut` boolean operation between the two
    given shapes.

    Parameters
    ----------
    shape1 : Any
        The first shape of the `cut` operation.
    shape2 : Any
        The second shape of the `cut` operation.

    Returns
    -------
    Any
        A shape object resulting from the `cut` operation.
    """
    return geompy.MakeCut(shape1, shape2)


def make_edge(vertex1: Any, vertex2: Any) -> Any:
    """
    Function that returns an edge object, given the vertex objects being
    its start-end points.

    Parameters
    ----------
    vertex1 : Any
        The vertex object being the edge start point.
    vertex2 : Any
        The vertex object being the edge end point.

    Returns
    -------
    Any
        The edge object built from the given start-end points.
    """
    return geompy.MakeEdge(vertex1, vertex2)


def make_face(borders: List[Any]) -> Any:
    """
    Function that returns a 2D face object, given the list of its edge
    objects being the face borders.

    Parameters
    ----------
    borders : List[Any]
        The list of edge objects being the face borders.

    Returns
    -------
    Any
        The face object built on the given borders.
    """
    return geompy.MakeFaceWires(borders, True)


def make_fuse(shapes: List[Any]) -> Any:
    """
    Function that performs a `fuse` boolean operation on the given list
    of shapes.

    Parameters
    ----------
    shapes : List[Any]
        The list of geometrical shapes to be fused into a single one.

    Returns
    -------
    Any
        The face object resulting from fusing all the given shapes.
    """
    return make_face([geompy.MakeFuseList(shapes, True, True)])


def make_line(point1: Any, point2: Any) -> Any:
    """
    Function that returns a line object (i.e. a straight edge), given
    the point objects being its start-end points.

    Parameters
    ----------
    point1 : Any
        The point object being the line start point.
    point2 : Any
        The point object being the line end point.

    Returns
    -------
    Any
        The line object built from the given start-end points.
    """
    return geompy.MakeLineTwoPnt(point1, point2)


def make_partition(
        shapes: List[Any], tools: List[Any], shape_type: ShapeType) -> Any:
    """
    Function that performs a partition operation on the given list of shapes
    by means of the tool shapes intersecting the first ones.
    The result is a shape made by the intersection of all the provided ones
    with type given as input.

    Parameters
    ----------
    shapes : List[Any]
        The list of shapes to be intersected.
    tools : List[Any]
        The list of shapes intersecting.
    shape_type : ShapeType
        The type of the shape resulting from the partition operation.

    Returns
    -------
    Any
        A shape made by the intersection of all the provided ones with the
        type specified as input.
    """
    return geompy.MakePartition(ListShapes=shapes,
                                ListTools=tools,
                                Limit=shape_type.value)


def make_rotation(shape: Any, axis: Any, angle: float) -> Any:
    """
    Function that rotates the given shape by the given angle in radians.

    Parameters
    ----------
    shape : Any
        The generic shape object the rotation should be applied to.
    axis : Any
        The vector object representing the rotation axis.
    angle : float
        The rotation angle in radians.

    Returns
    -------
    Any
        The rotated shape object.
    """
    return geompy.MakeRotation(shape, axis, angle)


def make_translation(shape: Any, vector: Any) -> Any:
    """
    Function that translates the given shape along the vector object.

    Parameters
    ----------
    shape : Any
        The generic shape object the translation should be applied to.
    vector : Any
        The translation vector object.

    Returns
    -------
    Any
        The translated shape object.
    """
    return geompy.MakeTranslationVector(shape, vector)


def make_vector(vect_elem: Tuple[float, float, float]) -> Any:
    """
    Function that returns a vector object, given its XYZ components.

    Parameters
    ----------
    vect_elem : Tuple[float, float, float]
        The XYZ components of the vector to build.

    Returns
    -------
    Any
        The vector object built on the given XYZ components.
    """
    return geompy.MakeVectorDXDYDZ(*vect_elem)


def make_vector_from_points(point1: Any, point2: Any) -> Any:
    """
    Function that returns a vector object, given its start-end points.

    Parameters
    ----------
    point1 : Any
        The start point object of the vector to build.
    point2 : Any
        The end point object of the vector to build.

    Returns
    -------
    Any
        The vector object built on the given XYZ components.
    """
    return geompy.MakeVector(point1, point2)


def make_vertex(coords: Tuple[float, float, float]) -> Any:
    """
    Function that returns a vertex object, given its XYZ coordinates.

    Parameters
    ----------
    coords : Tuple[float, float, float]
        The XYZ coordinates of the vertex to build.

    Returns
    -------
    Any
        The vertex object positioned at the given XYZ coordinates.
    """
    return geompy.MakeVertex(*coords)


def make_vertex_inside_face(face: Any) -> Any:
    """
    Function that returns a vertex object which lays on the given face, in
    an arbitrary position with the only condition of having a non-zero
    distance to the face boundary.

    Parameters
    ----------
    face : Any
        The reference face object to put a vertex into.

    Returns
    -------
    Any
        The vertex object positioned within the given face object.
    """
    return geompy.MakeVertexInsideFace(face)


def make_vertex_on_curve(curve: Any, u_param: float) -> Any:
    """
    Function that returns a vertex object built on the given edge object
    (being a circle, an arc or a segment) and placed at a position identified
    by the adimensional parameter `u_param`, which expresses the scaled
    position wrt the edge length.

    Parameters
    ----------
    curve : Any
        The reference edge object along which the vertex should be placed.
    u_param : float
        Identifying the position of the vertex along the edge, expressed
        in the range [0-1] as it is scaled to the length of the edge.

    Returns
    -------
    Any
        The vertex object positioned along the edge length.
    """
    return geompy.MakeVertexOnCurve(curve, u_param, True)


def make_vertex_on_lines_intersection(line1: Any, line2: Any) -> Any:
    """
    Function that returns a vertex object built on the intersection of
    the two given line objects.

    Parameters
    ----------
    line1 : Any
        The first reference line object.
    line2 : Any
        The second reference line object.

    Returns
    -------
    Any
        The vertex object being the intersection point of the two given lines.
    """
    return geompy.MakeVertexOnLinesIntersection(line1, line2)


def make_wire(edges: List[Any]) -> Any:
    """
    Function that builds a wire object from the given list of connected
    edge objects.

    Parameters
    ----------
    edges : List[Any]
        The list of edge objects to build a wire from.

    Returns
    -------
    Any
        A wire object build from the given edges.
    """
    return geompy.MakeWire(edges)


def remove_from_study(entry_id: str) -> None:
    """
    Function that removes the geometrical object whose ID is provided as
    input.

    Parameters
    ----------
    entry_id : str
        The ID of the geometrical object to be removed from the study.
    """
    gst.removeFromStudy(entry_id)


def set_color_face(face: Any, color: Tuple[int, int, int]) -> None:
    """
    Function that assigns a color to the given face object. When displayed
    in the 3D viewer, the face will be shown with the assigned color.

    Parameters
    ----------
    face : Any
        The face object to which a color is assigned when displayed in the
        3D viewer.
    color : Tuple[int, int, int]
        The RGB code of the color to assign to the face.
    """
    face.SetColor(SALOMEDS.Color(*(rgb/255 for rgb in color)))


def set_shape_name(shape: Any, name: str) -> None:
    """
    Function that assigns a value to the ``name`` attribute for the given
    shape.

    Parameters
    ----------
    shape : Any
        The geometrical shape whose name to set.
    name : str
        The name to assign to the shape.
    """
    shape.SetName(name)


def update_salome_study() -> None:
    """
    Function that updates the current SALOME study by updating the 3D
    viewer and the object browser with the shapes currently added to
    the study.
    """
    # Show everything on the SALOME application
    if is_gui_available():
        # Update the object browser with the objects added to the study
        salome.sg.updateObjBrowser()
        # Set the view to top
        salome.sg.ViewTop()
        # Fit all content in the viewer
        salome.sg.FitAll()
        # Update the view by showing the built geometry objects
        salome.sg.UpdateView()
