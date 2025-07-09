"""
Module containing classes providing the means for creating basic geometries
in SALOME.
"""
import math

from abc import ABC, abstractmethod
from typing import Any, List, Tuple, Union

from glow.geometry_layouts.utility import check_shape_expected_types, \
    update_relative_pos
from glow.interface.geom_interface import ShapeType, add_to_study, \
    add_to_study_in_father, extract_sub_shapes, get_basic_properties, \
    get_bounding_box, get_min_distance, get_point_coordinates, \
    get_shape_name, make_arc_edge, make_cdg, make_circle, make_edge, \
    make_face, make_partition, make_rotation, make_translation, make_vector, \
    make_vector_from_points, make_vertex, make_vertex_on_curve, \
    remove_from_study, update_salome_study


class Surface(ABC):
    """
    Abstract class for representing any geometric surface.

    Parameters
    ----------
    center        : Union[Tuple[float, float, float], None] = None
                    The X-Y-Z coordinates of the generic geometric surface.

    Attributes
    ----------
    o             : Any
                    A vertex object representing the surface center
    borders       : List[Any]
                    A list of edge objects representing the border edges
                    of the surface
    face          : Any
                    A face object representing the geometric surface
    name          : str
                    The name of the surface when displayed in the SALOME study
    face_entry_id : Union[str, None]
                    An ID associated to the surface in the SALOME study
    vertices      : List[Any]
                    A list of vertex objects representing the vertices of the
                    generic surface
    out_circle    : Any
                    An edge object representing the construction circle which
                    the generic surface is inscribed into
    rotation      : float
                    The rotation angle of the geometric surface in radians
    lx            : float
                    The characteristic dimension of the generic geometric
                    surface along the X-axis
    ly            : float
                    The characteristic dimension of the generic geometric
                    surface along the Y-axis
    """
    def __init__(
            self,
            center: Union[Tuple[float, float, float], None] = None) -> None:
        super().__init__()
        # Initialize the instance attributes
        if not center:
            center = (0.0, 0.0, 0.0)
        self.o: Any = make_vertex(center)
        self.borders: List[Any] = []
        self.face: Any | None = None
        self.face_entry_id: Union[str, None] = None
        self.name: str = ""
        self.vertices: List[Any] = []
        self.out_circle: Union[Any, None] = None
        self.rotation: float = 0.0
        self.lx: float = 0.0
        self.ly: float = 0.0

    @abstractmethod
    def _build_borders(self) -> List[Any]:
        """
        Method that builds the border edges of the generic surface this
        class is representing.
        Being abstract, a specific implementation must be provided by
        subclasses.

        Returns
        -------
        A list of edge objects representing the generic surface borders.
        """

    def rotate(self, angle: float) -> None:
        """
        Method that rotates the surface geometry. The face, its border
        edges, vertices and construction circle are rotated by the given
        angle expressed in degrees around the Z-axis.

        Parameters
        ----------
        angle : float
                The rotation angle in degrees
        """
        # Get the figure center coordinates
        center = get_point_coordinates(self.o)
        # Build the Z-axis of rotation positioned in the figure center
        z_axis = make_vector_from_points(
            self.o, make_vertex((center[0], center[1], 1)))
        # Rotate the surface elements
        self.rotate_from_axis(angle, z_axis)

    def rotate_from_axis(self, angle: float, axis: Any) -> None:
        """
        Method that rotates the surface geometry. The face, its border
        edges, vertices and construction circle are rotated by the given
        angle expressed in degrees around the given axis.

        Parameters
        ----------
        angle : float
                The rotation angle in degrees
        axis  : Any
                The axis object around which the rotation takes place
        """
        # Convert the rotation angle in radians
        self.rotation = math.radians(angle)
        # Rotate the geometric elements of the surface
        self.face = make_rotation(self.face, axis, self.rotation)
        for i, _ in enumerate(self.vertices):
            self.vertices[i] = make_rotation(
                self.vertices[i], axis, self.rotation)
        # Re-build the borders
        self.borders = extract_sub_shapes(self.face, ShapeType.EDGE)
        # Rotate the construction circle
        self.out_circle = make_rotation(self.out_circle, axis, self.rotation)

    def translate(self, new_pos: Tuple[float, float, float]) -> None:
        """
        Method that translates the surface geometry to the given position.
        All the geometrical elements (center, face, its borders, vertices
        and construction circle) needed to describe the shape are moved
        accordingly.

        Parameters
        ----------
        new_pos : Tuple[float, float, float]
                  The XYZ coordinates of the new center of the shape
        """
        # Build a vector from the current center to the new one
        transl_vect = make_vector_from_points(self.o, make_vertex(new_pos))

        # Store previous figure center
        pre_center = self.o
        # Translate the characteristic geometrical elements of the shape
        self.o = make_vertex(new_pos)
        self.face = make_translation(self.face, transl_vect)
        self.out_circle = make_translation(self.out_circle, transl_vect)
        # Re-build the vertices in the updated position relative to the
        # moved center
        for i, _ in enumerate(self.vertices):
            self.vertices[i] = make_vertex(
                update_relative_pos(self.vertices[i], pre_center, new_pos))
        # Re-build the borders
        self.borders = self._build_borders()

    def show_borders(self) -> None:
        """
        Method that adds the surface border edges to the current SALOME study.
        These edges are added to the Object Browser under the face they refer
        to only if any face has already been shown.

        Raises
        ------
        RuntimeError: if no face has been displayed yet in the SALOME study.
        """
        if not self.face_entry_id:
            raise RuntimeError("Bordes cannot be shown as no corresponding " \
            "face has been displayed yet.")
        # Show the face border edges as sub-elements of the face in the
        # current SALOME study
        for i, border in enumerate(self.borders):
            add_to_study_in_father(self.face, border, f'Border {i}')
        # Update the SALOME view
        update_salome_study()

    def show_edges_and_vertices(self) -> None:
        """
        Method that adds the edges and the vertices of the geometric surface
        to the current SALOME study.
        """
        self.show_borders()
        for i, v in enumerate(self.vertices):
            add_to_study_in_father(self.face, v, f"Vertex_{i}")
        # Update the SALOME view
        update_salome_study()

    def show_face(self) -> None:
        """
        Method that adds the face object corresponding to the geometric
        surface this class is representing to the current SALOME study.
        """
        # If no custom name has been provided, use the default one of the
        # object.
        if self.name == "":
            self.name = get_shape_name(self.face)
        # Delete the face from the study if already present
        if self.face_entry_id:
            remove_from_study(self.face_entry_id)
        # Add the surface to the current SALOME study
        self.face_entry_id = add_to_study(self.face, self.name)
        # Update the SALOME view
        update_salome_study()

    @abstractmethod
    def update_from_face(self, face: Any) -> None:
        """
        Method for updating the geometric characteristics of the surface
        from the given face object.

        Parameters
        ----------
        face  : Any
                The new face object to substitute the current face with
        """


class Circle(Surface):
    """
    Class representing a circle surface in SALOME. It is built by
    providing its center, normal vector and radius: if no center
    center (X-Y-Z coordinates) or normal vector (X-Y-Z components)
    are provided, the circle is built in the origin, normal to the
    Z plane.

    Parameters
    ----------
    center      : Union[Tuple[float, float, float], None] = None
                  The X-Y-Z coordinates of the circle.
    normal_vect : Union[Tuple[float, float, float], None] = None
                  The normal vector to the surface of the circle
    radius      : float = 1.0
                  The radius of the circle
    name        : str = "Circle"
                  The name to be assigned to the face when displayed
                  in the current SALOME study

    Attributes
    ----------
    o             : Any
                    A vertex object representing the surface center
    borders       : List[Any]
                    A list of edge objects representing the border edges
                    of the geometric surface
    face          : Any
                    A GEOM object representing the face of the geometry
    name          : str
                    The name of the geometry surface when displayed in
                    the SALOME study
    face_entry_id : Union[str, None]
                    An ID associated to the geometry surface in the SALOME
                    study
    vertices      : List[Any]
                    A list of GEOM vertex objects representing the points
                    of the generic surface (only the circle center)
    out_circle    : Any
                    A GEOM object representing the construction circle the
                    generic surface is inscribed into (the circle itself)
    ON            : Any
                    A GEOM vector object representing the normal vector
                    to the surface of the circle
    radius        : float
                    The radius of the circle
    lx            : float
                    The characteristic dimension of the circle along the
                    X-axis
    ly            : float
                    The characteristic dimension of the circle along the
                    Y-axis
    """
    def __init__(self,
                 center: Union[Tuple[float, float, float], None] = None,
                 normal_vect: Union[Tuple[float, float, float], None] = None,
                 radius: float = 1.0,
                 name: str = "Circle"):
        super().__init__(center)
        self.radius: float = radius
        # ----------------------------------------------
        # Build the GEOM objects representing the circle
        # ----------------------------------------------
        if not normal_vect:
            # The Z unit vector is selected
            normal_vect = (0.0, 0.0, 1.0)
        self.on: Any = make_vector(normal_vect)
        self.borders = self._build_borders()
        self.face = make_face(self.borders)
        self.vertices.append(self.o)
        self.out_circle = self.borders[0]
        # Store the characteristic dimensions of the circle
        self.lx = radius
        self.ly = radius

        # Store the circle name
        self.name = name

    def _build_borders(self) -> List[Any]:
        """
        Method that builds the border edge of the circle surface, being
        its circumference.

        Returns
        -------
        A list of GEOM edge objects representing the surface borders,
        i.e. the circle circumference.
        """
        return [make_circle(self.o, self.on, self.radius)]

    def update_from_face(self, face: Any) -> None:
        """
        Method for updating the geometrical characteristics of the surface
        from the given GEOM face object.

        Parameters
        ----------
        face  : Any
                The new GEOM face object to substitute
        """
        # Check whether the received argument is a FACE-type object
        check_shape_expected_types(face, [ShapeType.FACE])
        # Update the GEOM face object
        self.face = face
        # Re-evaluate all the geometrical characteristics from the face
        self.o = make_cdg(face)
        self.borders = extract_sub_shapes(face, ShapeType.EDGE)
        self.vertices = [self.o]
        self.out_circle = self.borders[0]
        self.radius = get_min_distance(self.o, self.borders[0])
        self.lx = self.radius
        self.ly = self.radius


class Rectangle(Surface):
    """
    Class representing a 2D rectangle surface in SALOME. It is built by
    providing its center, height and width: if no center (X-Y-Z coordinates)
    is provided, the rectangle is built centered in the XYZ origin.
    By default, its height and width have the same value, thus allowing to
    represent a square surface.
    In addition, it supports the presence of rounded corners by providing a
    list of tuples, where the first element is the ID of the corner, while
    the second element the radius of the rounded corner.

    The convention for numbering the vertices and edges of the rectangle is
    the following:

    .. code-block:: text

        3    (2)    2
        + <———————— +
        |           ^
        | (3)       |
        |           | (1)
        v           |
        + ————————> +
        0    (0)    1

    Parameters
    ----------
    center          : Union[Tuple[float, float, float], None] = None
                      The X-Y-Z coordinates of the surface center
    height          : float = 1.0
                      The rectangle height
    width           : float = 1.0
                      The rectangle width
    rounded_corners : List[Tuple[int, float]] | None = None
                      List of tuples containing the rectangle corner ID and
                      the corresponding radius for building a rounded corner
    name            : str = "Rectangle"
                      The name to be assigned to the face when displayed
                      in the current SALOME study

    Attributes
    ----------
    O             : Any
                    A GEOM vertex object representing the surface center
    borders       : List[Any]
                    A list of GEOM objects representing the border edges
                    of the geometry surface
    face          : Any
                    A GEOM object representing the face of the geometry
    name          : str
                    The name of the geometry surface when displayed in
                    the SALOME study
    face_entry_id : Union[str, None]
                    An ID associated to the geometry surface in the SALOME
                    study
    vertices      : List[Any]
                    A list of GEOM vertex objects collecting all the surface
                    construction points (corners + arc points)
    out_circle    : Any
                    A GEOM object representing the construction circle the
                    generic surface is inscribed into
    lx            : float
                    The characteristic dimension of the rectangle along the
                    X-axis
    ly            : float
                    The characteristic dimension of the rectangle along the
                    Y-axis
    """
    def __init__(self,
                 center: Union[Tuple[float, float, float], None] = None,
                 height: float = 1.0,
                 width: float = 1.0,
                 rounded_corners: Union[List[Tuple[int, float]], None] = None,
                 name: str = "Rectangle") -> None:
        super().__init__(center)
        # Store the characteristic dimensions of the rectangle
        self.lx = width
        self.ly = height
        # -------------------------------------------------
        # Build the GEOM objects representing the rectangle
        # -------------------------------------------------
        # Build the rectangle corner vertices
        o_xyz = get_point_coordinates(self.o)
        v0 = make_vertex((o_xyz[0] - width/2, o_xyz[1] - height/2, 0.0))
        v1 = make_vertex((o_xyz[0] + width/2, o_xyz[1] - height/2, 0.0))
        v2 = make_vertex((o_xyz[0] + width/2, o_xyz[1] + height/2, 0.0))
        v3 = make_vertex((o_xyz[0] - width/2, o_xyz[1] + height/2, 0.0))
        self.vertices: List[Any] = [v0, v1, v2, v3]
        # Build the rectangle border edges
        self.borders: List[Any] = self._build_borders()
        # Build the rectangle face
        self.face = make_face(self.borders)
        # Build the arc edges, if any rounded corner is provided and update
        # the rectangular shape
        if rounded_corners is not None:
            rc = self.__build_borders_with_rounded_corners(
                rounded_corners, o_xyz, height, width)
            face = make_partition([self.face], rc, ShapeType.FACE)
            self.face = sorted(
                extract_sub_shapes(face, ShapeType.FACE),
                key=lambda subface: get_min_distance(self.o, subface)
            )[0]
            self.borders = extract_sub_shapes(self.face, ShapeType.EDGE)
        # Build the construction circle the rectangle is inscribed into
        diag = 0.5 * math.sqrt(height*height + width*width)
        self.out_circle: Any = make_circle(self.o, None, diag)
        self.name: str = name

    def _build_borders(self):
        """
        Method that builds the borders of the rectangle from its vertices.

        Returns
        -------
        A list of edge objects representing the rectangle borders.
        """
        l0 = make_edge(self.vertices[0], self.vertices[1])
        l1 = make_edge(self.vertices[1], self.vertices[2])
        l2 = make_edge(self.vertices[2], self.vertices[3])
        l3 = make_edge(self.vertices[3], self.vertices[0])
        return [l0, l1, l2, l3]

    def __build_borders_with_rounded_corners(
            self,
            rounded_corners: List[Tuple[int, float]],
            center: Tuple[float, float, float],
            height: float,
            width: float) -> List[Any]:
        """
        Build the surface borders as GEOM edge objects in the case the
        rectangle has rounded corners. The information about which corners
        are rounded and the corresponding radius is provided by an input
        list.
        This method builds the arc edge elements of the surface, stored
        within the corresponding instance attribute.

        Parameters
        ----------
        rounded_corners : List[Tuple[int, float]]
                          List of tuples containing the rectangle corner ID
                          and the corresponding radius for building a rounded
                          corner
        center          : Tuple[float, float, float]
                          The X-Y-Z coordinates of the circle center
        height          : float
                          The rectangle height
        width           : float
                          The rectangle width

        Returns
        -------
        A list of GEOM edge objects representing the arcs of circle for the
        rounded corners on the borders of the rectangle.
        """
        arcs = []
        max_radius = min(self.lx/2, self.ly/2)
        for rc in rounded_corners:
            # Check the radius of the corner does not exceed the minimum
            # among the characteristic dimensions
            if rc[1] > max_radius:
                raise RuntimeError(
                    f"The corner no. {rc[0]} has a radius of {rc[1]} that "
                    f"exceeds the maximum allowed value of {max_radius}.")
            # Build the XY coordinates of the arc's start/end points and its
            # center
            match rc[0]:
                case 0:
                    # Build the GEOM objects to create an arc
                    xy1 = center[0] - width/2, center[1] - height/2 + rc[1]
                    xy2 = center[0] - width/2 + rc[1], center[1] - height/2
                    cxy = (center[0] - width/2 + rc[1],
                           center[1] - height/2 + rc[1])
                case 1:
                    xy1 = center[0] + width/2 - rc[1], center[1] - height/2
                    xy2 = center[0] + width/2, center[1] - height/2 + rc[1]
                    cxy = (center[0] + width/2 - rc[1],
                           center[1] - height/2 + rc[1])
                case 2:
                    xy1 = center[0] + width/2, center[1] + height/2 - rc[1]
                    xy2 = center[0] + width/2 - rc[1], center[1] + height/2
                    cxy = (center[0] + width/2 - rc[1],
                           center[1] + height/2 - rc[1])
                case 3:
                    xy1 = center[0] - width/2 + rc[1], center[1] + height/2
                    xy2 = center[0] - width/2, center[1] + height/2 - rc[1]
                    cxy = (center[0] - width/2 + rc[1],
                           center[1] + height/2 - rc[1])
            v1 = make_vertex((*xy1, 0.0))
            v2 = make_vertex((*xy2, 0.0))
            cv = make_vertex((*cxy, 0.0))
            # Add the arc construction points to the list of GEOM points
            self.vertices.append(v1)
            self.vertices.append(v2)
            self.vertices.append(cv)
            arcs.append(make_arc_edge(cv, v1, v2))
        # Return the built arcs
        return arcs

    def show_edges_and_vertices(self) -> None:
        """
        Method that adds the rectangle edges and the corners vertex to the
        current SALOME study.
        """
        self.show_borders()
        for i, v in enumerate(self.vertices):
            add_to_study_in_father(self.face, v, f"Vertex_{i}")
        # Update the SALOME view
        update_salome_study()

    def update_from_face(self, face: Any) -> None:
        """
        Method for updating the geometrical characteristics of the surface
        from the given GEOM face object.

        Parameters
        ----------
        face  : Any
                The new GEOM face object to substitute
        """
        # Check whether the received argument is a FACE-type object
        check_shape_expected_types(face, [ShapeType.FACE])
        # Update the GEOM face object
        self.face = face
        # Re-evaluate all the geometrical characteristics from the face
        self.o = make_cdg(face)
        self.borders = extract_sub_shapes(face, ShapeType.EDGE)
        # Chech if the borders represent a rectangular shape
        if len(self.borders) != 4 or len({
            get_basic_properties(b)[0] for b in self.borders}) != 2:
            raise RuntimeError(
                "The provided face does not represent a rectangular shape.")
        # Get the buonding box dimension
        b_box = get_bounding_box(self.face)
        # Update the characteristic dimensions of the rectangle
        self.lx = b_box[1] - b_box[0]
        self.ly = b_box[3] - b_box[2]

        # Update the construction circle the rectangle is inscribed into
        self.out_circle = make_circle(
            self.o, None, 0.5*math.sqrt(self.lx*self.lx + self.ly*self.ly))
        # Update the list of vertices representing the rectangle corners
        self.vertices = extract_sub_shapes(self.face, ShapeType.VERTEX)


class Hexagon(Surface):
    r"""
    Class representing a 2D hexagon surface in SALOME. It is built by
    providing its center and the length of edge: if no center (X-Y-Z
    coordinates) is provided, the hexagon is built centered in the XYZ
    origin.

    The hexagon is built by means of a construction circle the hexagon
    is inscribed into. By default, the hexagon orientation is zero, i.e.
    the (1)-(4) vertices are horizontally aligned.

    .. code-block:: text

            (3)     (2)
             + ————— +
            /         \
           /           \
      (4) +             + (1)
           \           /
            \         /
             + ————— +
            (5)     (6)


    Parameters
    ----------
    center          : Union[Tuple[float, float, float], None] = None
                      The X-Y-Z coordinates of the hexagon center
    edge_length     : float = 1.0
                      The lenght of the hexagon edge
    name            : str = "Hexagon"
                      The name to be assigned to the face when displayed
                      in the current SALOME study

    Attributes
    ----------
    O             : Any
                    A GEOM vertex object representing the surface center
    borders       : List[Any]
                    A list of GEOM objects representing the border edges
                    of the geometry surface
    face          : Any
                    A GEOM object representing the face of the geometry
    name          : str
                    The name of the geometry surface when displayed in
                    the SALOME study
    face_entry_id : Union[str, None]
                    An ID associated to the geometry surface in the SALOME
                    study
    vertices      : List[Any]
                    A list of GEOM vertex objects collecting all the surface
                    construction points (corners)
    radius        : float
                    The radius of the construction circle the hexagon is
                    inscribed into
    apothem       : float
                    The apothem of the hexagon
    out_circle    : float
                    The construction circle the hexagon is inscribed into
    lx            : float
                    The characteristic dimension of the hexagon along the
                    X-axis
    ly            : float
                    The characteristic dimension of the hexagon along the
                    Y-axis
    """
    def __init__(self, center: Union[Tuple[float, float, float], None] = None,
                 edge_length: float = 1.0,
                 name: str = "Hexagon") -> None:
        super().__init__(center)
        # Calculate the radius of the construction circle the hexagon is
        # inscribed into
        self.radius: float = edge_length
        # Calculate the apothem of the hexagon
        self.apothem: float = edge_length * math.sin(math.pi/3)
        # Build the construction circle the hexagon is inscribed into
        self.out_circle: Any = make_circle(self.o, None, self.radius)
        # Build the list of vertices representing the hexagon corners
        self.vertices: List[Any] = [
            make_vertex_on_curve(self.out_circle, i/6) for i in range(6)]
        # Build the list of edges connecting successive vertices
        self.borders: List[Any] = self._build_borders()
        # Build the hexagon face
        self.face = make_face(self.borders)
        # Store the characteristic dimensions of the hexagon
        self.lx = edge_length
        self.ly = self.apothem

        self.name: str = name

    def _build_borders(self) -> List[Any]:
        """
        Method that builds the borders of the hexagon.

        Returns
        -------
        A list of GEOM edge objects representing the hexagon surface borders.
        """
        return [make_edge(
            self.vertices[i], self.vertices[(i+1) % 6]) for i in range(6)]

    def update_from_face(self, face: Any) -> None:
        """
        Method for updating the geometrical characteristics of the surface
        from the given GEOM face object.

        Parameters
        ----------
        face  : Any
                The new GEOM face object to substitute
        """
        # Check whether the received argument is a FACE-type object
        check_shape_expected_types(face, [ShapeType.FACE])
        # Update the GEOM face object
        self.face = face
        # Re-evaluate all the geometrical characteristics from the face
        self.o = make_cdg(self.face)
        self.borders = extract_sub_shapes(self.face, ShapeType.EDGE)
        # Chech if the borders represent a hexagonal shape
        if len(self.borders) != 6 or len({
            round(get_basic_properties(b)[0], 6) for b in self.borders}) != 1:
            raise RuntimeError(
                "The provided face does not represent a hexagonal shape.")
        self.radius = get_basic_properties(self.borders[0])[0]
        self.apothem = 0.5 * self.radius / math.tan(math.radians(30))
        # Update the characteristic dimensions of the hexagon
        self.lx = self.radius
        self.ly = self.apothem
        # Update the construction circle the hexagon is inscribed into
        self.out_circle = make_circle(self.o, None, self.radius)
        # Update the list of vertices representing the hexagon corners
        self.vertices = extract_sub_shapes(self.face, ShapeType.VERTEX)


class GenericSurface(Surface):
    """
    Class representing a 2D generic surface defined starting from a given
    face object.
    Borders and vertices are directly extracted from the shape, whereas no
    construction circle is declared: this is because the shape might not have
    a regular geometric shape, hence there could be no circle within which
    the shape is perfectly inscribed.
    The characteristic dimensions of the surface are taken from the bounding
    box enclosing the given surface.

    Parameters
    ----------
    face : Any
        The face object representing the generic surface

    Attributes
    ----------
    o             : Any
                    A vertex object representing the surface center
    borders       : List[Any]
                    A list of edge objects representing the border edges
                    of the surface
    face          : Any
                    A face object representing the geometric surface
    name          : str
                    The name of the surface when displayed in the SALOME study
    face_entry_id : Union[str, None]
                    An ID associated to the surface in the SALOME study
    vertices      : List[Any]
                    A list of vertex objects representing the vertices of the
                    generic surface
    out_circle    : Any
                    An edge object representing the construction circle which
                    the generic surface is inscribed into
    rotation      : float
                    The rotation angle of the geometric surface in radians
    lx            : float
                    The characteristic dimension of the generic geometric
                    surface along the X-axis
    ly            : float
                    The characteristic dimension of the generic geometric
                    surface along the Y-axis
    """
    def __init__(self, face: Any):
        super().__init__(get_point_coordinates(make_cdg(face)))
        self.face = face
        self.name = get_shape_name(face)
        # Build the list of vertices representing the hexagon corners
        self.vertices = extract_sub_shapes(self.face, ShapeType.VERTEX)
        # Build the list of edges connecting successive vertices
        self.borders = self._build_borders()
        # Set the characteristic dimensions of the shape from the bounding box
        # extension
        b_box = get_bounding_box(self.face)
        self.lx = (b_box[1] - b_box[0]) / 2
        self.ly = (b_box[3] - b_box[2]) / 2

        self.rotation: float = 0.0

    def _build_borders(self) -> List[Any]:
        """
        Method that extracts all the edge objects of the generic surface
        this class is representing.

        Returns
        -------
        A list of the edge objects of the surface this class is built on.
        """
        return extract_sub_shapes(self.face, ShapeType.EDGE)

    def update_from_face(self, face: Any) -> None:
        """
        Method for updating the geometric characteristics of the surface
        from the given face object.

        Parameters
        ----------
        face  : Any
                The new face object to substitute the current face with
        """
        # Check whether the received argument is a FACE or COMPOUND-type
        # object
        check_shape_expected_types(face, [ShapeType.FACE, ShapeType.COMPOUND])
        # Update the GEOM face object
        self.face = face
        # Re-evaluate all the geometrical characteristics from the face
        self.o = make_cdg(face)
        self.vertices = extract_sub_shapes(self.face, ShapeType.VERTEX)
        self.borders = self._build_borders()
        # Store the characteristic dimensions of the hexagon
        b_box = get_bounding_box(self.face)
        self.lx = b_box[1] - b_box[0]
        self.ly = b_box[3] - b_box[2]


def build_hexagon(
        apothem: float,
        center: Tuple[float, float, float] | None = None) -> Hexagon:
    """
    Function that allows to build an hexagon from its apothem. The resulting
    face is placed at the indicated center, if any, otherwise it is centered
    in the XYZ space origin.

    Parameters
    ----------
    apothem : float
        The value of the hexagon apothem
    center  : Tuple[float, float, float] | None = None
        The XYZ coordinates of the hexagon center, if any

    Returns
    -------
    A 'Hexagon' object with dimensions and center as indicated.
    """
    # Calculate the value of the hexagon side
    hex_side = 2 * apothem / math.tan(math.pi/3)
    # Build and return the 'Hexagon' object
    return Hexagon(center, hex_side)


if __name__ == "__main__":
    # ------------------------------------------------
    # Testing the geometry surface classes for SALOME.
    # ------------------------------------------------
    # Build an hexagon and show it in the SALOME study
    shape = Hexagon(edge_length=1.5)
    shape.rotate(-90)
    shape.show_face()
    shape.show_edges_and_vertices()

    r = Rectangle(height=1, width=2, rounded_corners=[
        (0, 0.1), (1, 0.1), (2, 0.1), (3, 0.1)])
    r.show_face()

    # Show everything on the SALOME application
    update_salome_study()
