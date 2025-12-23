"""
Module containing classes providing the means for creating basic geometric
shapes in SALOME.
"""
import math

from typing import Any, List, Tuple

from glow.geometry_layouts.layouts import Layout
from glow.interface.geom_entities import Edge, Face, wrap_shape
from glow.interface.geom_interface import ShapeType, add_to_study, clear_view, \
    display_shape, extract_sub_shapes, get_angle_between_shapes, \
    get_basic_properties, get_bounding_box, get_min_distance, \
    get_object_from_id, get_point_coordinates, get_shape_name, get_shape_type, \
    make_cdg, make_circle, make_edge, make_face, make_partition, make_rotation, \
    make_scale, make_translation,make_vector_from_points, make_vertex, \
    make_vertex_on_curve, remove_from_study, update_salome_study
from glow.support.utility import build_arcs_for_rounded_corners, \
    build_contiguous_edges, build_z_axis_from_vertex, \
    check_shape_expected_types, sort_shapes_from_vertex


class Surface(Face, Layout):
    """
    Abstract class for representing any geometric surface.

    Parameters
    ----------
    geom_obj: Face | None
        The ``Face`` object (or the GEOM face directly) the surface
        corresponds to.
    center : Union[Tuple[float, float, float], None] = None
        The XYZ coordinates of the generic geometric surface center.

    Attributes
    ----------
    borders : List[Any]
        A list of edge objects representing the border edges of the surface.
    dimensions : Tuple[float, float]
        The X-Y characteristic dimensions of the geometric surface.
    entry_id : str | None
        The ID attributed by SALOME when the geometric surface is added
        to the study.
    name : str
        The name assigned to the geometric surface when added to the study.
    o : Vertex
        The ``Vertex`` object representing the centre of the geometric
        surface.
    out_circle : Edge | None
        The ``Edge`` object representing the construction circle which the
        geometric surface is inscribed into.
    rot_angle : float
        The rotation angle (in degrees) of the geometric surface wrt the
        X-axis.
    """
    def __init__(
            self,
            geom_obj: Face | None,
            center: Tuple[float, float, float] | None = None) -> None:
        super().__init__(geom_obj)
        if geom_obj:
            self._initialize_geom_object(geom_obj)
        # Initialize the instance attributes
        if not center:
            center = (0.0, 0.0, 0.0)
        self.borders: List[Edge] = []
        self.o = wrap_shape(make_vertex(center))
        self.name = "Surface"
        self.out_circle: Edge | None = None

    def rotate(self, angle: float, axis: Edge | None = None) -> None:
        """
        Method that rotates the surface geometry. The face, its border
        edges, vertices and construction circle are rotated by the given
        angle expressed in degrees around the given axis, if any is provided,
        otherwise the Z-axis.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        axis : Edge | None = None
            The ``Edge`` object representing the rotation axis, if any.
        """
        # Return immediately if the angle is zero
        if math.isclose(angle, 0.0, abs_tol=1e-6):
            return
        # Build a Z-axis, if none is provided
        if not axis:
            # Build the Z-axis of rotation positioned in the figure center
            axis = build_z_axis_from_vertex(self.o)
        # Rotate the surface elements
        self._rotate_from_axis(angle, axis)

    def scale(self, factor: float) -> None:
        """
        Method for scaling the surface by the given factor.

        Parameters
        ----------
        factor : float
            The scaling factor.
        """
        # Apply the scaling and re-build the borders
        self.geom_obj = make_scale(self, self.o, factor)
        self.borders = extract_sub_shapes(self.geom_obj, ShapeType.EDGE)

    def show(self, *args: Any) -> None:
        """
        Method for displaying the layout in the 3D viewer of SALOME
        according to the given settings.

        Parameters
        ----------
        *args : Any
            Positional arguments. Must be empty.

        Raises
        ------
        ValueError
            If any positional arguments are provided.
        """
        # Check that no args are provided
        if args:
            raise ValueError(f"No arguments are accepted for 'show()'.")
        # Erase all objects from the current view
        clear_view()
        # Delete the face from the study if already present
        if self.entry_id and get_object_from_id(self.entry_id):
            remove_from_study(self.entry_id)
        # Add the surface to the current SALOME study and display it
        self.entry_id = add_to_study(self, self.name)
        display_shape(self.entry_id)
        # Update the SALOME view
        update_salome_study()

    def translate(self, new_pos: Tuple[float, float, float]) -> None:
        """
        Method that translates the surface geometry to the given position.
        All the geometrical elements (center, face, its borders, vertices
        and construction circle) needed to describe the shape are moved
        accordingly.

        Parameters
        ----------
        new_pos : Tuple[float, float, float]
            The XYZ coordinates of the new center of the shape.
        """
        # Build a vector from the current center to the new one
        transl_vect = make_vector_from_points(self.o, make_vertex(new_pos))
        # Translate the characteristic geometrical elements of the shape
        self.o = wrap_shape(make_vertex(new_pos))
        self.geom_obj = make_translation(self, transl_vect)
        # Translate the construction circle
        if self.out_circle:
            self.out_circle = make_translation(self.out_circle, transl_vect)
        # Re-build the borders
        self.borders = extract_sub_shapes(self.geom_obj, ShapeType.EDGE)

    def _initialize_geom_object(self, geom_obj: Face) -> None:
        """
        Method that sets the GEOM object the instance refers to. It checks
        if the given object is a face; if not, an exception is raised.

        Parameters
        ----------
        geom_obj : Face
            The GEOM face object this instance refers to.

        Raises
        ------
        RuntimeError
            If the given GEOM object is not a face.
        """
        # Check the given GEOM object is a face
        if get_shape_type(geom_obj) != ShapeType.FACE:
            raise RuntimeError("A 'Surface' object must be of type 'FACE'")
        self.geom_obj = geom_obj

    def _rotate_from_axis(self, angle: float, axis: Edge) -> None:
        """
        Method that rotates the surface geometry. The face, its border
        edges, vertices and construction circle are rotated by the given
        angle expressed in degrees around the given axis.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees.
        axis : Edge
            The ``Edge`` object representing the rotation axis.
        """
        # Convert the rotation angle in radians
        angle_rad = math.radians(angle)
        # Rotate the geometric elements of the surface
        self.geom_obj = make_rotation(self, axis, angle_rad)
        # Re-build the borders
        self.borders = extract_sub_shapes(self.geom_obj, ShapeType.EDGE)
        # Rotate the construction circle
        if self.out_circle:
            self.out_circle = make_rotation(self.out_circle, axis, angle_rad)
        # Update the rotation angle of the surface wrt X-axis
        self.rot_angle += angle


class Circle(Surface):
    """
    Class representing a circle surface in SALOME. It is built by providing
    its center and radius: if no center (XYZ coordinates) is provided, the
    circle is built in the origin, normal to the XY plane.

    Parameters
    ----------
    center : Tuple[float, float, float] | None = None
        The XYZ coordinates of the circle center.
    radius : float = 1.0
        The radius of the circle.
    name : str = "Circle"
        The name to be assigned to the face when displayed in the current
        SALOME study.

    Attributes
    ----------
    borders : List[Edge]
        The list of ``Edge`` objects representing the border edges of the
    dimensions : Tuple[float, float]
        The characteristic dimensions of the circle along the X-Y axes.
    entry_id : str | None
        The ID associated to the surface in the SALOME study.
    name : str
        The name of the surface when displayed in the SALOME study.
    o : Vertex
        The ``Vertex`` object representing the surface center.
    out_circle : Edge | None
        The ``Edge`` object representing the construction circle which the
        circle surface is inscribed into (i.e. the circle itself).
    radius : float
        The radius of the circle.
    rot_angle : float
        The rotation angle (in degrees) of the geometric surface wrt the
        X-axis.
    """
    def __init__(self,
                 center: Tuple[float, float, float] | None = None,
                 radius: float = 1.0,
                 name: str = "Circle") -> None:
        super().__init__(None, center)
        # Build the Z-axis of rotation positioned in the figure center
        axis = build_z_axis_from_vertex(self.o)
        # Initialize instance attribute
        self.radius: float = radius
        self.borders = [wrap_shape(make_circle(self.o, axis, self.radius))]
        self.out_circle = self.borders[0]
        self.dimensions = (radius, radius)
        self._initialize_geom_object(make_face(self.borders))
        self.name = name

    def update(self, layout: Face) -> None:
        """
        Method for updating the geometrical characteristics of the surface
        from the given GEOM face object.

        Parameters
        ----------
        face : Face
            The new ``Face`` object to substitute the current with.

        Raises
        ------
        RuntimeError
            If the provided shape type is not ``ShapeType.FACE``.
        """
        # Check whether the received argument is a FACE-type object
        check_shape_expected_types(layout, [ShapeType.FACE])
        # Update the GEOM face object
        self.geom_obj = layout
        # Re-evaluate all the geometrical characteristics from the face
        self.o = wrap_shape(make_cdg(layout))
        self.borders = extract_sub_shapes(layout, ShapeType.EDGE)
        self.out_circle = self.borders[0]
        self.radius = get_min_distance(self.o, self.borders[0])
        self.dimensions = (self.radius, self.radius)


class Rectangle(Surface):
    """
    Class representing a 2D rectangle surface in SALOME. It is built by
    providing its center, height and width: if no center (XYZ coordinates)
    is provided, the rectangle is built centered in the XYZ origin.
    By default, its height and width have the same value, thus allowing to
    represent a square surface.
    In addition, it supports the construction of rounded corners by the input
    list of tuples. Each tuple has as first element the ID of the corner and
    as second element the radius of the rounded corner.

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
    center : Union[Tuple[float, float, float], None] = None
        The XYZ coordinates of the surface center.
    height : float = 1.0
        The height of the rectangle.
    width : float = 1.0
        The width of the rectangle.
    rounded_corners : List[Tuple[int, float]] | None = None
        List of tuples with corner ID and corresponding curvature radius.
        The numbering convention for ID is presented above.
    name : str = "Rectangle"
        The name to be assigned to the face when displayed in the current
        SALOME study.

    Attributes
    ----------
    borders : List[Edge]
        The list of ``Edge`` objects representing the border edges of the
        surface.
    dimensions : Tuple[float, float]
        The characteristic dimensions of the rectangle along the X-Y axes (
        i.e. the width and the height).
    entry_id : str | None
        The ID associated to the surface in the SALOME study.
    name : str
        The name of the surface when displayed in the SALOME study.
    o : Vertex
        The ``Vertex`` object representing the surface center.
    out_circle : Edge | None
        The ``Edge`` object representing the construction circle the rectangle
        is inscribed into.
    rot_angle : float
        The rotation angle (in degrees) of the geometric surface.
    """
    def __init__(
            self,
            center: Tuple[float, float, float] | None = None,
            height: float = 1.0,
            width: float = 1.0,
            rounded_corners: List[Tuple[int, float]] | None = None,
            name: str = "Rectangle"
        ) -> None:
        super().__init__(None, center)
        # Store the characteristic dimensions of the rectangle
        self.dimensions = (width, height)
        # Build the rectangle corner vertices
        o_xyz = get_point_coordinates(self.o)
        v0 = make_vertex((o_xyz[0] - width/2, o_xyz[1] - height/2, 0.0))
        v1 = make_vertex((o_xyz[0] + width/2, o_xyz[1] - height/2, 0.0))
        v2 = make_vertex((o_xyz[0] + width/2, o_xyz[1] + height/2, 0.0))
        v3 = make_vertex((o_xyz[0] - width/2, o_xyz[1] + height/2, 0.0))
        # Build the rectangle border edges and face
        borders = build_contiguous_edges([v0, v1, v2, v3])
        rect_face = make_face(borders)

        # Build the arc edges, if any rounded corner is provided and update
        # the rectangular shape
        if rounded_corners is not None:
            arcs = build_arcs_for_rounded_corners(
                rounded_corners, o_xyz, height, width)
            # Partition the rectangle with the arcs and extract the cut shape
            p_face = make_partition(
                [make_face(borders)], arcs, ShapeType.FACE
            )
            rect_face = sort_shapes_from_vertex(
                extract_sub_shapes(p_face, ShapeType.FACE), self.o)[0]
        # Initialize borders and GEOM face attributes
        self.borders = [
            wrap_shape(e)
            for e in extract_sub_shapes(rect_face, ShapeType.EDGE)
        ]
        self._initialize_geom_object(rect_face)

        # Build the construction circle the rectangle is inscribed into
        diag = 0.5 * math.sqrt(height*height + width*width)
        self.out_circle = make_circle(
            self.o, build_z_axis_from_vertex(self.o), diag
        )
        # Set the surface name
        self.name = name

    def update(self, layout: Face) -> None:
        """
        Method for updating the geometrical characteristics of the surface
        from the given GEOM face object.

        Parameters
        ----------
        layout : Face
            The new face object to substitute the current face with.

        Raises
        ------
        RuntimeError
            If the provided shape type is not ``ShapeType.FACE``.
        RuntimeError
            If the provided shape does not represent a rectangle.
        """
        # Check whether the received argument is a FACE-type object
        check_shape_expected_types(layout, [ShapeType.FACE])
        # Check if the borders represent a rectangular shape
        borders = extract_sub_shapes(layout, ShapeType.EDGE)
        borders_lengths = {
            round(get_basic_properties(b)[0], 6): b for b in borders}
        if len(borders) != 4 or len(borders_lengths) != 2:
            raise RuntimeError(
                "The provided face does not represent a rectangular shape.")

        # Update the GEOM face object
        self.geom_obj = layout
        # Re-evaluate all the geometrical characteristics from the face
        self.o = wrap_shape(make_cdg(self.geom_obj))
        self.borders = borders

        # Update the construction circle the rectangle is inscribed into by
        # calculating the diagonal of the rectangle
        l1, l2 = tuple(borders_lengths.keys())
        self.out_circle = make_circle(
            self.o,
            build_z_axis_from_vertex(self.o),
            math.sqrt(l1*l1 + l2*l2)/2
        )
        # Update the characteristic dimensions by checking if any border is
        # parallel to the X-axis
        o_x, o_y, _ = get_point_coordinates(self.o)
        x_axis = make_vector_from_points(
            self.o, make_vertex((o_x+1, o_y, 0.0)))
        for l, b in borders_lengths.items():
            if math.isclose(
                get_angle_between_shapes(b, x_axis), 0.0, abs_tol=1e-6
            ):
                self.dimensions[0] = l
                break
            self.dimensions[1] = l
        else:
            self.dimensions = list(borders_lengths.keys())


class Hexagon(Surface):
    r"""
    Class representing a 2D hexagon surface in SALOME. It is built by
    providing its center and the length of edge: if no center (XYZ
    coordinates) is provided, the hexagon is built centered in the XYZ
    origin.

    The hexagon is built by means of a construction circle the hexagon
    is inscribed into.
    The convention for numbering the vertices of the hexagon is the
    following:

    .. code-block:: text

             3       2
             + ————— +
            /         \
           /           \
        4 +             + 1
           \           /
            \         /
             + ————— +
             5       6

    By default, the hexagon orientation is zero, i.e. the vertices `1` and `4`
    are horizontally aligned.

    Parameters
    ----------
    center : Union[Tuple[float, float, float], None] = None
        The XYZ coordinates of the center of the hexagon.
    edge_length : float = 1.0
        The lenght of the edge of the hexagon.
    name : str = "Hexagon"
        The name to be assigned to the face when displayed in the current
        SALOME study.

    Attributes
    ----------
    borders : List[Any]
        The list of edge objects representing the border edges of the surface.
    dimensions : Tuple[float, float]
        The characteristic dimensions of the hexagon along the X-Y axes (i.e.
        the edge and the apothem).
    entry_id : str | None
        The ID associated to the surface in the SALOME study.
    name : str
        The name of the surface when displayed in the SALOME study.
    o : Vertex
        The ``Vertex`` object representing the surface center.
    out_circle : Edge | None
        The ``Edge`` object representing the construction circle the hexagon
        is inscribed into.
    rot_angle : float
        The rotation angle (in degrees) of the geometric surface.
    """
    def __init__(
            self,
            center: Tuple[float, float, float] | None = None,
            edge_length: float = 1.0,
            name: str = "Hexagon"
        ) -> None:
        super().__init__(None, center)
        # Build the construction circle the hexagon is inscribed into
        self.out_circle = make_circle(
            self.o, build_z_axis_from_vertex(self.o), edge_length
        )
        # Build the list of vertices representing the hexagon corners
        vertices: List[Any] = [
            make_vertex_on_curve(self.out_circle, i/6) for i in range(6)]
        # Build the list of edges connecting successive vertices
        self.borders = [
            wrap_shape(
                make_edge(vertices[i], vertices[(i+1) % 6]) for i in range(6)
            )
        ]
        # Build the hexagon face
        self._initialize_geom_object(make_face(self.borders))
        # Store the characteristic dimensions of the hexagon
        self.dimensions = (edge_length, edge_length * math.sin(math.pi/3))
        # Set the surface name
        self.name = name

    def update(self, layout: Face) -> None:
        """
        Method for updating the geometrical characteristics of the surface
        from the given GEOM face object.

        Parameters
        ----------
        layout : Face
            The new face object to substitute the current face with.

        Raises
        ------
        RuntimeError
            If the provided shape type is not ``ShapeType.FACE``.
        RuntimeError
            If the provided shape does not represent a hexagon.
        """
        # Check whether the received argument is a FACE-type object
        check_shape_expected_types(layout, [ShapeType.FACE])
        # Check if the borders represent a hexagonal shape
        borders = extract_sub_shapes(layout, ShapeType.EDGE)
        borders_lengths = {
            round(get_basic_properties(b)[0], 6) for b in borders}
        if len(borders) != 6 or len(borders_lengths) != 1:
            raise RuntimeError(
                "The provided face does not represent a hexagonal shape.")

        # Update the GEOM face object
        self.geom_obj = layout
        # Re-evaluate all the geometrical characteristics from the face
        self.o = wrap_shape(make_cdg(self.geom_obj))
        self.borders = borders

        # Update the construction circle the hexagon is inscribed into by
        # considering the length of the border as its radius
        radius = get_basic_properties(self.borders[0])[0]
        self.out_circle = make_circle(self.o, None, radius)
        # Update the characteristic dimensions of the hexagon
        self.dimensions = (
            radius, 0.5 * radius / math.tan(math.radians(30))
        )


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
        The face object representing the generic surface.
    name : str | None = None
        The name to be assigned to the face when displayed in the current
        SALOME study.

    Attributes
    ----------
    borders : List[Any]
        The list of ``Edge`` objects representing the border edges of the
        surface.
    dimensions : Tuple[float, float]
        The characteristic dimensions of the generic surface along the X-Y
        axes (i.e. the dimensions of the bounding box).
    entry_id : str | None
        The ID associated to the surface in the SALOME study.
    name : str
        The name of the surface when displayed in the SALOME study.
    o : Vertex
        The ``Vertex`` object representing the surface center.
    out_circle : Edge | None
        The ``Edge`` object representing the construction circle the generic
        surface is inscribed into. In this case, no circle is declared.
    rot_angle : float
        The rotation angle (in degrees) of the geometric surface.
    """
    def __init__(self, face: Any, name: str | None = None) -> None:
        super().__init__(face, get_point_coordinates(make_cdg(face)))
        if not name:
            name = get_shape_name(face)
        self.name = name
        # Build the list of edges connecting successive vertices
        self.borders = extract_sub_shapes(self.geom_obj, ShapeType.EDGE)
        # Set the characteristic dimensions of the shape from the bounding box
        # extension
        self._set_surface_dimensions()

    def update(self, layout: Face) -> None:
        """
        Method for updating the geometric characteristics of the surface
        from the given GEOM face object.

        Parameters
        ----------
        layout : Face
            The new face object to substitute the current face with.

        Raises
        ------
        RuntimeError
            If the provided shape type is not ``ShapeType.FACE`` or
            ``ShapeType.COMPOUND``.
        """
        # Check whether the received argument is a FACE or COMPOUND-type
        # object
        check_shape_expected_types(
            layout, [ShapeType.FACE, ShapeType.COMPOUND]
        )
        # Update the GEOM face object
        self.geom_obj = layout
        # Re-evaluate all the geometrical characteristics from the face
        self.o = wrap_shape(make_cdg(layout))
        self.borders = extract_sub_shapes(layout, ShapeType.EDGE)
        # Store the characteristic dimensions of the generic surface
        self._set_surface_dimensions()

    def _set_surface_dimensions(self) -> None:
        """
        Method that sets the characteristic dimensions along the XY axes by
        calculating the bounding box of the GEOM face this instance refers
        to.
        """
        b_box = get_bounding_box(self.geom_obj)
        self.dimensions = (
            (b_box[1] - b_box[0]) / 2, (b_box[3] - b_box[2]) / 2)


# -------------------------------------------------------------------------- #
#                                FUNCTIONS                                   #
# -------------------------------------------------------------------------- #

def build_hexagon_from_apothem(
        apothem: float,
        center: Tuple[float, float, float] | None = None) -> Hexagon:
    """
    Function that builds an instance of the ``Hexagon`` class from the value
    of the apothem.
    The resulting face is placed at the indicated center, if any, otherwise
    it is centered in the XYZ space origin.

    Parameters
    ----------
    apothem : float
        The value of the hexagon apothem.
    center : Tuple[float, float, float] | None = None
        The XYZ coordinates of the hexagon center, if any.

    Returns
    -------
    Hexagon
        A ``Hexagon`` object with dimensions and center as indicated.
    """
    # Calculate the value of the hexagon side
    hex_side = 2 * apothem / math.tan(math.pi/3)
    # Build and return the 'Hexagon' object
    return Hexagon(center, hex_side)


def build_right_triangle(
        hypotenuse: float,
        cathetus: float,
        left_corner: Tuple[float, float, float] | None = None
    ) -> GenericSurface:
    """
    Function that builds a ``GenericSurface`` instance representing a right
    triangle from the given values for the hypotenuse and its left cathetus.
    The resulting face is placed with its left corner that coincides with the
    given coordinates, if any, otherwise the corner is the XYZ space origin.
    The hypotenuse of the right triangle is considered to be parallel to the
    X-axis.

    Parameters
    ----------
    hypotenuse : float
        The length of the hypotenuse.
    cathetus : float
        The length of the left cathetus.
    left_corner : Tuple[float, float, float] | None = None
        The XYZ coordinates of the left corner, if any.

    Returns
    -------
    GenericSurface
        The ``GenericSurface`` instance representing a right triangle.
    """
    # Calculate the length of the smallest cathetus
    cat_2_len = math.sqrt(hypotenuse*hypotenuse - cathetus*cathetus)
    # Calculate the lenght of the projection of the smallest cathetus on the
    # hypotenuse
    proj = (cat_2_len*cat_2_len) / hypotenuse
    # Calculate the height relative to the hypotenuse
    height = math.sqrt(cat_2_len*cat_2_len - proj*proj)
    # Build the vertices of the right triangle, considering the coordinates
    # of the left corner, if any
    if not left_corner:
        left_corner = (0.0, 0.0, 0.0)
    lc_x, lc_y = left_corner[0], left_corner[1]
    vertices = [
        make_vertex(left_corner),
        make_vertex((lc_x + hypotenuse, lc_y, 0.0)),
        make_vertex((lc_x + hypotenuse - proj, lc_y + height, 0.0))
    ]
    # Build the face object
    return GenericSurface(make_face(build_contiguous_edges(vertices)))


def build_regular_triangle(
        side_length: float,
        left_corner: Tuple[float, float, float] | None = None
    ) -> GenericSurface:
    """
    Function that builds a ``GenericSurface`` instance representing a regular
    triangle from the given length of its side.
    The resulting face is placed with its left corner that coincides with the
    given coordinates, if any, otherwise the corner is the XYZ space origin.

    Parameters
    ----------
    side_length : float
        The length of the side of the regular triangle.
    left_corner : Tuple[float, float, float] | None = None
        The XYZ coordinates of the left corner, if any.

    Returns
    -------
    GenericSurface
        The ``GenericSurface`` instance representing a regular triangle.
    """
    # Calculate the height of the triangle
    height = math.sqrt(3) / 2 * side_length
    # Build the vertices of the right triangle, considering the coordinates
    # of the left corner, if any
    if not left_corner:
        left_corner = (0.0, 0.0, 0.0)
    lc_x, lc_y = left_corner[0], left_corner[1]
    vertices = [
        make_vertex(left_corner),
        make_vertex((lc_x + side_length, lc_y, 0.0)),
        make_vertex((lc_x + side_length/2, lc_y + height, 0.0))
    ]
    # Build the face object
    return GenericSurface(make_face(build_contiguous_edges(vertices)))
