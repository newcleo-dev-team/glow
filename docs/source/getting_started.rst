===============
Getting started
===============

This chapter presents all the information the user needs to make the best
use of |TOOL|.
Section :ref:`geom-def` provides details about the available functionalities
for setting up a geometry layout and displaying it in the 3D viewer of *SALOME*.
Section :ref:`lattice-export` provides information about the process that
|TOOL| automatically performs to generate the output *.dat* file containing
the description of the geometry layout.
In both sections, the available functionalities for building and exporting a
geometry layout are described with code snippets showing their usage. Images
are also provided to graphically present the results displayed in the 3D viewer
of *SALOME*.
Lastly, section :ref:`usage` gives indications about how to include the functionalities
of |TOOL| in a Python script and how to run the script in the *SALOME*
environment.

.. _geom-def:

Geometry Definition
-------------------

From a topological point of view, the *GEOM* module of *SALOME* enables the
creation of the following entities:

  - *vertex*, which is a point in the XYZ space;
  - *edge*, made by two vertices, can be classified as *segment*, *arc of circle*
    or *circle*;
  - *wire*, a closed set of edges;
  - *face*, a 2D area made from one or two *wires*;
  - *shell*, made from a group of *faces*;
  - *solid*, a 3D object made from a group of *shells*;
  - *compound*, a container grouping together several *GEOM* entities.

|TOOL| relies on most of the above-mentioned topological entities to assemble
the geometry layouts and visualize them in the 3D viewer of *SALOME*.
In addition, the specific *GEOM* functions enabling the operations for building
topological entities and applying operations among them are provided in module
:py:mod:`geom_interface<glow.interface.geom_interface>`.

In |TOOL|, a *cell*, identified by any of the subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>`,
is the base unit of the geometry layout.
*Cells* are built from a characteristic *surface*, i.e. a subclass of
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>`, which represents
a *GEOM face*. Cells can either have a rectangular or an hexagonal shape.
Several *surfaces* can be juxtaposed or overlapped to determine the final layout
of a *cell* by assembling the corresponding *GEOM faces* together.
The construction of a *cell* geometry layout relies on boolean operations, in
particular the *partition* (function :py:func:`make_partition()<glow.interface.geom_interface.make_partition>`)
one. This operation cuts one *GEOM face* with another, and viceversa, grouping
both cut faces with the intersected area. The result is a *GEOM compound*, which
represents the geometry layout of a *cell*, as made of *GEOM faces* identifying
the different areas.

In |TOOL|, the geometry layout of a *cell* is described in terms of:

  - the **technological geometry**, which is the one delimiting the different
    cell's *regions* in terms of the materials;
  - the **sectorized geometry**, which further subdivides the cell's *regions*
    of the technological geometry into sectors.

*Lattices*, identified by :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>`
instances, are made by the repetition of adjacent *cells* in the 2D space.
As for *cells*, also *lattices* can be displayed in terms of one of these two
types of geometry by accessing the corresponding geometry layout of its *cells*.

All the areas (i.e. the *GEOM faces*) of a *cell*, or a *lattice*, can be
associated with properties, such as the *material*. In |TOOL|, these areas are
referred to as *regions*.
To enable the visualization of the single *GEOM faces* (i.e. the *regions*)
constituing a *cell*, or a *lattice*, with a specific colorset for a type of
property, |TOOL| relies on instances of the dataclass
:py:class:`Region<glow.geometry_layouts.cells.Region>`. Each instance associates
any of the *GEOM face* of the geometry layout (either technologial or sectorized)
with a color corresponding to the value of the property type to be visualized
in the *SALOME* 3D viewer.

The three main classes that are involved in the definition of the geometry layout
of a *cell* or a *lattice* are :py:class:`Surface<glow.geometry_layouts.geometries.Surface>`,
:py:class:`Cell<glow.geometry_layouts.cells.Cell>` and
:py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>`.
A :py:class:`Surface<glow.geometry_layouts.geometries.Surface>` represents a
single *GEOM face*, while :py:class:`Cell<glow.geometry_layouts.cells.Cell>` and
:py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` are identified by
a *GEOM compound*.
*Surfaces* are intended to facilitate the definition of the *GEOM compound*
that defines the technological geometry of the *cell*. Conversely, the sectorized
geometry and the associated *GEOM compound* are defined by acting on the
technological geometry using methods from the
:py:class:`Cell<glow.geometry_layouts.cells.Cell>` and
:py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` classes.
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>`,
:py:class:`Cell<glow.geometry_layouts.cells.Cell>` and
:py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` classes rely on the
same functionalities that apply **transformation** and **visualization**
operations.
Here, transformation refers to rotating and translating the geometric object in
the XYZ space, and visualization refers to displaying the object in the 3D viewer
of *SALOME*. Each of the aforementioned classes provides its own implementation
of these operations.

In the following, the three main classes describing *surfaces*, *cells* and
*lattices* are discussed and their public methods are detailed.

Surfaces Definition
^^^^^^^^^^^^^^^^^^^

|TOOL| comes with classes to quickly build specific *surfaces*, identified by
*GEOM faces*, in *SALOME*.
In the *Object-Oriented Programming* (*OOP*) view, the types of *surface* that
are available in |TOOL| inherit from the same superclass
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>`, which is an
abstract class characterised by both abstract and concrete methods.
The classes to build specific surfaces are the following:

  - class :py:class:`Circle<glow.geometry_layouts.geometries.Circle>` that
    addresses circular surfaces.
  - class :py:class:`Hexagon<glow.geometry_layouts.geometries.Hexagon>` that
    addresses hexagonal surfaces.
  - class :py:class:`Rectangle<glow.geometry_layouts.geometries.Rectangle>`
    that addresses rectangular surfaces.
  - class :py:class:`GenericSurface<glow.geometry_layouts.geometries.GenericSurface>`
    that addresses any 2D surface created in *SALOME*.

Depending on the specific type of surface, the instantiation requires to specify
the centre of the surface and its characteristic dimensions (i.e. radius for a
circle, width and height for a rectangle, the edge length for a hexagon).
Classes are implemented with default values for the characteristic dimensions.
When an object of the :py:class:`Surface<glow.geometry_layouts.geometries.Surface>`
subclasses is instantiated, the *GEOM* objects for the vertices, the edges and
the corresponding face are automatically built. In this way, the surface can
be shown in the *SALOME* 3D viewer right after the initialization by means of the
method :py:meth:`show_face()<glow.geometry_layouts.geometries.Surface.show_face>`.

The following code snippet shows how to instantiate and display a geometric
surface for a hexagonal case.

.. code-block:: python

  from glow.geometry_layouts.geometries import Hexagon

  surface = Hexagon(center=(1.0, 1.0, 0.0), edge_length=2.0)
  surface.show_face()

:numref:`hex-shape` shows the graphical result obtained by running the code
above in a Python script or directly from the Python console of *SALOME*.

.. _hex-shape:
.. figure:: images/hexagon.png
   :alt: Hexagon in SALOME
   :width: 400px
   :align: center

   Hexagon displayed in the *SALOME* 3D viewer.

Transformation operations can be applied by calling the methods for rotating
or translating the surface, which are declared in the base class
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>`, and are accessible
for any of its subclasses.
The method :py:meth:`rotate()<glow.geometry_layouts.geometries.Surface.rotate>`
performs a rotation of the *GEOM* elements of :py:class:`Surface<glow.geometry_layouts.geometries.Surface>`
around the centre of the corresponding *GEOM face* by a given rotation angle, in
degrees. The direction of the rotation follows the standard *right-hand* rule.
The method :py:meth:`translate()<glow.geometry_layouts.geometries.Surface.translate>`
moves the *GEOM* elements of :py:class:`Surface<glow.geometry_layouts.geometries.Surface>`
so that the center of the corresponding *GEOM face* coincides with the given
XYZ coordinates.
For the hexagonal surface declared above, the code instructions are the
following:

.. code-block:: python

  surface.rotate(90)
  surface.translate((0.0, 0.0, 0.0))
  surface.show_face()

By applying these methods, the resulting *GEOM face* is shown in :numref:`hex-transf`.

.. _hex-transf:
.. figure:: images/hexagon_rot_transl.png
   :alt: Hexagon rotated and translated in SALOME
   :width: 400px
   :align: center

   Hexagon after applying rotation and traslation operations, as shown in the
   *SALOME* 3D viewer.

The *GEOM face* that is specific of the subclass of
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>` can be directly
modified within *SALOME* and the modified *GEOM face* applied to the
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>` object by calling
the method :py:meth:`update_from_face()<glow.geometry_layouts.geometries.Surface.update_from_face>`.
The implementation of this method is specific for each of the subclasses of
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>`. In general, the
method receives as parameter a *GEOM face* and updates the instance attributes
of :py:class:`Surface<glow.geometry_layouts.geometries.Surface>` accordingly.
A check ensures that only *GEOM faces* are provided, and that the given *GEOM
face* corresponds to the characteristic surface the
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>` class refers to.

Cell Definition
^^^^^^^^^^^^^^^

|TOOL| comes with classes to build cells having either a hexagonal or a
rectangular characteristic *surface*.
The module :py:mod:`glow.geometry_layouts.cells` provides the abstract base class
:py:class:`Cell<glow.geometry_layouts.cells.Cell>`, which represents a cell
described in terms of an instance from the subclasses of the
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>` base class.
The subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>` are the
following ones:

  - class :py:class:`RectCell<glow.geometry_layouts.cells.RectCell>` that
    represents rectangular cells.
  - class :py:class:`HexCell<glow.geometry_layouts.cells.HexCell>` that
    represents hexagonal cells.
  - class :py:class:`GenericCell<glow.geometry_layouts.cells.GenericCell>`
    that represents cells characterised by any *GEOM face* or *GEOM compound*
    created in *SALOME*.

When instantiating any of the aforementioned sublasses, the corresponding instance
of the :py:class:`Surface<glow.geometry_layouts.geometries.Surface>` subclasses
is built. Depending on the cell type, the instantiation requires either the
characteristic dimensions of the *surface* (e.g. width and height for a rectangle
or edge length for a hexagon) or the *GEOM face*, or *GEOM compound*, directly
(:py:class:`GenericCell<glow.geometry_layouts.cells.GenericCell>` case).

The following code snippet shows how to instantiate the different type of cells
available in |TOOL|.

.. code-block:: python

  from glow.geometry_layouts.cells import GenericCell, HexCell, RectCell

  hex_cell = HexCell(
      center=(0.0, 0.0, 0.0),
      edge_length=1.0,
      name='HexCell')

  rect_cell = RectCell(
      center=(0.0, 0.0, 0.0),
      height_x_width=(1.0, 2.0),
      rounded_corners=[(1, 0.1), (3, 0.1)],
      name='RectCell')

  gnrc_cell = GenericCell(shape=surface)

For a rectangular cell, the ``rounded_corners`` parameter indicates the index
of the corner of the rectangle and the associated curvature radius to generate
a rectangle with rounded corners.
For a :py:class:`GenericCell<glow.geometry_layouts.cells.GenericCell>`, the
instantiation requires to specify a *GEOM face* or a *GEOM compound* representing
its geometric *surface*.

The class :py:class:`Cell<glow.geometry_layouts.cells.Cell>` declares attributes
and methods common to all its subclasses. Public methods cover the following
functionalities:

  - displaying the cell's geometry layout in the *SALOME* 3D viewer;
  - adding and removing circles within the cell's boundaries;
  - applying transformation operations for rotating and translating the cell's
    characteristic *GEOM* elements;
  - applying a sectorization operation of the cell's geometry layout;
  - setting up values for the available property types associated to one or all
    the *regions* of the cell's technological geometry;
  - inspecting the information (name and properties value) related to a specific
    *region* of the cell that has been selected in the *SALOME* 3D viewer;
  - updating the cell's geometry layout with a *GEOM face* or a *GEOM compound*;
  - restoring the cell to its original state in terms of both geometry and the
    properties associated with its *regions*.

In the following, all these methods are detailed.

.. _cell-show:

Displaying the Cell's Geometry Layout
"""""""""""""""""""""""""""""""""""""

The cell's geometry layout can be displayed in the *SALOME* 3D viewer by calling
the method :py:meth:`show()<glow.geometry_layouts.cells.Cell.show>`.
The method has two parameters, each associated with a default value:

  - ``property_type_to_show``, an item of the enumeration :py:class:`PropertyType<glow.support.types.PropertyType>`,
    it identifies the property type (e.g. the *material*) according to which
    the cell's *regions* (i.e. the *GEOM faces*) are displayed with a color.
    Each *region* has a colour associated with the value of the indicated
    property type. If no property type is provided, the *regions* are displayed
    with a default colour.
  - ``geometry_type_to_show``, an item of the enumeration :py:class:`GeometryType<glow.support.types.GeometryType>`,
    it identifies the type of geometry to show, i.e. either the technological
    or the sectorized one. The cell's *regions*, identified by a list of objects
    of the dataclass :py:class:`Region<glow.geometry_layouts.cells.Region>`,
    are build from the *GEOM compound* associated with the technological or
    sectorized layout. By default, the method displays the *regions* of the
    technological geometry.

Users should note that the method :py:meth:`show()<glow.geometry_layouts.cells.Cell.show>`
will raise an exception if they request to display the *regions* according to
a property type for which a *region* has no corresponding value.

The following code snippet shows how to display the *regions* of the cell's
technological geometry (indicated by the :py:attr:`TECHNOLOGICAL<glow.support.types.GeometryType.TECHNOLOGICAL>`
type of geometry) with a colorset in terms of the property type
:py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`.

.. code-block:: python

  hex_cell.show(
      property_type_to_show=PropertyType.MATERIAL,
      geometry_type_to_show=GeometryType.TECHNOLOGICAL
  )

*Regions* are added to the *Object Browser* in *SALOME* as children of the cell
they belong to. If not displayed automatically (it can happen when running a
new *SALOME* instance with a script), they can be shown by selecting the
"*Show Only Children*" item in the contextual menu for the cell (see
:numref:`show-children`).

.. _show-children:
.. figure:: images/cell_show_children.png
   :alt: How to display the cell's regions in SALOME
   :width: 400px
   :align: center

   How to display the *regions* associated to a cell in *SALOME*.

The geometry layout resulting from the aforementioned code is shown in
:numref:`cell-mat`.

.. _cell-mat:
.. figure:: images/cell_show_col.png
   :alt: Cell's technological geometry with MATERIAL colorset
   :width: 400px
   :align: center

   Hexagonal cell's technological geometry with the :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`
   colorset.

Circles Addition and Removal
""""""""""""""""""""""""""""

Typically, fuel pin cells, having either a cartesian or a hexagonal geometry,
are characterised by several concentric circles to represent the different
regions of a cell, each having its own properties.
In general, circles can be placed either in the cell's centre or in any other
point within its boundaries.

In |TOOL|, the method :py:meth:`add_circle()<glow.geometry_layouts.cells.Cell.add_circle>`
allows to position a circle, with a specified radius, inside the cell. The
addition is performed only if the circle's radius does not exceeds the
characteristic dimensions (e.g. the apothem for a hexagon) of the *surface* (
:py:class:`Surface<glow.geometry_layouts.geometries.Surface>` subclasses) the
cell is based on.
Given the circle's radius, a *GEOM face* object is built in the given position,
if any is specified, otherwise the circle is added in the cell's centre.
In any case, a *partition* operation between the *GEOM compound* representing
the current technological geometry of the cell and the *GEOM face* of the new
circle is performed, resulting in a *GEOM compound* that comprises both.

The following code snippet shows how to add circles in specific positions within
a hexagonal cell.

.. code-block:: python

  hex_cell.add_circle(radius=0.5)
  hex_cell.add_circle(radius=0.1, position=(0.2, 0.2, 0.0))
  hex_cell.show()

:numref:`cell-circles` shows the result of adding two circles, the first in the
cell's centre, the second in a specific position. The resulting updated
technological geometry is shown in the *SALOME* 3D viewer after calling the
method :py:meth:`show()<glow.geometry_layouts.cells.Cell.show>`.

.. _cell-circles:
.. figure:: images/cell_add_circle.png
   :alt: Hexagonal cell with two circular regions in SALOME
   :width: 400px
   :align: center

   Hexagonal cell's geometry layout after adding two circles to its
   technological geometry.

Calling the method :py:meth:`add_circle()<glow.geometry_layouts.cells.Cell.add_circle>`
updates the technological geometry of the cell. The same goes for the method
:py:meth:`remove_circle()<glow.geometry_layouts.cells.Cell.remove_circle>`.

When any property type (e.g. a *material*) has been assigned to the cell's *region*
where the circle is added, the *regions* resulting from partitioning the cell
with the circle inherit the properties of the overlapped *regions* (see
:numref:`prop-regions`).

.. _prop-regions:
.. figure:: images/cell_prop_regions.png
   :alt: Hexagonal cell with property colorset in SALOME
   :width: 400px
   :align: center

   Hexagonal cell's technological geometry shown with a properties colorset;
   the new circular *regions* have the same property type value as the *region*
   they overlap.

If the added circle is cell-centred, then it also inherits the sectorization
options of the overlapped centred *region* (see :numref:`sect-regions`).

.. _sect-regions:
.. figure:: images/cell_sect_regions.png
   :alt: Hexagonal cell with sectorization visualization in SALOME
   :width: 400px
   :align: center

   Hexagonal cell's sectorized geometry; only the cell-centred circle is
   subdivided in six regions as the the overlapped *region*.

When removing a circular *region* having any property type or sectorization option
associated, the *region* resulting from its removal keeps the same values of the
*region* that surrounded the removed circular *region*.

Transformation Operations
"""""""""""""""""""""""""

Transformation operations can be applied by calling the methods for rotating
or translating the cell's geometric elements, i.e. the *GEOM compounds*
representing the cell's technological and sectorized layouts, as well as the
:py:class:`Region<glow.geometry_layouts.cells.Region>` objects corresponding
to the layout currently displayed.
The method :py:meth:`rotate()<glow.geometry_layouts.cells.Cell.rotate>`
requires the rotation angle, in degrees, and assumes that the rotation is
performed around the Z-axis. The direction of the rotation follows the standard
*right-hand* rule.
The method :py:meth:`translate()<glow.geometry_layouts.cells.Cell.translate>`
needs the XYZ coordinates of the new centre of the cell.
While the former operates on the same instance, the latter returns a deep copy
of the original instance positioned in the new centre.
For a hexagonal cell, the code instructions for rotating and translating the
cell are the following:

.. code-block:: python

  hex_cell.rotate(90)
  new_cell = hex_cell.translate((1.0, 1.0, 0.0))
  new_cell.show()

Sectorization Operation
"""""""""""""""""""""""

Other than the technological geometry, cells can be displayed also in terms of
the sectorized one.
This type of geometry consists in subdividing the cell's *regions* of the
technological geometry in a number of angular regions (the *sectors*) which is
specific for the type of cell. Subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>` declares
the available number of sectors a *region* of the technological geometry can have,
as well as the starting angle from which the subdivision starts.
We can have the following values:

  - :py:class:`HexCell<glow.geometry_layouts.cells.HexCell>` - admitted number
    of sectors are either `1` or `6`, while `0` or `30` for the starting angle.
  - :py:class:`RectCell<glow.geometry_layouts.cells.RectCell>` - admitted number
    of sectors are `1`, `4`, `8` and `16`, while the corresponding angles are
    `0` and `45.0` for a subdivision in four sectors, `0` and `22.5` for a
    subdivision in eight sectors, `0` for a subdivision in one or sixteen
    sectors.

Rectangular cells also have the option of applying a *windmill* sectorization
to the region farthest from the cell's center, provided that the *region* is
subdivided into eight sectors.

Each of the subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>`
provide their own configuration for applying the sectorization. In particular,
for a :py:class:`RectCell<glow.geometry_layouts.cells.RectCell>` the ``windmill``
parameter can be provided to apply a *windmill* sectorization, while for
:py:class:`HexCell<glow.geometry_layouts.cells.HexCell>` and
:py:class:`GenericCell<glow.geometry_layouts.cells.GenericCell>` this parameter
is absent. In any case, the logic for subdividing the *regions* in sectors is
common to all subclasses.

The following code snippet shows how to apply a sectorization, with ``windmill``
option enabled, for a cartesian cell having two cell-centred circles.

.. code-block:: python

  rect_cell.sectorize([1, 4, 8], [0, 45, 22.5], windmill=True)
  rect_cell.show(geometry_type_to_show=GeometryType.SECTORIZED)

Elements in the two lists provided to the method
:py:meth:`sectorize()<glow.geometry_layouts.cells.RectCell.sectorize>` are
associated to the *regions* from the closest to the farthest one from the cell's
centre.
:numref:`cart-cell-sect` shows the result after applying the indicated sectorization.

.. _cart-cell-sect:
.. figure:: images/cell_sectorize.png
   :alt: Cartesian cell after its sectorization
   :width: 400px
   :align: center

   Cartesian cell after applying the sectorization operation. The number of
   subdivisions of the cell's *regions* matches the order in which sectorization
   numbers are provided to the method.

.. _set-cell-prop:

Setting Up the Cell's Regions Properties
""""""""""""""""""""""""""""""""""""""""

Cells' *regions* can be displayed by applying a colorset that depends on the type
of property to show, as item of the :py:class:`PropertyType<glow.support.types.PropertyType>`
enumeration. An example of property type is the *material* constituing each
*region*, identified by the item :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`.
To set values for a specific property type, users can rely on two methods:

  - :py:meth:`set_properties()<glow.geometry_layouts.cells.Cell.set_properties>`,
    which allows users to set values for different types of properties for *all*
    the regions of the cell's technological geometry.
    The convention for declaring the values of a property is from the closest
    to the farthest *region* with respect to the cell's centre.
  - :py:meth:`set_region_property()<glow.geometry_layouts.cells.Cell.set_region_property>`,
    which allows to set a value for the indicated type of property of a *single*
    region of a cell; this can be either the *GEOM face* currently selected in
    the *SALOME* 3D viewer or the one provided as parameter to the method.

The following code snippet shows how to apply values for the
:py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>` type of property,
which is the only one currently implemented.

.. code-block:: python

  rect_cell.set_properties(
      {PropertyType.MATERIAL: ['GAP', 'FUEL', 'COOLANT']}
  )
  rect_cell.add_circle(0.1)
  rect_cell.set_region_property(
      PropertyType.MATERIAL,
      'MAT',
      Circle(radius=0.1).face
  )
  rect_cell.show(PropertyType.MATERIAL)

In particular, given a cartesian cell with two cell-centred circles, the first
method enables all the material values to be set at the same time.
A new circular *region* is added, and the corresponding *GEOM face* is used to
identify the *region* within the cell to which the property should be assigned.
From within the *SALOME* 3D viewer, the *region* can be provided by simply
selecting the corresponding *GEOM face* and calling the method from the
integrated Python console.
In any case, the cell's geometry layout with the :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`
colorset is shown in :numref:`cell-after-props`.

.. _cell-after-props:
.. figure:: images/cell_properties.png
   :alt: Cartesian cell after setting up the properties
   :width: 400px
   :align: center

   Cartesian cell after setting up values for the :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`
   property type for each region. It is shown with a colorset highlighting the
   different values assigned to the cell's *regions*.

Inspection of Regions
"""""""""""""""""""""

When *regions* of a cell are displayed in the *SALOME* 3D viewer, users can
obtain information about an individual *region*, including its assigned
properties. This is done by calling the method :py:meth:`get_regions_info()<glow.geometry_layouts.cells.Cell.get_regions_info>`
directly in the Python console of *SALOME* from an object
of any of the subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>`.
If no *region* (as *GEOM* face), or more than one, is selected when calling the
method, an exception is raised.
The available information, which is printed in the Python console, includes the
name of the cell's *region* and the value for each of the assigned type of
properties (see :numref:`reg-info`).

.. _reg-info:
.. figure:: images/region_info.png
   :alt: Information about a selected region of the cell
   :width: 400px
   :align: center

   Information about a selected *region* of the cell; its name and values for
   the assigned properties are printed.

Updating the Cell's Geometry Layout
"""""""""""""""""""""""""""""""""""

The methods of the class :py:class:`Cell<glow.geometry_layouts.cells.Cell>`
enable the cell's technological and sectorized geometries to be customized
by means of circles and lines, where the latter must follow the rules tied to
the sectorization operation (i.e. lines subdivides *regions* of the technological
geometry in fixed numbers of angular sectors).
To support any customization of the cell's geometry layout, while keeping the
base *surface* (subclass of :py:class:`Surface<glow.geometry_layouts.geometries.Surface>`)
the same, two methods are provided:

  - :py:meth:`update_geometry()<glow.geometry_layouts.cells.Cell.update_geometry>`,
    which enables to update the *GEOM compound*, representing either the
    technological or the sectorized geometry, that is displayed in the *SALOME*
    3D viewer with the *GEOM face* or *GEOM compound* currently selected.
  - :py:meth:`update_geometry_from_face()<glow.geometry_layouts.cells.Cell.update_geometry_from_face>`,
    which enables to update the *GEOM compound* corresponding to the indicated
    :py:class:`GeometryType<glow.support.types.GeometryType>` with the given
    *GEOM face* or *GEOM compound*.

In both cases, the result is a new layout for the technological or the sectorized
geometry where the new *regions* inherit the already assigned properties, if
any; the same goes for the sectorization options.

The following code snippet shows how the cell's technological geometry could
be updated with a non-standard geometry built by overlapping two hexagonal
*surfaces* with different dimensions.

.. code-block:: python

  hex_1 = Hexagon(edge_length=1)
  hex_2 = Hexagon(edge_length=1.5)

  shape = make_partition([hex_2.face], [hex_1.face], ShapeType.COMPOUND)

  hex_cell = HexCell()
  hex_cell.update_geometry_from_face(GeometryType.TECHNOLOGICAL, shape)
  hex_cell.show()

The function :py:func:`make_partition()<glow.interface.geom_interface.make_partition>`
cuts a list of *GEOM faces* (in the first argument) with those provided in the
list as second argument; the resulting type of shape is indicated as third argument.
After applying the built geometry to the cell, the result can be displayed in
the *SALOME* 3D viewer (see :numref:`updated-cell`).

.. _updated-cell:
.. figure:: images/updated_cell.png
   :alt: Cell's geometry after update
   :width: 400px
   :align: center

   Hexagonal cell's layout after updating its technological geometry.

Restoring Cell's State
""""""""""""""""""""""

There could be cases where users need to reset the cell's geometry layout and
the properties associated to its regions (see :ref:`tutorial-overlap`).
The method :py:meth:`restore()<glow.geometry_layouts.cells.Cell.restore>`
satisfies this need by restoring the *GEOM compound* of the cell's technological
layout to its base *surface* (e.g. a *GEOM face* identifying a rectangle) without
any inner circle.
The sectorized layout, as well as the properties and sectorization options, are
completely removed.

Lattice Definition
^^^^^^^^^^^^^^^^^^

|TOOL| comes with classes to build lattices characterised by either hexagonal
or cartesian cells.
The module :py:mod:`glow.geometry_layouts.lattices` provides the class
:py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` to describe any
kind of lattice of cells.
The type of lattice is determined by the type of the cells, either cartesian or
hexagonal. All the cells in the lattice must be of the same type, identified by
an item of the enumeration :py:class:`CellType<glow.support.types.CellType>`.
This is automatically set at instantiation time or when adding cells to the
lattice.

The :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` class can be
instantiated either without any cell or by providing a list of objects of the
subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>`.

In |TOOL|, the logic behind the construction of a lattice relies on the *layer*
concept: when a new cell, or a group of cells is added to the lattice, the cells
are associated to a layer (either a new layer or an existing one already
containing some cells). The layer to which the cells are added depends on the
specific method used to add them.
By adopting this logic, |TOOL| can easily handle the construction of the *GEOM
compound* that identifies the lattice geometry layout, especially in the case
of lattices made by superimposing cells with different dimensions.

The following code snippet shows how to instantiate a lattice with the cartesian
or hexagonal cells available in |TOOL|.

.. code-block:: python

  from glow.geometry_layouts.cells import HexCell, RectCell
  from glow.geometry_layouts.lattices import Lattice

  hex_cell = HexCell()
  rect_cell = RectCell()

  # Lattice instantiation by providing all the cartesian cells at once
  cart_lattice = Lattice(
      cells=[
          rect_cell.translate((0.5, 0.5, 0.0)),
          rect_cell.translate((-0.5, 0.5, 0.0)),
          rect_cell.translate((-0.5, -0.5, 0.0)),
          rect_cell.translate((0.5, -0.5, 0.0)),
      ],
      name="Cartesian Lattice",
      center=(0.0, 0.0, 0.0),
      boxes_thick=[0.075, 0.075]
  )
  # Lattice instantiation without any cell
  lattice = Lattice()
  # Lattice instantiation with a hexagonal central cell
  hex_lattice = Lattice([hex_cell])

The three examples show different instantiations; in particular, we have:

  - a cartesian lattice built from a list of cells positioned to recreate a
    2x2 pattern; by specifying the ``boxes_thick`` parameter, the built lattice
    is enclosed within a rectangular box made by two layers of given thicknesses.
  - a lattice built without any cell. The lattice's methods for adding cells
    need to be called to define its geometry layout (see :ref:`add-cells`).
  - a hexagonal lattice built from a single cell positioned in the centre of
    the lattice.

Similarly to the cells, the two types of geometry layout, the technological and
the sectorized ones, apply to the lattice (see ::ref:`geom-def`).

The :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` public methods
cover the following functionalities:

  - building the lattice's *regions*, as elements of the dataclass
    :py:class:`Region<glow.geometry_layouts.cells.Region>`, according to either
    the technological or the sectorized type of geometry of the cells in the
    lattice;
  - displaying the lattice's geometry layout in the *SALOME* 3D viewer;
  - adding a single cell or a group of the same cell organised in one or more
    rings around the lattice's centre;
  - transformation operations for rotating or translating the lattice's cells;
  - enclosing the lattice in a box declared from the thicknesses of its layers
    or by means of an instance of the subclasses of
    :py:class:`Cell<glow.geometry_layouts.cells.Cell>`;
  - setting up the properties associated to one *region* of the lattice or to
    the ones of the box;
  - applying a specific type of symmetry in accordance with the type of lattice;
  - setting the type of geometry in accordance with the type of lattice and of
    applied symmetry;
  - inspecting the information related to a specific *region* of the lattice
    that has been selected in the *SALOME* 3D viewer;
  - restoring a list of cells of the lattice to their original state, both in
    terms of geometry and properties.

Building Lattice's Regions
""""""""""""""""""""""""""

To facilitate displaying and exporting the lattice's geometry layout, the method
:py:meth:`build_regions()<glow.geometry_layouts.lattices.Lattice.build_regions>`
is provided. It builds a list of :py:class:`Region<glow.geometry_layouts.cells.Region>`
objects that are representative of the *regions* in which the lattice is subdivided
when assembling all the cells together with the box, if present.
Cells can be associated with different layers of cells in the lattice. When the
lattice's regions are built to be displayed in the *SALOME* 3D viewer, a process
is carried out. This can be imagined as if all the layers were collapsed into a
single layer of cells. The layers are traversed from top to bottom, and any cells
that are found to be overlapped by those of a higher layer are either cut or
removed from the lattice.
:numref:`overlap` shows the result of overlapping a cell with others.

.. _overlap:
.. figure:: images/lattice_overlap_cells.png
   :alt: Lattice with a cell overlapping other cells
   :width: 400px
   :align: center

   Hexagonal lattice where a cell overlaps other cells of an inferior layer.

If any symmetry is applied or the lattice is enclosed in a box, the *GEOM compound*
of the assembled cells is either cut to extract the portion that replicates the
symmetry or assembled with the geometry layout of the box.
Given the final *GEOM compound*, the contained *GEOM faces* are extracted and
a :py:class:`Region<glow.geometry_layouts.cells.Region>` object is built for
each one.
In any case, the property assignment involves identifying the corresponding
*region* among the ones of the technological geometry of the lattice's cells.

According to the type of geometry of the cells that is provided to the method
:py:meth:`build_regions()<glow.geometry_layouts.lattices.Lattice.build_regions>`,
the resulting regions describe either the technological or the sectorized
geometry of the lattice.

Displaying the Lattice's Geometry Layout
""""""""""""""""""""""""""""""""""""""""

The lattice's geometry layout can be displayed in the *SALOME* 3D viewer by
calling the method :py:meth:`show()<glow.geometry_layouts.lattices.Lattice.show>`.
Depending on its parameters, it builds and displays the corresponding *regions*
(i.e. the *GEOM faces*) of the lattice.

Regions are built and shown according to either the technological or the
sectorized geometry by specifying it as parameter of the method.
The same considerations on the parameters done for the method
:py:meth:`show()<glow.geometry_layouts.cells.Cell.show>` of the subclasses of
:py:class:`Cell<glow.geometry_layouts.cells.Cell>` are valid for the lattice
as well (see :ref:`cell-show`).
It is important to note that when displaying the lattice's *regions* with a
colorset according to the indicated :py:class:`PropertyType<glow.support.types.PropertyType>`,
regions with the same property type value are coloured the same.

In *SALOME*, regions are added to the *Object Browser* as children of the
lattice they belong to, similarly to what happens for cells (see
:numref:`show-children`).

The following code snippet shows how to display the regions of the lattice's
technological geometry (indicated by the :py:attr:`TECHNOLOGICAL<glow.support.types.GeometryType.TECHNOLOGICAL>`
type of geometry) with a colorset in terms of the property type
:py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`.

.. code-block:: python

  cart_lattice.show(
      property_type_to_show=PropertyType.MATERIAL,
      geometry_type_to_show=GeometryType.TECHNOLOGICAL
  )

:numref:`lattice-show` shows the resulting geometry layout of the lattice after
running the above code.

.. _lattice-show:
.. figure:: images/lattice_show_col.png
   :alt: Lattice's technological geometry with the MATERIAL colorset
   :width: 400px
   :align: center

   Cartesian lattice's technological geometry with the :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`
   colorset.

.. _add-cells:

Adding cell(s)
""""""""""""""

A lattice can be built by instantianting a :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>`
object, providing a list of :py:class:`Cell<glow.geometry_layouts.cells.Cell>`
subclasses. In addition to this approach, it is often useful to contruct a
lattice by adding a cell or a ring of cells with simple methods. For this reason,
the following methods have been introduced:

  - :py:meth:`add_cell()<glow.geometry_layouts.lattices.Lattice.add_cell>`,
    which allows to add a single cell at an indicated position;
  - :py:meth:`add_ring_of_cells()<glow.geometry_layouts.lattices.Lattice.add_ring_of_cells>`,
    which allows to add a ring of the same cell at the indicated ring index;
  - :py:meth:`add_rings_of_cells()<glow.geometry_layouts.lattices.Lattice.add_rings_of_cells>`,
    which allows to add the indicated number of rings of the same cell, starting
    from the current ring index occupied by cells.

The method :py:meth:`add_cell()<glow.geometry_layouts.lattices.Lattice.add_cell>`
adds the cell to the specified position, if any is provided, otherwise the cell
is placed at the position indicated by the cell's centre. It is important to
note that any cell added with this method is included in a new *layer*, i.e. a
new sub-list is created for the attribute :py:attr:`layers<glow.geometry_layouts.lattices.Lattice.layers>`
containing the cell itself.

The layout of a lattice can be considered as consisting of several rings, each
occupied by an increasing number of cells as the ring index increases. The two
methods :py:meth:`add_ring_of_cells()<glow.geometry_layouts.lattices.Lattice.add_ring_of_cells>`
and :py:meth:`add_rings_of_cells()<glow.geometry_layouts.lattices.Lattice.add_rings_of_cells>`
provide a quick way for adding one or more rings of cells. The former adds the
cells at the given ring index while the latter adds the indicated number of
rings of cells starting from the maximum value of ring index currently present
in the lattice.
Users should also note that, while the former method enables them to specify
the *layer* to which the ring of cells is added (by providing its index), the
latter always adds the rings of cells to a new *layer*.

All the aforementioned methods do not allow to mix cells with different types
(i.e. having different item of the enumeration :py:class:`CellType<glow.support.types.CellType>`);
this ensures that all the cells have either a cartesian or a hexagonal type.

The following code snippet shows the different ways to add cells to a lattice.

.. code-block:: python

  cell = HexCell()
  lattice = Lattice([cell])

  lattice.add_ring_of_cells(cell, 1)
  lattice.add_rings_of_cells(cell, 2)
  lattice.add_cell(cell, (1.5, 1.5, 0.0))
  lattice.show()

The lattice's geometry layout resulting from adding hexagonal cells using the
three methods is shown in :numref:`lattice-add`.

.. _lattice-add:
.. figure:: images/lattice_add_cells.png
   :alt: Lattice after adding cells
   :width: 400px
   :align: center

   Hexagonal lattice built by applying the three methods for adding cells.

Lattice's Transformation Operations
"""""""""""""""""""""""""""""""""""

Transformation operations can be applied by calling the methods for rotating
or translating the lattice's geometric elements, i.e. the *GEOM compound* objects
representing its full and partial (if any symmetry is applied) geometry layouts,
the contained cells, including the box (if present), and all the *regions*.
The method :py:meth:`rotate()<glow.geometry_layouts.lattices.Lattice.rotate>`
requires the rotation angle, in degrees, and assumes that the rotation is
performed around the Z-axis. The direction of the rotation follows the standard
*right-hand* rule.
The method :py:meth:`translate()<glow.geometry_layouts.lattices.Lattice.translate>`
needs the new XYZ coordinates of the centre of the lattice.
Users should note that both methods operate on the same instance and the result
of the transformation is directly shown in the *SALOME* 3D viewer.

Enclosing the Lattice in a Box
""""""""""""""""""""""""""""""

In nuclear reactors, fuel assemblies are typically framed in a metallic container.
To replicate exactly the same kind of layouts, |TOOL| allows to insert a lattice
within a box.
A box is an instance of the subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>`
which can be built either from the thickness of its layers or by instantiating
the corresponding :py:class:`Cell<glow.geometry_layouts.cells.Cell>` object
directly.
The former case relies on the method :py:meth:`build_lattice_box()<glow.geometry_layouts.lattices.Lattice.build_lattice_box>`,
which, given the type of lattice (i.e. hexagonal or cartesian), automatically
instantiates a :py:class:`Cell<glow.geometry_layouts.cells.Cell>` object built
by overlapping as many rectangles or hexagons as the number of the indicated
thicknesses of the layers.
If all the values provided to the :py:meth:`build_lattice_box()<glow.geometry_layouts.lattices.Lattice.build_lattice_box>`
method are positive (independently from the value), the borders of the layer
closest to the centre of the lattice touch the outermost ring of cells without
overlapping it (see :numref:`box-pos`).
The method also allows the first thickness value in the list to be negative,
which handles a situation where the layer closest to the centre cuts the
farthest ring of cells (see :numref:`box-neg`).

The following code snippet shows how to build a box for the lattice using the
method :py:meth:`build_lattice_box()<glow.geometry_layouts.lattices.Lattice.build_lattice_box>`
with the thickness of the first layer either being positive or negative.

.. code-block:: python

  lattice.build_lattice_box([0.1, 0.1])
  lattice.show()

  lattice.build_lattice_box([-0.1, 0.1])
  lattice.show()

The result of applying both method calls separately, for a hexagonal lattice,
is shown in :numref:`box-pos` and in :numref:`box-neg` respectively.

.. _box-pos:
.. figure:: images/lattice_box_pos.png
   :alt: Lattice within a box with positive thicknesses
   :width: 400px
   :align: center

   Hexagonal lattice framed in a box with all positive thicknesses for the
   layers.

.. _box-neg:
.. figure:: images/lattice_box_neg.png
   :alt: Lattice within a box with negative first thickness
   :width: 400px
   :align: center

   Hexagonal lattice framed in a box with a negative thickness for the first
   layer. The box cuts the farthest ring of cells.

The lattice's box can also be declared by setting the corresponding property
:py:attr:`lattice_box<glow.geometry_layouts.lattices.Lattice.lattice_box>` with
an object of the subclasses of :py:class:`Cell<glow.geometry_layouts.cells.Cell>`.
The setter of the property requires the cell's centre to coincide with that of
the lattice, otherwise an exception is raised.
Both :py:class:`Cell<glow.geometry_layouts.cells.Cell>` objects or ``None`` are
valid inputs for the setter. The latter can be used to remove any box previously
set.

Both approaches to setting a box lead to the same result: the *GEOM compound*
representing the geometry layout of the lattice is updated by assembling the
*GEOM compound* of each cell with that of the box, which can potentially cut the
*GEOM compound* of the cells of the farthest ring.

Setting Up Properties
"""""""""""""""""""""

Just like for cells, the *regions* of a lattice can be displayed with a colorset
according to the type of property to display, as item of the
:py:class:`PropertyType<glow.support.types.PropertyType>` enumeration.

There are different ways for users to set values for a specific property type
of a *region* of the lattice.
If the *region* belongs to any cell, the methods previously described (see
:ref:`set-cell-prop`) for a :py:class:`Cell<glow.geometry_layouts.cells.Cell>`
object remain valid, provided they are applied to the correct instance stored
in the attribute :py:attr:`layers<glow.geometry_layouts.lattices.Lattice.layers>`.

In addition, users can rely on the following methods of the class
:py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>`:

  - :py:meth:`set_region_property()<glow.geometry_layouts.lattices.Lattice.set_region_property>`,
    which allows to set a value for the indicated type of property of a single
    lattice's *region* (i.e. a *GEOM face*); this can be either the *GEOM face*
    currently selected in the *SALOME* 3D viewer or the one provided as parameter
    to the method.
  - :py:meth:`set_lattice_box_properties()<glow.geometry_layouts.lattices.Lattice.set_lattice_box_properties>`,
    which allows users to set values for different types of properties for all
    the regions of the :py:class:`Cell<glow.geometry_layouts.cells.Cell>`
    instance, which is the box that encloses the lattice.
    The convention for declaring the values of a property is always the same,
    i.e. from the *region* closest to the center to the farthest *region*.
    Users should note that for hexagonal boxes, the number of values to provide
    is always equal to that of the layers plus one. The reason is that the first
    value in the list is associated with the *regions* between the cells and the
    first layer of the box. These regions all share the same property type value.

The following code snippet shows the different ways to apply values for the
:py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>` property type,
i.e. either to all the cells or to an indicated *region* or to the regions of
the lattice's box.

.. code-block:: python

  # Build the lattice geometry layout
  hex_cell = HexCell()
  hex_cell.rotate(90)
  lattice = Lattice([hex_cell])
  lattice.add_ring_of_cells(hex_cell, 1)
  lattice.build_lattice_box([0.1])
  # The same value for the 'MATERIAL' property is assigned to all the cells
  for layer in lattice.layers:
      for cell in layer:
          cell.set_properties(
              {PropertyType.MATERIAL: ['COOLANT']}
          )
  # A different value for the 'MATERIAL' property is assigned to the central
  # cell
  lattice.set_region_property(PropertyType.MATERIAL, 'GAP', hex_cell.face)
  # Values for the 'MATERIAL' property are assigned to the box's regions
  lattice.set_lattice_box_properties(
      {PropertyType.MATERIAL: ['COOLANT', 'METAL']}
  )
  lattice.show(PropertyType.MATERIAL)

The resulting lattice's geometry layout with the :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`
colorset is shown in :numref:`lattice-set-props`.

.. _lattice-set-props:
.. figure:: images/lattice_properties.png
   :alt: Lattice after setting up the properties
   :width: 400px
   :align: center

   Lattice after setting up the values for a type of property. It is shown
   with the corresponding colorset.

Applying Symmetries
"""""""""""""""""""

Solving the Boltzmann transport equation on the full geometry layout of a fuel
assembly can be computationally expensive, in particular if the geometry contains
many rings of cells.
To speed up the calculations, users can rely on cuts to extract parts out of
the existing layout, thereby isolating the minimum portion of the geometry
required to describe the entire pattern.
|TOOL| supports the application of specific types of symmetries to the lattice.
According to the type of cells in the lattice, we can have:

  - Half, quarter, and eighth symmetries for a cartesian lattice. Half and quarter
    symmetries cut out the corresponding rectangular portion of the lattice,
    while the eighth symmetry cuts out a right triangular portion with a centre
    angle of 45°.
  - Third, sixth and twelfth symmetries for a hexagonal lattice framed in a box.
    The third symmetry cuts out a parallelogram of the lattice, the sixth symmetry
    a regular triangle and the twelfth a right triangle with a centre angle of
    30°.

The method :py:meth:`apply_symmetry()<glow.geometry_layouts.lattices.Lattice.apply_symmetry>`
allows users to apply the indicated type of symmetry as item of the enumeration
:py:class:`SymmetryType<glow.support.types.SymmetryType>`.
Since |TOOL| considers that only specific types of symmetry are allowed for
each type of lattice, an exception is raised if the user tries to apply an
invalid symmetry.
Independently from the type of symmetry, the method
:py:meth:`apply_symmetry()<glow.geometry_layouts.lattices.Lattice.apply_symmetry>`
automatically performs *cut* operations on the *GEOM compound* of the lattice
so that the remaining part describes the requested symmetry.

For cartesian lattices, the operation of applying a symmetry is performed
independently of the presence of a box. However, for hexagonal lattices, |TOOL|
requires the lattice to be framed in a box. This is because the *SALT* module
of *DRAGON5* cannot track the resulting geometry layout if the shape is not
triangular or quadrilateral.

The following code snippet shows different applications of a symmetry type
for a cartesian and a hexagonal lattice.

.. code-block:: python

  rect_lattice.apply_symmetry(SymmetryType.QUARTER)
  hex_lattice.apply_symmetry(SymmetryType.TWELFTH)

When calling the method :py:meth:`apply_symmetry()<glow.geometry_layouts.lattices.Lattice.apply_symmetry>`,
the geometry layout of the lattice is automatically updated and displayed in
the *SALOME* 3D viewer (if the method is called from its Python console).
If the :py:attr:`FULL<glow.support.types.SymmetryType.FULL>` is provided to the
method, any previously applied symmetry is removed and the entire geometry layout
of the lattice is displayed.

:numref:`quarter-symm` and :numref:`twelfth-symm` show the results of applying
a :py:attr:`QUARTER<glow.support.types.SymmetryType.QUARTER>` and a
:py:attr:`TWELFTH<glow.support.types.SymmetryType.TWELFTH>` symmetry to a
cartesian and a hexagonal lattice, respectively.

.. _quarter-symm:
.. figure:: images/lattice_qsym.png
   :alt: Cartesian lattice after applying a quarter symmetry
   :width: 400px
   :align: center

   Cartesian lattice after applying the :py:attr:`QUARTER<glow.support.types.SymmetryType.QUARTER>`
   type of symmetry.

.. _twelfth-symm:
.. figure:: images/lattice_twsym.png
   :alt: Hexagonal lattice after applying a twelfth symmetry
   :width: 400px
   :align: center

   Hexagonal lattice after applying the :py:attr:`TWELFTH<glow.support.types.SymmetryType.TWELFTH>`
   type of symmetry.

Users should note that |TOOL| does not recognize whether the layout of cells
replicates the full layout when any valid symmetry is applied.
It is up to the user to apply a symmetry that can be representative for the
specific layout of the lattice.

Setting the Lattice's Type of Geometry
""""""""""""""""""""""""""""""""""""""

The *SALT* module of *DRAGON5* identifies each type of geometry layout of the
lattice with a specific index value. In the *TDT* file, this is identified by
the *typgeo* value which is representative of the geometry layout (either full
or partial, if any symmetry is applied) and the type of BCs on the lattice's
borders.
User should note that specific values of *typgeo* are also associated to the
two different types of tracking allowed by the *SALT* module of *DRAGON5*
:cite:`dragon5-ug`. In particular, we have that:

  - values of `0`, `1` and `2` for *typgeo* are associated with a *TISO* tracking
    type, which produces non-cycling tracks distributed uniformally over the
    domain.
  - values greater that `2` for *typgeo* are associated with a *TSPC* tracking
    type, which indicates a cyclic tracking over a closed domain.

The items of the enumeration :py:class:`LatticeGeometryType<glow.support.types.LatticeGeometryType>`
identify the different *typgeo* values available in |TOOL|. In particular, we have:

  - :py:attr:`ISOTROPIC<glow.support.types.LatticeGeometryType.ISOTROPIC>` to
    represent a layout having an isotropic reflection on its boundaries. It is
    associated with a *TISO* tracking.
  - :py:attr:`SYMMETRIES_TWO<glow.support.types.LatticeGeometryType.SYMMETRIES_TWO>`
    to represent a layout having symmetries of two axis of angle ``pi/n`` (
    :math:`n>0`) on its boundaries. It is associated with a *TISO* tracking.
  - :py:attr:`ROTATION<glow.support.types.LatticeGeometryType.ROTATION>` to
    represent a layout with a rotation of angle ``2*pi/n`` (:math:`n>1`) for
    its boundaries. It is associated with a *TISO* tracking.
  - :py:attr:`RECTANGLE_TRAN<glow.support.types.LatticeGeometryType.RECTANGLE_TRAN>`
    to represent a cartesian layout having a translation BC to its boundaries.
    It is associated with a *TSPC* tracking.
  - :py:attr:`RECTANGLE_SYM<glow.support.types.LatticeGeometryType.RECTANGLE_SYM>`
    to represent a full, half and quarter symmetry for a cartesian layout.
    It is associated with a *TSPC* tracking.
  - :py:attr:`RECTANGLE_EIGHT<glow.support.types.LatticeGeometryType.RECTANGLE_EIGHT>`
    to represent a layout with an eighth symmetry. It is associated with a *TSPC*
    tracking.
  - :py:attr:`SA60<glow.support.types.LatticeGeometryType.SA60>` to represent a
    layout with an sixth symmetry. It is associated with a *TSPC* tracking.
  - :py:attr:`HEXAGON_TRAN<glow.support.types.LatticeGeometryType.HEXAGON_TRAN>`
    to represent a full hexagonal layout having a translation BC to its boundaries.
    It is associated with a *TSPC* tracking.
  - :py:attr:`RA60<glow.support.types.LatticeGeometryType.RA60>` to represent a
    layout with an sixth symmetry with both rotation and translation BCs to its
    boundaries. It is associated with a *TSPC* tracking.
  - :py:attr:`R120<glow.support.types.LatticeGeometryType.R120>` to represent a
    layout with an third symmetry with both rotation and translation BCs to its
    boundaries. It is associated with a *TSPC* tracking.
  - :py:attr:`S30<glow.support.types.LatticeGeometryType.S30>` to represent a
    layout with a twelfth symmetry. It is associated with a *TSPC* tracking.

When a :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` class is
instantiated, a default value for the property :py:attr:`type_geo<glow.geometry_layouts.lattices.Lattice.type_geo>`
is assigned according to the number and the type of cells.
Users can assign a value to this property directly, provided it is valid for
the lattice's geometry layout. This means that values specific for a type of
lattice and symmetry cannot be applied if not matching the current state of the
lattice.
For any values of *typgeo* involving BCs of type *translation*, the assignement
is performed only if the lattice is either made by a single cell or if enclosed
in a box.

|TOOL| provides also the method :py:meth:`set_type_geo()<glow.geometry_layouts.lattices.Lattice.set_type_geo>`
to set the item of the enumeration
:py:class:`LatticeGeometryType<glow.support.types.LatticeGeometryType>`.

The following code snippet shows different applications of the property
:py:attr:`type_geo<glow.geometry_layouts.lattices.Lattice.type_geo>`.

.. code-block:: python

  rect_lattice.type_geo = LatticeGeometryType.RECTANGLE_TRAN
  hex_lattice.set_type_geo(LatticeGeometryType.SA60)

Setting the value for the property does not result in any change in the lattice's
geometry layout. It influences the information written in the output *TDT* file
in terms of the BCs section, as this is strictly related to the *typgeo*.

Lattice's Regions Inspection
""""""""""""""""""""""""""""

When the regions of the lattice's technological or sectorized geometry are
displayed in the *SALOME* 3D viewer, information about a selected *region*,
including the assigned properties, can be inspected.
The method :py:meth:`get_regions_info()<glow.geometry_layouts.lattices.Lattice.get_regions_info>`
can be called directly in the Python console of *SALOME* from an object
of :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>`.
If no *region* (i.e. a *GEOM face*), or more than one, is selected when calling
the method, an exception is raised.
The available information, that is printed in the Python console, includes the
name of the lattice's *region* and the value for each of the assigned type of
properties.

Restoring Lattice's Cells
"""""""""""""""""""""""""

Similarly to the class :py:class:`Cell<glow.geometry_layouts.cells.Cell>`, also
the class :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` offers
a *restore* functionality.
The method :py:meth:`restore_cells()<glow.geometry_layouts.lattices.Lattice.restore_cells>`
allows users to restore the geometry layout of a group of cells of the lattice
by calling the method :py:meth:`restore()<glow.geometry_layouts.cells.Cell.restore>`
for each cell. The result is that any circular *region* of the cells is removed,
while also setting the cell's properties accordingly with the ones passed as
input to the method.
If any cells have no centered circular regions, the *restore* operation is not
performed for those specific cells.
In addition, users can specify whether the operation should be ignored for cells
whose circular regions (being part of the technological geometry) have not been
cut when overlapping with another cell (see :numref:`overlap`).

This method can be combined with the function :py:func:`get_changed_cells()<glow.geometry_layouts.lattices.get_changed_cells>`
to retrieve any cells whose geometry layout has been modified, making it easy
to restore them.

The following code snippet shows the case of a hexagonal lattice where a
central cell overlaps those of the layer below it. The *restore* operation
is applied to all the overlapped cells resulting in the lattice's geometry
layout of :numref:`restored-cells`.

.. code-block:: python

  # Build the lattice geometry layout
  cell = HexCell()
  cell.add_circle(0.2)
  cell.add_circle(0.3)
  cell.add_circle(0.4)
  cell.rotate(90)
  cell.set_properties({PropertyType.MATERIAL: ['MAT_1', 'MAT_2', 'MAT_3', 'MAT_4']})
  lattice = Lattice([])
  lattice.add_ring_of_cells(cell, 2)
  # A cell with greater dimensions is added in the lattice centre, overlapping
  # those of the layer below
  central_cell = HexCell(edge_length=1.5)
  central_cell.rotate(90)
  central_cell.set_properties({PropertyType.MATERIAL: ['MAT_4']})
  lattice.add_cell(central_cell, ())
  # Assemble all the layers
  lattice.build_regions()
  # Restore the overlapped cells
  lattice.restore_cells(
      get_changed_cells(lattice),
      {PropertyType.MATERIAL: 'MAT_4'},
      ignore_not_cut=False
  )
  lattice.show(PropertyType.MATERIAL)

.. _restored-cells:
.. figure:: images/lattice_restore.png
   :alt: Lattice's after restoring overlapped cells shown with MATERIAL colorset
   :width: 400px
   :align: center

   Hexagonal lattice's technological geometry showing the result of restoring
   the overlapped cells. The geometry layout is displayed with the
   :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>` colorset.

.. _lattice-export:

Lattice Analysis and Export
---------------------------

The aim of |TOOL| is to provide neutronics code users with a tool that allows
them to create geometry layouts and export the surface geometry representation
to a file. This file can then be used to perform a tracking with the *SALT*
module of *DRAGON5*.
The generated file is in the format *APOLLO2* requires for its *TDT* solver.

To meet this requirement, |TOOL| comes with a functionality for extracting the
necessary information about the geometry and generate the output file in the
required format.

Once the lattice geometry has been created using a
:py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>` instance, users can
run the export process by calling the function
:py:func:`analyse_and_generate_tdt()<glow.main.analyse_and_generate_tdt>`.
This function first analyses the lattice, then generates the output file
containing the extracted information.

This function operates on the provided :py:class:`Lattice<glow.geometry_layouts.lattices.Lattice>`
instance on the basis of specific configuration options defined in the dataclass
:py:class:`TdtSetup<glow.main.TdtSetup>`. Values for these options influence
the data about the surface geometry representation of the layout contained in
the output *TDT* file.
The available settings in the :py:class:`TdtSetup<glow.main.TdtSetup>` instance
include:

  - the type of geometry layout of the lattice's cells, as item of the enumeration
    :py:class:`GeometryType<glow.support.types.GeometryType>`. A value different
    from that used to display the lattice in the *SALOME* 3D viewer can be
    specified.
  - the type of property associated to the lattice's *regions*, as item of the
    enumeration :py:class:`PropertyType<glow.support.types.PropertyType>`.
    A value different to that used to apply the colorset to the *regions* can
    be specified.
  - the value of the *albedo*, indicating how much reflective the BCs are,
    i.e. the ratio of exiting to entering neutrons. This attribute can assume
    values between `0.0` (no reflection) and `1.0` (full reflection) for a
    :py:attr:`ISOTROPIC<glow.support.types.LatticeGeometryType.ISOTROPIC>`
    type of geometry of the lattice. If nothing is provided, a default value
    that corresponds to the lattice's geometry type is adopted (i.e. `1.0` for
    :py:attr:`ISOTROPIC<glow.support.types.LatticeGeometryType.ISOTROPIC>`
    geometry layouts, `0.0` for the others). An exception is raised if users
    provide a value different from `0.0` for a geometry type other than
    :py:attr:`ISOTROPIC<glow.support.types.LatticeGeometryType.ISOTROPIC>`,
    as this would not make sense.

The values set in the given :py:class:`TdtSetup<glow.main.TdtSetup>` instance
drives the analysis step in which the needed geometric data is extracted from
the lattice.
The first step consists in determining the lattice's *GEOM compound* to analyse;
this is selected on the basis of the :py:class:`GeometryType<glow.support.types.GeometryType>`
and on the applied :py:class:`SymmetryType<glow.support.types.SymmetryType>`.
Each :py:class:`Region<glow.geometry_layouts.cells.Region>` object, which
corresponds to the *regions* of the lattice compound, is associated with its
value of property type (:py:class:`PropertyType<glow.support.types.PropertyType>`)
for which the lattice is analysed. In addition, an index is assigned to each
*region* to ensure its identification.
The *GEOM edge* objects are then extracted and associated to the corresponding
regions. This means that each edge, identified with another index, has one or
two regions associated with it. Those associated with two regions are internal
edges, shared by two adjacent regions, while those associated with only one
*region* are border edges.
Lastly, the indices of the border edges are associated to a *boundary*, whose
type (as item of the enumeration :py:class:`BoundaryType<glow.support.types.BoundaryType>`)
and geometric data are determined on the basis of the
:py:class:`LatticeGeometryType<glow.support.types.LatticeGeometryType>` and the
applied :py:class:`SymmetryType<glow.support.types.SymmetryType>`.

.. only:: html

   :numref:`tdt-types` provides the association between
   :py:class:`LatticeGeometryType<glow.support.types.LatticeGeometryType>` and
   :py:class:`BoundaryType<glow.support.types.BoundaryType>` for the two type of
   cells with the symmetries available in |TOOL|.
   The first group of coloumns *LatticeGeometryType*-*BoundaryType* indicates the
   values for which a uniform tracking (i.e. *TISO*) should be performed in *SALT*;
   the second group refers to values which correspond to a cyclic tracking (i.e.
   *TSPC*).
   An :py:attr:`ISOTROPIC<glow.support.types.LatticeGeometryType.ISOTROPIC>` type
   of geometry does not correspond to any BC, whereas those having two types of
   BCs applies a :py:attr:`ROTATION<glow.support.types.BoundaryType.ROTATION>`
   on the internal boundaries and a :py:attr:`TRANSLATION<glow.support.types.BoundaryType.TRANSLATION>`
   on the external ones (see :numref:`tran-rota`).

.. only:: latex

   The following tables provides the association between
   :py:class:`LatticeGeometryType<glow.support.types.LatticeGeometryType>` and
   :py:class:`BoundaryType<glow.support.types.BoundaryType>` for the two type of
   cells with the symmetries available in |TOOL|.
   The first table indicates the values for which a uniform tracking (i.e.
   *TISO*) should be performed in *SALT*; the second table refers to values
   which correspond to a cyclic tracking (i.e. *TSPC*).
   An :py:attr:`ISOTROPIC<glow.support.types.LatticeGeometryType.ISOTROPIC>` type
   of geometry does not correspond to any BC, whereas those having two types of
   BCs applies a :py:attr:`ROTATION<glow.support.types.BoundaryType.ROTATION>`
   on the internal boundaries and a :py:attr:`TRANSLATION<glow.support.types.BoundaryType.TRANSLATION>`
   on the external ones (see :numref:`tran-rota`).

.. only:: html

   .. _tdt-types:

   .. table:: Available combinations for *TISO* and *TSPC* cases.
      :widths: auto
      :align: center

      +----------+--------------+---------------------+----------------------+---------------------+----------------------+
      | CellType | SymmetryType | LatticeGeometryType | BoundaryType         | LatticeGeometryType | BoundaryType         |
      +==========+==============+=====================+======================+=====================+======================+
      |          | FULL         | ISOTROPIC           |          /           | HEXAGON_TRAN        | TRANSLATION          |
      |          +--------------+---------------------+----------------------+---------------------+----------------------+
      |          | THIRD        | ROTATION            | TRANSLATION/ROTATION | R120                | TRANSLATION/ROTATION |
      |          +--------------+---------------------+----------------------+---------------------+----------------------+
      |  HEX     |              | SYMMETRIES_TWO      | AXIAL_SYMMETRY       | SA60                | AXIAL_SYMMETRY       |
      |          | SIXTH        +---------------------+----------------------+---------------------+----------------------+
      |          |              | ROTATION            | TRANSLATION/ROTATION | RA60                | TRANSLATION/ROTATION |
      |          +--------------+---------------------+----------------------+---------------------+----------------------+
      |          | TWELFTH      | SYMMETRIES_TWO      | AXIAL_SYMMETRY       | S30                 | AXIAL_SYMMETRY       |
      +----------+--------------+---------------------+----------------------+---------------------+----------------------+
      |          |              |                     |                      | RECTANGLE_TRAN      | TRANSLATION          |
      |          | FULL         | ISOTROPIC           |          /           +---------------------+----------------------+
      |          |              |                     |                      | RECTANGLE_SYM       | AXIAL_SYMMETRY       |
      |          +--------------+---------------------+----------------------+---------------------+----------------------+
      |  RECT    | HALF         | SYMMETRIES_TWO      | AXIAL_SYMMETRY       | RECTANGLE_SYM       | AXIAL_SYMMETRY       |
      |          +--------------+---------------------+----------------------+---------------------+----------------------+
      |          | QUARTER      | SYMMETRIES_TWO      | AXIAL_SYMMETRY       | RECTANGLE_SYM       | AXIAL_SYMMETRY       |
      |          +--------------+---------------------+----------------------+---------------------+----------------------+
      |          | EIGHTH       | SYMMETRIES_TWO      | AXIAL_SYMMETRY       | RECTANGLE_EIGHTH    | AXIAL_SYMMETRY       |
      +----------+--------------+---------------------+----------------------+---------------------+----------------------+

.. only:: latex

   .. raw:: latex

      {\small
      \begin{table}[ht]
      \centering
      \begin{tabularx}{.95\textwidth}{|X|X|X|X|}\hline
        \textbf{CellType} & \textbf{SymmetryType} & \textbf{LatticeGeometryType} & \textbf{BoundaryType} \\ \hline
        HEX & FULL & ISOTROPIC & N.D.\\ \hline
        HEX & THIRD & ROTATION & TRANSLATION\-ROTATION\\ \hline
        HEX & SIXTH & SYMMETRIES\_TWO & AXIAL\_SYMMETRY\\ \hline
        HEX & SIXTH & ROTATION & TRANSLATION\-ROTATION\\ \hline
        HEX & TWELFTH & SYMMETRIES\_TWO & AXIAL\_SYMMETRY\\ \hline
        RECT & FULL & ISOTROPIC & N.D.\\ \hline
        RECT & HALF & SYMMETRIES\_TWO & AXIAL\_SYMMETRY\\ \hline
        RECT & QUARTER & SYMMETRIES\_TWO & AXIAL\_SYMMETRY\\ \hline
        RECT & EIGHTH & SYMMETRIES\_TWO & AXIAL\_SYMMETRY\\ \hline
      \end{tabularx}
      \caption{Available combinations for \textit{TISO} case}
      \end{table}
      }

      {\small
      \begin{table}[ht]
      \centering
      \begin{tabularx}{.95\textwidth}{|X|X|X|X|}\hline
        \textbf{CellType} & \textbf{SymmetryType} & \textbf{LatticeGeometryType} & \textbf{BoundaryType} \\ \hline
        HEX & FULL & HEXAGON\_TRAN & TRANSLATION\\ \hline
        HEX & THIRD & R120 & TRANSLATION\-ROTATION\\ \hline
        HEX & SIXTH & SA60 & AXIAL\_SYMMETRY\\ \hline
        HEX & SIXTH & RA60 & TRANSLATION\-ROTATION\\ \hline
        HEX & TWELFTH & S30 & AXIAL\_SYMMETRY\\ \hline
        RECT & FULL & RECTANGLE\_TRAN & TRANSLATION\\ \hline
        RECT & FULL & RECTANGLE\_SYM & AXIAL\_SYMMETRY\\ \hline
        RECT & HALF & RECTANGLE\_SYM & AXIAL\_SYMMETRY\\ \hline
        RECT & QUARTER & RECTANGLE\_SYM & AXIAL\_SYMMETRY\\ \hline
        RECT & EIGHTH & RECTANGLE\_SYM & AXIAL\_SYMMETRY\\ \hline
      \end{tabularx}
      \caption{Available combinations for \textit{TSPC} case}
      \end{table}
      }

The different values of BCs that are automatically applied by |TOOL| to the
boundaries of the lattice's geometry layout are identified by the items of the
enumeration :py:class:`BoundaryType<glow.support.types.BoundaryType>`. Their
meaning and usage is the same as specified in :cite:`dragon5-ug`:

  - :py:attr:`VOID<glow.support.types.BoundaryType.VOID>`, indicating that
    boundaries have zero re-entrant angular flux;
  - :py:attr:`REFL<glow.support.types.BoundaryType.REFL>`, indicating a
    reflective boundary condition;
  - :py:attr:`TRANSLATION<glow.support.types.BoundaryType.TRANSLATION>`,
    indicating that the analysed layout is connected to another one for all its
    boundaries, thus treating an infinite geometry with translation symmetry;
  - :py:attr:`ROTATION<glow.support.types.BoundaryType.ROTATION>`, indicating
    a rotation symmetry;
  - :py:attr:`AXIAL_SYMMETRY<glow.support.types.BoundaryType.AXIAL_SYMMETRY>`,
    indicating a reflection symmetry;
  - :py:attr:`CENTRAL_SYMMETRY<glow.support.types.BoundaryType.CENTRAL_SYMMETRY>`,
    indicating a mirror reflective boundary condition.

.. _tran-rota:
.. figure:: images/lattice_tran_rota.png
   :alt: Assignment of ROTATION and TRANSLATION BC types to boundaries
   :width: 400px
   :align: center

   Showing to which boundaries the :py:attr:`ROTATION<glow.support.types.BoundaryType.ROTATION>`
   and :py:attr:`TRANSLATION<glow.support.types.BoundaryType.TRANSLATION>` BC
   types are assigned to (third symmetry case).

Given all the geometric data extracted from the lattice, the output file is
generated. Its structure consists of five sections, that are:

  - the *header* section, providing information about the type of geometry
    (*typgeo* value), the number of *folds* (*nbfold* value), which is
    consistent with the *typgeo*, the number of *nodes* (i.e. the regions),
    the number of *elements* (i.e. the edges).
  - the *regions* section, providing a list of indices attributed to the
    *regions* in the lattice. It also contains the definition of the *macros*
    to indicate subvolumes of the assembly.
  - the *edges* section, providing the geometric information about all the edges
    in the geometry layout, as well as the indices of the regions they belong
    to.
  - the *boundary conditions* section, providing information about the BC types
    and the indices of the edges belonging to each boundary.
  - the *property* section, indicating the index of each value of the considered
    property type (e.g. the :py:attr:`MATERIAL<glow.support.types.PropertyType.MATERIAL>`
    one). The order in which values are present matches the indices of the
    regions.

.. _usage:

Usage
-----

|TOOL| can be used directly by writing down a Python script where the single
needed modules can be imported; alternatively, users can import all the modules
at once to have them available by setting the following import instruction:

.. code-block:: python

  from glow import *

Given that, classes and methods are directly accessible and users can exploit
them to:

- assemble the geometry;
- assign properties to regions;
- visualize the result in the *SALOME* 3D viewer;
- perform the geometry analysis and the output file generation.

To run this script, users can:

- provide it as argument when running *SALOME*;

    .. code-block:: bash

      salome my_script.py

- load it directly from within the *SALOME* application.

In addition, since *SALOME* comes with an embedded Python console, users can
import the |TOOL| modules and exploit its functionalities directly.

To see some of the |TOOL| functionalities in action, please refer to the script
files present in the ``tutorials`` folder: they are intended to show a few
case studies and how they are managed in |TOOL|.
For further information about the available classes and methods, please refer
to the :doc:`api_guide` section.
