.. only:: latex

  .. raw:: latex

     \setcounter{secnumdepth}{3}

============
Introduction
============

|TOOL| (*Geometry Layout Oriented Workflow*) is a Python package
providing 2D unstructured geometries to the `DRAGON5 <http://merlin.polymtl.ca/version5.htm>`_
lattice transport computer code for tracking calculations.

The *DRAGON5* computer code provides several models, in terms of FORTRAN modules,
to simulate the neutronic behaviour of the fuel elements
in a nuclear reactor :cite:`dragon5-ug`.
Among the available modules, some focus on solution techniques for the
neutron transport equation, which describes the statistical behaviour of a
large population of particles :cite:`hebert2009`.
In particular, the *SALT* module can process a geometry mesh to compute the
tracking information requested to solve the neutron transport equation.
The geometry processed by *DRAGON5* is of the surface type where cell mesh
boundaries are described by equations.

Among the many methods for the solution of the neutron transport equation, the
following two are used with *DRAGON5*:

  - the *Method of Characteristics* (*MOC*), an approach used to solve the
    neutron transport equation by tracing the paths (*characteristics*) of
    neutrons through a discretized geometry. It evaluates the neutron flux
    distribution by integrating the transport equation along the characteristics,
    which are straight lines in neutral particle transport, accounting for
    interactions with the material.
  - the *Collision Probability Method* (*CPM*), an approach that computes the
    probability of a neutron emitted in one region to have the first collision
    in another region. It generally assumes isotropic neutron emission.

Having a proper description of the geometry layout is essential to perform
tracking calculations in *DRAGON5*.
Geometries for fuel assemblies used in nuclear reactors are typically made of
several fuel pins arranged according to regular patterns to form Cartesian or
hexagonal lattices.
Even though *DRAGON5* can track any convex surface-type geometry without empty
concavities, it can only natively handle the generation of *structured geometries*.
These rely on regular basic surface elements, such as polygons and circles,
repeated to create uniform patterns, which in *DRAGON5* are built using a small
set of predefined operations.
This is done with the *GEO* module where the *G2S* module converts the geometry
layout to a *.dat* description file used by the *SALT* tracking module.
More complex cases deal with *unstructured geometries* which are made by general
elements having different shapes and dimensions. The combination of surface
elements can be obtained as the result of applying boolean operations, such
as cuts and unions.
Unstructured geometries require external tools based on *Constructive Solid
Geometry* (*CSG*), a method for creating 2D/3D geometries from combinations of
simpler primitive shapes.

So far, several applications have been developed to address the generation of
unstructured surface geometries, among the main ones we can cite *SALOMON* and
*ALAMOS* :cite:`tomatis2022`.
*SALOMON* is a Python2-based application developed in 2001 by EDF that exploits
the *GEOM* module of the `SALOME <https://www.salome-platform.org/>`_ platform,
an open-source environment offering 2D/3D CAD modelling capabilities
:cite:`salome2007`.
*SALOMON* relies on an input file with a specific structure containing the
description of the lattice's geometry layout in terms of cells and materials
associated to each region of the lattice. It allows users to construct a geometric
description file in the format of the *TDT* solver of *APOLLO2*.
Currently, SALOMON's development seems no longer active, resulting in the tool
becoming obsolete.
*ALAMOS* is GUI-based application developed by CEA for *APOLLO3* that exploits
the Python libraries offered by the *SALOME* platform :cite:`tomatis2022`. In
particular, it makes use of the *MEDCoupling* module of *SALOME* to build the
geometries and apply meshes to them. The regions of the geometry can also be
characterised in terms of material property maps.
While *ALAMOS* can handle any kind of complex geometry layout, *SALOMON*
can only deal with lattices made of cartesian-type cells, thus greatly limiting
its application.
In addition, *SALOMON* does not provide any support to the construction of
surface geometries by exploiting the most common boolean operations, which,
instead, are supported by *ALAMOS*.
Unfortunately, *ALAMOS*, despite being a more complete tool with respect to
*SALOMON*, does not have an open-source distribution.

|TOOL| was developed to offer an open-source alternative for generating
unstructured surface geometries, exploiting the *Constructive Solid Geometry*
functionalities provided by the APIs of the *GEOM* module of the *SALOME*
platform.
In |TOOL|, complex surface geometry layouts are constituted by 2D areas bounded
by closed sets of edges, which are referred to as *regions*.
A *cell* is a base geometry layout built from a characteristic surface, such as
a rectangle or a hexagon, subdivided in different areas, each representing a
*region*.
The repetition of adjacent *cells* in the 2D space constitutes a *lattice*.
|TOOL| can handle *lattices* made of cartesian or hexagonal *cells*, including
the option to frame the lattice in a box.
Different types of properties, such as the *material*, can be assigned to each
*region*.
|TOOL| supports Euclidean geometric **transformations**, such as translation and
rotation, and **boolean** operations, such as union, intersection, cut and
partitioning, among geometric shapes.

Since calculations on smaller geometries are expected to be computationally
cheaper, symmetries should be considered whenever possible.
|TOOL| can perform cuts to extract parts from an existing layout, thereby
isolating the minimum portion of the geometry required to describe the entire
pattern. As an example for a Cartesian lattice with a 180 degree reflection
symmetry, the cut part would be half of the lattice.
Special symmetries according to the main lattice types are implemented:

  - *Cartesian* - half, quarter and eighth symmetries.
  - *Hexagonal* - third, sixth and twelfth symmetries.

|TOOL| does not only exploit the geometric functionalities offered by the *GEOM*
module of *SALOME*; it also displays the built geometries in *SALOME* 3D viewer
through its graphical user interface (*GUI* module).
Additionally, if the basic functionalities offered by |TOOL| are not enough to
construct specific layouts for cells and lattices, users can rely on the
geometry-building functionalities offered by the *GEOM* module directly. The
resulting custom layouts can still be used in |TOOL|.

In *DRAGON5*, tracking relies on a description file (*.dat*) of the geometry
layout. |TOOL| allows the collection of geometric and property information from
the built lattice layout, producing this *.dat* file in the format of the *TDT*
solver of *APOLLO2*.

|TOOL| is developed by the **Codes & Methods** Department of |newcleo| and it
is released under the |LICENSE| **License**.

.. admonition:: DISCLAIMER

   The *DRAGON5* lattice code is not distributed with |TOOL|. Please, refer
   to `<http://merlin.polymtl.ca/version5.htm>`_ for downloading the latest
   version and for the installation instructions.
