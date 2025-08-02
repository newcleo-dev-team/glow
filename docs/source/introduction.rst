.. only:: latex

  .. raw:: latex

     \setcounter{secnumdepth}{3}

============
Introduction
============

|TOOL| (*Geometry Layout Oriented Workflows*) is a Python application
intended for providing the `DRAGON5 <http://merlin.polymtl.ca/version5.htm>`_
lattice code with a tool for building 2D unstructured geometries for performing
tracking analyses.

In the context of enabling neutronics experts to perform analyses on the
different parts of the core of a nuclear reactor, a detailed geometry layout,
describing how fuel pins are positioned, the subdivision of the layout into
specific regions of calculus and the assignment of properties to the single
regions are required.

Among the different available analyses, we have the following:

  - the *Method of Characteristics* (*MoC*), an approach used to solve the
    neutron transport equation by tracing the paths (*characteristics*) of
    neutrons through a discretized geometry. It evaluates the neutron flux
    distribution by integrating the transport equation along straight-line
    trajectories, accounting for interactions with the material.
  - the *Collision Probability Method* (*CPM*), an approach that computes the
    probability of a neutron emitted in one region to have a collision in
    another region. It assumes isotropic neutron emission and transport in
    one-dimensional or radially symmetric geometries. It solves for the flux
    by balancing sources and collisions based on precomputed collision
    probabilities.

The `DRAGON5 <http://merlin.polymtl.ca/version5.htm>`_ lattice code offers
several calculation modules focusing on solution techniques for the neutron
transport equation :cite:`dragon5-ug`.
In particular, the *SALT* module enables the tracking of neutrons by solving
the transport equation on the basis of a 2D representation of the geometry
layout to study (i.e. a surface geometry). The two type of analyses mentioned
above are both supported by *SALT*.

Having a proper description of the geometry layout to analyse is essential to
perform tracking calculations in *DRAGON5*.
Geometries for fuel assemblies are typically structured in terms of repeated
fuel pin cells to form specific patterns, or lattices, on a Cartesian grid or
a hexagonal one (as in fast reactors).
*DRAGON5* can natively handle the generation of simple geometry layouts by means
of the *GEO* module; to perform a tracking they need to be translated into a
description file by the *G2S* module.
However, complex layouts in reactors, require to rely on boolean operations
among geometrical shapes that are offered by *Constructive Solid Geometry*
(*CSG*) tools.
So far, several applications have been developed to address the generation of
surface geometries, among the main ones we can cite *SALOMON* and *ALAMOS*.

*SALOMON* is a Python2-based application that exploits the *GEOM* module of
*SALOME* to construct a geometric description file in the *APOLLO2* format
for its TDT solver. It relies on input files with a specific structure
containing the description of the lattice's geometry layout in terms of cells
and materials associated to each region of the lattice.

*ALAMOS* is GUI-based application developed for APOLLO3 that exploits the
Python libraries offered by the *SALOME* platform :cite:`tomatis2022`. In
particular, it makes use of the *MEDCoupling* module of *SALOME* to build the
geometries and apply meshes to them; it also allows users to assign material
property maps to the regions of the geometry.

While the *ALAMOS* can handle any kind of complex geometry layout, *SALOMON*
can only deal with lattices made by cartesian-type cells, thus greatly limiting
its application.
In addition, *SALOMON* does not provide any support to the construction of
surface geometries by exploiting the most common boolean operations, which,
instead, are supported by *ALAMOS*.
On the other hand, *ALAMOS*, despite being a more complete tool with respect
to *SALOMON*, does not have an open-source distribution.

|TOOL| was developed with the purpose of offering an open-source alternative to
*ALAMOS*; the idea was to include all the functionalities needed to assemble
complex geometry layouts made of *regions*, intended as any geometrical surface
defined by a closed set of edges. Different types of properties can be assigned
to each region. An example of a property type are the *materials*.

Similarly to *SALOMON*, |TOOL| geometry building functionalities relies on the
APIs of the *GEOM* module of the `SALOME <https://www.salome-platform.org/>`_
platform, an open-source environment offering 2D/3D CAD modelling capabilities.
However, contrary to *SALOMON*, |TOOL| exploit these APIs to support users even
in the construction of lattices made of hexagonal cells and not only of
cartesian ones.
Operations among geometrical shapes, such as **transformation** (e.g.
translation, rotation, etc.) and **boolean** ones (e.g. cuts, partitioning,
etc.) are also supported; the same goes for applying specific **symmetry**
operations to the lattice according to the type of its cells.
The following symmetries, associated to each lattice type, are handled:

  - *Cartesian* - half, quarter and eighth symmetries.
  - *Hexagonal* - third, sixth and twelfth symmetries.

|TOOL| does not only exploits the geometrical functionalities offered by the
*GEOM* module of *SALOME*; it also relies on its graphical user interface
(*GUI* module) to visualize each geometrical object built with |TOOL| in the
3D viewer of *SALOME*. In addition, users can exploit the geometry building
functionalities offered by the *GEOM* module directly by assembling custom
layouts for cells and lattices and providing them to the instances of the
|TOOL| classes.

Tracking analyses in *DRAGON5* relies on a the description file of the geometry
to analyse: |TOOL| allows the collection of the required information from the
built lattice layout, so that its surface geometry representation can be
exported to a *.dat* file in the format *APOLLO2* requires for its *TDT* solver.

|TOOL| is developed by the **Codes & Methods** group of |newcleo| and it is
released under the |LICENSE| **License**.

**DISCLAIMER:**

  - The *DRAGON5* lattice code is not distributed with |TOOL|. Please, refer
    to `<http://merlin.polymtl.ca/version5.htm>`_ for downloading the latest
    version and for the installation instructions.
