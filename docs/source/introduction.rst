.. only:: latex

  .. raw:: latex

     \setcounter{secnumdepth}{3}

============
Introduction
============

|TOOL| (*Geometry Layout Oriented Workflows*) is a Python application
intended for providing the `DRAGON5 <http://merlin.polymtl.ca/version5.htm>`_
lattice code with a tool for building **non-native** geometries for the
*Method of Characteristics* (*MoC*) calculations, which requires a 2D detailed
geometry.

Natively, *DRAGON5* can handle the generation of simple geometry layouts, in
terms of lattices made up of cells, which represent the base unit of an
assembly. This is performed by the *GEO* module: however, this cannot exploit
any of the boolean operations among geometrical shapes that are offered by
*Constructive Solid Geometry* (*CSG*) tools. So far, several tools have been
developed to address the generation of surface geometries, among the main ones
we can cite *SALOMON* and *ALAMOS*.
While the latter can handle any kind of complex geometry layout, the former can
only deal with cartesian-type cells, thus greatly limiting its application.
In addition, *SALOMON* does not provide any support to the construction of
surface geometries by exploiting the most common boolean operations, which,
instead, are supported by *ALAMOS*.
On the other hand, *ALAMOS*, despite being a more complete tool with respect
to *SALOMON*, does not have an open-source distribution.

|TOOL| was developed with the purpose of offering an open-source alternative to
*ALAMOS*; the idea was to include all the functionalities needed to assemble
complex geometry layouts made of:

  - regions, intended as any geometrical surface defined by a closed set of
    edges;
  - the properties (e.g. materials, etc.) assigned to each region.

Similarly to *SALOMON*, |TOOL| geometry building functionalities relies on the
APIs of the *GEOM* module of the `SALOME <https://www.salome-platform.org/>`_
platform, an open-source environment offering 2D/3D CAD modelling capabilities.
However, contrary to *SALOMON*, |TOOL| exploit these APIs to support users even
in the construction of lattices made of hexagonal cells and not only of
cartesian ones.
Operations among geometrical shapes, such as **transformation** (e.g.
translation, rotation, etc.) and **boolean** ones (e.g. fuse, intersection,
etc.) are also supported; the same goes for applying specific **symmetry**
operations to the lattice according to the type of its cells.
Depending on the cell types, the following symmetries are handled:

  - *Cartesian* - half, quarter and eighth sections.
  - *Hexagonal* - SA60 and S30 symmetries.

|TOOL| does not only exploits the geometrical functionalities offered by the
*GEOM* module of *SALOME*; it also relies on its graphical user interface
(*GUI* module) to visualize each geometrical object built with |TOOL| in the
3D viewer of *SALOME*.

*MoC* calculations in *DRAGON5* relies on a proper description of the geometry
to analyse: |TOOL| allows the collection of the required information from the
built lattice layout, so that its surface geometry representation can be
exported to a *.dat* file in the *TDT* *APOLLO2* format.

|TOOL| is developed by the **Codes & Methods** group of |newcleo| and it is
released under the |LICENSE| **License**.

**DISCLAIMER:**

  - The *DRAGON5* lattice code is not distributed with |TOOL|. Please, refer
    to `<http://merlin.polymtl.ca/version5.htm>`_ for downloading the latest
    version and for the installation instructions.
  - The present version of |TOOL| is distributed as a **beta** version, meaning
    that the implementation is still in progress. Should users find any issue,
    please refer to the :doc:`contacts` section to know how to report it.
