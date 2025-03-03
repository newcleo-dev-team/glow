|TOOL|
======

.. image:: resources/icons/newcleo_logo.png
   :width: 100
   :alt: newcleo logo

.. |newcleo| replace:: *new*\cleo
.. |TOOL| replace:: **GLOW**
.. |LICENSE| replace:: **LGPL-2.1**
.. _newcleo: https://www.newcleo.com/


:Authors: Davide Manzione, Daniele Tomatis
:Contributors: Gianluca Cirillo
:Date: 03/03/2025

Introduction
------------

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
    please refer to the :doc:`How to Contribute` section to know how to report
    it.

Project Structure
-----------------

The project is organized according to the following folder structure:

.. code:: text

  <glow parent folder>
    ├── docs/
    ├── glow/
    ├── resources/
    ├── tests/
    ├── LICENSE
    └── README.rst


- ``docs``: contains the files for generating the project documentation with
  Sphinx;
- ``glow``: contains all the modules that provide the functionalities
  implemented by |TOOL|;
- ``resources``: contains files that support the configuration and operation
  of |TOOL|;
- ``tests``: contains Python scripts that provide examples showing how to use
  the |TOOL| functionalities for building simple geometries and exporting them
  in the *TDT* format.

Dependencies
------------

To run the code, the following dependencies must be satisfied:

- ``Python`` :math:`>= 3.11+`
- ``SALOME`` :math:`>= 9.12`
- ``typing-extensions`` :math:`>= 4.12.2`
- ``PyQt5`` :math:`>= 5.15.11`
- ``psutil`` :math:`>= 7.0.0`

To build the documentation in both *html* and *LaTeX* formats, the following
dependencies must be satisfied:

- ``sphinx`` :math:`>= 7.3.7`
- ``sphinx-rtd-theme`` :math:`>= 3.0.2`
- ``myst-parser`` :math:`>= 3.0.1`
- ``sphinxcontrib-bibtex`` :math:`>= 2.6.3`

How to Install
--------------

To install the |TOOL| application, please check that all the dependencies
are met, and then clone the repository at
https://github.com/newcleo-dev-team/glow using the following command:

  .. code-block:: bash

    git clone https://github.com/newcleo-dev-team/glow

Since |TOOL| exploits the *GEOM* module of *SALOME*, a correct installation
of *SALOME* is required. Please, refer to the *Building and installing* section
of the *SALOME* `FAQ <https://www.salome-platform.org/?page_id=428>`_ page for
the installation instructions according to the user's specific OS.

Please, note that the |TOOL| usage is limited to the OSs supported by *SALOME*
itself.

How to Use
----------

|TOOL| can be used directly by writing down a Python script that exploits the
provided classes and methods to:

- assemble the geometry;
- assign properties to regions;
- visualize the result in the *SALOME* 3D viewer;
- perform the geometry analysis and the output *TDT* file generation.

To run this script, users can:

- provide it as argument when running *SALOME*;

    .. code-block:: bash

      salome my_script.py

- load it directly from within the *SALOME* application.

In addition, since *SALOME* comes with an embedded Python console, users can
import the |TOOL| modules and exploit its functionalities directly.

To see some of the |TOOL| functionalities in action, please refer to the script
files present in the `test/examples` folder. They are intended to show few case
studies and how they are managed in |TOOL|.

Documentation
-------------

The Sphinx documentation can be built both in *html* and *LaTeX* formats by
executing the following command in the folder ``docs/``:

  .. code-block:: bash

      make html

  .. code-block:: bash

      make latexpdf

To see the available templates for generating the documentation in *PDF*
format and to choose among them, please look at the ``docs/conf.py`` file.

How to Contribute
-----------------

For anyone wishing to contribute to the development of the |TOOL| project,
report issues or problems with the software, or request support, please refer
to this
`web page <https://github.com/newcleo-dev-team/glow/blob/master/CONTRIBUTIONS.rst>`_.

Acknowledgements
----------------

|newcleo| is thankful to prof. Alain Hébert and the whole *DRAGON5* development
team of the **Polytechnique of Montreal** for their constant support.
