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
:Contributors: Gianluca Cirillo, Camilla De Santis
:Date: 01/08/2025

Introduction
------------

|TOOL| (*Geometry Layout Oriented Workflow*) is a Python package
providing 2D unstructured geometries to the `DRAGON5 <http://merlin.polymtl.ca/version5.htm>`_
lattice transport computer code.

Natively, *DRAGON5* handles the generation of simple geometry layouts in
terms of lattices made up of cells, which represent the base unit of a fuel
assembly. This is performed by the *GEO* module. However, this cannot exploit
any of the boolean operations among geometrical shapes that are offered by
*Constructive Solid Geometry* (*CSG*) tools.

|TOOL| was developed with the purpose of offering a tool that includes all the
functionalities needed to assemble complex geometry layouts. It relies on the
APIs of the *GEOM* module of the `SALOME <https://www.salome-platform.org/>`_
platform, an open-source environment offering 2D/3D CAD modelling capabilities.
|TOOL| supports euclidean geometric **transformation**, such as translation and
rotation, and **boolean** operations.

|TOOL| can operate cuts to extract parts out of an existing layout, thus
focusing on the minimum portion of the geometry that can reproduce the whole
motif by unfolding symmetries. Calculations on smaller geometries are
expected to be computationally cheaper. The application of special symmetries
to the main lattice types are implemented:

  - *Cartesian* lattice - half, quarter and eighth symmetries.
  - *Hexagonal* lattice - third, sixth and twelfth symmetries.

The geometries built with |TOOL| can be visualized in the 3D viewer of *SALOME*
through the graphical user interface (*GUI* module).

|TOOL| can export geometry layouts to *.dat* files using the *TDT* format of
*APOLLO2* where the cells mesh boundaries are given by surface equations.

|TOOL| is developed by the **Codes & Methods** Department of |newcleo| and it is
released under the |LICENSE| **License**.

.. admonition:: DISCLAIMER

   The *DRAGON5* lattice code is not distributed with |TOOL|. Please, refer
   to `<http://merlin.polymtl.ca/version5.htm>`_ for downloading the latest
   version and for the installation instructions.

Project Structure
-----------------

The project is organized according to the following folder structure:

.. code:: text

  <glow parent folder>
    ├── docs/
    ├── glow/
    ├── resources/
    ├── tests/
    ├── tutorials/
    ├── LICENSE
    └── README.rst


- ``docs``: contains the files for generating the project documentation with
  Sphinx;
- ``glow``: contains all the modules that provide the functionalities
  implemented by |TOOL|;
- ``resources``: contains files that support the configuration and operation
  of |TOOL|;
- ``tests``: contains both unit tests and functional tests, the former to ensure
  the correctness of individual code units, the latter to ensure the overall
  functionalities behave as expected;
- ``tutorials``: contains Python scripts that provide use cases showing how to
  use the |TOOL| functionalities for building different geometry layouts.

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

- ``sphinx`` :math:`>= 8.2.3`
- ``sphinx-rtd-theme`` :math:`>= 3.0.2`
- ``myst-parser`` :math:`>= 4.0.1`
- ``sphinxcontrib-bibtex`` :math:`>= 2.6.3`

How to Install
--------------

To install the |TOOL| application, clone the repository at
https://github.com/newcleo-dev-team/glow using the following command:

  .. code-block:: bash

    git clone https://github.com/newcleo-dev-team/glow

Now, from the root folder of |TOOL| the following command can be run to
automatically install all the needed dependencies:

  .. code-block:: bash

      pip install .

To upgrade the |TOOL| package, please type the ``install`` command along with
the ``--upgrade`` or ``-U`` flag:

  .. code-block:: bash

      pip install --upgrade .

Since |TOOL| exploits the *GEOM* module of *SALOME*, a correct installation
of *SALOME* is required. Please, refer to the *Building and installing* section
of the *SALOME* `FAQ <https://www.salome-platform.org/?page_id=428>`_ page for
the installation instructions according to the user's specific OS.

Please, note that the |TOOL| usage is limited to the OSs supported by *SALOME*
itself.

How to Use
----------

|TOOL| can be used directly by writing down a Python script where the single
needed modules can be imported; alternatively, users can import all the modules
at once to have them available by setting the following import instruction:

.. code-block:: python

  from glow import *

Given that, classes and methods are directly accessible and users can exploit
them to:

- assemble the geometry;
- assign properties to regions;
- visualize the result in the *SALOME* 3D viewer (see image below);
- perform the geometry analysis and generate the output *TDT* file.

![Cartesian lattice after applying a one-eighth symmetry](resources/example_glow_geometry.png)

The above image shows the resulting geometry layout obtained by applying a
symmetry that extracts a eighth of the entire Cartesian lattice.

To run any script using the |TOOL| functionalities, users can:

- provide it as argument when running *SALOME*;

    .. code-block:: bash

      salome my_script.py

- load it directly from within the *SALOME* application.

In addition, since *SALOME* comes with an embedded Python console, users can
import the |TOOL| modules and exploit its functionalities directly.

For a detailed description of the functionalities provided by |TOOL|, please
refer to the *Getting Started* chapter of the documentation.
Python scripts are also provided in the ``tutorials`` folder. They are
intended to show some case studies and how they are managed in |TOOL|.

Documentation
-------------

The Sphinx documentation can be built both in *html* and *LaTeX* formats by
executing the following command in the folder ``docs``:

  .. code-block:: bash

      make html

  .. code-block:: bash

      make latexpdf

To see the available templates for generating the documentation in *PDF*
format and to choose among them, please look at the ``docs/conf.py`` file.

.. _How to Contribute:

How to Contribute
-----------------

For anyone wishing to contribute to the development of the |TOOL| project,
report issues or problems with the software, or request support, please refer
to this
`web page <https://github.com/newcleo-dev-team/glow/blob/master/CONTRIBUTIONS.rst>`_.

Acknowledgements
----------------

|newcleo| is thankful to Professor Alain Hébert and the entire *DRAGON5*
development team at the **Polytechnique of Montreal**, Canada, for their
constant support.
