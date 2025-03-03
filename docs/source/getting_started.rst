===============
Getting started
===============

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
files present in the `test/examples` folder: they are intended to show few case
studies and how they are managed in |TOOL|.
For further information about the available classes and methods, please refer
to the :doc:`api_guide` section.

