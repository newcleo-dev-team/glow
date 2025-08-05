=========================
Application configuration
=========================

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
of the *SALOME* platform is required.
Please, refer to the *Building and installing* section of the *SALOME*
`FAQ <https://www.salome-platform.org/?page_id=428>`_ page for the
installation instructions according to the user's specific OS.

Please, note that the |TOOL| usage is limited to the OSs supported by *SALOME*
itself.
