import salome
salome.salome_init()
import SALOMEDS

import GEOM
# Different import handling of the 'geomBuilder' and 'geomtools' modules.
# The first import instructions only works within VSCode to enable docstrings
# and autocompletion, not when running the code. In the second case, the
# 'ImportError' exception is raised and the import instructions, as specified
# by the SALOME user guide, are used.
try:
    from geom import geomBuilder, geomtools
    geompy = geomBuilder.New()
except ImportError:
    from salome.geom import geomBuilder, geomtools
    geompy = geomBuilder.New()

from salome.kernel.studyedit import getStudyEditor
studyEditor = getStudyEditor()
gst = geomtools.GeomStudyTools(studyEditor)
