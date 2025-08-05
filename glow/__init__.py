"""
GLOW (Geometry Layout Oriented Workflows) is a Python application intended
for providing the DRAGON5 lattice code with a tool for building non-native
geometries for performing tracking analyses with the SALT module.
"""
from glow.interface.geom_interface import *
from glow.geometry_layouts.geometries import *
from glow.geometry_layouts.cells import *
from glow.geometry_layouts.lattices import *
from glow.main import *
from glow.support.types import *
from glow.support.utility import *

__version__ = "1.0.0"
__author__ = "Davide Manzione, Daniele Tomatis"
__company__ = "newcleo"
__date__ = "01 August 2025"
