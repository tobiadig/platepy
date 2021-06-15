"""
platepy
=====

Provides\n
* Methods to easily model structural components for plates via the gmsh API (see https://gmsh.info/). \n
* Mesh generation via gmsh API. \n
* FEM calculation of the system. \n
* Results visualization. \n
* Analythical solutions according to Kirchhoff plate theory.
"""

from .analyticPlateSolutions import *
from .createModel import *
from .generateMesh import *
from .solveModel import *
from .displayModel import *

