'''
Provides

* Methods to easily model structural components for plates via the gmsh API (see https://gmsh.info/).
* Mesh generation via gmsh API.
* FEM calculation of displacements, rotations, bending moments, shear forces and required bending moment resistance according to SIA262.
* Results visualization.
* Analythical solutions according to Kirchhoff plate theory.

*Copywrite Tobia Diggelmann (ETH Zurich) 23.06.2021.*
'''

from .analyticPlateSolutions import *
from .createModel import *
from .generateMesh import *
from .solveModel import *
from .displayModel import *

