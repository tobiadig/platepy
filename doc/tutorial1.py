# ------------------------------------------------------------------------------
#
#  platepy tutorial 1
#
# Creating a simple slab geometry and computing the FEM solution
#
# ------------------------------------------------------------------------------

# The program is entirely defined in the `platepy.py' module (which contains the
# full documentation of all the function). The easier way is to directely import
# all function with the * command:
from platepy import *
import numpy as np

# The model will now be initialized
tutorialModel = PlateModel()

# To define the structural components of the model a dictionary has to be created.
# The required entries of the dictionary are described in the documentation of the 
# relative class.

# let's define the concrete by defining the dictionary and the required entries:
# (attention to the units! by default all lengths are considered to be in meters,
# forces in kN).
concreteDict = {}
concreteDict["eModule"] = 32.1*1e6 #kN/m2
concreteDict["gModule"] =  14.36*1e6 #kN/m2
concreteDict["nu"] = 0.17
C25_30 = Concrete(concreteDict)

# alternatively, a standard concrete type can be used:
C30_37 = StandardConcrete("C30_37")


# Let's now define the structural components with the same procedure.
# First a square plate with a 10m side:
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [10,0],[10,10], [0,10], [0, 0]])
plateDict["thickness"] = 0.3
plateDict["body"]=C30_37
plate = Plate(plateDict)

# and let's add the plate to the model:
tutorialModel.addPlate(plate)

# now we define two walls on the left and right sides of the plate:
# the wall dictionary can be re-used by only changing the outline.
wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0],[0,10]])
wallDict["thickness"] = 0.5 # m
# The entries in the support array define is the relative DOF is free (0)
# or blocked (1). The first entry is the vertical displacement, the second
# is the rotation on the axis in the direction of the wall, the third entry
# is the rotation on the axis perpendicular to the wall. For simply supported
# hard boundary condition:
wallDict["support"] = np.array([1, 0, 1])
wall1 = Wall(wallDict)

# the second wall:
wallDict["outlineCoords"] = np.array([[10,0],[10,10]])
wall2 = Wall(wallDict)

# and add to the model:
tutorialModel.addWall(wall1)
tutorialModel.addWall(wall2)

# add a column in the center:
columnDict = {}
columnDict["outlineCoords"] = np.array([[5,5]])
columnDict["support"] = np.array([1, 0, 0])
columnDict["width"] = 0.5
col1 = Column(columnDict, isInPlate=True) # if the plate is inside the plate, set "isInPlate = True"
tutorialModel.addColumn(col1)

# Lastly, a constant load distributed over the plate:
load = Load('area', np.array([-10,0,0]))
tutorialModel.addLoad(load)

# The model is ready! Let's check the geomety and then 
# generate the mesh and see if the finess is
# appropriate or has the be adjusted
plotInputGeometry(tutorialModel)
plt.show()
generateMesh(tutorialModel, showGmshMesh=True)

# The mesh has been generated correctly, let's solve the model.
# The scale for the displacements has been changed to 1e3, in
# order to transform m to mm
solveModel(tutorialModel, resultsScales=(1e3,1,1))
#
# Now the results can be plotted, for example
# Vertical displacements, bending moments and shear in x-direction
plotResults(tutorialModel, ['vDisp', 'Mx', 'Vx'])
plt.show()