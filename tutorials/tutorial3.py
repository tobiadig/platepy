# ------------------------------------------------------------------------------
#
# platepy tutorial 3
#
# Manually create a simple mesh
#
# ------------------------------------------------------------------------------


#The present tutorial shows how to manually define nodes, element, boundary conditions
# and nodal loads. This can be useful when a simple, yet precise defined mesh is required
# For exampl for a patch test.

# Let's prepare the model and the material
import numpy as np
from platepy import *

tutorialModel = PlateModel()
C25_30 = StandardConcrete("C25_30")

# the plate has still to be assigned to the model in order to specify the material properties
# of the elements
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [10,0]]) # the outline is irrelevant, since 
                                                        # the shape will be manually defined later
plateDict["thickness"] = 0.1
plateDict["body"]=C25_30
plate = Plate(plateDict)
tutorialModel.addPlate(plate)

# now, instead of defining loads, structural components and generate the mesh, all steps are
# performed manually.

# Firstly, the array of nodes coordinates has to be defined.
nodesArray = np.array([[0,0],    # Node 1
                        [10, 0], # Node 2
                        [10,10], # Node 3
                        [0,10],  # Node 4
                        [2, 2],  # Node 5
                        [8,3],   # Node 6
                        [8,7],   # Node 7
                        [4,7]])  # Node 8

#Then the elements connectivity
elements = np.array([[1,2,6,5], # Element 1
                    [2,3,7,6],  # Element 2
                    [7,3,4,8],  # Element 3
                    [1,5,8,4],  # Element 4
                    [5,6,7,8]]) # Element 5

#Then, the boundary condition. The first column is the node of referenc, the other columns
# define if ther relative DOF is blocked or free (vDisp/xRot/yRot). 
# Let's clamp the left-hand side and block all rotations.
BCs = np.array([[1, 1, 1, 1],
                [4, 1, 1, 1],
                [2, 0, 1, 1],
                [3, 0, 1, 1],
                [5, 0, 1, 1],
                [6, 0, 1, 1],
                [7, 0, 1, 1],
                [8, 0, 1, 1]])

# The nodal forces have to be defined. Let's put two vertical loads on the right-hand side.
forces = np.array([[2, -5, 0, 0],
                    [3, -5, 0, 0]])

# A load object is created, then the "nodePattern" attribute is assigned with 
# the above defined array of nodal forces
forcePattern = Load('nodes', np.array([0,0,0])) # the magnitude is not used
forcePattern.nodePattern=forces
tutorialModel.addLoad(forcePattern)

# nodes, elements and boundary conditions are assigned to the model with the setMesh function
# differently from the generateMesh function, there is no default element type
setMesh(tutorialModel, nodesArray, elements, BCs,  elementDefinition='MITC-4-N')

solveModel(tutorialModel, resultsScales=(1e6,1,1))

# Let's now look at the displacement in the individual nodes with the "tect+mesh" plot type
plotResults(tutorialModel,plotType='text+mesh',valuesToPlotList=['vDisp'])
plt.show()