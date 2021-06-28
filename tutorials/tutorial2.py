# ------------------------------------------------------------------------------
#
# platepy tutorial 2
#
# Downstand beams and plotting results over a line ('Schnitte')
#
# ------------------------------------------------------------------------------

# Let's create a model with a square plate supported by two walls on opposite
# Edges, but this time with a downstand beam spanning between the walls.
# Important is to add coordinates in the outline of walls and plates for
#  where the downstand beam will come!

from platepy import *
import numpy as np
tutorialModel = PlateModel()

C30_37 = StandardConcrete("C25_30")

plateDict = {}
# points [10,5] and [0,5] have been added
plateDict["outlineCoords"]=np.array([[0,0], [10,0],[10,5],[10,10], [0,10],[0,5], [0, 0]])
plateDict["thickness"] = 0.3
plateDict["body"]=C30_37
plate = Plate(plateDict)
tutorialModel.addPlate(plate)

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0],[0,5],[0,10]])
wallDict["thickness"] = 0.5 # m
wallDict["support"] = np.array([1, 0, 1])
wall1 = Wall(wallDict)
wallDict["outlineCoords"] = np.array([[10,0],[10,5],[10,10]])
wall2 = Wall(wallDict)
tutorialModel.addWall(wall1)
tutorialModel.addWall(wall2)

# Let's now define a downstand beam spanning between the two walls
dsbDic = {}
dsbDic['outlineCoords'] = np.array([[0,5],[10,5]])
dsbDic['body'] = C30_37
# The cross section object stores information about the dsb:
# Area, Iy, Iz, width and high (m and mm^4)
dsbDic['crossSection'] = CrossSection(0.3*0.3, 0.3*4/12, 0, 0.3, 0.3)
dsb = DownStandBeam(dsbDic)
tutorialModel.addDownStandBeam(dsb)

# Add the load:
load = Load('area', np.array([-10,0,0]))
tutorialModel.addLoad(load)

# As in tutorial 1, let's check the geometry, generate the mesh and solve
plotInputGeometry(tutorialModel)
plt.show()
generateMesh(tutorialModel)
solveModel(tutorialModel, resultsScales=(1e3,1,1))

# Instead of plotting the results with isolines, let's inspect a vertical cut
# for Vertical displacements, bending moments and shear in y-direction

# The cut is defined by starting and end coordinates, the number of points to be evaluated and 
# the scale of the result
computeBeamComponents(tutorialModel, (5,0), (5,10), 150, resultsScales=(1e3,1,1))

# the results can be now plotted with the function plotBeamComponent
# the suffix '_line' at the end of the desired plot has to be added
plotBeamComponent(tutorialModel, ['vDisp_line', 'My_line', 'Vy_line'])
plt.show()