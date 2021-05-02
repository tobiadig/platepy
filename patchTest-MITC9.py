#%%
import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress

ConcreteDict = {}
# ConcreteDict["eModule"] = 30e6 #kN/m2
# ConcreteDict["gModule"] =  15e6 #kN/m2
# ConcreteDict["density"] = 2.5 # t/m3
# ConcreteDict["nu"] = 0.0

ConcreteDict["eModule"] = 10920 #kN/m2
ConcreteDict["gModule"] =  10920/(2*1.3) #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.0
C25_30 = Concrete(ConcreteDict)

a=10
b=10
h=0.001
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('line',np.array([-1,0,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
# wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 0]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

patchTestModel.addLoad(lineLoad)

import matplotlib.pyplot as plt
import pandas as pd
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
# generateMesh(patchTestModel, showGmshMesh=True, elementType='QUAD', meshSize=5e-1, order = 'quadratic')
generateMesh(patchTestModel, showGmshMesh=False, elementType='QUAD', nEdgeNodes=3, order='quadratic')

#manually define elements and nodes

# nodesArray = np.array([[0,0], #1
#                         [10, 0], #2
#                         [10,10], #3
#                         [0,10], #4
#                         [2, 2], #5
#                         [8,3], #6
#                         [8,7], #7
#                         [4,7], #8
#                         [5,0], #9
#                         [1,1], #10
#                         [5,1.25], #11
#                         [9,1.5], #12
#                         [5,2.5], #13
#                         [0,5], #14
#                         [1.5,4.75], #15
#                         [3,4.5], #16
#                         [5.5,4.75], #17
#                         [8,5], #18
#                         [9,5], #19
#                         [10,5], #20
#                         [2,8.5], #21
#                         [6,7], #22
#                         [9,8.5], #23
#                         [5.5,8.5], #24
#                         [5,10]]) #25

# elements = np.array([[1,2,6,5,9,12,13,10,11],
#                     [6,2,3,7,12,20,23,18,19],
#                     [8,7,3,4,22,23,25,21,24],
#                     [1,5,8,4,10,16,21,14,15],
#                     [5,6,7,8,13,16,22,18,17]])


# BCs = np.array([[1, 1, 1, 1],
#                 [4, 1, 1, 1],
#                 [14, 1, 1, 1],
#                 [2,0,1,1],
#                 [3,0,1,1],
#                 [5,0,1,1],
#                 [6,0,1,1],
#                 [7,0,1,1],
#                 [8,0,1,1]])



# forces = np.array([[2, -2.5, 0, 0],
#                     [3, -2.5, 0, 0],
#                     [20,-5, 0, 0]])

# forcePattern = Load('nodes', np.array([0,0,0]))
# forcePattern.nodePattern=forces
# patchTestModel.addLoad(forcePattern)

# setMesh(patchTestModel, nodesArray, elements, BCs)

# plotMesh(patchTestModel)
# plt.show()
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'center')

# # print('MITC')
# solveModel(patchTestModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, elemType=elemTypes[1], internalForcePosition = 'nodes')
# display results
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=['x', 'y'])

plt.show()