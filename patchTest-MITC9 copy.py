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
ConcreteDict["nu"] = 0.3
C25_30 = Concrete(ConcreteDict)

a=10
b=10
h=0.001
h=1
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('area',np.array([-1,0,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
# wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 0, 1]))
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
generateMesh(patchTestModel, showGmshMesh=False, elementType='QUAD', nEdgeNodes=31, order='quadratic')
#manually define elements and nodes

# nodesArray = np.array([[0,0], #1
#                         [10, 0], #2
#                         [10,10], #3
#                         [0,10], #4
#                         [7, 0], #5
#                         [10,5], #6
#                         [5,10], #7
#                         [0,5], #8
#                         [5,5]]) #25

# elements = np.array([[1,2,3,4,5,6,7,8,9]])

# BCs = np.array([[1, 1, 1, 1],
#                 [4, 1, 1, 1],
#                 [8, 1, 1, 1]])

# forces = np.array([[2, 0, 0, 2.5],
#                     [3, 0, 0, 2.5],
#                     [6,0, 0, 5]])

# forcePattern = Load('nodes', np.array([0,0,0]))
# forcePattern.nodePattern=forces
# patchTestModel.addLoad(forcePattern)

# setMesh(patchTestModel, nodesArray, elements, BCs)

# plotMesh(patchTestModel)
# # plt.show()
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e2, elementDefinition=elemTypes[3], internalForcePosition = 'nodes')

# # print('MITC')
# solveModel(patchTestModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, elemType=elemTypes[1], internalForcePosition = 'nodes')
# display results
plotResults(patchTestModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])

plt.show()