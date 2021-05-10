import numpy as np
from createModel import *
ConcreteDict = {}
ConcreteDict["eModule"] = 32.1*1e6 #kN/m2
ConcreteDict["gModule"] =  14.36*1e6 #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.17
C25_30 = Concrete(ConcreteDict)

distributedLoad = Load('area',np.array([-1, 0, 0]))
a=10
b=10
b1=4.75
b2=b1+0.2
h=0.2
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b1], [0,b1],[0,0]])
plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

plateDict["thickness"] = h*4
plateDict["outlineCoords"]=np.array([[0,b1], [a,b1], [a,b2], [0,b2],[0,b1]])
unterZug = Plate(plateDict)

plateDict["thickness"] = h
plateDict["outlineCoords"]=np.array([[0,b2], [a,b2], [a,b], [0,b],[0,b2]])

plate2 = Plate(plateDict)

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b1],[0,b2],[0,b]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 1]))
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[a,0], [a,b1],[a,b2],[a,b]])
wall2 = Wall(wallDict)

firstModel = PlateModel("plateModel1")
firstModel.addPlate(plate1)
firstModel.addWall(wall1)

firstModel.addPlate(unterZug, isUnterZug=True)
firstModel.addWall(wall2)

firstModel.addPlate(plate2)

firstModel.addLoad(distributedLoad)

#%%
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
import cubusGeometry
from displayModel import *

# create mesh
from generateMesh import *
generateMesh(firstModel, showGmshMesh=False, elementType='QUAD', meshSize=2e-1)
# generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

# compute
from solveModel import *
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3,elementDefinition=elemTypes[1], internalForcePosition = 'nodes')

# display results
plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
plt.show()