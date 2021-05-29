import numpy as np
from createModel import *
E=10920
nu=0.3
ConcreteDict = {}
ConcreteDict["eModule"] = E #kN/m2
ConcreteDict["gModule"] =   E/(2*(1+nu)) #kN/m2
ConcreteDict["density"] = 2.5 # t/m3

ConcreteDict["nu"] = nu
C25_30 = Concrete(ConcreteDict)

distributedLoad = Load('area',np.array([-1, 0, 0]))
a=10
b=10
h=0.01
bUZ = a/20

plateDict = {}
plateDict["outlineCoords"]=np.array([[0,bUZ], [a,bUZ],[a,b-bUZ], [0,b-bUZ],[0,bUZ]])
plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

plateDict["thickness"] = h*3
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,bUZ], [0,bUZ],[0,0]])
unterZug1 = Plate(plateDict, isUnterZug=True, t=h)

plateDict["outlineCoords"]=np.array([[0,b-bUZ], [a,b-bUZ], [a,b], [0,b],[0,b-bUZ]])
unterZug2 = Plate(plateDict, isUnterZug=True, t=h)

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,bUZ],[0,b-bUZ],[0,b]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = np.array([1, 1, 0])
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[a,0], [a,bUZ],[a,b-bUZ],[a,b]])
wall2 = Wall(wallDict)


firstModel = PlateModel()
firstModel.addPlate(plate1)
firstModel.addPlate(unterZug1)
firstModel.addPlate(unterZug2)
firstModel.addWall(wall1)
firstModel.addWall(wall2)


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
elemDefinitions = ['DB-4-R', 'MITC-4-N', 'DB-9-R', 'MITC-9-N']
generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1],meshSize=8e-1, order ='linear')

# generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

# compute
from solveModel import *
solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e7*h**3/a**4, internalForcePosition = 'center', solveMethod = 'sparse', computeMoments=False)

# display results
plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
plt.show()