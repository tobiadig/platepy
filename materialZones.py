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
b2=b1+0.5
h=0.2
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b1], [0,b1],[0,0]])
plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

plateDict["thickness"] = h*3.7
plateDict["outlineCoords"]=np.array([[0,b1], [a,b1], [a,b2], [0,b2],[0,b1]])
unterZug = Plate(plateDict, isUnterZug=True, t=h)

plateDict["thickness"] = h
plateDict["outlineCoords"]=np.array([[0,b2], [a,b2], [a,b], [0,b],[0,b2]])
plate2 = Plate(plateDict)

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b1],[0,b2],[0,b]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 0, 1]))
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[a,0], [a,b1],[a,b2],[a,b]])
wall2 = Wall(wallDict)

firstModel = PlateModel("plateModel1")
firstModel.addPlate(plate1)

firstModel.addPlate(unterZug)
firstModel.addPlate(plate2)
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
generateMesh(firstModel, showGmshMesh=True,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], meshSize=4e-1, order ='linear')
# generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

# compute
from solveModel import *
solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, internalForcePosition = 'center', solveMethod = 'sparse', computeMoments=True)

# display results
plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
plt.show()