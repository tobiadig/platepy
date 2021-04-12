#%%
import numpy as np
from createModel import *


ConcreteDict = {}
ConcreteDict["eModule"] = 32.1*1e6 #kN/m2
ConcreteDict["gModule"] =  14.36*1e6 #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0
C25_30 = Concrete(ConcreteDict)

lineLoad = Load('line',np.array([-1,0,0]))
lineLoad.outlineCoords=np.array([[1,0], [1,1]])

a=1
b=1
h=0.01
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])
plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 1]))
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

patchTestModel.addLoad(lineLoad)


#%%
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it

from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
generateMesh(patchTestModel, showGmshMesh=True, elementType='QUAD', meshSize=5e-1)
# generateMesh(sampleModel, showGmshMesh=False, elementType='TRI', nEdgeNodes=51)

# compute
from solveModel import *
solveModel(patchTestModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3)

# display results
plotResults(patchTestModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
plt.show()