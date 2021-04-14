#%%
import numpy as np
from createModel import *


ConcreteDict = {}
ConcreteDict["eModule"] = 30*1e9 #kN/m2
ConcreteDict["gModule"] =  15*1e9 #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0
C25_30 = Concrete(ConcreteDict)



a=10
b=10
h=0.01
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
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
# generateMesh(patchTestModel, showGmshMesh=False, elementType='QUAD', meshSize=5e-1)
generateMesh(patchTestModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=3)

# compute
from solveModel import *
solveModel(patchTestModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e1)

# display results
plotResults(patchTestModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
plt.show()