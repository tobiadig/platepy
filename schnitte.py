#%%
import numpy as np
from createModel import *
from beamComponents import*
ConcreteDict = {}
ConcreteDict["eModule"] = 32.1*1e6 #kN/m2
ConcreteDict["gModule"] =  14.36*1e6 #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.17
C25_30 = Concrete(ConcreteDict)

distributedLoad = Load('area',np.array([-1, 0, 0]))
a=10
b=10
h=0.1

plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0],[a,b], [0,b],[0,0]])
plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)


wallDict = {}
wallDict["outlineCoords"]=np.array([[0,0], [a,0],[a,b], [0,b],[0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 0, 0]))
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)


firstModel = PlateModel("plateModel1")
firstModel.addPlate(plate1)


firstModel.addWall(wall1)



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
generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], nEdgeNodes=4, order ='linear')
# generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

# compute
from solveModel import *
solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, internalForcePosition = 'center', solveMethod = 'sparse', computeMoments=True)

# display results
plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=False, bendingMomentsToPlot=[],shearForcesToPlot=[])


#%%
startCoord = (1.67,5)
endCoord = (8.33,5)
nEvaluationPoints = 10
bendingMoments, shearForces, arrayEvaluationPoints  = beamComponents(firstModel,'line1', startCoord, endCoord,nEvaluationPoints)

plotBeamComponent(firstModel,'line1', verticalDisplacement = False, bendingMomentsToPlot = ['x'], shearForcesToPlot = [])

plt.show()