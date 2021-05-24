#%%
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
h=0.2

plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0],[a,0.5*b], [0,0.5*b],[0,0]])
plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

plateDict["outlineCoords"]=np.array([[a,0.5*b], [a,b], [0,b],[0,0.5*b],[a,0.5*b]])
plate2 = Plate(plateDict)


wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0],[0,0.5*b], [0,b]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 0]))
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[a,0],[a,0.5*b], [a,b]])
wall2 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[0,0.5*b], [a,0.5*b]])
wall3 = Wall(wallDict)

bUZ = 0.5
hUZ = 0.7*h
uzCrossSection = CrossSection(bUZ*hUZ, bUZ*hUZ**3/12,0, bUZ)
unterZugDict = {}
unterZugDict["outlineCoords"] = np.array([[0,0.5*b], [a,0.5*b]])
unterZugDict["body"] = C25_30
unterZugDict["crossSection"] = uzCrossSection
unterZugDict["thickness"] = 0.5
unterZug = downStandBeam(unterZugDict)

firstModel = PlateModel("plateModel1")
firstModel.addPlate(plate1)
firstModel.addPlate(plate2)

firstModel.addWall(wall1)
firstModel.addWall(wall2)
# firstModel.addWall(wall3)
firstModel.addDownStandBeam(unterZug)

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
generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], meshSize=5e-1, order ='linear')
# generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

# compute
from solveModel import *
solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, internalForcePosition = 'center', solveMethod = 'cho', computeMoments=True)

# display results
plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])

#%%

M = firstModel.results.bendingMomentsDSB
coords = firstModel.results.internalForcesPositionsDSB
V = firstModel.results.shearForcesDSB
N=firstModel.results.normalForcesDSB

ax1=plt.axes()
ax1.plot(coords[:,0], M)
ax1.set_title('M')
plt.show()

ax2=plt.axes()
ax2.plot(coords[:,0],V)
ax2.set_title('V')
plt.show()

ax2=plt.axes()
ax2.plot(coords[:,0],N)
ax2.set_title('N')
plt.show()

