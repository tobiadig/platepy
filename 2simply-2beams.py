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
a=1
b=1
h=0.001

plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0],[a,b], [0,b],[0,0]])
plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)


wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = np.array([1, 0, 0])
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[a,0], [a,b]])
wall2 = Wall(wallDict)


beamCoeff = 2
bUZ = 0.01


hUZ = h*(beamCoeff*a/((1-nu**2)*bUZ))**(1/3)
# hUZ = h*((beamCoeff*a/((1-nu**2)*bUZ))**(1/3)-1)

# print(hUZ)
uzCrossSection = CrossSection(bUZ*hUZ*0, bUZ*hUZ**3/12,0, bUZ, hUZ)
unterZugDict = {}
# E = 32e6
# ConcreteDict["eModule"] = E #kN/m2
# ConcreteDict["gModule"] =   E/(2*(1+nu)) #kN/m2
# C25_30 = Concrete(ConcreteDict)

unterZugDict["body"] = C25_30
unterZugDict["crossSection"] = uzCrossSection
unterZugDict["thickness"] = hUZ

unterZugDict["outlineCoords"] = np.array([[0,0], [a,0]])
unterZug1 = DownStandBeam(unterZugDict)

unterZugDict["outlineCoords"] = np.array([[0,b], [a,b]])
unterZug2 = DownStandBeam(unterZugDict)


firstModel = PlateModel()
firstModel.addPlate(plate1)
firstModel.addWall(wall1)
firstModel.addWall(wall2)

firstModel.addDownStandBeam(unterZug1)
firstModel.addDownStandBeam(unterZug2)

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
generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], \
    nEdgeNodes=21, order ='linear',deactivateRotation=True)



# compute
from solveModel import *
solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e6*h**3/a**4, internalForcePosition = 'center', solveMethod = 'cho', computeMoments=False)

# display results
plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
plt.show()