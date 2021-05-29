import numpy as np
from createModel import *
from beamComponents import *
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
wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 0, 0]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

columnDict = {}
columnDict["outlineCoords"] = np.array([[0,0]])
columnDict["high"] = 3
columnDict["body"] = C25_30
columnDict["support"] = Support(np.array([1, 0, 0]))
columnDict["crossSection"] = None
columnDict["width"] = 0.05
col1 = Column(columnDict)

columnDict["outlineCoords"] = np.array([[a,0]])
col2 = Column(columnDict)

columnDict["outlineCoords"] = np.array([[a,b]])
col3 = Column(columnDict)

columnDict["outlineCoords"] = np.array([[0,b]])
col4 = Column(columnDict)

beamCoeff = 5
bUZ = 0.05

hUZ = h*(beamCoeff*a/((1-nu**2)*bUZ))**(1/3)

uzCrossSection = CrossSection(bUZ*hUZ*0, bUZ*hUZ**3/12,0, bUZ)
unterZugDict = {}
unterZugDict["body"] = C25_30
unterZugDict["crossSection"] = uzCrossSection
unterZugDict["thickness"] = hUZ

unterZugDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
unterZug1 = downStandBeam(unterZugDict)

# unterZugDict["outlineCoords"] = np.array([[a,b], [a,b]])
# unterZug2 = downStandBeam(unterZugDict)

# unterZugDict["outlineCoords"] = np.array([[0,b], [a,b]])
# unterZug3 = downStandBeam(unterZugDict)

# unterZugDict["outlineCoords"] = np.array([[0,b], [a,b]])
# unterZug4 = downStandBeam(unterZugDict)


firstModel = PlateModel("plateModel1")
firstModel.addPlate(plate1)
firstModel.addColumn(col1)
firstModel.addColumn(col2)
firstModel.addColumn(col3)
firstModel.addColumn(col4)

firstModel.addDownStandBeam(unterZug1)
# firstModel.addWall(wall1)


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
    nEdgeNodes=21, order ='linear')

# generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

# compute
from solveModel import *
solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e6*h**3/a**4, internalForcePosition = 'center', solveMethod = 'cho', computeMoments=False)

# display results
# plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])

startCoord = (0,0.5)
endCoord = (1,0.5)

nEvaluationPoints = 100
bendingMoments, shearForces, arrayEvaluationPoints  = beamComponents(firstModel,'line1', startCoord, endCoord,nEvaluationPoints, integrationWidth = 0, nIntegrationPoints=10)

plotBeamComponent(firstModel,'line1', verticalDisplacement = False, bendingMomentsToPlot = [], shearForcesToPlot = ['x', 'y'], plotOnMesh = True)
plt.show()

