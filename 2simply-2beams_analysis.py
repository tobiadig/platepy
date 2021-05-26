#%%
import numpy as np
from createModel import *

import matplotlib.pyplot as plt
import pandas as pd

pd.set_option("display.max_rows", None, "display.max_columns", None)
import cubusGeometry
from displayModel import *
from generateMesh import *
from solveModel import *  
from evaluateAtPoints import *

def twoSimpy2Beams(beamCoeff):
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
    wallDict["support"] = Support(np.array([1, 0, 0]))
    wallDict["thickness"] = 0.5 # m
    wall1 = Wall(wallDict)

    wallDict["outlineCoords"] = np.array([[a,0], [a,b]])
    wall2 = Wall(wallDict)

    bUZ = 0.01

    hUZ = h*(beamCoeff*a/((1-nu**2)*bUZ))**(1/3)
    # hUZ = h*((beamCoeff*a/((1-nu**2)*bUZ))**(1/3)-1)

    uzCrossSection = CrossSection(bUZ*hUZ, bUZ*hUZ**3/12,0, bUZ)
    unterZugDict = {}
    # E = 32e6
    # ConcreteDict["eModule"] = E #kN/m2
    # ConcreteDict["gModule"] =   E/(2*(1+nu)) #kN/m2
    # C25_30 = Concrete(ConcreteDict)

    unterZugDict["body"] = C25_30
    unterZugDict["crossSection"] = uzCrossSection
    unterZugDict["thickness"] = hUZ

    unterZugDict["outlineCoords"] = np.array([[0,0], [a,0]])
    unterZug1 = downStandBeam(unterZugDict)

    unterZugDict["outlineCoords"] = np.array([[0,b], [a,b]])
    unterZug2 = downStandBeam(unterZugDict)


    firstModel = PlateModel("plateModel1")
    firstModel.addPlate(plate1)
    firstModel.addWall(wall1)
    firstModel.addWall(wall2)

    firstModel.addDownStandBeam(unterZug1)
    firstModel.addDownStandBeam(unterZug2)

    firstModel.addLoad(distributedLoad)

    #%%

    # import sample plate and show it

    # create mesh

    elemDefinitions = ['DB-4-R', 'MITC-4-N', 'DB-9-R', 'MITC-9-N']
    generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], nEdgeNodes=21, order ='linear')

    # generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

    # compute

    solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e6*h**3/a**4, internalForcePosition = 'center', solveMethod = 'cho', computeMoments=False)

    # display results
    # plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
    myCoords=np.array([[a/2,b/2]])
    verticalDisplacements,bendingMoments, shearForces  = evaluateAtPoints(firstModel, myCoords, displayPoints = False)

    outRes = verticalDisplacements
    return outRes

#%%
nPoints = 50
beamCoeffMin = 0
beamCoeffMax =10

beamCoeffarrayBeam = np.linspace(beamCoeffMin,beamCoeffMax,nPoints)
wMaxArrayBeam = np.zeros(len(beamCoeffarrayBeam))
wMaxArrayAnalythical = np.zeros(len(beamCoeffarrayBeam))
#%%
for i,value in enumerate(beamCoeffarrayBeam):
    wMaxArrayBeam[i] = twoSimpy2Beams(value)

#%%  import analytical values
from AnalyticPlateSolutions import *

pOpts = POpts()
pOpts.shape="rectangular"
pOpts.depth = "thin"
pOpts.support = "sSupportwBeams"
pOpts.geometry = (1,1)
pOpts.material = Material(10920, 0.3, 0.1) #E, nu and h


lOpts = LOpts()
lOpts.type = "distributed"
lOpts.magnitude = 1

sOpts = SOpts()
sOpts.nTerms = 40

inPos=np.array([[0.5,0]])
for i,value in enumerate(beamCoeffarrayBeam):
    pOpts.material.myLambda = value
    quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
    wMaxArrayAnalythical[i] = -values[0,0]*1000
# print('Maximum displacement, analythical: {}'.format(values[0,0]*1000))

#%%
ax=plt.axes()
# ax1 = ax

ax.plot(beamCoeffarrayBeam,wMaxArrayBeam, marker ='s',markerfacecolor = "None", color='k', label = 'numerical')
ax.plot(beamCoeffarrayBeam,wMaxArrayAnalythical, marker ='^',markerfacecolor = "None", color='k', label = 'analythical')
# ax.loglog(beamCoeffarrayBeam[3:],-wMaxArrayBeam[3:], marker ='s',markerfacecolor = "None", color='k', label = 'numerical')
# ax.loglog(beamCoeffarrayBeam[3:],-wMaxArrayAnalythical[3:], marker ='^',markerfacecolor = "None", color='k', label = 'analythical')

plt.xlabel('lambda [-]')
plt.ylabel('w max [mm]')
ax.legend()

# E=np.abs(wMaxArrayBeam-wMaxArrayAnalythical)/wMaxArrayAnalythical*-1


# ax.plot(beamCoeffarrayBeam,E,marker = '>',markerfacecolor = "None", color='k', label = 'Error')
# plt.xlabel('lambda [-]')
# plt.ylabel('Error [-]')
# ax.legend()



plt.show()

#%%
print(wMaxArrayBeam[-1])
print(wMaxArrayAnalythical[-1])