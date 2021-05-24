#%%
import numpy as np
from createModel import *
import matplotlib.pyplot as plt
import pandas as pd
from displayModel import *
from generateMesh import *
from solveModel import *
from tqdm import tqdm
def uz_high_sensitivity(hUZ, bUZ):
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
    # wallDict["outlineCoords"] = np.array([[0,0],[a,0],[a,0.5*b], [a,b],[0,b],[0,0.5*b],[0,0]])
    wallDict["high"] = 3 # m
    wallDict["body"] = C25_30
    wallDict["support"] = Support(np.array([1, 0, 0]))
    wallDict["thickness"] = 0.5 # m
    wall1 = Wall(wallDict)

    wallDict["outlineCoords"] = np.array([[a,0],[a,0.5*b], [a,b]])
    wall2 = Wall(wallDict)

    wallDict["outlineCoords"] = np.array([[0,0.5*b], [a,0.5*b]])
    wall3 = Wall(wallDict)


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

    pd.set_option("display.max_rows", None, "display.max_columns", None)

    # import sample plate and show it


    # create mesh

    elemDefinitions = ['DB-4-R', 'MITC-4-N', 'DB-9-R', 'MITC-9-N']
    generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], meshSize=6e-1, order ='linear')

    # generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')
    # compute

    solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, internalForcePosition = 'center', solveMethod = 'cho', computeMoments=False)

    # display results
    # plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])

    outRes = firstModel.results.wMax[2]
    return outRes

def uz_materialZones_high_sensitivity(hUZ,bUZ):

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

    plateDict["thickness"] = hUZ
    plateDict["outlineCoords"]=np.array([[0,b1], [a,b1], [a,b2], [0,b2],[0,b1]])
    unterZug = Plate(plateDict, isUnterZug=True, t=h)

    plateDict["thickness"] = h
    plateDict["outlineCoords"]=np.array([[0,b2], [a,b2], [a,b], [0,b],[0,b2]])

    plate2 = Plate(plateDict)

    wallDict = {}
    wallDict["outlineCoords"] = np.array([[0,0], [0,b1],[0,b2],[0,b]])
    wallDict["high"] = 3 # m
    wallDict["body"] = C25_30
    wallDict["support"] = Support(np.array([1, 1, 0]))
    wallDict["thickness"] = 0.5 # m
    wall1 = Wall(wallDict)

    wallDict["outlineCoords"] = np.array([[a,0], [a,b1],[a,b2],[a,b]])
    wall2 = Wall(wallDict)

    firstModel = PlateModel("plateModel1")
    firstModel.addPlate(plate1)
    firstModel.addWall(wall1)

    firstModel.addPlate(unterZug)
    firstModel.addWall(wall2)
    firstModel.addPlate(plate2)
    firstModel.addLoad(distributedLoad)


    elemDefinitions = ['DB-4-R', 'MITC-4-N', 'DB-9-R', 'MITC-9-N']
    generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], meshSize=4e-1, order ='linear')
    solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, internalForcePosition = 'center', solveMethod = 'sparse', computeMoments=True)
    outRes = firstModel.results.wMax[2]
    return outRes

#%%
nPoints = 10
hMin = 0
hMax = 1.4
bUZ = 0.5
h=0.2
hUZarrayBeam = np.linspace(hMin,hMax,nPoints)
wMaxArrayBeam = np.zeros(len(hUZarrayBeam))
#%%
for i,hUZ in tqdm(enumerate(hUZarrayBeam)):
    wMaxArrayBeam[i] = uz_high_sensitivity(hUZ, bUZ)

#%%
wMaxArrayPlate = np.zeros(nPoints)
hUZarray = np.linspace(hMin+h,hMax+h,nPoints)
for i,hUZ in tqdm(enumerate(hUZarray)):
    wMaxArrayPlate[i] = uz_materialZones_high_sensitivity(hUZ, bUZ)
effectiveH = hUZarray - np.ones(len(hUZarray))*0.2


#%%  import cubus analysis
import pandas as pd
myPath = 'C:/Users/Diggelmann/Desktop/comparing uzs/'
data = np.loadtxt(myPath+'cubus uz analysis.csv', delimiter=',')
xCubus = data[:,1]
yCubus = data[:,2]



#%%
ax=plt.axes()
ax.plot(hUZarrayBeam,wMaxArrayBeam, marker ='s',markerfacecolor = "None", color='k', label = 'beam')
ax.plot(effectiveH,wMaxArrayPlate, marker ='^',markerfacecolor = "None", color='k', label = 'plate')
ax.plot(xCubus, yCubus, marker='^',color='g', label = 'cubus' )
plt.xlabel('Effective high [m]')
plt.ylabel('w max [mm]')
ax.legend()

plt.show()

#%%





