#%%
import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress
from displayModel import *
import os

def convergenceAnalysis(self, analyticalValue,nEdgeNodes, elemTypes, meshDistortion = False, distVal=100 ):
    myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\convergence\\"
    nVal = len(nEdgeNodes)*len(elemTypes)
    outRes = np.zeros((nVal,2))
    k=0
    for n in nEdgeNodes:
        for elemType in elemTypes:

            self.clearMesh()
            t = elemType.split('-')[0]
            if t =='L' or t=='MITC4':
                myOrder = 'linear'
            elif t == 'Q'or t=='MITC9':
                myOrder = 'quadratic'
            else:
                raise TypeError('something went wrong')
            generateMesh(self, showGmshMesh=False, elementType='QUAD', nEdgeNodes=n, order=myOrder, meshDistortion=meshDistortion, distVal=distVal)
            # if meshDistortion and (n==9):
            #     plotMesh(self, plotStrucElements=False, plotNodes=False, plotPoints=True)
            #     figRatio = 1
            #     figSize = 5
            #     currentFig = plt.gcf()
            #     currentAx = plt.gca()
            #     currentFig.set_size_inches(figRatio*figSize,1*figSize)
            #     # currentAx.set_xlim(xlim1, xlim2)
            #     # currentAx.set_ylim(ylim1, ylim2)
            #     plt.axis('off')
            #     imgName = 'mesh'+t+ '.svg'
            #     savingPath = os.path.join(myPath, imgName)
            #     plt.savefig(savingPath)
            
            # if n<30:
            #     plotMesh(self)
            #     plt.show()
            solveModel(self, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3/1000**4, elementDefinition=elemType, internalForcePosition = 'center')
            # plotResults(patchTestModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
            # plt.show()
            outRes[k,1] = np.abs(patchTestModel.results.wMax[2]-analyticalValue)/np.abs(analyticalValue)
            print('E: ', outRes[k,1])
            print('anal: ', analyticalValue)
            print(patchTestModel.results.wMax[2])
            outRes[k,0] = 1/n
            k+=1

    return outRes

ConcreteDict = {}
# ConcreteDict["eModule"] = 30e6 #kN/m2
# ConcreteDict["gModule"] =  15e6 #kN/m2
# ConcreteDict["density"] = 2.5 # t/m3
# ConcreteDict["nu"] = 0.0

ConcreteDict["eModule"] = 10920 #kN/m2
ConcreteDict["gModule"] =  10920/(2*1.3) #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.3
C25_30 = Concrete(ConcreteDict)

a=1000
b=1000
h=0.1
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('area',np.array([-1,-1*0,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
# wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 0, 1]))
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

patchTestModel.addLoad(lineLoad)

from AnalyticPlateSolutions import *

pOpts = POpts()
pOpts.shape="rectangular"
pOpts.depth = "thin"
pOpts.support = "simplySupported"
pOpts.geometry = (1,1)
pOpts.material = Material(10920, 0.3, 0.1) #E, nu and h

lOpts = LOpts()
lOpts.type = "distributed"
lOpts.magnitude = 1

sOpts = SOpts()
sOpts.nTerms = 40

inPos=np.array([[0,0]])

quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
analyticalValue = values[0,0]*1000


# elemTypes = ['MITC4-N', 'MITC9-N']
elemTypes = ['L-R', 'Q-R']
markerSymbols = ['s', '^']




myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\convergence\\"
figSize=5
figRatio = 1/1.7
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1



#%% L-Q regular
elemTypes = ['L-R', 'Q-R']
patchTestModel.walls[0].support=Support(np.array([1, 1, 0]))
ax=plt.axes()
outRes = convergenceAnalysis(patchTestModel,-analyticalValue, nEdgeNodes = [3, 5, 9, 17, 33, 65, 129], elemTypes = elemTypes, meshDistortion=False, distVal = 120 )

nAnalysis = len(elemTypes)
for i in range(0,nAnalysis):
    x=np.log10(outRes[i::nAnalysis,0])
    y=np.log10(outRes[i::nAnalysis,1])
    ax.plot(x, y, marker =markerSymbols[i],markerfacecolor = "None", color='k', label = elemTypes[i])
    reg = linregress(x, y)
    print('slope ',elemTypes[i],' : ', reg.slope)
plt.xlabel("LOG(h)")
plt.ylabel("LOG(E)")
ax.legend()

currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
imgName = 'L-Q regular.svg'
currentFig.subplots_adjust(left = 0.15)
plt.tight_layout()

savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

tableName =imgName[0:-3]+'csv'
savingPath = os.path.join(myPath, tableName)
myTable = np.array([x,y])
np.savetxt(savingPath, myTable, delimiter = ',')
#%%L-Q distorted
patchTestModel.walls[0].support=Support(np.array([1, 1, 0]))
elemTypes = ['L-R', 'Q-R']
ax=plt.axes()
outRes = convergenceAnalysis(patchTestModel,-analyticalValue, nEdgeNodes = [3, 5, 9, 17, 33, 65, 129], elemTypes = elemTypes, meshDistortion=True, distVal = 120 )

nAnalysis = len(elemTypes)
for i in range(0,nAnalysis):
    x=np.log10(outRes[i::nAnalysis,0])
    y=np.log10(outRes[i::nAnalysis,1])
    ax.plot(x, y, marker =markerSymbols[i],markerfacecolor = "None", color='k', label = elemTypes[i])
    reg = linregress(x, y)
    print('slope ',elemTypes[i],' : ', reg.slope)
plt.xlabel("LOG(h)")
plt.ylabel("LOG(E)")
ax.legend()

currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
imgName = 'L-Q distorted.svg'
plt.tight_layout()
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
tableName =imgName[0:-3]+'csv'
savingPath = os.path.join(myPath, tableName)
myTable = np.array([x,y])
np.savetxt(savingPath, myTable, delimiter = ',')
#%% MITC4-9 regular
elemTypes = ['MITC4-N', 'MITC9-N']
patchTestModel.walls[0].support=Support(np.array([1, 0, 1]))
ax=plt.axes()
outRes = convergenceAnalysis(patchTestModel,-analyticalValue, nEdgeNodes = [3, 5, 9, 17, 33, 65], elemTypes = elemTypes, meshDistortion=False, distVal = 120 )

nAnalysis = len(elemTypes)
for i in range(0,nAnalysis):
    x=np.log10(outRes[i::nAnalysis,0])
    y=np.log10(outRes[i::nAnalysis,1])
    ax.plot(x, y, marker =markerSymbols[i],markerfacecolor = "None", color='k', label = elemTypes[i])
    reg = linregress(x, y)
    print('slope ',elemTypes[i],' : ', reg.slope)
plt.xlabel("LOG(h)")
plt.ylabel("LOG(E)")
ax.legend()

currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
imgName = 'MITC4-9 regular.svg'
plt.tight_layout()
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
tableName =imgName[0:-3]+'csv'
savingPath = os.path.join(myPath, tableName)
myTable = np.array([x,y])
np.savetxt(savingPath, myTable, delimiter = ',')
#%%MITC4-9 distorted
elemTypes = ['MITC4-N', 'MITC9-N']
patchTestModel.walls[0].support=Support(np.array([1, 0, 1]))
ax=plt.axes()
outRes = convergenceAnalysis(patchTestModel,-analyticalValue, nEdgeNodes = [3, 5, 9, 17, 33, 65], elemTypes = elemTypes, meshDistortion=True, distVal = 120 )

nAnalysis = len(elemTypes)
for i in range(0,nAnalysis):
    x=np.log10(outRes[i::nAnalysis,0])
    y=np.log10(outRes[i::nAnalysis,1])
    ax.plot(x, y, marker =markerSymbols[i],markerfacecolor = "None", color='k', label = elemTypes[i])
    reg = linregress(x, y)
    print('slope ',elemTypes[i],' : ', reg.slope)
plt.xlabel("LOG(h)")
plt.ylabel("LOG(E)")
ax.legend()

currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
imgName = 'MITC4-9 distorted.svg'
plt.tight_layout()
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
tableName =imgName[0:-3]+'csv'
savingPath = os.path.join(myPath, tableName)
myTable = np.array([x,y])
np.savetxt(savingPath, myTable, delimiter = ',')