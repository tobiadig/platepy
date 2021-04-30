#%%
import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress


def convergenceAnalysis(self, analyticalValue,nEdgeNodes, elemTypes, meshDistortion = False ):
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
            generateMesh(self, showGmshMesh=False, elementType='QUAD', nEdgeNodes=n, order=myOrder, meshDistortion=meshDistortion)
            solveModel(self, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3/1000**4, elemType=elemType, internalForcePosition = 'center')
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

a=10
b=10
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

# analythical vetical displacement, rectangular, simply supported distributed
# all inputs in kN and m
    # pOpts : plate options (e.g. dataclass or named tuple)
    #     pOpts.shape    = "rectangular" | "circular"
    #     pOpts.depth    = "thick" | "thin"
    #     pOpts.support  = "simplySupported" | "clamped"
    #     pOpts.geometry = list of pertinent parameters sufficient to describe geometry
    #     pOpts.material = list of pertinent parameters sufficient to describe material
    # lOpts : load options (e.g. dataclass or named tuple)
    #     lOpts.type      = "concentrated" | "distributed"
    #     lOpts.position  = list of pertinent parameters sufficient to describe position
    #     lOpts.magnitude = magnitude of vertical force
    # sOpts : solution options (e.g. dataclass or named tuple)
    #     sOpts.nTerms = list (e.g. describing the amount of series expansion terms requested)
    # inPos : positions at which output quantities are requested (i.e. array of 2D points)
#%%##################
#ANALYTHICAL SOLUTION
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
print('Maximum displacement, analythical: {}'.format(values[0,0]*1000))

elemTypes = ['MITC4-N', 'MITC9-N']
markerSymbols = ['s', '^']
ax=plt.axes()
outRes = convergenceAnalysis(patchTestModel,-analyticalValue, nEdgeNodes = [3, 5, 9, 17, 33], elemTypes = elemTypes, meshDistortion=True )

nAnalysis = len(elemTypes)
for i in range(0,nAnalysis):
    x=np.log10(outRes[i::nAnalysis,0])
    y=np.log(outRes[i::nAnalysis,1])
    ax.plot(x, y, marker =markerSymbols[i],markerfacecolor = "None", color='k', label = elemTypes[i])
    reg = linregress(x, y)
    print('slope ',elemTypes[i],' : ', reg.slope)
plt.xlabel("LOG(h)")
plt.ylabel("LOG(E)")
# outRes = convergenceAnalysis(patchTestModel,-analyticalValue, nEdgeNodes = [3, 5, 9, 17, 33, 65], elemTypes = ['Q-R'] )
# x=np.log10(outRes[:,0])
# y=np.log(outRes[:,1])
# ax.plot(x, y, marker ="^",markerfacecolor = "None", color='k', label = "Q-R")
# reg = linregress(x, y)
# print('coeff2: ', reg)
ax.legend()
plt.show()