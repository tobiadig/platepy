#%%
import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress
from displayModel import *

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
            solveModel(self, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3/1000**4, elementDefinition=elemType, internalForcePosition = 'center')
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

# from AnalyticPlateSolutions import *

# pOpts = POpts()
# pOpts.shape="rectangular"
# pOpts.depth = "thin"
# pOpts.support = "simplySupported"
# pOpts.geometry = (1,1)
# pOpts.material = Material(10920, 0.3, 0.1) #E, nu and h

# lOpts = LOpts()
# lOpts.type = "distributed"
# lOpts.magnitude = 1

# sOpts = SOpts()
# sOpts.nTerms = 40

# inPos=np.array([[0,0]])

# quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
# analyticalValue = values[0,0]*1000
# print('Maximum displacement, analythical: {}'.format(values[0,0]*1000))

# elemTypes = ['MITC4-N', 'MITC9-N']
# # elemTypes = ['L-R', 'Q-R']
# markerSymbols = ['s', '^']
# ax=plt.axes()
# outRes = convergenceAnalysis(patchTestModel,-analyticalValue, nEdgeNodes = [3, 5, 9, 17, 33, 65], elemTypes = elemTypes, meshDistortion=False )

# nAnalysis = len(elemTypes)
# for i in range(0,nAnalysis):
#     x=np.log10(outRes[i::nAnalysis,0])
#     y=np.log(outRes[i::nAnalysis,1])
#     ax.plot(x, y, marker =markerSymbols[i],markerfacecolor = "None", color='k', label = elemTypes[i])
#     reg = linregress(x, y)
#     print('slope ',elemTypes[i],' : ', reg.slope)
# plt.xlabel("LOG(h)")
# plt.ylabel("LOG(E)")

# ax.legend()
# plt.show()


#%% mesh distortion

generateMesh(patchTestModel, showGmshMesh=False, elementType='QUAD', nEdgeNodes=17, order='linear', meshDistortion=False)


# nodesArray = patchTestModel.mesh.nodesArray

# def distortMesh(nodesArray, alpha)
#     myIndex = nodesArray.index.to_numpy()

#     nodesArrayNumpy = nodesArray.to_numpy()
#     v1=np.ones(nodesArrayNumpy.shape[0])
#     x0 = nodesArrayNumpy[:,0]
#     y0 = nodesArrayNumpy[:,1]
#     a=np.max(x0)


#     newNodes = np.zeros(nodesArray.shape)
#     newNodes[:,0] = x0+2*(v1-x0/a)*(2*y0/a-v1)*alpha
#     newNodes[:,1] = y0+2*(v1-y0/a)*(2*x0/a-v1)*alpha

#     xMask = np.logical_or(x0==0, x0==a)
#     yMask = np.logical_or(y0==0, y0==a)
#     newNodes[xMask,0] = x0[xMask]
#     newNodes[yMask,1] = y0[yMask]
#     newNodesArray = pd.DataFrame(newNodes, index = myIndex)
#     return newNodesArray
# print(newNodesArray)
# # patchTestModel.mesh.nodesArray = newNodesArray
plotMesh(patchTestModel)



plt.show()


