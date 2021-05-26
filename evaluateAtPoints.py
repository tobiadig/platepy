from copy import Error

import gmsh
from shapeFunctions import *
from scipy.integrate import trapezoid
import numpy as np
from createModel import *
from beamComponents import*
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option("display.max_rows", None, "display.max_columns", None)
import cubusGeometry
from displayModel import * 
from generateMesh import *
from solveModel import *   
def evaluateAtPoints(self, coords, displayPoints = False):

    uGlob = self.results.uGlobPlate
    # print(uGlob)
    nEvaluationPoints = coords.shape[0]
    elementsList = self.mesh.elementsList
    arrayEvaluationPoints = np.zeros((nEvaluationPoints,2))

    bendingMoments = np.zeros((nEvaluationPoints,3))
    shearForces = np.zeros((nEvaluationPoints,2))
    verticalDisplacements = np.zeros(nEvaluationPoints)
    resultsScaleVertDisp = self.results.resultsScaleVertDisp

    for k,evaluationPoint in enumerate(coords):
        elementTaggetEl, elementTypegetEl, nodeTagsgetEl, ugetEl, vgetEl, wgetEl = gmsh.model.mesh.getElementByCoordinates(evaluationPoint[0],evaluationPoint[1],0, dim=2)
        # print('element: ', elementTaggetEl)
        element = self.mesh.getElementByTagDictionary[elementTaggetEl]
        plateOfTheElement = element.whichPlate
        elementType = element.type
        elemNodes = element.connectivity
        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        nNodes=element.nNodes
        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        elementShape = len(xi)
        kCoeff = np.zeros((3*nNodes),dtype=int)
        for i in range(0,3):
            kCoeff[0+i::3]=coherentElemNodes*3+i
        vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
        # print('vLoc: ', vLoc)
        Df = self.plates[plateOfTheElement].Df
        Dc = self.plates[plateOfTheElement].Dc

        ri =ugetEl
        si = vgetEl

        N, Bb,Bs, detJ =getShapeFunctionForElementType(elementType,ri, si, xi, yi)

        tempDispl = N@vLoc

        verticalDisplacements[k] = tempDispl[0]*resultsScaleVertDisp
        bendingMoments[k,0:3] = np.matmul(Df,np.matmul(Bb, vLoc))[:,0]*1
        shearForces[k,0:2] = np.matmul(Dc, np.matmul(Bs, vLoc))[:,0]*1

    if displayPoints:
        fig,outAx = plotInputGeometry(self)
        for k,evaluationPoint in enumerate(coords):
            outAx.scatter(evaluationPoint[0],evaluationPoint[1], facecolor='r', marker='.')



    return verticalDisplacements,bendingMoments, shearForces

        
if __name__ == "__main__":

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
    wallDict["support"] = Support(np.array([1, 0, 1]))
    wallDict["thickness"] = 0.5 # m
    wall1 = Wall(wallDict)


    firstModel = PlateModel("plateModel1")
    firstModel.addPlate(plate1)


    firstModel.addWall(wall1)



    firstModel.addLoad(distributedLoad)

    #%%

    # import sample plate and show it

    # create mesh

    elemDefinitions = ['DB-4-R', 'MITC-4-N', 'DB-9-R', 'MITC-9-N']
    generateMesh(firstModel, showGmshMesh=False,showGmshGeometryBeforeMeshing=False, elementDefinition=elemDefinitions[1], nEdgeNodes=6, order ='linear', deactivateRotation=False)
    # generateMesh(sampleModel, showGmshMesh=True, elementType='QUAD', nEdgeNodes=11, order ='linear')

    # compute
    
    solveModel(firstModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3, internalForcePosition = 'center', solveMethod = 'sparse', computeMoments=True)
    plotResults(firstModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
    coords=np.array([[5,5]])
    verticalDisplacements,bendingMoments, shearForces  = evaluateAtPoints(firstModel, coords, displayPoints = True)
    print(verticalDisplacements.shape)
    plt.show()