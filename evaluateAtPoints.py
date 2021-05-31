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