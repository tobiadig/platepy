from copy import Error
import numpy as np
import gmsh
from shapeFunctions import *

def beamComponents(self,lineName, startCoord, endCoord, nEvaluationPoints):

    uGlob = self.results.uGlobPlate
    # print(uGlob)

    elementsList = self.mesh.elementsList
    arrayEvaluationPoints = np.zeros((nEvaluationPoints,2))

    arrayEvaluationPoints[:,0] = np.linspace(startCoord[0], endCoord[0], num=nEvaluationPoints)
    arrayEvaluationPoints[:,1] = np.linspace(startCoord[1], endCoord[1], num=nEvaluationPoints)
    bendingMoments = np.zeros((nEvaluationPoints,3))
    shearForces = np.zeros((nEvaluationPoints,2))

    for k,evaluationPoint in enumerate(arrayEvaluationPoints):
        
        elementTaggetEl, elementTypegetEl, nodeTagsgetEl, ugetEl, vgetEl, wgetEl = gmsh.model.mesh.getElementByCoordinates(evaluationPoint[0],evaluationPoint[1],0)
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
        bendingMoments[k,0:3] = np.matmul(Df,np.matmul(Bb, vLoc))[:,0]*1
        shearForces[k,0:2] = np.matmul(Dc, np.matmul(Bs, vLoc))[:,0]*1

        self.results.schnittList[lineName] = Schnitt(bendingMoments, shearForces, arrayEvaluationPoints)

    return bendingMoments, shearForces, arrayEvaluationPoints

class Schnitt:
    def __init__(self, bendingMoments, shearForces, arrayEvaluationPoints):
        self.bendingMoments = bendingMoments
        self.shearForces = shearForces
        self.arrayEvaluationPoints = arrayEvaluationPoints
        
