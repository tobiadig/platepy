import numpy as np
import gmsh
from .shapeFunctions import *
from scipy.integrate import trapezoid

def beamComponents(self,lineName, startCoord, endCoord, nEvaluationPoints,integrationWidth = 0, nIntegrationPoints =0):

    uGlob = self.results.uGlobPlate
    # print(uGlob)

    elementsList = self.mesh.elementsList
    arrayEvaluationPoints = np.zeros((nEvaluationPoints,2))

    arrayEvaluationPoints[:,0] = np.linspace(startCoord[0], endCoord[0], num=nEvaluationPoints)
    arrayEvaluationPoints[:,1] = np.linspace(startCoord[1], endCoord[1], num=nEvaluationPoints)
    
    bendingMoments = np.zeros((nEvaluationPoints,3))
    shearForces = np.zeros((nEvaluationPoints,2))
    verticalDisplacements = np.zeros(nEvaluationPoints)
    resultsScaleVertDisp = self.results.resultsScaleVertDisp

    if integrationWidth ==0:
        for k,evaluationPoint in enumerate(arrayEvaluationPoints):
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

    elif integrationWidth>0:
        verticalDisplacements = None
        sampleWidth = integrationWidth/nIntegrationPoints
        lineDir = np.array([endCoord[0]-startCoord[0],endCoord[1]-startCoord[1]])
        lineDir = lineDir/np.sqrt(lineDir[0]**2+lineDir[1]**2)
        normLineDir = np.array([lineDir[1], -lineDir[0]])

        for k,evaluationPoint in enumerate(arrayEvaluationPoints):
            pCord = np.array([evaluationPoint[0],evaluationPoint[1]])
            startIntegrationLine = pCord-normLineDir*integrationWidth/2
            endIntegrationLine = pCord + normLineDir*integrationWidth/2

            arrayIntegrationPoints = np.zeros((nIntegrationPoints,2))
            arrayIntegrationPoints[:,0] = np.linspace(startIntegrationLine[0], endIntegrationLine[0], num=nIntegrationPoints)

            arrayIntegrationPoints[:,1] = np.linspace(startIntegrationLine[1], endIntegrationLine[1], num=nIntegrationPoints)
            integrationBendingMoments = np.zeros((nIntegrationPoints,3))
            integrationShearForces = np.zeros((nIntegrationPoints,2))

            for j, integrationPoint in enumerate(arrayIntegrationPoints):
                elementTaggetEl, elementTypegetEl, nodeTagsgetEl, ugetEl, vgetEl, wgetEl = gmsh.model.mesh.getElementByCoordinates(integrationPoint[0],integrationPoint[1],0, dim=2)
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

                integrationBendingMoments[j,0:3] = np.matmul(Df,np.matmul(Bb, vLoc))[:,0]*1
                integrationShearForces[j,0:2] = np.matmul(Dc, np.matmul(Bs, vLoc))[:,0]*1
            
            for j in range(0,3):
                arrayToIntegrate = integrationBendingMoments[:,j]
                bendingMoments[k,j]=trapezoid(arrayToIntegrate,dx=sampleWidth)

            for j in range(0,2):
                arrayToIntegrate = integrationShearForces[:,j]
                shearForces[k,j]=trapezoid(arrayToIntegrate,dx=sampleWidth)

    self.results.schnittList[lineName] = Schnitt(bendingMoments, shearForces,verticalDisplacements, arrayEvaluationPoints)

    return bendingMoments, shearForces, arrayEvaluationPoints

class Schnitt:
    def __init__(self, bendingMoments, shearForces,verticalDisplacements, arrayEvaluationPoints):
        self.bendingMoments = bendingMoments
        self.shearForces = shearForces
        self.verticalDisplacements = verticalDisplacements
        self.arrayEvaluationPoints = arrayEvaluationPoints
        
