import numpy as np 
import pandas as pd
from ._shapeFunctions import *
from ._slicingFunctions import *

def getInternalForces(elementsList,uGlob,internalForcePosition, nodesArray, smoothedValues):
    ''' Redirects to the right subroutine, depending on the position where the internal forces should be computed. \n
        Input: \n
            * elementsList: List of element objects. \n
            * uGlob: numpy array with all displacements and rotations computed by solveModel(). \n
            * internalForcePosition = 'nodes': String defining the desired position where the internal forces should be computed.
            Possible values are "nodes" (default), "center" and "intPoints" for the positions used in the Gauss quadratur.\n
            * nodesArray: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates.\n
            * smoothedValues = False: Experimental. If True, the values of the shear forces by displacement based elements are smoothed according to
            the values at the Gauss points. \n
        Return: \n
            * bendingMoments: Numpy array of shape (nInternalForceValues x 3), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively mx, my and mxy.\n
            * shearForces: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively vx and vy.\n
            * internalForcesPositions: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point, the column are respectively the x and the y coordinates.\n
    '''
    if internalForcePosition == 'center':
        bendingMoments, shearForces, internalForcesPositions = getInternalForcesCenter(elementsList,uGlob)
    elif internalForcePosition == 'nodes':
        bendingMoments, shearForces, internalForcesPositions = getInternalForcesNodes(elementsList,uGlob, nodesArray, smoothedValues)
    elif internalForcePosition == 'intPoints':
        bendingMoments, shearForces, internalForcesPositions = getInternalForcesIntPoints(elementsList,uGlob)
    return bendingMoments, shearForces, internalForcesPositions

def getInternalForcesCenter(elementsList,uGlob):
    ''' Computes internal forces at the geometrical center of each element. \n
        Input: \n
            * elementsList: List of element objects. \n
            * uGlob: numpy array with all displacements and rotations computed by solveModel(). \n
        Return: \n
            * bendingMoments: Numpy array of shape (nInternalForceValues x 3), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively mx, my and mxy.\n
            * shearForces: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively vx and vy.\n
            * internalForcesPositions: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point, the column are respectively the x and the y coordinates.\n
    '''
    bendingMoments = np.zeros((len(elementsList),3))
    internalForcesPositions = np.zeros((len(elementsList),2))
    shearForces = np.zeros((len(elementsList),2))
    for k,element in enumerate(elementsList):
        Df = element.Db
        Dc = element.Ds
        elemNodes = element.connectivity
        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        elementType = element.type
        nNodes=element.nNodes
        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        xm = np.average(xi)
        ym = np.average(yi)
        internalForcesPositions[k,0]=xm
        internalForcesPositions[k,1]=ym
        if elementType == 'DB' and nNodes ==4:
            ri = np.array([0])
            si = np.array([0])
        else:
            ri = 0
            si = 0
        _, Bf,Bc, _ = getShapeFunctionForElementType(elementType,ri, si, xi, yi)
        if elementType == 'DB' and nNodes ==4:
            Bf = Bf[0,:,:]
            Bc = Bc[0,:,:]
        kCoeff, _ = getKCoeff(elementType, coherentElemNodes)
        vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
        if elementType == 'DB' and nNodes ==9 or elementType=='MITC':
            changeSign = -1
        else:
            changeSign = 1
        bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*changeSign
        shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]
    return bendingMoments, shearForces, internalForcesPositions

def getInternalForcesNodes(elementsList,uGlob, nodesArray, smoothedValues):
    ''' Computes internal forces at the nodes of each element. The values coming from multiple elements are averaged. \n
        Input: \n
            * elementsList: List of element objects. \n
            * uGlob: numpy array with all displacements and rotations computed by solveModel(). \n
            * internalForcePosition = 'nodes': String defining the desired position where the internal forces should be computed.
            Possible values are "nodes" (default), "center" and "intPoints" for the positions used in the Gauss quadratur.\n
            * nodesArray: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates.\n
            * smoothedValues = False: Experimental. If True, the values of the shear forces by displacement based elements are smoothed according to
            the values at the Gauss points. \n
        Return: \n
            * bendingMoments: Numpy array of shape (nInternalForceValues x 3), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively mx, my and mxy.\n
            * shearForces: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively vx and vy.\n
            * internalForcesPositions: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point, the column are respectively the x and the y coordinates.\n
    '''
    bendingMomentsSum = np.zeros((nodesArray.shape[0],3+1))
    shearForcesSum = np.zeros((nodesArray.shape[0],2+1))
    internalForcesPositions = nodesArray.to_numpy()[:,0:2]
    for k,element in enumerate(elementsList):
        Df = element.Db
        Dc = element.Ds
        elementType = element.type
        elemNodes = element.connectivity
        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        nNodes=element.nNodes
        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        kCoeff, discartedDOF = getKCoeff(elementType, coherentElemNodes)
        vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
        if elementType =='DB' and nNodes ==4 : 
            ri = np.array([-1, 1, 1, -1])
            si = np.array([-1,-1,1,1])
            sigma = np.zeros((2,4))
            for i in range(0, len(ri)):
                N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,np.array([ri[i]]), np.array([si[i]]), xi, yi)
                Bf = Bf[0,:,:]
                Bc = Bc[0,:,:]
                bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                bendingMomentsSum[coherentElemNodes[i],3] += 1
                sigma[:,i] = np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]
            if smoothedValues == False:
                shearForcesSum[coherentElemNodes,2] += 1
                shearForcesSum[coherentElemNodes,0:2] += sigma.transpose()
            elif smoothedValues == True:
                gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',9)
                Ksmooth = np.zeros((4,4))
                fSmooth = np.zeros(4)
                for i in range(0, gaussPoints.shape[0]):
                    ri = gaussPoints[i,0]
                    si = gaussPoints[i,1]
                    wi = gaussWeights[i]
                    N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,np.array([ri]), np.array([si]), xi, yi)
                    N=N[0,0, 0::3]
                    Bf = Bf[0,:,:]
                    Bc = Bc[0,:,:]
                    m1 = np.expand_dims(N,axis=0)
                    m2 = np.expand_dims(N,axis=1)
                    sigmax = np.expand_dims(sigma[0,:], axis = 1)
                    Ksmooth += wi*m2@m1*detJ
                    deb = wi*m2@m1@sigmax*detJ
                    fSmooth += deb[:,0]
                fSmooth = np.expand_dims(fSmooth, 1)
                sigmaSmooth = np.linalg.inv(Ksmooth)@fSmooth
                sigma[0,:] = sigmaSmooth[:,0]
                shearForcesSum[coherentElemNodes,2] += 1
                shearForcesSum[coherentElemNodes,0:2] += sigma.transpose()
        elif elementType =='MITC' and nNodes ==4:
            ri = np.array([1, -1, -1, 1])
            si = np.array([1,1,-1,-1])
            for i in range(0, len(ri)):
                N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,ri[i], si[i], xi, yi)
                bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
                shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 
                bendingMomentsSum[coherentElemNodes[i],3] += 1
                shearForcesSum[coherentElemNodes[i],2] += 1 
        elif elementType=='DB' and nNodes ==9:
            ri = np.array([-1, 1, 1, -1, 0, 1, 0, -1, 0])
            si = np.array([-1, -1, 1, 1, -1, 0, 1, 0, 0])
            for i in range(0, len(ri)):
                N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,ri[i], si[i], xi, yi)
                bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
                shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1
                bendingMomentsSum[coherentElemNodes[i],3] += 1
                shearForcesSum[coherentElemNodes[i],2] += 1 
        elif elementType == 'MITC' and nNodes == 9:
            ri = np.array([1, -1, -1, 1, 0, -1, 0, 1, 0])
            si = np.array([1, 1, -1, -1, 1, 0, -1, 0, 0])
            for i in range(0, len(ri)):
                N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,ri[i], si[i], xi, yi)
                bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
                shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1
                bendingMomentsSum[coherentElemNodes[i],3] += 1
                shearForcesSum[coherentElemNodes[i],2] += 1 
        else:
            raise TypeError('Node calculation of internal forces not implemented yet for this element type')
    bendingMoments = bendingMomentsSum[:,0:3]/bendingMomentsSum[:,3, None]
    shearForces = shearForcesSum[:,0:2]/shearForcesSum[:,2, None]
    return bendingMoments, shearForces, internalForcesPositions

def getInternalForcesIntPoints(elementsList,uGlob):
    ''' Computes internal forces at the integration points of each element. \n
        Input: \n
            * elementsList: List of element objects. \n
            * uGlob: numpy array with all displacements and rotations computed by solveModel(). \n
        Return: \n
            * bendingMoments: Numpy array of shape (nInternalForceValues x 3), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively mx, my and mxy.\n
            * shearForces: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point which coordinates are in internalForcesPositions, the column are respectively vx and vy.\n
            * internalForcesPositions: Numpy array of shape (nInternalForceValues x 2), each row is an evaluation point, the column are respectively the x and the y coordinates.\n
    '''
    nGaussPoints = 4
    bendingMoments = np.zeros((len(elementsList)*nGaussPoints,3))
    internalForcesPositions = np.zeros((len(elementsList)*nGaussPoints,2))
    shearForces = np.zeros((len(elementsList)*nGaussPoints,2))
    for k,element in enumerate(elementsList):
        Df = element.Db
        Dc = element.Ds
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
        if elementType == 'DB' and nNodes == 4:
            ri =np.linspace(-1, 1, num=nGaussPoints)
            si = np.zeros(nGaussPoints)
            allN, allBf,allBc, alldetJ =getLinearVectorizedShapeFunctions(ri, si, xi, yi)
            for i in range(0, len(ri)):
                N=allN[i,:,:]
                Bf = allBf[i,:,:]
                Bc = allBc[i,:,:]
                Nval = N[0, 0::3]
                bendingMoments[k*nGaussPoints+i,0:3] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                shearForces[k*nGaussPoints+i,0:2] = np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*1 
                xr = Nval@xi
                yr = Nval@yi
                internalForcesPositions[k*nGaussPoints+i,0]=xr
                internalForcesPositions[k*nGaussPoints+i,1]=yr
        elif elementType == 'DB' and nNodes == 9:
            gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',4)
            for i in range(0,nGaussPoints):
                ri = gaussPoints[i,0]
                si = gaussPoints[i,1]
                wi = gaussWeights[i]
                N, Bf,Bc, detJ = getQuadraticShapeFunctions(ri, si, xi, yi)
                Nval=N[0, 0::3]
                bendingMoments[k*nGaussPoints+i,0:3] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                shearForces[k*nGaussPoints+i,0:2] = np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 
                xr = Nval@xi
                yr = Nval@yi
                internalForcesPositions[k*nGaussPoints+i,0]=xr
                internalForcesPositions[k*nGaussPoints+i,1]=yr
    return bendingMoments, shearForces, internalForcesPositions

def  getInternalForcesDSB(elementsList,uDownStandBeam,internalForcePosition, downStandBeam):
    ''' Computes internal forces at the center of each downstand beam element.\n
        Input: \n
            * elementsList: List of element objects from the downstandbeam. \n
            * uDownStandBeam: numpy array with all displacements and rotations of the downstand beam computed by solveModel(). \n
            * internalForcePosition = 'nodes': String defining the desired position where the internal forces should be computed.
            Possible values are "nodes" (default), "center" and "intPoints" for the positions used in the Gauss quadratur.\n
            * downStandBeam: DownStandBeam object.\n
        Return: \n
            * Nforces: Numpy vector of shape (number downstand beam elements,) containing normal forces. \n
            * Vforces: Numpy vector of shape (number downstand beam elements,) containing shear forces.\n
            * Mforces: Numpy vector of shape (number downstand beam elements,) containing bending moments.\n
            * internalForcesPositions: Numpy array of shape (number downstand beam elements x 2), each row is an evaluation point, the column are respectively the x and the y coordinates.\n
    '''
    if internalForcePosition == 'center':
        Nforces = np.zeros((len(elementsList)))
        Vforces =  np.zeros((len(elementsList)))
        Mforces = np.zeros((len(elementsList)))
        internalForcesPositions = np.zeros((len(elementsList),2))
        for k,element in enumerate(elementsList):
            elemNodes = element.connectivity
            coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
            nNodes=element.nNodes
            xi=element.coordinates[:,0]
            yi=element.coordinates[:,1]
            xm = np.average(xi)
            ym = np.average(yi)
            internalForcesPositions[k,0]=xm
            internalForcesPositions[k,1]=ym
            ri = np.array([0])
            si = np.array([0])
            N, Bc,Bb,Bs, detJ = getLinearVectorizedShapeFunctions(ri, si, xi, yi)

            kCoeff, discartedDOF = getKCoeff('timo', coherentElemNodes)

            Emod = downStandBeam.body.eModule
            Gmod = downStandBeam.body.gModule
            crossA = downStandBeam.crossSection.A
            crossI = downStandBeam.crossSection.Iy
            Dc =Emod*crossA
            Db = Emod*crossI
            Ds = 5/6*Gmod*crossA

            vLoc = uDownStandBeam[kCoeff]
            Nforces[k] = Dc*Bc@vLoc
            Mforces[k] = Db*Bb@vLoc
            Vforces[k]= Ds*Bs@vLoc
    return Nforces, Vforces, Mforces, internalForcesPositions