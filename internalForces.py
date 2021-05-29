import numpy as np 
import pandas as pd
pd.set_option('display.max_rows', None)
from shapeFunctions import *
from slicingFunctions import *

def  getInternalForcesDSB(elementsList,uDownStandBeam,internalForcePosition, downStandBeam):
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

            vLoc = np.matmul(element.rotationMatrix, uDownStandBeam[kCoeff])
            Nforces[k] = Dc*Bc@vLoc
            Mforces[k] = Db*Bb@vLoc
            Vforces[k]= Ds*Bs@vLoc

    return Nforces, Vforces, Mforces, internalForcesPositions

def getInternalForces(elementsList,uGlob,internalForcePosition, nodesArray, smoothedValues):

    if internalForcePosition == 'center':
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
            N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,ri, si, xi, yi)
            if elementType =='L':
                Bf = Bf[0,:,:]
                Bc = Bc[0,:,:]

            kCoeff, discartedDOF = getKCoeff(elementType, coherentElemNodes)

            vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
            bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]
            shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1

    elif internalForcePosition == 'nodes':
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
            # kCoeff = np.zeros((3*nNodes),dtype=int)
            # for i in range(0,3):
            #     kCoeff[0+i::3]=coherentElemNodes*3+i
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
        
    elif internalForcePosition == 'intPoints':
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
                # ri =np.linspace(-1, 1, num=nGaussPoints)
                # si = np.ones(nGaussPoints)/np.sqrt(3)
                gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',4)
                for i in range(0,nGaussPoints):
                    ri = gaussPoints[i,0]
                    si = gaussPoints[i,1]
                    wi = gaussWeights[i]
                    N, Bf,Bc, detJ = getQuadraticShapeFunctions(ri, si, xi, yi)
                    # N, Bf,Bc, detJ = getQuadraticShapeFunctions(ri[i], si[i], xi, yi)
                    Nval=N[0, 0::3]

                    bendingMoments[k*nGaussPoints+i,0:3] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                    shearForces[k*nGaussPoints+i,0:2] = np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 
                    xr = Nval@xi
                    yr = Nval@yi
                    internalForcesPositions[k*nGaussPoints+i,0]=xr
                    internalForcesPositions[k*nGaussPoints+i,1]=yr
    
    return bendingMoments, shearForces, internalForcesPositions