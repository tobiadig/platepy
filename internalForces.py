import numpy as np 
from shapeFunctions import *
from slicingFunctions import *

def getInternalForces(elementType,elementsList,uGlob,internalForcePosition, Df, Dc, nodesArray):
    if internalForcePosition == 'center':
        bendingMoments = np.zeros((len(elementsList),3))
        internalForcesPositions = np.zeros((len(elementsList),2))
        shearForces = np.zeros((len(elementsList),2))
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
            if elementType == 'L':
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
            bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
            shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1

    elif internalForcePosition == 'nodes':
        bendingMomentsSum = np.zeros((nodesArray.shape[0],3+1))
        shearForcesSum = np.zeros((nodesArray.shape[0],2+1))
        internalForcesPositions = nodesArray.to_numpy()[:,0:2]
        for k,element in enumerate(elementsList):
            elemNodes = element.connectivity
            coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
            nNodes=element.nNodes
            xi=element.coordinates[:,0]
            yi=element.coordinates[:,1]
            kCoeff = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeff[0+i::3]=coherentElemNodes*3+i
            vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
            if elementType =='L' :
                ri = np.array([-1, 1, 1, -1])
                si = np.array([-1,-1,1,1])
                for i in range(0, len(ri)):
                    N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,np.array([ri[i]]), np.array([si[i]]), xi, yi)
                    Bf = Bf[0,:,:]
                    Bc = Bc[0,:,:]
                    bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
                    shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 
                    bendingMomentsSum[coherentElemNodes[i],3] += 1
                    shearForcesSum[coherentElemNodes[i],2] += 1
            elif elementType =='MITC4':
                ri = np.array([1, -1, -1, 1])
                si = np.array([1,1,-1,-1])
                for i in range(0, len(ri)):
                    N, Bf,Bc, detJ = getShapeFunctionForElementType(elementType,ri[i], si[i], xi, yi)
                    bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                    shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*1 
                    bendingMomentsSum[coherentElemNodes[i],3] += 1
                    shearForcesSum[coherentElemNodes[i],2] += 1 
            else:
                raise TypeError('Node calculation of internal forces not implemented yet for this element type')
        bendingMoments = bendingMomentsSum[:,0:3]/bendingMomentsSum[:,3, None]
        shearForces = shearForcesSum[:,0:2]/shearForcesSum[:,2, None]
        
    elif internalForcePosition == 'intPoints':
        nGaussPoints = 7
        bendingMoments = np.zeros((len(elementsList)*nGaussPoints,3))
        internalForcesPositions = np.zeros((len(elementsList)*nGaussPoints,2))
        shearForces = np.zeros((len(elementsList)*nGaussPoints,2))
        for k,element in enumerate(elementsList):
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
            if elementType == 'L':
                ri =np.linspace(-1, 1, num=nGaussPoints)
                si = np.zeros(nGaussPoints)
                allN, allBf,allBc, alldetJ, allNval =getLinearVectorizedShapeFunctions(ri, si, xi, yi)
                for i in range(0, len(ri)):
                    N=allN[i,:,:]
                    Bf = allBf[i,:,:]
                    Bc = allBc[i,:,:]
                    Nval = allNval[i,:]
                    bendingMoments[k*nGaussPoints+i,0:3] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                    shearForces[k*nGaussPoints+i,0:2] = np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*1 
                    xr = Nval@xi
                    yr = Nval@yi
                    internalForcesPositions[k*nGaussPoints+i,0]=xr
                    internalForcesPositions[k*nGaussPoints+i,1]=yr
    
    return bendingMoments, shearForces, internalForcesPositions