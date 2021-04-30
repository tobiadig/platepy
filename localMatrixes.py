import numpy as np
from shapeFunctions import *
from getGaussQuadrature import *

def GetLocalMatrix(xi, yi, Df,Dc, p,nNodes, elementDefinition):
    ''' Input/Output descriptions
    ElemType: Quadrangluar or Triangular or Beam + Linear or Quadratic or MITC + Reduced or Normal Integration
        xi, yi: element's nodes coordinates
        Df, Dc: material of the plate
        p: load vector
        reducedIntegration: bool to manage the number of Gauss integration points 

        return:
        kLocal, fLocal
    '''
    temp = elementDefinition.split('-')
    elementType = temp[0]
    elementIntegration = temp[1]

    if elementType == 'L':
        kLocal, fLocal = getLinearMatrix(xi, yi, Df, Dc, p, nNodes, elementIntegration)
    elif elementType == 'MITC4':
        kLocal, fLocal = getMITC4Matrix(xi, yi, Df, Dc, p, nNodes)
    elif elementType == 'Q':
        kLocal, fLocal = getQuadraticMatrix(xi, yi, Df, Dc, p, nNodes, elementIntegration)
    elif elementType == 'MITC9':
        kLocal, fLocal = getMITC9Matrix(xi, yi, Df, Dc, p, nNodes)

    return kLocal, fLocal

def getLinearMatrix(xi, yi, Df, Dc, p, nNodes, elementIntegration):
    if nNodes==2:
        gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('linear',1)      
        gaussPoints, gaussWeights =  getGaussQuadrature('linear',2)
    if nNodes == 3:
        gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('triangular',1)          
        gaussPoints, gaussWeights =  getGaussQuadrature('triangular',3)
    elif nNodes == 4:
        gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('rectangular',1)
        gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',4)

    fLocal = np.zeros(3*nNodes)
    ri = gaussPoints[:,0]
    si = gaussPoints[:,1]
    wi = gaussWeights
    N, Bf,Bc, detJ = getLinearVectorizedShapeFunctions(ri, si, xi, yi)

    m11 = np.matmul(Bf.transpose((0,2,1)),Df)
    m12 = np.matmul(m11,Bf)

    m12 = np.moveaxis(m12,0,-1)
    kLocal = np.dot(m12,wi*detJ)

    kBDEbug = np.dot(m12,wi*detJ)
    if elementIntegration == 'R': # reduces number of integration points, else keep full integration for shear
        ri = gaussPointsRed[:,0]
        si = gaussPointsRed[:,1]
        wi = gaussWeightsRed[:]
        N, Bf,Bc, detJ = getLinearVectorizedShapeFunctions(ri, si, xi, yi)

    m21 = np.matmul(Bc.transpose((0,2,1)),Dc)
    m22 = np.matmul(m21,Bc)
    m22 = np.moveaxis(m22,0,-1)
    kLocal += np.dot(m22,wi*detJ)
    ksDEbug = np.dot(m22,wi*detJ)

    ri = gaussPoints[:,0]
    si = gaussPoints[:,1]
    wi = gaussWeights[:]
    N, Bf,Bc, detJ = getLinearVectorizedShapeFunctions(ri, si, xi, yi)
    mTemp = np.matmul(N.transpose((0,2,1)),p.magnitude)
    mTemp=np.expand_dims(mTemp,2)
    mTemp = np.moveaxis(mTemp,0,-1)
    fLocal = np.dot(mTemp,wi*detJ)

    return kLocal, fLocal

def getMITC4Matrix(xi, yi, Df, Dc, p, nNodes):
    if nNodes==2:        
        gaussPoints, gaussWeights =  getGaussQuadrature('linear',2)
    if nNodes == 3:         
        gaussPoints, gaussWeights =  getGaussQuadrature('triangular',3)
    elif nNodes == 4:
        gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',4)
    
    fLocal = np.zeros(3*4)
    kLocal = np.zeros((3*4,3*4))
    kBDEbug = np.zeros((3*4,3*4))
    ksDEbug = np.zeros((3*4,3*4))
    for i in range(0,gaussPoints.shape[0]):
        ri = gaussPoints[i,0]
        si = gaussPoints[i,1]
        wi = gaussWeights[i]
        N, Bf,Bc, detJ = getMITCShapefunctions(ri, si, xi, yi)
        m11 = np.matmul(Bf.transpose(),Df)
        m12 = np.matmul(m11,Bf)
        m21 = np.matmul(Bc.transpose(),Dc) 
        m22 = np.matmul(m21,Bc)
        ksDEbug += wi*m22*detJ
        kBDEbug += wi*m12*detJ

        kLocal += wi*m12*detJ + wi*m22*detJ
        fLocal += wi*np.matmul(N.transpose(),p.magnitude)*detJ

    return kLocal, fLocal

def getQuadraticMatrix(xi, yi, Df, Dc, p, nNodes, elementIntegration):
    elementShape = len(xi)
    if elementShape==2:
        gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('linear',1)
        gaussPoints, gaussWeights =  getGaussQuadrature('linear',2)
    if elementShape == 3:
        gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('triangular',1)          
        gaussPoints, gaussWeights =  getGaussQuadrature('triangular',3)
    elif elementShape == 9:
        gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('rectangular',4)
        gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',9)

    fLocal = np.zeros(3*9)
    kLocal = np.zeros((3*9,3*9))

    kBDEbug = np.zeros((3*9,3*9))
    ksDEbug = np.zeros((3*9,3*9))
    for i in range(0,gaussPoints.shape[0]):
        ri = gaussPoints[i,0]
        si = gaussPoints[i,1]
        wi = gaussWeights[i]
        N, Bf,Bc, detJ = getQuadraticShapeFunctions(ri, si, xi, yi)
        m11 = np.matmul(Bf.transpose(),Df)
        m12 = np.matmul(m11,Bf)
        kBDEbug += wi*m12*detJ
        kLocal += wi*m12*detJ
        fLocal += wi*np.matmul(N.transpose(),p.magnitude)*detJ

    if elementIntegration == 'N':
        gaussPointsRed = gaussPoints
        gaussWeightsRed = gaussWeights
    elif elementIntegration != 'R':
        raise TypeError('Integration not recognised')

    for i in range(0, gaussPointsRed.shape[0]):
        ri = gaussPointsRed[i,0]
        si = gaussPointsRed[i,1]
        wi = gaussWeightsRed[i]
        N, Bf,Bc, detJ = getQuadraticShapeFunctions(ri, si, xi, yi)
        m21 = np.matmul(Bc.transpose(),Dc) 
        m22 = np.matmul(m21,Bc)
        kLocal += wi*m22*detJ
    return kLocal, fLocal

def getMITC9Matrix(xi, yi, Df, Dc, p, nNodes):
    elementShape = len(xi)
    if elementShape==2:     
        gaussPoints, gaussWeights =  getGaussQuadrature('linear',2)
    if elementShape == 3:
        gaussPoints, gaussWeights =  getGaussQuadrature('triangular',3)
    elif elementShape == 9:
        gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',9)

    nDOFs = 8+2*9
    fLocal = np.zeros(nDOFs)
    kLocal = np.zeros((nDOFs,nDOFs))

    kBDEbug = np.zeros((nDOFs,nDOFs))
    ksDEbug = np.zeros((nDOFs,nDOFs))
    for i in range(0,gaussPoints.shape[0]):
        ri = gaussPoints[i,0]
        si = gaussPoints[i,1]
        wi = gaussWeights[i]
        N, Bf,Bc, detJ = getMITC9Shapefunctions(ri, si, xi, yi)

        m21 = np.matmul(Bc.transpose(),Dc) 
        m22 = np.matmul(m21,Bc)
        ksDEbug += wi*Bc.transpose() @ Dc @ Bc *detJ
        kBDEbug += wi*Bf.transpose() @ Df @ Bf * detJ
        temp1 = np.expand_dims(p.magnitude, 1)
        temp2 = wi* N.transpose() @ temp1 *detJ
        kLocal += wi*Bf.transpose() @ Df @ Bf *detJ + wi*Bc.transpose() @ Dc @ Bc*detJ
        fLocal += temp2[:,0]
    return kLocal, fLocal
