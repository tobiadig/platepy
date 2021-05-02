import numpy as np
from getGaussQuadrature import *
def getShapeFunctionForElementType(elementType,ri, si, xi, yi):
    if elementType=='L':
        N, Bb,Bs, detJ = getLinearVectorizedShapeFunctions(ri, si, xi, yi)
    elif elementType == 'Q':
        N, Bb,Bs, detJ = getQuadraticShapeFunctions(ri, si, xi, yi)
    elif elementType == 'MITC4':
        N, Bb,Bs, detJ = getMITCShapefunctions(ri, si, xi, yi)
    elif elementType == 'MITC9':
        N, Bb,Bs, detJ = getMITC9Shapefunctions(ri, si, xi, yi)

    return N, Bb,Bs, detJ

def getLinearVectorizedShapeFunctions(ri, si, xi, yi):
    '''
    INPUT-->    ri: calculation point along the r-axis [float]
                si: calculation point along the s-axis [float]
                xi: coordinates of the element's nodes on the x-axis [1x3 array]
                yi: coordinates of the element's nodes on the y-axis [1x3 array]
    OUTPUT-->   N: values of the shape function matrix at ri, si [1x3 array]
                B: values of the deformation matrix at ri, si [1x3 array]
                detJ: Determinant of the Jacobian matrix [float]
    '''
    nPoints = ri.size
    elemType = len(xi)
    nodeCoordinates = np.zeros((elemType, 2))
    nodeCoordinates[:,0]=xi
    nodeCoordinates[:,1]=yi

    if elemType==2:
        # Define shape functions
        N1 = lambda r, s: 0.5*(1-r)
        N2 = lambda r, s: 0.5*(1+r)

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s)]
        Nval=np.array(Nfun(ri, si))
        Nval=np.moveaxis(Nval, -1, 0)

        N=np.zeros((nPoints,3, 3*elemType))
        N[:,0, 0::3]=Nval
        N[:,1, 1::3]=Nval
        N[:,2, 2::3]=Nval

        # Define shape function derivatives, derive deformation matrix
        N1r = lambda r, s: -0.25*(1-s)
        N2r = lambda r, s: 0.25*(1-s)
        N3r = lambda r, s: 0.25*(1+s)
        N4r = lambda r, s: -0.25*(1+s)

        N1s = lambda r, s: -0.25*(1-r)
        N2s = lambda r, s: -0.25*(1+r)
        N3s = lambda r, s: 0.25*(1+r)
        N4s = lambda r, s: 0.25*(1-r)

        NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)],[N4r(r, s), N4s(r, s)]])
        NrsVal=np.array(NrsFun(ri,si))
        NrsVal = np.moveaxis(NrsVal,-1,0)
        # matmul treat NrsVal as stack of matrixes residing in the LAST 2 indexes

        J=np.matmul(nodeCoordinates.transpose(), NrsVal)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)

        NrsVal = np.matmul(NrsVal,invJ)
        Bf=np.zeros((nPoints,3,3*elemType))
        Bf[:,0,1::3]=NrsVal[:,:,0]
        Bf[:,1,2::3]=NrsVal[:,:,1]
        Bf[:,2,1::3]=NrsVal[:,:,1]
        Bf[:,2,2::3]=NrsVal[:,:,0]

        Bc=np.zeros((nPoints,2,3*elemType))
        Bc[:,0,0::3]=NrsVal[:,:,0]
        Bc[:,0,1::3]=Nval
        Bc[:,1,0::3]=NrsVal[:,:,1]
        Bc[:,1,2::3]=Nval

    elif elemType==3:
        # Define shape functions
        N1 = lambda r, s: r
        N2 = lambda r, s: s
        N3 = lambda r, s: 1-r-s

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s)]
        Nval=np.array(Nfun(ri, si))
        Nval=np.moveaxis(Nval, -1, 0)
        N=np.zeros((nPoints,3, 3*elemType))
        N[:,0, 0::3]=Nval
        N[:,1, 1::3]=Nval
        N[:,2, 2::3]=Nval

        # Define shape function derivatives, derive deformation matrix
        N1r = lambda r, s: 1*np.ones(len(r))
        N2r = lambda r, s: 0*np.ones(len(r))
        N3r = lambda r, s: -1*np.ones(len(r))

        N1s = lambda r, s: 0*np.ones(len(r))
        N2s = lambda r, s: 1*np.ones(len(r))
        N3s = lambda r, s: -1*np.ones(len(r))

        NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)]])
        NrsVal=np.array(NrsFun(ri,si))
        NrsVal=np.moveaxis(NrsVal, -1, 0)

        # Jacobian matrix
        J=np.matmul(nodeCoordinates.transpose(), NrsVal)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)

        NrsVal = np.matmul(NrsVal,invJ)
        Bf=np.zeros((nPoints,3,3*elemType))
        
        Bf[:,0,1::3]=NrsVal[:,:,0]
        Bf[:,1,2::3]=NrsVal[:,:,1]
        Bf[:,2,1::3]=NrsVal[:,:,1]
        Bf[:,2,2::3]=NrsVal[:,:,0]

        Bc=np.zeros((nPoints,2,3*elemType))
        Bc[:,0,0::3]=NrsVal[:,:,0]
        Bc[:,0,1::3]=Nval[:]
        Bc[:,1,0::3]=NrsVal[:,:,1]
        Bc[:,1,2::3]=Nval[:]

    elif elemType==4:
        # Define shape functions
        N1 = lambda r, s: 0.25*(1-r)*(1-s)
        N2 = lambda r, s: 0.25*(1+r)*(1-s)
        N3 = lambda r, s: 0.25*(1+r)*(1+s)
        N4 = lambda r, s: 0.25*(1-r)*(1+s)

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s), N4(r,s)]
        Nval=np.array(Nfun(ri, si))
        # print('N: ', pd.DataFrame(Nval))
        Nval=np.moveaxis(Nval, -1, 0)

        N=np.zeros((nPoints,3, 3*elemType))
        N[:,0, 0::3]=Nval
        N[:,1, 1::3]=Nval
        N[:,2, 2::3]=Nval
 
        # Define shape function derivatives, derive deformation matrix
        N1r = lambda r, s: -0.25*(1-s)
        N2r = lambda r, s: 0.25*(1-s)
        N3r = lambda r, s: 0.25*(1+s)
        N4r = lambda r, s: -0.25*(1+s)

        N1s = lambda r, s: -0.25*(1-r)
        N2s = lambda r, s: -0.25*(1+r)
        N3s = lambda r, s: 0.25*(1+r)
        N4s = lambda r, s: 0.25*(1-r)

        NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)],[N4r(r, s), N4s(r, s)]])
        NrsVal=np.array(NrsFun(ri,si))
        NrsVal = np.moveaxis(NrsVal,-1,0)
        # matmul treat NrsVal as stack of matrixes residing in the LAST 2 indexes

        J=np.matmul(nodeCoordinates.transpose(), NrsVal)
        
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)
        # print('Jac: ', invJ)
        NrsVal = np.matmul(NrsVal,invJ)
        # print('XY der: ', NrsVal)
        Bf=np.zeros((nPoints,3,3*elemType))
        Bf[:,0,1::3]=NrsVal[:,:,0]
        Bf[:,1,2::3]=NrsVal[:,:,1]
        Bf[:,2,1::3]=NrsVal[:,:,1]
        Bf[:,2,2::3]=NrsVal[:,:,0]

        Bc=np.zeros((nPoints,2,3*elemType))
        Bc[:,0,0::3]=NrsVal[:,:,0]
        Bc[:,0,1::3]=-Nval
        Bc[:,1,0::3]=NrsVal[:,:,1]
        Bc[:,1,2::3]=-Nval
        # print('Nval: ', Nval)
        # Bc[:,0,0:elemType]=NrsVal[:,:,0]
        # Bc[:,1,0:elemType]=NrsVal[:,:,1]
        # Bc[:,0,elemType:2*elemType]=Nval
        
        # Bc[:,1,2*elemType:3*elemType]=Nval
    return N, Bf,Bc, detJ

def getQuadraticShapeFunctions(ri, si, xi, yi):
    '''
    INPUT-->    ri: calculation point along the r-axis [float]
                si: calculation point along the s-axis [float]
                xi: coordinates of the element's nodes on the x-axis [1x3 array]
                yi: coordinates of the element's nodes on the y-axis [1x3 array]
    OUTPUT-->   N: values of the shape function matrix at ri, si [1x3 array]
                B: values of the deformation matrix at ri, si [1x3 array]
                detJ: Determinant of the Jacobian matrix [float]
    '''

    nodeCoordinates = np.zeros((9, 2))
    nodeCoordinates[:,0]=xi
    nodeCoordinates[:,1]=yi
    nCoord = len(xi)
    elemType = len(xi)

    if elemType==3:
        # Define shape functions
        N1 = lambda r, s: 0.5*r*(r-1)
        N2 = lambda r, s: 1-r**2
        N3 = lambda r,s: 0.5*r*(1+r)

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s)]
        Nval=np.array(Nfun(ri, si))

        N=np.zeros((3, 3*elemType))
        N[0, 0::3]=Nval
        N[1, 1::3]=Nval
        N[2, 2::3]=Nval

        # Define shape function derivatives, derive deformation matrix
        N1r = lambda r, s: 0.5*(2*r-1)
        N2r = lambda r, s: 1-2*r
        N3r = lambda r, s: 0.5*(2*r+1)

        NrsFun = lambda r,s: np.array([N1r(r, s), N2r(r, s), N3r(r,s)])
        NrsVal=np.array(NrsFun(ri,si))

        # matmul treat NrsVal as stack of matrixes residing in the LAST 2 indexes
        L = np.sqrt((xi[1]-xi[0])**2+(yi[1]-yi[0])**2)
        J= L/2
        detJ = L/2
        invJ = 2/L


        NrsVal = np.matmul(NrsVal,invJ)
        Bf=np.zeros((3,3*elemType))
        Bf[0,1::3]=NrsVal[:,0]
        Bf[1,2::3]=NrsVal[:,1]
        Bf[2,1::3]=NrsVal[:,1]
        Bf[2,2::3]=NrsVal[:,0]

        Bc=np.zeros((nPoints,2,3*elemType))
        Bc[0,0::3]=NrsVal[:,0]
        Bc[0,1::3]=Nval
        Bc[1,0::3]=NrsVal[:,1]
        Bc[1,2::3]=Nval

    elif elemType==7:
        # Define shape functions
        N1 = lambda r, s: r
        N2 = lambda r, s: s
        N3 = lambda r, s: 1-r-s

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s)]
        Nval=np.array(Nfun(ri, si))

        N=np.zeros((3, 3*elemType))
        N[0, 0::3]=Nval
        N[1, 1::3]=Nval
        N[2, 2::3]=Nval

        # Define shape function derivatives, derive deformation matrix
        N1r = lambda r, s: 1*np.ones(len(r))
        N2r = lambda r, s: 0*np.ones(len(r))
        N3r = lambda r, s: -1*np.ones(len(r))

        N1s = lambda r, s: 0*np.ones(len(r))
        N2s = lambda r, s: 1*np.ones(len(r))
        N3s = lambda r, s: -1*np.ones(len(r))

        NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)]])
        NrsVal=np.array(NrsFun(ri,si))

        # Jacobian matrix
        J=np.matmul(nodeCoordinates.transpose(), NrsVal)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)

        NrsVal = np.matmul(NrsVal,invJ)
        Bf=np.zeros((3,3*elemType))

        Bf[0,1::3]=NrsVal[:,0]
        Bf[1,2::3]=NrsVal[:,1]
        Bf[2,1::3]=NrsVal[:,1]
        Bf[2,2::3]=NrsVal[:,0]

        Bc=np.zeros((2,3*elemType))
        Bc[0,0::3]=NrsVal[:,0]
        Bc[0,1::3]=Nval[:]
        Bc[1,0::3]=NrsVal[:,1]
        Bc[1,2::3]=Nval[:]

    elif nCoord==9:
        # Define shape functions
        N1 = lambda r, s: 0.25*(r**2-r)*(s**2-s)
        N2 = lambda r, s: 0.25*(r**2+r)*(s**2-s)
        N3 = lambda r, s: 0.25*(r**2+r)*(s**2+s)
        N4 = lambda r, s: 0.25*(r**2-r)*(s**2+s)
        N5 = lambda r, s: 0.5*(s**2-s)*(1-r**2)
        N6 = lambda r, s: 0.5*(r**2+r)*(1-s**2)
        N7 = lambda r, s: 0.5*(s**2+s)*(1-r**2)
        N8 = lambda r, s: 0.5*(r**2-r)*(1-s**2)
        N9 = lambda r, s: (1-r**2)*(1-s**2)

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s), N4(r,s), N5(r,s), N6(r,s), N7(r,s), N8(r,s), N9(r,s)]
        Nval=np.array(Nfun(ri, si))

        N=np.zeros((3, 3*9))
        N[0, 0::3]=Nval
        N[1, 1::3]=Nval
        N[2, 2::3]=Nval

        # Define shape function derivatives, derive deformation matrix
        N1r = lambda r, s: 0.25*(2*r-1)*(s**2-s)
        N2r = lambda r, s: 0.25*(2*r+1)*(s**2-s)
        N3r = lambda r, s: 0.25*(r*2+1)*(s**2+s)
        N4r = lambda r, s: 0.25*(r*2-1)*(s**2+s)
        N5r = lambda r, s: 0.5*(s**2-s)*(-r*2)
        N6r = lambda r, s: 0.5*(r*2+1)*(1-s**2)
        N7r = lambda r, s: 0.5*(s**2+s)*(-r*2)
        N8r = lambda r, s: 0.5*(r*2-1)*(1-s**2)
        N9r = lambda r, s: (-r*2)*(1-s**2)

        N1s = lambda r, s: 0.25*(r**2-r)*(s*2-1)
        N2s = lambda r, s: 0.25*(r**2+r)*(s*2-1)
        N3s = lambda r, s: 0.25*(r**2+r)*(s*2+1)
        N4s = lambda r, s: 0.25*(r**2-r)*(s*2+1)
        N5s = lambda r, s: 0.5*(s*2-1)*(1-r**2)
        N6s = lambda r, s: 0.5*(r**2+r)*(-s*2)
        N7s = lambda r, s: 0.5*(s*2+1)*(1-r**2)
        N8s = lambda r, s: 0.5*(r**2-r)*(-s*2)
        N9s = lambda r, s: (1-r**2)*(-s*2)

        NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], 
                                        [N2r(r, s), N2s(r, s)], 
                                        [N3r(r, s), N3s(r, s)],
                                        [N4r(r, s), N4s(r, s)],
                                        [N5r(r, s), N5s(r, s)],
                                        [N6r(r, s), N6s(r, s)],
                                        [N7r(r, s), N7s(r, s)],
                                        [N8r(r, s), N8s(r, s)],
                                        [N9r(r, s), N9s(r, s)]])
        NrsVal=np.array(NrsFun(ri,si))

        # matmul treat NrsVal as stack of matrixes residing in the LAST 2 indexes

        J=np.matmul(nodeCoordinates.transpose(), NrsVal)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)
        NrsVal = np.matmul(NrsVal,invJ)
        # print('J of Q elem: ', J)
        # print('det J: ', detJ)
        Bf=np.zeros((3,3*9))
        Bf[0,1::3]=NrsVal[:,0]
        Bf[1,2::3]=NrsVal[:,1]
        Bf[2,1::3]=NrsVal[:,1]
        Bf[2,2::3]=NrsVal[:,0]

        Bc=np.zeros((2,3*9))
        Bc[0,0::3]=NrsVal[:,0]
        Bc[0,1::3]=Nval
        Bc[1,0::3]=NrsVal[:,1]
        Bc[1,2::3]=Nval
    return N, Bf,Bc, detJ

def getMITCShapefunctions(ri, si, xi, yi):
    nodeCoordinates = np.zeros((4, 2))
    nodeCoordinates[:,0]=xi
    nodeCoordinates[:,1]=yi
    # ACCORDING TO BATHES
    N1 = lambda r, s: 0.25*(1+r)*(1+s)
    N2 = lambda r, s: 0.25*(1-r)*(1+s)
    N3 = lambda r, s: 0.25*(1-r)*(1-s)
    N4 = lambda r, s: 0.25*(1+r)*(1-s)

    # Define shape function derivatives, derive deformation matrix
    N1r = lambda r, s: 0.25*(1+s)
    N2r = lambda r, s: -0.25*(1+s)
    N3r = lambda r, s: -0.25*(1-s)
    N4r = lambda r, s: 0.25*(1-s)

    N1s = lambda r, s: 0.25*(1+r)
    N2s = lambda r, s: 0.25*(1-r)
    N3s = lambda r, s: -0.25*(1-r)
    N4s = lambda r, s: -0.25*(1+r)

    # Form the shape function matrix
    Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s), N4(r,s)]
    Nval=Nfun(ri, si)
    N=np.zeros((3, 3*4))
    N[0, 0::3]=Nval
    N[1, 1::3]=Nval
    N[2, 2::3]=Nval

    NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)],[N4r(r, s), N4s(r, s)]])
    NrsVal=NrsFun(ri,si)

    J=np.matmul(nodeCoordinates.transpose(), NrsVal)
    detJ = np.linalg.det(J)

    invJ = np.linalg.inv(J)
    NrsVal = np.matmul(NrsVal,invJ)

    Bb=np.zeros((3,3*4))
    #qui cambia un po tutto wtf, prima:    Bb[0,1::3]=  NrsVal[:,0] Bb[1,2::3]=  -NrsVal[:,1] Bb[2,1::3]=  NrsVal[:,1]   Bb[2,2::3]=  -NrsVal[:,0]
    Bb[0,2::3]=  NrsVal[:,0]
    Bb[1,1::3]=  -NrsVal[:,1]
    Bb[2,2::3]=  NrsVal[:,1]  
    Bb[2,1::3]=  -NrsVal[:,0]

    alpha, beta = naturalToCartesian(xi,yi)
    # print('xi: ', xi)
    # print('yi: ', yi)
    ROTab = np.array([[np.sin(beta), -np.sin(alpha)], 
                        [-np.cos(beta), np.cos(alpha)]])


    Ax = xi[0] - xi[1] - xi[2] + xi[3]
    Ay = yi[0] - yi[1] - yi[2] + yi[3]

    Bx = xi[0] - xi[1] + xi[2] - xi[3]
    By = yi[0] - yi[1] + yi[2] - yi[3]

    Cx = xi[0] + xi[1] - xi[2] - xi[3]
    Cy = yi[0] + yi[1] - yi[2] - yi[3]

    Dx = np.sqrt((Cx+ri*Bx)**2 + (Cy + ri*By)**2)/(8*detJ)
    Dy = np.sqrt((Ax+si*Bx)**2 + (Ay + si*By)**2)/(8*detJ)
    # print('Dx: ', Dx)
    # print('Dy: ', Dy)
    a1 = 0.5
    a2 = (-yi[0] + yi[1])*0.25
    a3 = (xi[0]-xi[1])*0.25
    a4 = (xi[3]-xi[2])*0.25
    a5 = -0.25*(yi[3]-yi[2])
    Arz1 = np.array([a1, a2, a3, -a1, a2, a3, 0, 0, 0, 0, 0, 0])
    Arz2 = np.array([0,0,0,0,0,0, -a1, a5, a4, a1, a5, a4])
    grz = Dx*((1+si)*Arz1+(1-si)*Arz2)

    b1 = 0.5
    b2 = -0.25*(yi[0]-yi[3])
    b3 = 0.25*(xi[0]-xi[3])
    b4 = 0.25*(xi[1]-xi[2])
    b5 = -0.25*(yi[1]-yi[2])
    Bsz1 = np.array([b1, b2, b3, 0, 0, 0, 0, 0, 0, -b1, b2, b3])
    Bsz2 = np.array([0,0,0,b1,b5,b4,-b1, b5, b4, 0, 0, 0])
    gsz = Dy*((1+ri)*Bsz1+(1-ri)*Bsz2)

    gz = np.array([grz, gsz])
    Bs = np.matmul(ROTab, gz)

    return N, Bb,Bs, detJ 

def shapeFun8(GPr, GPs):
    N=np.zeros(8)
    Nr=np.zeros(8)
    Ns=np.zeros(8)

    N[0] = 0.25*(1+GPr)*(1+GPs) - 0.25*(1-GPr**2)*(1+GPs) - 0.25*(1-GPs**2)*(1+GPr)
    N[1] = 0.25*(1-GPr)*(1+GPs) - 0.25*(1-GPr**2)*(1+GPs) - 0.25*(1-GPs**2)*(1-GPr)
    N[2] = 0.25*(1-GPr)*(1-GPs) - 0.25*(1-GPs**2)*(1-GPr) - 0.25*(1-GPr**2)*(1-GPs)
    N[3] = 0.25*(1+GPr)*(1-GPs) - 0.25*(1-GPr**2)*(1-GPs) - 0.25*(1-GPs**2)*(1+GPr)
    N[4] = 0.5*(1-GPr**2)*(1+GPs)
    N[5] = 0.5*(1-GPs**2)*(1-GPr)
    N[6] = 0.5*(1-GPr**2)*(1-GPs)
    N[7] = 0.5*(1-GPs**2)*(1+GPr)

    Nr[0] =  0.25*(1+GPs) + 0.5*GPr*(1+GPs) - 0.25*(1-GPs**2)
    Nr[1] = -0.25*(1+GPs) + 0.5*GPr*(1+GPs) + 0.25*(1-GPs**2)
    Nr[2] = -0.25*(1-GPs) + 0.25*(1-GPs**2) + 0.5*GPr*(1-GPs)
    Nr[3] =  0.25*(1-GPs) + 0.5*GPr*(1-GPs) - 0.25*(1-GPs**2)
    Nr[4] = -GPr*(1+GPs)
    Nr[5] = -0.5*(1-GPs**2)
    Nr[6] = -GPr*(1-GPs)
    Nr[7] = 0.5*(1-GPs**2)

    Ns[0] =  0.25*(1+GPr) - 0.25*(1-GPr**2) + 0.5*GPs*(1+GPr)
    Ns[1] =  0.25*(1-GPr) - 0.25*(1-GPr**2) + 0.5*GPs*(1-GPr)
    Ns[2] = -0.25*(1-GPr) + 0.5*GPs*(1-GPr) + 0.25*(1-GPr**2)
    Ns[3] = -0.25*(1+GPr) + 0.25*(1-GPr**2) + 0.5*GPs*(1+GPr)
    Ns[4] = 0.5*(1-GPr**2)
    Ns[5] = -GPs*(1-GPr)
    Ns[6] = -0.5*(1-GPr**2)
    Ns[7] = -GPs*(1+GPr)

    return N, Nr, Ns

def shapeFun9(GPr, GPs):
    N=np.zeros(9)
    Nr=np.zeros(9)
    Ns=np.zeros(9)

    N[0] = 0.25*(1+GPr)*(1+GPs) - 0.25*(1-GPr**2)*(1+GPs) - 0.25*(1-GPs**2)*(1+GPr) + 0.25*(1-GPr**2)*(1-GPs**2)
    N[1] = 0.25*(1-GPr)*(1+GPs) - 0.25*(1-GPr**2)*(1+GPs) - 0.25*(1-GPs**2)*(1-GPr) + 0.25*(1-GPr**2)*(1-GPs**2)
    N[2] = 0.25*(1-GPr)*(1-GPs) - 0.25*(1-GPs**2)*(1-GPr) - 0.25*(1-GPr**2)*(1-GPs) + 0.25*(1-GPr**2)*(1-GPs**2)
    N[3] = 0.25*(1+GPr)*(1-GPs) - 0.25*(1-GPr**2)*(1-GPs) - 0.25*(1-GPs**2)*(1+GPr) + 0.25*(1-GPr**2)*(1-GPs**2)
    N[4] = 0.5*(1-GPr**2)*(1+GPs) - 0.5*(1-GPr**2)*(1-GPs**2)
    N[5] = 0.5*(1-GPs**2)*(1-GPr) - 0.5*(1-GPr**2)*(1-GPs**2)
    N[6] = 0.5*(1-GPr**2)*(1-GPs) - 0.5*(1-GPr**2)*(1-GPs**2)
    N[7] = 0.5*(1-GPs**2)*(1+GPr) - 0.5*(1-GPr**2)*(1-GPs**2)
    N[8] = (1-GPr**2)*(1-GPs**2)

    Nr[0] =  0.25*(1+GPs) + 0.5*GPr*(1+GPs) - 0.25*(1-GPs**2) - 0.5*GPr*(1-GPs**2)
    Nr[1] = -0.25*(1+GPs) + 0.5*GPr*(1+GPs) + 0.25*(1-GPs**2) - 0.5*GPr*(1-GPs**2)
    Nr[2] = -0.25*(1-GPs) + 0.25*(1-GPs**2) + 0.5*GPr*(1-GPs) - 0.5*GPr*(1-GPs**2)
    Nr[3] =  0.25*(1-GPs) + 0.5*GPr*(1-GPs) -0.25*(1-GPs**2) - 0.5*GPr*(1-GPs**2)
    Nr[4] = -GPr*(1+GPs) + GPr*(1-GPs**2)
    Nr[5] = -0.5*(1-GPs**2) + GPr*(1-GPs**2)
    Nr[6] = -GPr*(1-GPs) + GPr*(1-GPs**2)
    Nr[7] = 0.5*(1-GPs**2) + GPr*(1-GPs**2) 
    Nr[8] = -2*GPr*(1-GPs**2)

    Ns[0] =  0.25*(1+GPr) - 0.25*(1-GPr**2) + 0.5*GPs*(1+GPr) - 0.5*GPs*(1-GPr**2)
    Ns[1] =  0.25*(1-GPr) - 0.25*(1-GPr**2) + 0.5*GPs*(1-GPr) - 0.5*GPs*(1-GPr**2)
    Ns[2] = -0.25*(1-GPr) + 0.5*GPs*(1-GPr) + 0.25*(1-GPr**2) - 0.5*GPs*(1-GPr**2)
    Ns[3] = -0.25*(1+GPr) + 0.25*(1-GPr**2) + 0.5*GPs*(1+GPr) - 0.5*GPs*(1-GPr**2)
    Ns[4] = 0.5*(1-GPr**2) + GPs*(1-GPr**2)
    Ns[5] = -GPs*(1-GPr) + GPs*(1-GPr**2)
    Ns[6] = -0.5*(1-GPr**2) + GPs*(1-GPr**2)
    Ns[7] = -GPs*(1+GPr) + GPs*(1-GPr**2)
    Ns[8] = -2*GPs*(1-GPr**2)
    return N, Nr, Ns 

def getJac(Nr, Ns, xi, yi):
    xr = np.matmul(Nr,xi)
    xs = np.matmul(Ns,xi)
    yr = np.matmul(Nr,yi)
    ys=np.matmul(Ns,yi)
    Jac = np.array([[xr, yr],[xs, ys]])
    JacInv = np.linalg.inv(Jac)
    JacDet = np.linalg.det(Jac)
    return Jac, JacInv, JacDet

def naturalToCartesian(xi,yi):
    #according to bathes book, p. 425
    #alpha: 
    # xb= 0.5*(xi[1]+xi[2])
    # yb = 0.5*(yi[1]+yi[2])

    # xd = 0.5*(xi[0]+xi[3])
    # yd = 0.5*(yi[0]+yi[3])

    #my numeration:
    xb= 0.5*(xi[1]+xi[2])
    yb = 0.5*(yi[1]+yi[2])

    xd = 0.5*(xi[0]+xi[3])
    yd = 0.5*(yi[0]+yi[3])

    u = np.array([xd-xb, yd-yb, 0])
    v = np.array([1,0,0])
    # %using dot product a * b = |a| |b| * cos(angle)
    uLength = np.sqrt(np.dot(u,u))
    vLength = np.sqrt(np.dot(v,v))

    alpha = np.arccos(np.dot(u,v)/(uLength*vLength))

    #beta:  
    # xc= 0.5*(xi[2]+xi[3])
    # yc = 0.5*(yi[2]+yi[3])

    # xa = 0.5*(xi[0]+xi[1])
    # ya = 0.5*(yi[0]+yi[1])

    xc= 0.5*(xi[2]+xi[3])
    yc = 0.5*(yi[2]+yi[3])

    xa = 0.5*(xi[0]+xi[1])
    ya = 0.5*(yi[0]+yi[1])

    u = np.array([xa-xc, ya-yc, 0])
    v = np.array([1,0,0])
    # %using dot product a * b = |a| |b| * cos(angle)
    uLength = np.sqrt(np.dot(u,u))
    vLength = np.sqrt(np.dot(v,v))

    beta = np.arccos(np.dot(u,v)/(uLength*vLength))

    return alpha, beta 

def getMITC9Shapefunctions(ri, si, xi, yi):

    v1 = lambda r, s: np.array([1, r, s, r*s, s**2, 0, 0, 0, 0, 0])
    v2 = lambda r, s: np.array([0, 0, 0, 0, 0, 1, r, s, r*s, r**2])

    M1=np.zeros((2,10))
    M1[0,:] = v1(ri, si)
    M1[1,:] = v2(ri, si)

    a=1/np.sqrt(3)
    pA = (a, 1)  #0
    pB = (-a, 1) #1
    pC = (-1, a) #2
    pD = (-1, -a) #3
    pE = (-a, -1) #4 
    pF = (a, -1) #5
    pG = (1, -a) #6
    pH = (1, a) #7

    points = [pA, pB, pE, pF, pC, pD, pG, pH ]

    M2 = np.zeros((10,10))
    M3 = np.zeros((10,8+9*2))

    for i in range(0,8):
        rP=points[i][0]
        sP = points[i][1]

        N8, Nr8, Ns8 = shapeFun8(rP, sP)


        N9, Nr9, Ns9 = shapeFun9(rP, sP)

        if i < 4:
            M2[i,:] = v1(rP, sP)
            M3[i,0:8] = Nr8
            
            M3[i, 8:17] = -np.dot(Nr9, yi)*N9
            M3[i, 17:] = np.dot(Nr9, xi)*N9
        else:
            M2[i,:] = v2(rP, sP)
            M3[i,0:8] = Ns8
            M3[i, 8:17] = -np.dot(Ns9, yi)*N9
            M3[i, 17:] = np.dot(Ns9, xi)*N9

    nPoints=9
    gaussPoints, gaussWeights = getGaussQuadrature('rectangular', nPoints)
    for i in range(0,nPoints):
        wi = gaussWeights[i]
        N8, Nr8, Ns8 = shapeFun8(gaussPoints[i,0], gaussPoints[i,1])

        N9, Nr9, Ns9 = shapeFun9(gaussPoints[i,0], gaussPoints[i,1])

        M2[8,:] += wi*v1(gaussPoints[i,0], gaussPoints[i,1])
        M2[9, :] += wi*v2(gaussPoints[i,0], gaussPoints[i,1])

        M3[8,0:8] += wi*Nr8
        M3[8, 8:17] += -wi*np.dot(Nr9, yi)*N9
        M3[8, 17:] += wi*np.dot(Nr9, xi)*N9

        M3[9,0:8] += wi*Ns8
        M3[9, 8:17] += -wi*np.dot(Ns9, yi)*N9
        M3[9, 17:] += wi*np.dot(Ns9, xi)*N9

    Bs = M1@np.linalg.inv(M2)@M3

    # Transofrmation in global coordinates
    N8, Nr8, Ns8 = shapeFun8(ri, si)
    N9, Nr, Ns = shapeFun9(ri, si)
    Jac, JacInv, JacDet = getJac(Nr, Ns, xi, yi)
    # print('Jacobian: ', Jac)
    # print('JacDet: ', JacDet)

    endNTemp = np.zeros((3,3*9))
    endNTemp[0,0:-3:3] = N8
    endNTemp[1,1::3] = N9
    endNTemp[2,2::3] = N9

    endN = np.zeros((3, 26))
    endN[:,0:-2] = endNTemp[:,0:-3]
    endN[:,-2:] = endNTemp[:,-2:]


    # alpha, beta = naturalToCartesian(xi,yi)

    # Dx=np.sqrt((np.dot(Ns,xi))**2+(np.dot(Ns,yi))**2)/(JacDet)
    # Dy=np.sqrt((np.dot(Nr,xi))**2+(np.dot(Nr,yi))**2)/(JacDet)

    # transformMat = np.array([[np.sin(beta)*Dx, -np.sin(alpha)*Dy], [-np.cos(beta)*Dx, np.cos(alpha)*Dy]])
    transformMat = JacInv

    Bs = transformMat@Bs

    #Bending component!  #TODO: I dont know what im doing
    NrsVal = np.array([Nr, Ns])
    NrsVal = np.matmul(JacInv,NrsVal)
    Nr = NrsVal[0,:]
    Ns = NrsVal[1,:]

    Bb = np.zeros((3,26))
    Bb[0,17:] = Nr
    Bb[1, 8:17] = -Ns

    Bb[2,8:17] = -Nr
    Bb[2,17:] = Ns


    #reorder matrixes:
    BsTemp = np.zeros((2,27))
    BbTemp = np.zeros((3,27))
    BsToReturn = np.zeros((2,26))
    BbToReturn = np.zeros((3,26))

    BsTemp[:,0:-3:3] = Bs[:,0:8]
    BsTemp[:,1::3] = Bs[:,8:17]
    BsTemp[:,2::3] = Bs[:,17:]
    BsToReturn[:,0:-2] = BsTemp[:,0:-3]
    BsToReturn[:,-2:] = BsTemp[:,-2:]

    BbTemp[:,0:-3:3] = Bb[:,0:8]
    BbTemp[:,1::3] = Bb[:,8:17]
    BbTemp[:,2::3] = Bb[:,17:]
    BbToReturn[:,0:-2] = BbTemp[:,0:-3]
    BbToReturn[:,-2:] = BbTemp[:,-2:]

    return endN, BbToReturn,BsToReturn, JacDet