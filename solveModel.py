''' Module Information
-----------------------------------------------------------
Purpose of module: computes the equilibrium solution using finite elements, stores the result as attribute
                    of a plateModel object
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''
#%% Basic modules
import numpy as np

import gmsh # To create CAD model and mesh
from scipy import sparse    # To work with sparse matrix
from scipy.sparse import linalg # Linear sparse matrices operations
from scipy.linalg import block_diag # to create the rotation matrix

# for debug purposes
import time
from tqdm import tqdm

def solveModel(self, reducedIntegration = False, resultsScaleIntForces = 1, resultsScaleVertDisp = 1):
    ''' Input/Output descriptions
        self: PlateModel class, where the geometry and the mesh are initialized
        reducedIntegration: to manage the number of Gauss integration points
    '''
    # Loop over elements and assemble stiffness matrices
    nodes=self.mesh.nodesArray
    nodesRotations = self.mesh.nodesRotations
    elementsList = self.mesh.elementsList
    BCs = self.mesh.BCs
    nGDofs = nodes.shape[0]*3

    #initialize arrays for the creation of the sparse matrixes
    nSparseData = len(elementsList)*(4*3)**2
    rowsForStiffnessSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForStiffnessSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForStiffnessSparseMatrix = np.zeros(nSparseData)
    rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForForceSparseMatrix = np.zeros(nSparseData)

    p=self.loads[0].magnitude   #TODO: implement calculation for loads as functions of x and y
                                        
    startIndexStiffness = 0
    startIndexForce = 0
    k=0
    # for element in tqdm(elementsList):
    for element in elementsList:
        # iterate over each element, get stiffness and creates basis for the creation of the global stiffness matrix
        elemNodes = element.connectivity
        
        nNodes=element.nNodes
        elemNodesRotations = nodesRotations.loc[elemNodes].to_numpy()

        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]

        Df = self.plates[0].Df   #TODO: the elements should know to which plate it belongs and take the right material param.
        Dc = self.plates[0].Dc

        kLocalNotRotated,fLocal = GetLocalMatrix(xi, yi, Df,Dc,p, reducedIntegration)

        # create rotation matrix
        Ri = []
        RiInit = False
        for dofRot in elemNodesRotations:
            if not(RiInit):
                R=rotMatrix(dofRot)
                RiInit=True
            else:
                R = block_diag(R, rotMatrix(dofRot))

        element.rotationMatrix = R
        #rotate stiffness matrix
        kTemp = np.matmul(kLocalNotRotated, R)
        kLocal = np.matmul(R.transpose(), kTemp)
        # kLocal = kLocalNotRotated

        # coefficients of the DOFs and assignment of the stiffness matrix / force vector
        kCoeff = np.zeros((3*nNodes),dtype=int)
        for i in range(0,3):
            kCoeff[0+i::3]=elemNodes*3-3+i
        
        rows = np.zeros(kLocal.size,dtype=int)
        columns = np.zeros(kLocal.size,dtype=int)
        c=0
        for j in kCoeff:
            for i in kCoeff:
                rows[c] = i
                columns[c] = j
                c+=1

        # create vectors to assemble sparse matrixes
        rowsForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = rows
        columnsForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = columns
        dataForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = kLocal.flatten(order='F')
        startIndexStiffness += rows.size
        
        rowsForForceSparseMatrix[startIndexForce:startIndexForce+kCoeff.size] = kCoeff
        dataForForceSparseMatrix[startIndexForce:startIndexForce+kCoeff.size] = fLocal
        startIndexForce += kCoeff.size
        k+=1

    # create global matrixes
    sparseGlobalMatrix = sparse.csr_matrix((dataForStiffnessSparseMatrix,(rowsForStiffnessSparseMatrix,columnsForStiffnessSparseMatrix)))
    sparseForceGlobal = sparse.csr_matrix((dataForForceSparseMatrix,(rowsForForceSparseMatrix,columnsForForceSparseMatrix)))

    # apply boundary conditions
    rDofsBool = np.zeros((nGDofs),dtype=bool)
    for constraint in BCs:
        node=int(constraint[0])
        rDofsBool[node*3-3:node*3] = constraint[1:].astype(bool)
    allDofs =np.arange(0,nGDofs)
    fDofsBool = np.invert(rDofsBool)
    fDofsInt = allDofs[fDofsBool]

    kMatFree = sparseGlobalMatrix[fDofsInt,:]
    kMatFree = kMatFree[:,fDofsInt]
    fVecFree = sparseForceGlobal[fDofsInt]

    # SOLVE
    Uf = sparse.linalg.spsolve(kMatFree,fVecFree)
    Uf=Uf.reshape(-1,1)

    # global displacement and force vector
    uGlob=np.zeros((nGDofs,1))
    uGlob[fDofsInt]=Uf

    uGlobSparse = sparse.csr_matrix(uGlob)
    globalForce = np.matmul(sparseGlobalMatrix.toarray(),uGlob)

    # elaborate and store results
    reactionForces = globalForce[rDofsBool]
    nodes = self.mesh.nodesArray
    outPos = np.zeros((nodes.shape[0],2))
    values = np.zeros((nodes.shape[0],3))
    outPos[:,0:2]=nodes[:,0:2]
    values[:,0]=uGlob[0::3,0]
    values[:,1]=uGlob[1::3,0]
    values[:,2]=uGlob[2::3,0]
    self.results = Result(outPos,values[:,0], values[:,1], values[:,2],resultsScale=(resultsScaleVertDisp,resultsScaleIntForces))

    ###################################################################3!!!
    #compute MOMENTSSS
    bendingMoments = np.zeros((len(elementsList),3))
    internalForcesPositions = np.zeros((len(elementsList),2))

    shearForces = np.zeros((len(elementsList),2))

    k=0
    for element in elementsList:
        elemNodes = element.connectivity
        nNodes=element.nNodes
        # elemNodesRotations = nodesRotations[elemNodes-1]

        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        xm = np.average(xi)
        ym = np.average(yi)
        internalForcesPositions[k,0]=xm
        internalForcesPositions[k,1]=ym


        Df = self.plates[0].Df   #TODO: the elements should know to which plate it belongs and take the right material param.
        Dc = self.plates[0].Dc

        #Bb at the center for later moment calculation
        elemType = len(xi)
        N, Bf,Bc, detJ = GetShapeFunction(elemType,0, 0, xi, yi)
        element.BbMat = Bf

        kCoeff = np.zeros((3*nNodes),dtype=int)
        for i in range(0,3):
            kCoeff[0+i::3]=elemNodes*3-3+i
        # uLoc = R*uGlob
        # print('R shape: ', element.rotationMatrix.shape)
        # print('uGlob shape: ',uGlob[kCoeff].shape )

        vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])

        bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]
        shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]
        k+=1
    self.results.bendingMoments=bendingMoments*resultsScaleIntForces
    self.results.internalForcesPositions=internalForcesPositions
    
    self.results.shearForces = shearForces*resultsScaleIntForces

    return outPos, values

def rotMatrix(theta):
    '''
        creates rotation matrix of angle theta for a single node
    '''
    A = np.array([[1., 0, 0],
                  [0, np.cos(theta), -np.sin(theta)],
                  [0, np.sin(theta), np.cos(theta)]])
    return A

def GetLocalMatrix(xi, yi, Df,Dc, p, reducedIntegration):
    ''' Input/Output descriptions
        xi, yi: element's nodes coordinates
        Df, Dc: material of the plate
        p: load vector
        reducedIntegration: bool to manage the number of Gauss integration points 

        return:
        kLocal, fLocal
    '''
    elemType = len(xi)

    if elemType==3:
    # 1 point
        gaussPointsRed =  np.array([[1/3, 1/3]])
        gaussWeightsRed =  np.array([1/2])

    # 3 points            
        gaussPoints =  np.array([[1/6, 1/6],
                                [2/3,   1/6],
                                [1/6, 2/3]])
        gaussWeights =  np.array([1/6, 1/6, 1/6])
    # 4 points
        # gaussPoints =  np.array([[1/3, 1/3],
        #                         [1/5, 1/5],
        #                         [3/5, 1/5],
        #                         [1/5, 3/5]])
        # gaussWeights =  np.array([-27/96, 25/96, 25/96, 25/96])

    elif elemType == 4:
    # 1 point
        gaussPointsRed =  np.array([[0, 0]])
        gaussWeightsRed =  np.array([4])

    # 4 points
        tempVal = 1/np.sqrt(3)
        gaussPoints =  np.array([[-tempVal, -tempVal],
                                [tempVal,   -tempVal],
                                [tempVal,tempVal],
                                [-tempVal, tempVal]])
        gaussWeights =  np.array([1, 1, 1, 1])

    # 9 points
        # tempVal = np.sqrt(3/5)

        # gaussPoints =  np.array([[-tempVal, -tempVal],
        #                         [ 0      , -tempVal],
        #                         [+tempVal, -tempVal],
        #                         [-tempVal,    0    ],
        #                         [   0    ,    0    ] ,
        #                         [+tempVal,    0    ],
        #                         [-tempVal, +tempVal],
        #                         [   0    , +tempVal],
        #                         [+tempVal, +tempVal]])
        # gaussWeights =  np.array([25, 40, 25, 40, 64, 40, 25, 40, 25])/81

    kLocal = np.zeros((3*elemType,3*elemType))
    fLocal = np.zeros(3*elemType)

    #full integration for bending component
    for i in range(0,gaussPoints.shape[0]):
        ri = gaussPoints[i,0]
        si = gaussPoints[i,1]
        wi = gaussWeights[i]
        N, Bf,Bc, detJ = GetShapeFunction(elemType,ri, si, xi, yi)
        
        m11 = np.matmul(Bf.transpose(),Df)
        m12 = np.matmul(m11,Bf)

        kLocal += wi*m12*detJ

    #Reduced integration for shear contribution
    for i in range(0,gaussPointsRed.shape[0]):
        ri = gaussPointsRed[i,0]
        si = gaussPointsRed[i,1]
        wi = gaussWeightsRed[i]
        N, Bf,Bc, detJ = GetShapeFunction(elemType,ri, si, xi, yi)

        m21 = np.matmul(Bc.transpose(),Dc) 
        m22 = np.matmul(m21,Bc)

        kLocal += wi*m22*detJ

    for i in range(0,gaussPointsRed.shape[0]):
        ri = gaussPointsRed[i,0]
        si = gaussPointsRed[i,1]
        wi = gaussWeightsRed[i]
        N, Bf,Bc, detJ = GetShapeFunction(elemType,ri, si, xi, yi)

        fLocal += wi*np.matmul(N.transpose(),p)*detJ

    return kLocal, fLocal

def GetShapeFunction(elemType,ri, si, xi, yi):
    '''
    INPUT-->    ri: calculation point along the r-axis [float]
                si: calculation point along the s-axis [float]
                xi: coordinates of the element's nodes on the x-axis [1x3 array]
                yi: coordinates of the element's nodes on the y-axis [1x3 array]
    OUTPUT-->   N: values of the shape function matrix at ri, si [1x3 array]
                B: values of the deformation matrix at ri, si [1x3 array]
                detJ: Determinant of the Jacobian matrix [float]
    '''

    nodeCoordinates = np.zeros((elemType, 2))
    nodeCoordinates[:,0]=xi
    nodeCoordinates[:,1]=yi

    if elemType == 3:
        # Define shape functions
        N1 = lambda r, s: r
        N2 = lambda r, s: s
        N3 = lambda r, s: 1-r-s
        
        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s)]
        Nval=Nfun(ri, si)
        N=np.zeros((3, 3*3))
        N[0, 0:N.shape[1]:3]=Nval
        N[1, 1:N.shape[1]:3]=Nval
        N[2, 2:N.shape[1]:3]=Nval



        # Define shape function derivatives, derive deformation matrix
        N1r = lambda r, s: 1
        N2r = lambda r, s: 0
        N3r = lambda r, s: -1

        N1s = lambda r, s: 0
        N2s = lambda r, s: 1
        N3s = lambda r, s: -1

        NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)]])
        NrsVal=NrsFun(ri,si)

        # Jacobian matrix
        J=np.matmul(nodeCoordinates.transpose(), NrsVal)
        detJ = np.linalg.det(J)

        Bf=np.zeros((3,3*3))
        NrsVal = np.matmul(NrsVal,np.linalg.inv(J))
        Bf[0,1::3]=NrsVal[:,0]
        Bf[1,2::3]=NrsVal[:,1]
        Bf[2,1::3]=NrsVal[:,1]
        Bf[2,2::3]=NrsVal[:,0]

        Bc=np.zeros((2,3*3))
        Bc[0,0::3]=NrsVal[:,0]
        Bc[0,1::3]=Nval
        Bc[1,0::3]=NrsVal[:,1]
        Bc[1,2::3]=Nval
        return N, Bf,Bc, detJ
    
    elif elemType==4:
        # Define shape functions
        N1 = lambda r, s: 0.25*(1-r)*(1-s)
        N2 = lambda r, s: 0.25*(1+r)*(1-s)
        N3 = lambda r, s: 0.25*(1+r)*(1+s)
        N4 = lambda r, s: 0.25*(1-r)*(1+s)

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s), N4(r,s)]
        Nval=Nfun(ri, si)
        N=np.zeros((3, 3*elemType))
        N[0, 0::3]=Nval
        N[1, 1::3]=Nval
        N[2, 2::3]=Nval

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
        NrsVal=NrsFun(ri,si)

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
        Bc[0,1::3]=Nval
        Bc[1,0::3]=NrsVal[:,1]
        Bc[1,2::3]=Nval
        return N, Bf,Bc, detJ

class Result:
    '''
        class which stores all model output

    '''
    def __init__(self, outPos, wVert, xRot, yRot, resultsScale = (1,1)):
        self.outPos = outPos
        self.wVert = wVert*resultsScale[0]
        self.xRot = xRot
        self.yRot = yRot
        self.bendingMoments = None
        self.internalForcesPositions = None
        self.shearForces = None
        self.resultsScale = resultsScale

        z=np.abs(wVert)
        iMax = np.argmax(z)
        self.wMax = (outPos[iMax,0],outPos[iMax,1], self.wVert[iMax]*1000)

#unit testing
if __name__ == "__main__":
    # %matplotlib
    import matplotlib.pyplot as plt

    # import sample plate and show it
    import sampleGeometry
    from displayModel import *
    sampleModel=sampleGeometry.importModel()
    plotInputGeometry(sampleModel)

    # create mesh
    from generateMesh import *
    meshInput1=MeshInput(showGmshMesh=False, elementType='QUAD', nEdgeNodes=33)
    #meshInput1=MeshInput(showGmshMesh=True, elementType='QUAD', meshSize=5e-2)
    generateMesh(sampleModel, meshInput1)

    # compute
    from solveModel import *
    solveModel(sampleModel)
    # display results
    plotResults(sampleModel)

    coords = sampleModel.results.outPos
    z = sampleModel.results.wVert
    x = coords[:,0]
    y=coords[:,1]
    fig = plt.figure()

    ax = fig.gca(projection='3d')
    t=10
    ax.plot_trisurf(x*10,y*10,z,cmap=cm.jet)


    def axisEqual3D(ax):
        extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        sz = extents[:,1] - extents[:,0]
        centers = np.mean(extents, axis=1)
        # maxsize = max(abs(sz))
        r = sz[2]/2*2.6
        for ctr, dim in zip([centers[2]], 'z'):
            getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    # extents = np.array([getattr(ax, 'get_zlim'))
    # sz = extents[:,1] - extents[:,0]
    # centers = np.mean(extents, axis=1)
    # getattr(ax, 'set_zlim')()
    axisEqual3D(ax)
    # Hide grid lines
    ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    plt.show()

    # from mpl_toolkits.mplot3d import Axes3D # 3D plot
    # x= sampleModel.results.bendingMomentsPositions[:,0]
    # y = sampleModel.results.bendingMomentsPositions[:,1]
    # z= sampleModel.results.bendingMoments[:,0]
    # #z= sampleModel.results.shearForce[:,0]

    # # print('x, y: ',sampleModel.results.bendingMomentsPositions)
    # # print(z[0::10])

    # fig = plt.figure(

    # ax = fig.gca(projection='3d')
    # ax.plot_trisurf(x,y,z,cmap=cm.jet)
    # plt.show()


