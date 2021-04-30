''' Module Information
-----------------------------------------------------------
Purpose of module: computes the equilibrium solution using finite elements, stores the result as attribute
                    of a plateModel object
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''
#%% Basic modules
import numpy as np
import pandas as pd
import gmsh # To create CAD model and mesh
from scipy import sparse    # To work with sparse matrix
from scipy.sparse import linalg # Linear sparse matrices operations
from scipy.linalg import block_diag # to create the rotation matrix
from shapeFunctions import *
from localMatrixes import *
from internalForces import *
from slicingFunctions import *
from rotationMatrix import *

# for debug purposes
import time
from tqdm import tqdm

def solveModel(self, reducedIntegration = False, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1, elementDefinition=None, internalForcePosition = 'nodes'):
    ''' Input/Output descriptions
    ElemType:  or Quadratic or MITC + Reduced or Normal Integration
        self: PlateModel class, where the geometry and the mesh are initialized
        reducedIntegration: to manage the number of Gauss integration points
    '''
    temp = elementDefinition.split('-')
    elementType = temp[0]
    elementIntegration = temp[1]

    # Loop over elements and assemble stiffness matrices
    nodes=self.mesh.nodesArray
    nNodesTotal = nodes.shape[0]
    nodesRotations = self.mesh.nodesRotations # both dataframes
    # print('nodesRotation in solve model: ', nodesRotations)
    elementsList = self.mesh.elementsList
    nElements = len(elementsList)
    discartedDOFs = np.zeros(nElements, dtype=int)
    BCs = self.mesh.BCs
    nGDofs = nodes.shape[0]*3
    #initialize arrays for the creation of the sparse matrixes
    nSparseData = len(elementsList)*(9*3)**2
    rowsForStiffnessSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForStiffnessSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForStiffnessSparseMatrix = np.zeros(nSparseData)
    rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForForceSparseMatrix = np.zeros(nSparseData)
    p=self.loads[0]   

    startIndexStiffness = 0
    startIndexForce = 0
    k=0
    for element in tqdm(elementsList):
    # for element in elementsList:
        elemNodes = element.connectivity
        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        nNodes = element.nNodes
        elemNodesRotations = nodesRotations.loc[elemNodes].to_numpy()
        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        Df = self.plates[0].Df   
        Dc = self.plates[0].Dc

        kLocalNotRotated,fLocal = GetLocalMatrix(xi, yi, Df,Dc,p,nNodes , elementDefinition)

        # if the load is a line load IGNORE fLocal (set to zero), the force vector will be calculated in the next loop
        # bad solution, hopefully it works #TODO: adjust forces
        if p.case != "area":
            fLocal = np.zeros((3*nNodes,1))

        R = getRotationMatrix(elementType, elemNodesRotations)

        element.rotationMatrix = R

        # #rotate stiffness matrix
        kTemp = np.matmul(kLocalNotRotated, R)
        kLocal = np.matmul(R.transpose(), kTemp)
        nMatrixDofs = kLocal.size
        kCoeff, discartedDOF = getKCoeff(elementType, coherentElemNodes)
        if discartedDOF != None:
            discartedDOFs[k]=discartedDOF
        rows, columns = getRowsColumns(kCoeff, nMatrixDofs)

            # coefficients of the DOFs and assignment of the stiffness matrix / force vector
    
        # create vectors to assemble sparse matrixes
        rowsForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = rows
        columnsForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = columns
        dataForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = kLocal.flatten(order='F')
        startIndexStiffness += rows.size

        rowsForForceSparseMatrix[startIndexForce:startIndexForce+kCoeff.size] = kCoeff

        if fLocal.ndim == 1:
            dataForForceSparseMatrix[startIndexForce:startIndexForce+kCoeff.size] = fLocal[:]
        else:
            dataForForceSparseMatrix[startIndexForce:startIndexForce+kCoeff.size] = fLocal[:,0]

        startIndexForce += kCoeff.size
        k+=1

    #if line load, assemble HERE load vector
    if p.case == 'line':

        rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
        columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
        dataForForceSparseMatrix = np.zeros(nSparseData)
        startIndexForce = 0

        for element in p.elements1DList:
            elemNodes = element.connectivity

            coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
            nNodes=element.nNodes

            # elemNodesRotations = nodesRotations.loc[elemNodes].to_numpy()

            # xi=element.coordinates[:,0]
            # yi=element.coordinates[:,1]
            
            # Df = self.plates[0].Df   
            # Dc = self.plates[0].Dc

            # kLocalNotRotated,fLocal = GetLocalMatrix(xi, yi, Df,Dc,p, reducedIntegration)

            # # if the load is a line load IGNORE fLocal (set to zero), the force vector will be calculated in the next loop
            # # bad solution, hopefully it works
            # if p.case == "line":
            #     fLocal = np.zeros((3*nNodes,1))
            xi=element.coordinates[:,0]
            yi=element.coordinates[:,1]

            L = np.sqrt((xi[1]-xi[0])**2+(yi[1]-yi[0])**2)
            fLocal = np.zeros(3*nNodes)
            if nNodes ==2:
                fLocal[0:3] = p.magnitude*L/2
                fLocal[3:] = p.magnitude*L/2
            elif nNodes == 3:
                fLocal[0:3] = p.magnitude*L/4
                fLocal[3:6] = p.magnitude*L/2
                fLocal[6:] = p.magnitude*L/4
            # # create rotation matrix
            Ri = []
            RiInit = False
            for dofRot in elemNodesRotations:
                if not(RiInit):
                    R=rotMatrix(dofRot)
                    RiInit=True
                else:
                    R = block_diag(R, rotMatrix(dofRot))

            element.rotationMatrix = R
            # #rotate stiffness matrix
            # kTemp = np.matmul(kLocalNotRotated, R)
            # kLocal = np.matmul(R.transpose(), kTemp)
            # # kLocal = kLocalNotRotated
            # # coefficients of the DOFs and assignment of the stiffness matrix / force vector

            kCoeff = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeff[0+i::3]=coherentElemNodes*3+i

            rows = np.zeros(kLocal.size,dtype=int)
            columns = np.zeros(kLocal.size,dtype=int)
            c=0
            for j in kCoeff:
                for i in kCoeff:
                    rows[c] = i
                    columns[c] = j
                    c+=1

            # create vectors to assemble sparse matrixes
            rowsForForceSparseMatrix[startIndexForce:startIndexForce+kCoeff.size] = kCoeff
            dataForForceSparseMatrix[startIndexForce:startIndexForce+kCoeff.size] = fLocal[:]
            startIndexForce += kCoeff.size
            k+=1
    elif p.case == 'nodes':
        rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
        columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
        dataForForceSparseMatrix = np.zeros(nSparseData)
        startIndexForce = 0

        nNodes = p.nodePattern.shape[0]
        coherentNodes = p.nodePattern[:,0]

        for node in p.nodePattern:
            nodeTag = int(node[0]-1)
            rowsForForceSparseMatrix[startIndexForce:startIndexForce+3] = range(nodeTag*3,nodeTag*3+3)
            dataForForceSparseMatrix[startIndexForce:startIndexForce+3] = node[1:]
            startIndexForce += 3
            k+=1

    # create global matrixes
    sparseGlobalMatrix = sparse.csr_matrix((dataForStiffnessSparseMatrix,(rowsForStiffnessSparseMatrix,columnsForStiffnessSparseMatrix)))
    sparseForceGlobal = sparse.csr_matrix((dataForForceSparseMatrix,(rowsForForceSparseMatrix,columnsForForceSparseMatrix)), shape=(nNodesTotal*3,1))

    # np.savetxt('K.csv', sparseGlobalMatrix.toarray(), delimiter=",")
    # apply boundary conditions
    rDofsBool = np.zeros((nGDofs),dtype=bool)
    
    if elementType=='MITC9':
        keepedDisplacements = np.zeros((nGDofs),dtype=bool)
        keepedDisplacements[discartedDOFs] = np.ones((discartedDOFs.size), dtype=bool)
        keepedDisplacements = np.invert(keepedDisplacements)
        rDofsBool[discartedDOFs] = np.ones((discartedDOFs.size), dtype=bool)
    else:
        keepedDisplacements = np.zeros((nGDofs),dtype=bool)
        keepedDisplacements = np.invert(keepedDisplacements)

    for constraint in BCs:
        node=int(constraint[0])
        rDofsBool[node*3-3:node*3] = constraint[1:].astype(bool)

    #remove discarted nodes
    allDofs =np.arange(0,nGDofs)
    fDofsBool = np.invert(rDofsBool)
    fDofsInt = allDofs[fDofsBool]

    kMatFree = sparseGlobalMatrix[fDofsInt,:]
    kMatFree = kMatFree[:,fDofsInt]

    fVecFree = sparseForceGlobal[fDofsInt]
    # print('force vector: ', sparseForceGlobal.toarray())

    # SOLVE
    Uf = sparse.linalg.spsolve(kMatFree,fVecFree)
    Uf=Uf.reshape(-1,1)

    # global displacement and force vector
    
    uGlob=np.zeros((nGDofs,1))
    uGlob[fDofsInt]=Uf
    # if elementDegree == 'MITC9':
    # uGlob = uGlob[discartedDisplacements]
    
    uGlobSparse = sparse.csr_matrix(uGlob)

    # globalForce = np.matmul(sparseGlobalMatrix.toarray(),uGlob)
    # print(discartedDisplacements)
    # M1 = sparseGlobalMatrix[discartedDisplacements]
    # M1 = M1[:,discartedDisplacements]
    # print('m1size: ', M1.shape)
    # print('uGlob shape: ', uGlob.shape)
    globalForce = (sparseGlobalMatrix*uGlobSparse).toarray()

    # elaborate and store results
    reactionForces = globalForce[rDofsBool]
    nodes = self.mesh.nodesArray.to_numpy()

    outPos = np.zeros((nodes.shape[0],2))
    values = np.zeros((nodes.shape[0],3))
    outPos[:,0:2]=nodes[:,0:2]
    keepedDisplacements = keepedDisplacements[0::3]
    vDisp = uGlob[0::3,0]

    vDisp = vDisp[keepedDisplacements]
    outPos = outPos[keepedDisplacements, :]  #TODO: only works for vertical displacements!
    values[:,1]=uGlob[1::3,0]
    values[:,2]=uGlob[2::3,0]
    self.results = Result(outPos,vDisp, values[:,1], values[:,2],resultsScale=(resultsScaleVertDisp,resultsScaleIntForces))

    #compute MOMENTS
    Df = self.plates[0].Df 
    Dc = self.plates[0].Dc
    nodesArray = self.mesh.nodesArray
    bendingMoments, shearForces, internalForcesPositions = getInternalForces(elementType,elementsList,uGlob,internalForcePosition, Df, Dc, nodesArray)
    self.results.bendingMoments=bendingMoments*resultsScaleIntForces[0]
    self.results.internalForcesPositions=internalForcesPositions
    self.results.shearForces = shearForces*resultsScaleIntForces[1]
    return outPos, values

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
        self.wMax = (outPos[iMax,0],outPos[iMax,1], self.wVert[iMax])

if __name__ == "__main__":

    xi = np.array([5. , 0. , 0. , 5. , 2.5, 0. , 2.5, 5. , 2.5])
    yi = np.array([5. , 5. , 0. , 0. , 5. , 2.5, 0. , 2.5, 2.5])

    ri = -0.7745966692414834
    si = -0.7745966692414834

    N, Nr, Ns = shapeFun9(ri, si)
    N, Nr, Ns = shapeFun8(ri, si)

    nothing = 0