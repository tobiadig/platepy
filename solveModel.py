''' Module Information
-----------------------------------------------------------
Purpose of module: Computes the equilibrium solution of a plate model using finite elements and stores the results.
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''
#Basic modules
import numpy as np
# import pandas as pd
# from scipy import sparse
# from scipy.sparse import linalg
from scipy.linalg import ldl
from scipy.linalg import solve
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

# Import custom functions from other modules
from shapeFunctions import *
from localMatrixes import *
from internalForces import *
from slicingFunctions import *
from rotationMatrix import *

# for debug purposes
from tqdm import tqdm

def solveModel(self, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1,\
    internalForcePosition = 'nodes', smoothedValues = False, solveMethod = 'sparse', computeMoments=True, kBendingResistance = 1):
    ''' Given a plateModel object with an initialized mesh, this function computes displacements, rotations and internal forces at
        each node.\n
            Input: \n
            * self: plateModel object. \n
            * resultsScaleIntForces = (1, 1): Before being displayed, the computed beding moments are multiplied by resultsScaleIntForces[0]
                and the shear forces by resultsScaleIntForces[1]. \n
            * resultsScaleVertDisp = 1: Before being displayed, the computed displacements are multiplied by resultsScaleVertDisp. \n
            * internalForcePosition = 'nodes': String defining the desired position where the internal forces should be computed.
            Possible values are "nodes" (default), "center" and "intPoints" for the positions used in the Gauss quadratur.\n
            * smoothedValues = False: Experimental. If True, the values of the shear forces by displacement based elements are smoothed according to
            the values at the Gauss points. \n
            * solveMethod = 'sparse': select the algorithm used to solve the equilibrium equation. "cho" has to be used by 
            downstand beams (in order to solve non-positive definite matrix systems), "sparse" (default) in all other cases. \n
            * computeMoments = True: Deactivates the computation of internal forces. For debug purposes. \n
            * kBendingResistance = 1: 1/tan(alpha), to compute the plate bending resistance according to the SIA 262 swiss norm. \n
            Return: \n
            * outPos: (nNodes x 2) numpy matrix, with the x-y coordinates of the result points for the displacements. \n
            * values: (nNodes x 3) numpy matrix, with the values of (vertical displacement, rotation 1, rotation 2) at each position in outPos.
    '''
    # Loop over elements and assemble stiffness matrices
    nodes=self.mesh.nodesArray
    nNodesTotal = nodes.shape[0]
    nodesRotations = self.mesh.nodesRotations # both dataframes
    elementsList = self.mesh.elementsList
    BCs = self.mesh.BCs
    nGDofs = nodes.shape[0]*3
    p=self.loads[0]   

    platesList = self.plates
    downStandBeamsList = self.downStandBeams
    modelMesh = self.mesh

    sparseGlobalMatrix, sparseForceGlobal, discartedDOFs = getGlobalStiffnesAndForce(elementsList,platesList,downStandBeamsList, nodesRotations, modelMesh,p, nNodesTotal)
    # apply boundary conditions
    elementType = elementsList[0].type

    fDofsInt, rDofsBool,keepedDisplacements = getFreeDOFvector(BCs, nGDofs,elementType,discartedDOFs)

    kMatFree = sparseGlobalMatrix[fDofsInt,:]
    kMatFree = kMatFree[:,fDofsInt]
    fVecFree = sparseForceGlobal[fDofsInt]

    if solveMethod == 'cho':
        AmatList = self.mesh.AmatList
        A = assembleSystemMFCmatrix(AmatList)
        try:
            A=A[:,fDofsInt]
        except:
            raise ValueError('Please change solveMethod from "cho" to "sparse"')

        M, rightSide = getMmatrixAndRightSide(A, kMatFree,fVecFree)
        # M is the matrix used to solve the system using lagrangian MPC. It has form [[K, A.T],[A, 0]]

        lu, d, _ = ldl(M)
        y = solve(lu, rightSide)
        x=solve(d@lu.T, y)

        nConstraints = A.shape[0]
        Uf = x[0:-nConstraints]
        uGlob=np.zeros((nGDofs,1))
        uGlob[fDofsInt]=np.expand_dims(Uf, axis=1)

    elif solveMethod == 'sparse':
        Uf = spsolve(kMatFree,fVecFree)
        Uf=Uf.reshape(-1,1)
        uGlob=np.zeros((nGDofs,1))
        uGlob[fDofsInt]=Uf

    uGlobSparse = csr_matrix(uGlob)
    globalForce = (sparseGlobalMatrix*uGlobSparse).toarray()

    nodes = self.mesh.nodesArray.to_numpy()
    if solveMethod == 'cho':
        uDownStandBeam = uGlob
        uGlob = uGlob[:-nConstraints,:]
        nodes = nodes[:int(-nConstraints/3),:]

    outPos = np.zeros((nodes.shape[0],2))
    values = np.zeros((nodes.shape[0],3))
    outPos[:,0:2]=nodes[:,0:2]
    keepedDisplacements = keepedDisplacements[0::3]
    vDisp = uGlob[0::3,0]
    if solveMethod != 'cho':
        vDisp = vDisp[keepedDisplacements]
        outPos = outPos[keepedDisplacements, :]  #TODO: only works for vertical displacements!
    values[:,1]=uGlob[1::3,0]
    values[:,2]=uGlob[2::3,0]
    
    self.results = Result(outPos,vDisp, values[:,1], values[:,2],resultsScale=(resultsScaleVertDisp,resultsScaleIntForces))
    self.results.uGlobPlate = uGlob

    #compute MOMENTS
    if computeMoments:
        nodesArray = self.mesh.nodesArray
        mitcList = self.mesh.plateElementsList

        bendingMoments, shearForces, internalForcesPositions = getInternalForces(mitcList,uGlob,internalForcePosition,nodesArray, smoothedValues)

        self.results.bendingMoments=bendingMoments*resultsScaleIntForces[0]
        self.results.internalForcesPositions=internalForcesPositions
        self.results.shearForces = shearForces*resultsScaleIntForces[1]
        self.results.resultsScaleVertDisp = resultsScaleVertDisp
        self.results.bendingResistance = getBendingResistance(bendingMoments,kBendingResistance)

        if len(self.downStandBeams) > 0:
            uzList = self.downStandBeams[0].elementsList
            Nforces, Vforces, Mforces, internalForcesPositions = getInternalForcesDSB(uzList,uDownStandBeam,internalForcePosition, self.downStandBeams[0])
            self.results.bendingMomentsDSB=Mforces*resultsScaleIntForces[0]
            self.results.internalForcesPositionsDSB=internalForcesPositions
            self.results.shearForcesDSB = Vforces*resultsScaleIntForces[1]
            self.results.normalForcesDSB = Nforces*resultsScaleIntForces[1]

    return outPos, values

def getGlobalStiffnesAndForce(elementsList,platesList,downStandBeamsList, nodesRotations, modelMesh,p,nNodesTotal):
    nSparseData = len(elementsList)*(9*3)**2
    rowsForStiffnessSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForStiffnessSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForStiffnessSparseMatrix = np.zeros(nSparseData)
    rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForForceSparseMatrix = np.zeros(nSparseData)
    startIndexStiffness = 0
    startIndexForce = 0
    nElements = len(elementsList)
    discartedDOFs = np.zeros(nElements, dtype=int)

    print('Assembling stiffness matrix')
    for k,element in enumerate(tqdm(elementsList)):
        elementType=element.type
        elementIntegration = element.integration
        plateOfTheElement = element.whichPlate
        elemNodes = element.connectivity
        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        nNodes = element.nNodes
        elemNodesRotations = nodesRotations.loc[elemNodes].to_numpy()
        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        if elementType!='timo':
            Df = platesList[plateOfTheElement].Df
            Dc = platesList[plateOfTheElement].Dc
            kLocalNotRotated,fLocal = GetLocalMatrix(xi, yi, Df,Dc,p,nNodes , elementType, elementIntegration)
        else:
            Emod = downStandBeamsList[0].body.eModule
            Gmod = downStandBeamsList[0].body.gModule
            crossA = downStandBeamsList[0].crossSection.A
            crossI = downStandBeamsList[0].crossSection.Iy
            Dc =Emod*crossA
            Db = Emod*crossI
            Ds = 5/6*Gmod*crossA
            kLocalNotRotated, fLocal = gettimoBeamMatrix(xi, yi,Dc, Db, Ds, 0, nNodes)

        # if the load is a line load IGNORE fLocal (set to zero), the force vector will be calculated in the next loop
        # bad solution, hopefully it works #TODO: adjust forces
        if p.case != "area":
            fLocal = np.zeros((fLocal.size,1))

        R = getRotationMatrix(elementType, elemNodesRotations) 

        element.rotationMatrix = R
        if elementType!='timo':
            modelMesh.plateElementsList[k].rotationMatrix = R

        # #rotate stiffness matrix
        kTemp = np.matmul(kLocalNotRotated, R)
        kLocal = np.matmul(R.transpose(), kTemp)
        # kLocal = kLocalNotRotated


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

    #if line load, assemble HERE load vector
    if p.case == 'line':
        
        rowsForForceSparseMatrix,dataForForceSparseMatrix = getLineLoadForceVector(p,nSparseData,elemNodesRotations)

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
    sparseGlobalMatrix = csr_matrix((dataForStiffnessSparseMatrix,(rowsForStiffnessSparseMatrix,columnsForStiffnessSparseMatrix)))
    sparseForceGlobal = csr_matrix((dataForForceSparseMatrix,(rowsForForceSparseMatrix,columnsForForceSparseMatrix)), shape=(nNodesTotal*3,1))

    return sparseGlobalMatrix, sparseForceGlobal, discartedDOFs

def getLineLoadForceVector(p,nSparseData,elemNodesRotations):
    elements1DList = p.elements1DList
    rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForForceSparseMatrix = np.zeros(nSparseData)
    startIndexForce = 0

    for element in elements1DList:
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
        # gaussPoints, gaussWeights =  getGaussQuadrature('linear',2)
        # N, Bb,Bs, detJ = getShapeFunctionForElementType(elementType,ri, si, xi, yi)
        fLocal = np.zeros(3*nNodes)
        # for i in range(0, gaussPoints.shape[0]):
        #     if fLocal == None:
        #         ri = 1
        #         si= gaussPoints

        if nNodes ==2:
            fLocal[0:3] = p.magnitude*L/2
            fLocal[3:] = p.magnitude*L/2
        elif nNodes == 3:
            fLocal[0:3] = p.magnitude*L*(1/6)
            fLocal[3:6] = p.magnitude*L*(1/6)
            fLocal[6:] = p.magnitude*L*(2/3)
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

        rows = np.zeros((3*nNodes)**2,dtype=int)
        columns = np.zeros((3*nNodes)**2,dtype=int)
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
    return rowsForForceSparseMatrix,dataForForceSparseMatrix

def getFreeDOFvector(BCs, nGDofs,elementType,discartedDOFs):
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
    rDofsInt = np.array(range(0,nGDofs))[rDofsBool]

    #remove discarted nodes
    allDofs =np.arange(0,nGDofs)
    fDofsBool = np.invert(rDofsBool)
    fDofsInt = allDofs[fDofsBool]

    return fDofsInt, rDofsBool,keepedDisplacements

def assembleSystemMFCmatrix(AmatList):
    A = []
    for myA in AmatList:
        if len(A) == 0:
            A=myA
        else:
            A= np.concatenate((A, myA))
    return A

def getMmatrixAndRightSide(A, kMatFree,fVecFree):
    nConstraints = A.shape[0]
    nFreeDofs = A.shape[1]
    kMatFreeNp = kMatFree.toarray()
    rightSide = np.zeros(nFreeDofs+nConstraints)
    rightSide[0:nFreeDofs]=fVecFree.toarray()[:,0]
    M = block_diag(kMatFreeNp, np.zeros((nConstraints,nConstraints)))
    M[-nConstraints:,0:nFreeDofs] = A
    M[0:nFreeDofs, -nConstraints:] = A.T

    return M, rightSide

def getBendingResistance(bendingMoments,kBendingResistance):
    bendingResistance = np.zeros((bendingMoments.shape[0],2,3))
    bendingResistance[:,0,0] = bendingMoments[:,0] + kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,1,0] = bendingMoments[:,1] + 1/kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,0,1] = -bendingMoments[:,0] + kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,1,1] = -bendingMoments[:,1] + 1/kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,0,2] = np.abs(bendingMoments[:,0]) + kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,1,2] = np.abs(bendingMoments[:,1]) + 1/kBendingResistance*np.abs(bendingMoments[:,2])
    return bendingResistance

































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
        self.resultsScaleVertDisp = resultsScale[0]

        self.bendingResistance = None

        self.bendingMomentsDSB=None
        self.internalForcesPositionsDSB=None
        self.shearForcesDSB = None
        self.normalForcesDSB = None

        self.schnittList = {}

        self.uGlobPlate = None

        z=np.abs(wVert)
        iMax = np.argmax(z)
        self.wMax = (outPos[iMax,0],outPos[iMax,1], self.wVert[iMax])
