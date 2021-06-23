'''
Computes the equilibrium solution of a plate model using finite elements and stores the results.
'''
#Basic modules
import numpy as np
from scipy.linalg import ldl
from scipy.linalg import solve
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from scipy.integrate import trapezoid
# Import custom functions from other modules
from ._shapeFunctions import *
from ._localMatrixes import *
from ._internalForces import *
from ._slicingFunctions import *
from ._rotationMatrix import *
from .displayModel import plotInputGeometry
import gmsh
# for debug purposes
from tqdm import tqdm

def solveModel(self, resultsScales = (1, 1, 1),\
    internalForcePosition = 'center', smoothedValues = False, computeMoments=True, kBendingResistance = 1):
    ''' Given a plateModel object with an initialized mesh, this function computes displacements, rotations and internal forces at
        each node.
        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~
        * **self**: plateModel object. 
        * **resultsScales = (1e-3, 1, 1)**: Before being displayed the computed displacements are multiplied by resultsScales[0], the computed beding moments are multiplied by resultsScales[1]
            and the shear forces by resultsScales[2]. 
        * **internalForcePosition = 'center'**: String defining the desired position where the internal forces should be computed.
        Possible values are "center" (default), "nodes" and "intPoints" for the positions used in the Gauss quadrature.
        * **smoothedValues = False**: Experimental. If True, the values of the shear forces by displacement based elements are smoothed according to
        the values at the Gauss points. 
        * **computeMoments = True**: Deactivates the computation of internal forces. For debug purposes. 
        * **kBendingResistance = 1**: 1/tan(alpha), to compute the plate bending resistance according to the SIA 262 swiss norm. 
        ~~~~~~~~~~~~~~~~~~~
        RETURN
        ~~~~~~~~~~~~~~~~~~~
        * ResultsDictionary: Dictionary with Result objects.
    '''

    # chose the right solving algorithm based on the presence or not of downstand beams
    if len(self.downStandBeams)>0:
        solveMethod='cho'
    else:
        solveMethod = 'sparse'

    # re-import mesh information
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

    # for elastic vertical displacement of walls (experimental)
    elasticallySupportedNodes = self.mesh.elasticallySupportedNodes
    if len(self.walls)>0:
        if self.walls[0].isElasticallySupported:
            wallStiffness = self.walls[0].body.eModule / self.walls[0].high * self.walls[0].thickness *0.004 # what is the characteristic length?
        else:
            wallStiffness = None
    else:
        wallStiffness = None

    sparseGlobalMatrix, sparseForceGlobal, discartedDOFs = _getGlobalStiffnesAndForce(elementsList,platesList,downStandBeamsList, nodesRotations, modelMesh,p, nNodesTotal, elasticallySupportedNodes,wallStiffness)
    elementType = elementsList[0].type
    fDofsInt, rDofsBool,keepedDisplacements = _getFreeDOFvector(BCs, nGDofs,elementType,discartedDOFs)

    kMatFree = sparseGlobalMatrix[fDofsInt,:]
    kMatFree = kMatFree[:,fDofsInt]
    fVecFree = sparseForceGlobal[fDofsInt]

    if solveMethod == 'cho':
        AmatList = self.mesh.AmatList
        A = _assembleSystemMFCmatrix(AmatList)
        try:
            A=A[:,fDofsInt]
        except:
            raise ValueError('Please change solveMethod from "cho" to "sparse"')

        M, rightSide = _getMmatrixAndRightSide(A, kMatFree,fVecFree)
        # M is the matrix used to solve the system using lagrangian MPC. It has form [[K, A.T],[A, 0]]
        
        lu, d, _ = ldl(M)
        y = solve(lu, rightSide)
        x = solve(d@lu.T, y)
        
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
        outPos = outPos[keepedDisplacements, :] 
    values[:,1]=uGlob[1::3,0]
    values[:,2]=uGlob[2::3,0]

    resultsDictionary = {}
    resultsDictionary['vDisp'] = _Result(outPos[:,0], outPos[:,1], vDisp, resultsScales[0])
    resultsDictionary['xRot'] = _Result(outPos[:,0], outPos[:,1], values[:,1], resultsScales[0])
    resultsDictionary['yRot'] = _Result(outPos[:,0], outPos[:,1], values[:,2], resultsScales[0])
    self.resultsInformation = ResultsInformation()
    self.resultsInformation.uGlobPlate = uGlob

    #compute internal forces
    if computeMoments:
        nodesArray = self.mesh.nodesArray
        mitcList = self.mesh.plateElementsList

        bendingMoments, shearForces, internalForcesPositions = getInternalForces(mitcList,uGlob,internalForcePosition,nodesArray, smoothedValues)
        xPos = internalForcesPositions[:,0]
        yPos = internalForcesPositions[:,1]

        resultsDictionary['Mx'] = _Result(xPos,yPos,bendingMoments[:,0], resultsScales[1])
        resultsDictionary['My'] = _Result(xPos,yPos,bendingMoments[:,1], resultsScales[1])
        resultsDictionary['Mxy'] = _Result(xPos,yPos,bendingMoments[:,2], resultsScales[1])
        resultsDictionary['Vx'] = _Result(xPos,yPos,shearForces[:,0], resultsScales[2])
        resultsDictionary['Vy'] = _Result(xPos,yPos,shearForces[:,1], resultsScales[2])

        # store required bending resistance according to SIA262
        bendingResistance = _getBendingResistance(bendingMoments,kBendingResistance)
        resultsDictionary['Mx_Rd_+'] = _Result(xPos,yPos,bendingResistance[:,0,0], resultsScales[1])
        resultsDictionary['My_Rd_+'] = _Result(xPos,yPos,bendingResistance[:,1,0], resultsScales[1])

        resultsDictionary['Mx_Rd_-'] = _Result(xPos,yPos,bendingResistance[:,0,1], resultsScales[1])
        resultsDictionary['My_Rd_-'] = _Result(xPos,yPos,bendingResistance[:,1,1], resultsScales[1])

        resultsDictionary['Mx_Rd_max'] = _Result(xPos,yPos,bendingResistance[:,0,2], resultsScales[1])
        resultsDictionary['My_Rd_max'] = _Result(xPos,yPos,bendingResistance[:,1,2], resultsScales[1])

        if len(self.downStandBeams) > 0:
            uzList = self.downStandBeams[0].elementsList
            Nforces, Vforces, Mforces, internalForcesPositions = getInternalForcesDSB(uzList,uDownStandBeam,internalForcePosition, self.downStandBeams[0])
            resultsDictionary['N_DSB'] = _Result(internalForcesPositions[:,0], internalForcesPositions[:,1], Nforces, resultsScales[2])
            resultsDictionary['V_DSB'] = _Result(internalForcesPositions[:,0], internalForcesPositions[:,1], Vforces, resultsScales[2])
            resultsDictionary['N_DSB'] = _Result(internalForcesPositions[:,0], internalForcesPositions[:,1], Mforces, resultsScales[1])

    self.resultsInformation.resultsDictionary = resultsDictionary
    return resultsDictionary

def _getGlobalStiffnesAndForce(elementsList,platesList,downStandBeamsList, nodesRotations, modelMesh,p,nNodesTotal, elasticallySupportedNodes,wallStiffness):
    ''' Computes global stifness matrix and global force matrixes in sparse form.\n
        Input: \n
            * elementsList: List of element objects. \n
            * platesList: List of plate obkects. \n
            * downStandBeamsList: List of downStandBeam objects. \n
            * nodesRotations: Pandas dataframe where indexes are node tags (assigned by gmsh), values are the rotations in radians. \n
            * modelMesh: mesh object of the plateModel object. \n
            * p: load object.
            * nNodesTotal: nodesArray.shape[0], number of nodes in the entire model.\n
            * elasticallySupportedNodes: Array containing the nodes of the walls vertically elastic. \n
            * wallStiffness: value of the vertical stiffness of the wall.
        Return: \n
            * sparseGlobalMatrix: Scipy sparse matrix of the system's stiffness.
            * sparseForceGlobal: Scipy sparse matrix of the system's force vector.
            * discartedDOFs: List of the degrees of freedom removed (in the case of the MITC-9 element, since only 8 vertical displacements are used).
    '''
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

        # If the force is not on the entire plate's area, the local force vector has to be deleted and re-calculated
        if p.case != "area":
            fLocal = np.zeros((fLocal.size,1))

        R = getRotationMatrix(elementType, elemNodesRotations) 
        element.rotationMatrix = R
        if elementType!='timo':
            modelMesh.plateElementsList[k].rotationMatrix = R
        # rotate stiffness matrix
        kTemp = np.matmul(kLocalNotRotated, R)
        kLocal = np.matmul(R.transpose(), kTemp)
        if elementType == 'timo':
            kLocal = kLocalNotRotated

        nMatrixDofs = kLocal.size
        kCoeff, discartedDOF = getKCoeff(elementType, coherentElemNodes)
        if discartedDOF != None:
            discartedDOFs[k]=discartedDOF

        rows, columns = getRowsColumns(kCoeff, nMatrixDofs)

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

    # add elastically supported nodes:
    if len(elasticallySupportedNodes)>0:
        rows = elasticallySupportedNodes*3
        columns = elasticallySupportedNodes*3
        rowsForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = rows[:,0]
        columnsForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = columns[:,0]
        dataForStiffnessSparseMatrix[startIndexStiffness:startIndexStiffness+rows.size] = np.ones(rows.size)*wallStiffness

    #if line load, assemble HERE load vector
    if p.case == 'line':
        rowsForForceSparseMatrix,dataForForceSparseMatrix = _getLineLoadForceVector(p,nSparseData,elemNodesRotations)
    elif p.case == 'nodes':
        rowsForForceSparseMatrix,dataForForceSparseMatrix = _getNodesLoadForceVector(p, nSparseData)
    elif p.case == 'point':
        rowsForForceSparseMatrix,dataForForceSparseMatrix = _getPointLoadForceVector(p, nSparseData)

    # create global matrixes
    sparseGlobalMatrix = csr_matrix((dataForStiffnessSparseMatrix,(rowsForStiffnessSparseMatrix,columnsForStiffnessSparseMatrix)))
    sparseForceGlobal = csr_matrix((dataForForceSparseMatrix,(rowsForForceSparseMatrix,columnsForForceSparseMatrix)), shape=(nNodesTotal*3,1))
    return sparseGlobalMatrix, sparseForceGlobal, discartedDOFs

def _getLineLoadForceVector(p,nSparseData,elemNodesRotations):
    ''' Computes informations required to create the sparse force vector in case of line loads. \n
        Input: \n
            * p: Load object. \n
            * nSparseData: Total number of the non-zero entries in the global stiffness matrix. = len(elementsList)*(9*3)**2. \n
            * elemNodesRotations: Rotation of the node appartaineng at the element. \n
        Return: \n
            * rowsForForceSparseMatrix: numpy vector with the indexes of the entries in dataForForceSparseMatrix. \n
            * dataForForceSparseMatrix: numpy vector with the values of the entries in the global sparse force vector. \n
    '''
    elements1DList = p.elements1DList
    rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForForceSparseMatrix = np.zeros(nSparseData)
    startIndexForce = 0

    for element in elements1DList:
        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        nNodes=element.nNodes
        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        L = np.sqrt((xi[1]-xi[0])**2+(yi[1]-yi[0])**2)
        fLocal = np.zeros(3*nNodes)
        if nNodes ==2:
            fLocal[0:3] = p.magnitude*L/2
            fLocal[3:] = p.magnitude*L/2
        elif nNodes == 3:
            fLocal[0:3] = p.magnitude*L*(1/6)
            fLocal[3:6] = p.magnitude*L*(1/6)
            fLocal[6:] = p.magnitude*L*(2/3)
        # # create rotation matrix
        RiInit = False
        for dofRot in elemNodesRotations:
            if not(RiInit):
                R=rotMatrix(dofRot)
                RiInit=True
            else:
                R = block_diag(R, rotMatrix(dofRot))
        element.rotationMatrix = R
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

def _getNodesLoadForceVector(p, nSparseData):
    ''' Computes informations required to create the sparse force vector in case the user manually defines the node loads. \n
        Input: \n
            * p: Load object. \n
            * nSparseData: Total number of the non-zero entries in the global stiffness matrix. = len(elementsList)*(9*3)**2. \n
        Return: \n
            * rowsForForceSparseMatrix: numpy vector with the indexes of the entries in dataForForceSparseMatrix. \n
            * dataForForceSparseMatrix: numpy vector with the values of the entries in the global sparse force vector. \n
    '''
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
    return rowsForForceSparseMatrix,dataForForceSparseMatrix 

def _getPointLoadForceVector(p,nSparseData):
    ''' Computes informations required to create the sparse force vector in case of line loads. \n
        Input: \n
            * p: Load object. \n
            * nSparseData: Total number of the non-zero entries in the global stiffness matrix. = len(elementsList)*(9*3)**2. \n
        Return: \n
            * rowsForForceSparseMatrix: numpy vector with the indexes of the entries in dataForForceSparseMatrix. \n
            * dataForForceSparseMatrix: numpy vector with the values of the entries in the global sparse force vector. \n
    '''
    nodeWithLoad = p.pointLoadNode
    fLocal = p.magnitude
    rowsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    columnsForForceSparseMatrix = np.zeros(nSparseData, dtype=int)
    dataForForceSparseMatrix = np.zeros(nSparseData)
    startIndexForce = 0
    nNodes = 1

    coherentElemNodes = nodeWithLoad-1
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

    return rowsForForceSparseMatrix,dataForForceSparseMatrix

def _getFreeDOFvector(BCs, nGDofs,elementType,discartedDOFs):
    ''' Constructs the vector with the not-restrained degrees of freedom.\n
        Input: \n
            * BCs: n x 4 numpy array. n is the number of restrained nodes, the first column is the node tag, the other column represent the condition of the three DOFs (1 is blocked, 0 is free). \n
            * nGDofs: number of degrees of freedom in the system (=nNodes *3). \n
            * elementType: "DB" or "MITC".
            * discartedDOFs: numpy array with the tags of the DOFs discarted in the getGlobalStiffnesAndForce() function. \n 
        Return: \n
            * fDofsInt: Numpy array with the tags of the free DOFs. \n
            * rDofsBool: Numpy boolean array of length nGDofs. True if the relative DOF is restrained. \n
            * keepedDisplacements: DOFs which are not discarted in discartedDOFs.
    '''
    rDofsBool = np.zeros((nGDofs),dtype=bool)

    if discartedDOFs[0]!=0:
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

def _assembleSystemMFCmatrix(AmatList):
    ''' Constructs the MFC matrix for the system. \n
        Input: \n
            * AmatList: List of MFC matrix for the indivudual downstand beams. \n
        Return: \n
            * A: System MFC matrix of shape (nConstraints x nGDOFs)
    '''
    A = []
    for myA in AmatList:
        if len(A) == 0:
            A=myA
        else:
            A= np.concatenate((A, myA))
    return A

def _getMmatrixAndRightSide(A, kMatFree,fVecFree):
    ''' Returns left and right sides of the system of equations to solve. \n
            Input: \n
                * A: System MFC matrix of shape (nConstraints x nGDOFs). \n
                * kMatFree: Scipy sparse matrix of the system's stiffness. \n
                * fVecFree: Scipy sparse matrix of the system's force vector. \n
            Return: \n
                * M: Numpy matrix for the left side of the equation, has shape (nGDOFs + nConstraints,nGDOFs + nConstraints ). \n
                * rightSide: Numpy matrix for the right side of the equation, has shape (nGDOFs + nConstraints,). \n
    '''
    nConstraints = A.shape[0]
    nFreeDofs = A.shape[1]
    kMatFreeNp = kMatFree.toarray()
    rightSide = np.zeros(nFreeDofs+nConstraints)
    rightSide[0:nFreeDofs]=fVecFree.toarray()[:,0]
    M = block_diag(kMatFreeNp, np.zeros((nConstraints,nConstraints)))
    M[-nConstraints:,0:nFreeDofs] = A
    M[0:nFreeDofs, -nConstraints:] = A.T

    return M, rightSide

def _getBendingResistance(bendingMoments,kBendingResistance):
    '''Computes rquired bending resistance according to SIA 262.\n
        Input: \n
            * bendingMoments: Numpy array of shape (n evaluation points, 3). First column: mx, second column: my, third column: mxy. \n
            * kBendingResistance: 1/tan(alpha) according to SIA 262. \n
        Return: \n
            * bendingResistance: Numpy array of shape(n evaluation points,3,3). Dimension 2 are mx, my, mxy. Dimension 3 are positive resistance, negative resistance, maximal resistance.
    '''
    bendingResistance = np.zeros((bendingMoments.shape[0],2,3))
    bendingResistance[:,0,0] = bendingMoments[:,0] + kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,1,0] = bendingMoments[:,1] + 1/kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,0,1] = -bendingMoments[:,0] + kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,1,1] = -bendingMoments[:,1] + 1/kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,0,2] = np.abs(bendingMoments[:,0]) + kBendingResistance*np.abs(bendingMoments[:,2])
    bendingResistance[:,1,2] = np.abs(bendingMoments[:,1]) + 1/kBendingResistance*np.abs(bendingMoments[:,2])
    return bendingResistance

def computeBeamComponents(self, startCoord, endCoord, nEvaluationPoints,resultsScales = (1,1, 1),integrationWidth = 0, nIntegrationPoints =0):
    '''
    Computes the results over a custom-defined line.
    ~~~~~~~~~~~~~~~~~~~
    INPUT
    ~~~~~~~~~~~~~~~~~~~
    * **self**: PlateModel object.
    * **startCoord**: tuple with (x,y) coordinates defining the starting point of the line.
    * **endCoord**: tuple with (x,y) coordinates defining the ending point of the line. 
    * **resultsScales = (1e-3, 1, 1)**: Before being displayed the computed displacements are 
    multiplied by resultsScales[0], the computed beding moments are multiplied by resultsScales[1]
    and the shear forces by resultsScales[2]. 
    * **integrationWidth = 0**: if > 0, the values on a line normal to the main cut of length +-integrationWidth/2 are evaluated and 
    integrated to a single value, which will then be displayed. 
    * **nIntegrationPoint=0**: number of point for the integration normal to the direction of the main cut.
    ~~~~~~~~~~~~~~~~~~~
    RETURN
    ~~~~~~~~~~~~~~~~~~~
    * bendingMoments, shearForces, arrayEvaluationPoints
    '''
    uGlob = self.resultsInformation.uGlobPlate
    elementsList = self.mesh.elementsList
    arrayEvaluationPoints = np.zeros((nEvaluationPoints,2))
    arrayEvaluationPoints[:,0] = np.linspace(startCoord[0], endCoord[0], num=nEvaluationPoints)
    arrayEvaluationPoints[:,1] = np.linspace(startCoord[1], endCoord[1], num=nEvaluationPoints)
    bendingMoments = np.zeros((nEvaluationPoints,3))
    shearForces = np.zeros((nEvaluationPoints,2))
    verticalDisplacements = np.zeros(nEvaluationPoints)

    if integrationWidth ==0:
        for k,evaluationPoint in enumerate(arrayEvaluationPoints):
            elementTaggetEl, elementTypegetEl, nodeTagsgetEl, ugetEl, vgetEl, wgetEl = gmsh.model.mesh.getElementByCoordinates(evaluationPoint[0],evaluationPoint[1],0, dim=2)
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
            Df = self.plates[plateOfTheElement].Df
            Dc = self.plates[plateOfTheElement].Dc

            ri =-ugetEl
            si = -vgetEl

            N, Bb,Bs, detJ =getShapeFunctionForElementType(elementType,ri, si, xi, yi)

            tempDispl = N@vLoc
            verticalDisplacements[k] = tempDispl[0]
            bendingMoments[k,0:3] = np.matmul(Df,np.matmul(Bb, vLoc))[:,0]*-1
            shearForces[k,0:2] = np.matmul(Dc, np.matmul(Bs, vLoc))[:,0]*-1

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
                Df = self.plates[plateOfTheElement].Df
                Dc = self.plates[plateOfTheElement].Dc

                ri = -ugetEl
                si = -vgetEl
                N, Bb,Bs, detJ =getShapeFunctionForElementType(elementType,ri, si, xi, yi)
                integrationBendingMoments[j,0:3] = np.matmul(Df,np.matmul(Bb, vLoc))[:,0]*1
                integrationShearForces[j,0:2] = np.matmul(Dc, np.matmul(Bs, vLoc))[:,0]*1
            
            for j in range(0,3):
                arrayToIntegrate = integrationBendingMoments[:,j]
                bendingMoments[k,j]=trapezoid(arrayToIntegrate,dx=sampleWidth)

            for j in range(0,2):
                arrayToIntegrate = integrationShearForces[:,j]
                shearForces[k,j]=trapezoid(arrayToIntegrate,dx=sampleWidth)

    resultsDictionary = self.resultsInformation.resultsDictionary
    xPos = arrayEvaluationPoints[:,0]
    yPos = arrayEvaluationPoints[:,1]
    resultsDictionary['vDisp_line'] = _Result(xPos,yPos, verticalDisplacements, resultsScales[0])
    resultsDictionary['Mx_line'] = _Result(xPos,yPos,bendingMoments[:,0], resultsScales[1])
    resultsDictionary['My_line'] = _Result(xPos,yPos,bendingMoments[:,1], resultsScales[1])
    resultsDictionary['Mxy_line'] = _Result(xPos,yPos,bendingMoments[:,2], resultsScales[1])
    resultsDictionary['Vx_line'] = _Result(xPos,yPos,shearForces[:,0], resultsScales[2])
    resultsDictionary['Vy_line'] = _Result(xPos,yPos,shearForces[:,1], resultsScales[2])

    return bendingMoments, shearForces, arrayEvaluationPoints

def evaluateAtPoints(self, coords, displayPoints = False):
    '''
    Computes the displacement, rotation and internal forces at the location contained in coord.
    ~~~~~~~~~~~~~~~~~~~
    INPUT
    ~~~~~~~~~~~~~~~~~~~
    * **coords**: Numpy array of shape (nPoints, 2) with the x,y coordinates of the points to be evaluated. 
    * **displayPoints = False**: If True, plots the requested evaluation points on the inputGeometry. 
    ~~~~~~~~~~~~~~~~~~~
    RETURN
    ~~~~~~~~~~~~~~~~~~~
    * **verticalDisplacements,bendingMoments, shearForces**: length nPoints
    '''
    uGlob = self.resultsInformation.uGlobPlate
    nEvaluationPoints = coords.shape[0]
    elementsList = self.mesh.elementsList
    arrayEvaluationPoints = np.zeros((nEvaluationPoints,2))

    bendingMoments = np.zeros((nEvaluationPoints,3))
    shearForces = np.zeros((nEvaluationPoints,2))
    verticalDisplacements = np.zeros(nEvaluationPoints)
    resultsScaleVertDisp = self.resultsInformation.resultsDictionary['vDisp'].resultScale
    for k,evaluationPoint in enumerate(coords):
        elementTaggetEl, elementTypegetEl, nodeTagsgetEl, ugetEl, vgetEl, wgetEl = gmsh.model.mesh.getElementByCoordinates(evaluationPoint[0],evaluationPoint[1],0, dim=2)
        element = self.mesh.getElementByTagDictionary[elementTaggetEl]
        plateOfTheElement = element.whichPlate
        elementType = element.type
        elemNodes = element.connectivity
        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        nNodes=element.nNodes
        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        elementShape = len(xi)
        kCoeff, discartedDOF = getKCoeff(elementType, coherentElemNodes)

        vLoc = uGlob[kCoeff]
        Df = self.plates[plateOfTheElement].Df
        Dc = self.plates[plateOfTheElement].Dc
        if elementType == 'MITC':
            ri =-ugetEl
            si = -vgetEl
        else:
            ri =ugetEl
            si = vgetEl
        if elementType == 'DB' and nNodes ==4:
            ri = np.array([ri])
            si = np.array([si])
        N, Bb,Bs, detJ =getShapeFunctionForElementType(elementType,ri, si, xi, yi)
        tempDispl = N@vLoc
        if len(tempDispl.shape)>2:
            tempDispl = tempDispl[0,:,:]
        if (elementType == 'DB' and nNodes ==9) or (elementType == 'MITC' and nNodes==9):
            changeSign = -1
        else:
            changeSign = 1
        verticalDisplacements[k] = tempDispl[0]*resultsScaleVertDisp
        bendingMoments[k,0:3] = np.matmul(Df,np.matmul(Bb, vLoc))[:,0]*changeSign
        shearForces[k,0:2] = np.matmul(Dc, np.matmul(Bs, vLoc))[:,0]*changeSign

    if displayPoints:
        fig,outAx = plotInputGeometry(self)
        for k,evaluationPoint in enumerate(coords):
            outAx.scatter(evaluationPoint[0],evaluationPoint[1], facecolor='b', marker='.')
    return verticalDisplacements,bendingMoments, shearForces

class _Schnitt:
    def __init__(self, bendingMoments, shearForces,verticalDisplacements, arrayEvaluationPoints):
        ''' Class containing the results (displacements and internal forces) along a line cut.\n
        Input: \n
        * bendingMoments: Numpy array of shape (nEvaluationPoints, 3) \n
        * shearForces: Numpy array of shape (nEvaluationPoints, 2) \n
        * verticalDisplacements: Numpy array of shape (nEvaluationPoints, 1) \n
        * arrayEvaluationPoints: Numpy array of shape (nEvaluationPoints, 2) with the x y coordinates of the evaluation points. \n
        '''
        self.bendingMoments = bendingMoments
        self.shearForces = shearForces
        self.verticalDisplacements = verticalDisplacements
        self.arrayEvaluationPoints = arrayEvaluationPoints


class _Result:
    '''
    Attributes: \n
    * iMax, zMax
    * iMin, zMin
    * zMaxScaled, zMinScaled, zAbsScaled
    '''
    def __init__(self,x,y,z, resultScale):
        '''Class object to store results \n
        Inputs: \n
        * x \n
        * y \n
        * z \n
        * resultScale: Scaling factor for the display of the result. \n
        '''
        self.x = x
        self.y = y
        self.z = z
        iMax = np.argmax(z) 
        self.iMax = iMax
        self.zMax = z[iMax]
        iMin = np.argmin(z)
        self.iMin = iMin
        self.zMin = z[iMin]
        zAbs = np.abs(z)
        self.zAbsMax = np.max(zAbs)
        self.zMaxScaled = z[iMax]*resultScale
        self.zMinScaled = z[iMin]*resultScale
        self.zAbsScaled = np.max(zAbs)*resultScale
        self.resultScale = resultScale

class ResultsInformation:
    def __init__(self):
        '''Stores all kind of object regarding results.
        '''
        self.resultsDictionary = None
        self.bendingMoments = None
        self.internalForcesPositions = None
        self.shearForces = None
        self.bendingResistance = None
        self.bendingMomentsDSB=None
        self.internalForcesPositionsDSB=None
        self.shearForcesDSB = None
        self.normalForcesDSB = None
        self.schnittList = {}
        self.uGlobPlate = None



