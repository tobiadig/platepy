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

# for debug purposes
import time
from tqdm import tqdm

def solveModel(self, reducedIntegration = False, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1, elemType=None, internalForcePosition = 'nodes'):
    ''' Input/Output descriptions
    ElemType:  or Quadratic or MITC + Reduced or Normal Integration
        self: PlateModel class, where the geometry and the mesh are initialized
        reducedIntegration: to manage the number of Gauss integration points
    '''
    temp = elemType.split('-')
    elementDegree = temp[0]
    elementIntegration = temp[1]

    # Loop over elements and assemble stiffness matrices
    nodes=self.mesh.nodesArray

    nNodesTotal = nodes.shape[0]
    nodesRotations = self.mesh.nodesRotations # both dataframes
    elementsList = self.mesh.elementsList
    nElements = len(elementsList)
    discartedDOFs = np.zeros(nElements)
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
        # iterate over each element, get stiffness and creates basis for the creation of the global stiffness matrix
        elemNodes = element.connectivity
        # print('connectivity: ', elemNodes)

        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]

        nNodes = element.nNodes
        elemNodesRotations = nodesRotations.loc[elemNodes].to_numpy()

        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]
        elementShape = len(xi)
        # print('xi: ', xi)
        # print('yi: ', yi)
        # print('elementShape: ', elementShape)
        # print('elemType: ', elemType)

        Df = self.plates[0].Df   
        Dc = self.plates[0].Dc

        kLocalNotRotated,fLocal = GetLocalMatrix(xi, yi, Df,Dc,p,elementShape , elemType)

        # if the load is a line load IGNORE fLocal (set to zero), the force vector will be calculated in the next loop
        # bad solution, hopefully it works
        if p.case != "area":
            fLocal = np.zeros((3*nNodes,1))

        if elementDegree != 'MITC9':
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

            # #rotate stiffness matrix
            kTemp = np.matmul(kLocalNotRotated, R)
            kLocal = np.matmul(R.transpose(), kTemp)
            # print('kLocal: ', pd.DataFrame(kLocal))

            # kLocal = kLocalNotRotated
            # coefficients of the DOFs and assignment of the stiffness matrix / force vector

            kCoeff = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeff[0+i::3]=coherentElemNodes*3+i
            # print('kCoeff: ', kCoeff)

            rows = np.zeros(kLocal.size,dtype=int)
            columns = np.zeros(kLocal.size,dtype=int)
            c=0
            for j in kCoeff:
                for i in kCoeff:
                    rows[c] = i
                    columns[c] = j
                    c+=1
        else:
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

            # #rotate stiffness matrix
            # kTemp = np.matmul(kLocalNotRotated, R)
            # kLocal = np.matmul(R.transpose(), kTemp)

            kLocal = kLocalNotRotated #TODO: rotate the matrixes!

            # coefficients of the DOFs and assignment of the stiffness matrix / force vector

            kCoeffTemp = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeffTemp[0+i::3]=coherentElemNodes*3+i

            kCoeff = np.zeros((3*nNodes-1),dtype=int)
            kCoeff[:-2] = kCoeffTemp[:-3]
            kCoeff[-2:] = kCoeffTemp[-2:]
            discartedDOFs[k] = kCoeffTemp[-3]
            discartedDOFs = np.array(discartedDOFs, dtype = int)

            rows = np.zeros(kLocal.size,dtype=int)
            columns = np.zeros(kLocal.size,dtype=int)
            c=0
            for j in kCoeff:
                for i in kCoeff:
                    rows[c] = i
                    columns[c] = j
                    c+=1
            # print('kCoeffTemp: ', kCoeffTemp)
            # print('kCoeff: ', kCoeff)
            # print('disc: ', discartedDOFs)
    
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
    
    if elementDegree=='MITC9':
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
    # dear future me,
    # I0m sorry if you lost time with this line of code, I was too lazy to change
    values[:,1]=uGlob[1::3,0]
    values[:,2]=uGlob[2::3,0]
    self.results = Result(outPos,vDisp, values[:,1], values[:,2],resultsScale=(resultsScaleVertDisp,resultsScaleIntForces))
    #compute MOMENTS
    nodesArray = self.mesh.nodesArray

    bendingMomentsSum = np.zeros((nodesArray.shape[0],3+1))
    shearForcesSum = np.zeros((nodesArray.shape[0],2+1))
    internalForcesPositions = nodesArray.to_numpy()[:,0:2]
    k=0

    ##################### debug
    if internalForcePosition == 'center':
        bendingMoments = np.zeros((len(elementsList),3))
        internalForcesPositions = np.zeros((len(elementsList),2))
        shearForces = np.zeros((len(elementsList),2))

    elif internalForcePosition == 'intPoints':
        nGaussPoints = 7
        bendingMoments = np.zeros((len(elementsList)*nGaussPoints,3))
        internalForcesPositions = np.zeros((len(elementsList)*nGaussPoints,2))
        shearForces = np.zeros((len(elementsList)*nGaussPoints,2))

    ######################
    singlePrint = True
    for k,element in enumerate(elementsList):
        elemNodes = element.connectivity
        # print('conn: ', elemNodes)

        coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
        nNodes=element.nNodes

        xi=element.coordinates[:,0]
        yi=element.coordinates[:,1]

        # internalForcesPositions[k,0]=xm
        # internalForcesPositions[k,1]=ym

        Df = self.plates[0].Df   
        Dc = self.plates[0].Dc
        elementShape = len(xi)
        kCoeff = np.zeros((3*nNodes),dtype=int)
        for i in range(0,3):
            kCoeff[0+i::3]=coherentElemNodes*3+i

        vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])

        if elementDegree =='L' and internalForcePosition == 'nodes':
            ri = np.array([-1, 1, 1, -1])
            si = np.array([-1,-1,1,1])
            for i in range(0, len(ri)):
                # print('ri, si: ', ri[i], si[i])
                N, Bf,Bc, detJ, Nval = getLinearVectorizedShapeFunctions(elementShape,np.array([ri[i]]), np.array([si[i]]), xi, yi)
                Bf = Bf[0,:,:]
                Bc = Bc[0,:,:]

                bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
                shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1  # !!!! WHY -1?????????????????
                bendingMomentsSum[coherentElemNodes[i],3] += 1
                shearForcesSum[coherentElemNodes[i],2] += 1
                # if singlePrint:
                #     print('r, s: ', ri[i], si[i])
                #     print('node: ', coherentElemNodes[i])
                #     print('vloc: ', vLoc)
                #     print('shear strain matrix: ', Bc[0,:])
                #     print('Dc: ', Dc)
                #     print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc)) )
                #     singlePrint = True

        if elementDegree =='L' and internalForcePosition == 'center':
            
            elemNodes = element.connectivity

            coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
            nNodes=element.nNodes

            xi=element.coordinates[:,0]
            yi=element.coordinates[:,1]

            xm = np.average(xi)
            ym = np.average(yi)
            internalForcesPositions[k,0]=xm
            internalForcesPositions[k,1]=ym

            Df = self.plates[0].Df 
            Dc = self.plates[0].Dc

            elementShape = len(xi)
            N, Bf,Bc, detJ, Nval = getLinearVectorizedShapeFunctions(elementShape,np.array([0]), np.array([0]), xi, yi)
            Bf = Bf[0,:,:]
            Bc = Bc[0,:,:]

            kCoeff = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeff[0+i::3]=coherentElemNodes*3+i

            vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
            bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
            shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1
            # if singlePrint:
            #     print('r, s: ', 0, 0)
            #     print('node: ', coherentElemNodes[i])
            #     print('vloc: ', vLoc)
            #     print('shear strain matrix: ', Bc[0,:])
            #     print('Dc: ', Dc)
            #     print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc)) )
            #     singlePrint = True

        if elementDegree =='Q' and internalForcePosition == 'center':
            
            elemNodes = element.connectivity
            coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
            nNodes=element.nNodes

            xi=element.coordinates[:,0]
            yi=element.coordinates[:,1]

            xm = np.average(xi)
            ym = np.average(yi)
            internalForcesPositions[k,0]=xm
            internalForcesPositions[k,1]=ym

            Df = self.plates[0].Df 
            Dc = self.plates[0].Dc

            elementShape = len(xi)
            N, Bf,Bc, detJ = getQuadraticShapeFunctions(elementShape,0, 0, xi, yi)
            Bf = Bf[:,:]
            Bc = Bc[:,:]
                
            kCoeff = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeff[0+i::3]=coherentElemNodes*3+i

            vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
            bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
            shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1  

        elif elementDegree =='MITC4' and internalForcePosition == 'nodes':
            ri = np.array([1, -1, -1, 1])
            si = np.array([1,1,-1,-1])
            for i in range(0, len(ri)):
                # print('ri, si: ', ri[i], si[i])
                N, Bf,Bc, detJ =getMITCShapefunctions(ri[i], si[i], xi, yi)
                # print('node: ', coherentElemNodes[i])
                # print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 )
                bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*1 
                bendingMomentsSum[coherentElemNodes[i],3] += 1
                shearForcesSum[coherentElemNodes[i],2] += 1
                # if singlePrint:

                #     print('r, s: ', ri[i], si[i])
                #     print('node: ', coherentElemNodes[i])
                #     print('vloc: ', vLoc)
                #     print('shear strain matrix: ', Bc[0,:])
                #     print('Dc: ', Dc)
                #     print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc)) )
                #     singlePrint = True
                # print('node: ', coherentElemNodes[i])
                # print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 )
        elif elementDegree =='MITC4' and internalForcePosition == 'center' :     
            elemNodes = element.connectivity

            coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
            nNodes=element.nNodes

            xi=element.coordinates[:,0]
            yi=element.coordinates[:,1]

            xm = np.average(xi)
            ym = np.average(yi)
            internalForcesPositions[k,0]=xm
            internalForcesPositions[k,1]=ym

            Df = self.plates[0].Df 
            Dc = self.plates[0].Dc

            elementShape = len(xi)
            N, Bf,Bc, detJ = getMITCShapefunctions(0, 0, xi, yi)

            kCoeff = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeff[0+i::3]=coherentElemNodes*3+i

            vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
            bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]
            shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]
            # if singlePrint:
            #     print('r, s: ', 0, 0)
            #     print('node: ', coherentElemNodes[i])
            #     print('vloc: ', vLoc)
            #     print('shear strain matrix: ', Bc[0,:])
            #     print('Dc: ', Dc)
            #     print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc)) )
            #     singlePrint = True

        if elementDegree =='Q' and internalForcePosition == 'center':
            elemNodes = element.connectivity
            coherentElemNodes = element.coherentConnectivity.to_numpy()[:,0]
            nNodes=element.nNodes

            xi=element.coordinates[:,0]
            yi=element.coordinates[:,1]

            xm = np.average(xi)
            ym = np.average(yi)
            internalForcesPositions[k,0]=xm
            internalForcesPositions[k,1]=ym

            Df = self.plates[0].Df 
            Dc = self.plates[0].Dc

            elementShape = len(xi)
            N, Bf,Bc, detJ = getQuadraticShapeFunctions(elementShape,0, 0, xi, yi)
            Bf = Bf[:,:]
            Bc = Bc[:,:]
                
            kCoeff = np.zeros((3*nNodes),dtype=int)
            for i in range(0,3):
                kCoeff[0+i::3]=coherentElemNodes*3+i

            vLoc = np.matmul(element.rotationMatrix, uGlob[kCoeff])
            bendingMoments[k,:] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*-1
            shearForces[k,:]=np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1  

        elif elementDegree =='MITC4' and internalForcePosition == 'nodes':
            ri = np.array([1, -1, -1, 1])
            si = np.array([1,1,-1,-1])
            for i in range(0, len(ri)):
                # print('ri, si: ', ri[i], si[i])
                N, Bf,Bc, detJ =getMITCShapefunctions(ri[i], si[i], xi, yi)
                # print('node: ', coherentElemNodes[i])
                # print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 )
                bendingMomentsSum[coherentElemNodes[i],0:3] += np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                shearForcesSum[coherentElemNodes[i],0:2] += np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*1 
                bendingMomentsSum[coherentElemNodes[i],3] += 1
                shearForcesSum[coherentElemNodes[i],2] += 1 

        elif elementDegree == 'L' and internalForcePosition == 'intPoints':
            # nGaussPoints = 10
            # gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',nGaussPoints)
            # ri = gaussPoints[:,0]
            # si = gaussPoints[:,1]
            ri =np.linspace(-1, 1, num=nGaussPoints)
            # print('ri: ', ri)
            si = np.zeros(nGaussPoints)
            allN, allBf,allBc, alldetJ, allNval =getLinearVectorizedShapeFunctions(elementShape,ri, si, xi, yi)
            for i in range(0, len(ri)):
                # print('ri, si: ', ri[i], si[i])

                N=allN[i,:,:]
                Bf = allBf[i,:,:]
                Bc = allBc[i,:,:]

                Nval = allNval[i,:]

                # print('node: ', coherentElemNodes[i])
                # print('shear: ',np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*-1 )
                
                bendingMoments[k*nGaussPoints+i,0:3] = np.matmul(Df,np.matmul(Bf, vLoc))[:,0]*1
                shearForces[k*nGaussPoints+i,0:2] = np.matmul(Dc, np.matmul(Bc, vLoc))[:,0]*1 

                xr = Nval@xi
                yr = Nval@yi
                internalForcesPositions[k*nGaussPoints+i,0]=xr
                internalForcesPositions[k*nGaussPoints+i,1]=yr

    if internalForcePosition == 'nodes':
        # n = bendingMomentsSum.shape[0]
        # bendingMoments = np.zeros((n,3))
        # shearForces = np.zeros((n,2)) 
        # for i in range(0, n):
        #     bendingMoments[i,:] = bendingMomentsSum[i,0:3]/bendingMomentsSum[i,3]
        #     shearForces[i,:] = shearForcesSum[i,0:2]/shearForcesSum[i,2]

        bendingMoments = bendingMomentsSum[:,0:3]/bendingMomentsSum[:,3, None]
        shearForces = shearForcesSum[:,0:2]/shearForcesSum[:,2, None]
    elif internalForcePosition =='intPoints':
        pass



    self.results.bendingMoments=bendingMoments*resultsScaleIntForces[0]
    self.results.internalForcesPositions=internalForcesPositions
    
    self.results.shearForces = shearForces*resultsScaleIntForces[1]

    return outPos, values

def rotMatrix(theta):
    '''
        creates rotation matrix of angle theta for a single node
    '''
    A = np.array([[1., 0, 0],
                  [0, np.cos(theta), -np.sin(theta)],
                  [0, np.sin(theta), np.cos(theta)]], dtype=float)

    return A

def GetLocalMatrix(xi, yi, Df,Dc, p,elementShape, elemType):

    ''' Input/Output descriptions
    ElemType: Quadrangluar or Triangular or Beam + Linear or Quadratic or MITC + Reduced or Normal Integration
        xi, yi: element's nodes coordinates
        Df, Dc: material of the plate
        p: load vector
        reducedIntegration: bool to manage the number of Gauss integration points 

        return:
        kLocal, fLocal
    '''
    temp = elemType.split('-')
    elementDegree = temp[0]

    elementIntegration = temp[1]
    if elementDegree in ['L', 'MITC4']:
        if elementShape==2:
        # 1 point
            gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('linear',1)
        # 2 points            
            gaussPoints, gaussWeights =  getGaussQuadrature('linear',2)

        if elementShape == 3:
        # 1 point
            gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('triangular',1)

        # 3 points            
            gaussPoints, gaussWeights =  getGaussQuadrature('triangular',3)

        elif elementShape == 4:

        # 1 point
            gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('rectangular',1)

        # 4 points
            gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',4)
    elif elementDegree in ['Q', 'MITC9']:
        if elementShape==2:
        # 1 point
            gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('linear',1)
        # 2 points            
            gaussPoints, gaussWeights =  getGaussQuadrature('linear',2)

        if elementShape == 3:
        # 1 point
            gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('triangular',1)

        # 3 points            
            gaussPoints, gaussWeights =  getGaussQuadrature('triangular',3)

        elif elementShape == 9:
 
        # 4 point
            gaussPointsRed, gaussWeightsRed =  getGaussQuadrature('rectangular',4)

        # 9 points
            gaussPoints, gaussWeights =  getGaussQuadrature('rectangular',9)

    if elementDegree == 'L': #perform vectorized integration
        fLocal = np.zeros(3*elementShape)
        ri = gaussPoints[:,0]
        si = gaussPoints[:,1]
        wi = gaussWeights
        N, Bf,Bc, detJ, Nval = getLinearVectorizedShapeFunctions(elementShape,ri, si, xi, yi)

        m11 = np.matmul(Bf.transpose((0,2,1)),Df)
        m12 = np.matmul(m11,Bf)

        m12 = np.moveaxis(m12,0,-1)
        kLocal = np.dot(m12,wi*detJ)

        kBDEbug = np.dot(m12,wi*detJ)
        # CHANGING INTEGRATION
        ri = gaussPointsRed[:,0]
        si = gaussPointsRed[:,1]
        wi = gaussWeightsRed[:]
        N, Bf,Bc, detJ, Nval = getLinearVectorizedShapeFunctions(elementShape,ri, si, xi, yi)

        m21 = np.matmul(Bc.transpose((0,2,1)),Dc)
        m22 = np.matmul(m21,Bc)
        m22 = np.moveaxis(m22,0,-1)
        kLocal += np.dot(m22,wi*detJ)
        ksDEbug = np.dot(m22,wi*detJ)

        ri = gaussPoints[:,0]
        si = gaussPoints[:,1]
        wi = gaussWeights[:]
        N, Bf,Bc, detJ, Nval = getLinearVectorizedShapeFunctions(elementShape,ri, si, xi, yi)
        mTemp = np.matmul(N.transpose((0,2,1)),p.magnitude)
        mTemp=np.expand_dims(mTemp,2)
        mTemp = np.moveaxis(mTemp,0,-1)
        fLocal = np.dot(mTemp,wi*detJ)

    elif elementDegree == 'MITC4':
        #full integration for both components
        fLocal = np.zeros(3*4)
        kLocal = np.zeros((3*4,3*4))

        kBDEbug = np.zeros((3*4,3*4))
        ksDEbug = np.zeros((3*4,3*4))
        for i in range(0,gaussPoints.shape[0]):
            # print('NEW POINT')
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

    elif elementDegree == 'Q':
        fLocal = np.zeros(3*9)
        kLocal = np.zeros((3*9,3*9))

        kBDEbug = np.zeros((3*9,3*9))
        ksDEbug = np.zeros((3*9,3*9))
        for i in range(0,gaussPoints.shape[0]):
            # print('NEW POINT')
            ri = gaussPoints[i,0]
            si = gaussPoints[i,1]
            wi = gaussWeights[i]
            N, Bf,Bc, detJ = getQuadraticShapeFunctions(elemType,ri, si, xi, yi)

            m11 = np.matmul(Bf.transpose(),Df)
            m12 = np.matmul(m11,Bf)
            # m21 = np.matmul(Bc.transpose(),Dc) 
            # m22 = np.matmul(m21,Bc)
            # ksDEbug += wi*m22*detJ
            kBDEbug += wi*m12*detJ

            # kLocal += wi*m12*detJ + wi*m22*detJ
            kLocal += wi*m12*detJ
            fLocal += wi*np.matmul(N.transpose(),p.magnitude)*detJ

        for i in range(0, gaussPointsRed.shape[0]):
            ri = gaussPointsRed[i,0]
            si = gaussPointsRed[i,1]
            wi = gaussWeightsRed[i]
            N, Bf,Bc, detJ = getQuadraticShapeFunctions(elemType,ri, si, xi, yi)


            m21 = np.matmul(Bc.transpose(),Dc) 
            m22 = np.matmul(m21,Bc)


            kLocal += wi*m22*detJ


    elif elementDegree == 'MITC9':
        #full integration for both components
        nDOFs = 8+2*9
        fLocal = np.zeros(nDOFs)
        kLocal = np.zeros((nDOFs,nDOFs))

        kBDEbug = np.zeros((nDOFs,nDOFs))
        ksDEbug = np.zeros((nDOFs,nDOFs))
        for i in range(0,gaussPoints.shape[0]):
            # print('NEW POINT')
            ri = gaussPoints[i,0]
            si = gaussPoints[i,1]
            wi = gaussWeights[i]
            N, Bf,Bc, detJ = getMITC9Shapefunctions(ri, si, xi, yi)

            m21 = np.matmul(Bc.transpose(),Dc) 
            m22 = np.matmul(m21,Bc)
            ksDEbug += wi*Bc.transpose() @ Dc @ Bc *detJ #TODO: è giusto il Jacobian?
            kBDEbug += wi*Bf.transpose() @ Df @ Bf * detJ
            temp1 = np.expand_dims(p.magnitude, 1)
            temp2 = wi* N.transpose() @ temp1 *detJ
            kLocal += wi*Bf.transpose() @ Df @ Bf *detJ + wi*Bc.transpose() @ Dc @ Bc*detJ
            fLocal += temp2[:,0]


    return kLocal, fLocal

def getLinearVectorizedShapeFunctions(elemType,ri, si, xi, yi):
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
    return N, Bf,Bc, detJ, Nval

def getQuadraticShapeFunctions(elemType,ri, si, xi, yi):
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

    if elemType==2:
        # Define shape functions
        N1 = lambda r, s: 0.5*(1-r)
        N2 = lambda r, s: 0.5*(1+r)

        # Form the shape function matrix
        Nfun= lambda r, s: [N1(r,s), N2(r,s)]
        Nval=np.array(Nfun(ri, si))
        Nval=np.moveaxis(Nval, -1, 0)

        N=np.zeros((nPoints,3, 3*elemType))
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
        NrsVal=np.array(NrsFun(ri,si))

        # matmul treat NrsVal as stack of matrixes residing in the LAST 2 indexes

        J=np.matmul(nodeCoordinates.transpose(), NrsVal)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)


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

    elif elemType==3:
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
    Nr[3] =  0.25*(1-GPs) + 0.5*(1-GPs)     -0.25*(1-GPs**2)
    Nr[4] = -GPr*(1+GPs)
    Nr[5] = -0.5*(1-GPs**2)
    Nr[6] = -GPr*(1-GPs)
    Nr[7] = 0.5*(1-GPs**2)

    Ns[0] =  0.25*(1+GPr) - 0.25*(1-GPr**2) + 0.5*GPs*(1+GPr)
    Ns[1] =  0.25*(1-GPr) - 0.25*(1-GPr**2) + 0.5*GPs*(1-GPr)
    Ns[2] = -0.25*(1-GPr) + 0.5*GPs*(1-GPr) + 0.25*(1-GPr**2)
    Ns[3] = -0.25*(1+GPr) + 0.25*(1-GPr**2) + 0.5*(1+GPr)
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
    Nr[3] =  0.25*(1-GPs) + 0.5*GPr*(1-GPs)     -0.25*(1-GPs**2) - 0.5*GPr*(1-GPs**2)
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
    Jac = np.array([[xr, xs],[yr, ys]])
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

def getGaussQuadrature(shape, nPoints):
    gaussQuadrature={'linear': {1:{'points': np.array([[0,0]]),
                                    'weights': np.array([2])},

                                    2:{'points': np.array([[-1*np.sqrt(3),0],
                                                            [1*np.sqrt(3),0]]),
                                    'weights': np.array([1,1])}},
        
                    'triangular': {1:{'points': np.array([[1/3, 1/3]]),
                                    'weights': np.array([1/2])},

                                    3:{'points': np.array([[1/6, 1/6],
                                                            [2/3,   1/6],
                                                            [1/6, 2/3]]),
                                    'weights': np.array([1/6, 1/6, 1/6])},

                                    4:{'points': np.array([[1/3, 1/3],
                                                            [1/5, 1/5],
                                                            [3/5, 1/5],
                                                            [1/5, 3/5]]),
                                    'weights': np.array([-27/96, 25/96, 25/96, 25/96])}},

                    'rectangular':{1:{'points': np.array([[0, 0]]),
                                    'weights': np.array([4])},

                                    4:{'points': np.array([[-1/np.sqrt(3), -1/np.sqrt(3)],
                                                            [1/np.sqrt(3),   -1/np.sqrt(3)],
                                                            [1/np.sqrt(3),1/np.sqrt(3)],
                                                            [-1/np.sqrt(3), 1/np.sqrt(3)]]),
                                    'weights': np.array([1, 1, 1, 1])},

                                    9:{'points': np.array([[-np.sqrt(3/5), -np.sqrt(3/5)],
                                                            [ 0      , -np.sqrt(3/5)],
                                                            [+np.sqrt(3/5), -np.sqrt(3/5)],
                                                            [-np.sqrt(3/5),    0    ],
                                                            [   0    ,    0    ] ,
                                                            [+np.sqrt(3/5),    0    ],
                                                            [-np.sqrt(3/5), +np.sqrt(3/5)],
                                                            [   0    , +np.sqrt(3/5)],
                                                            [+np.sqrt(3/5), +np.sqrt(3/5)]]),
                                    'weights': np.array([25, 40, 25, 40, 64, 40, 25, 40, 25])/81}}}
    return gaussQuadrature[shape][nPoints]['points'], gaussQuadrature[shape][nPoints]['weights']

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
        xi8 = xi[0:-1]
        yi8 = yi[0:-1]
        # Jac, JacInv, JacDet = getJac(Nr, Ns, xi8, yi8)

        # T= np.matmul(JacInv, np.array([Nr, Ns]))
        
        # Nr = T[0,:]
        # Ns = T[1,:]
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
    
    #integration

    nPoints=9
    gaussPoints, gaussWeights = getGaussQuadrature('rectangular', nPoints)
    for i in range(0,nPoints):
        wi = gaussWeights[i]
        N8, Nr8, Ns8 = shapeFun8(gaussPoints[i,0], gaussPoints[i,1])
        xi8 = xi[0:-1]
        yi8 = yi[0:-1]

        N9, Nr9, Ns9 = shapeFun9(gaussPoints[i,0], gaussPoints[i,1])

        M2[8,:] += wi*v1(gaussPoints[i,0], gaussPoints[i,1])
        M2[9, :] += wi*v2(gaussPoints[i,0], gaussPoints[i,1])

        M3[8,0:8] += wi*Nr8
        M3[8, 8:17] += -wi*np.dot(Nr9, yi)*N9
        M3[8, 17:] += wi*np.dot(Nr9, xi)*N9

        M3[9,0:8] += wi*Ns8
        M3[9, 8:17] += -wi*np.dot(Ns9, yi)*N9
        M3[9, 17:] += wi*np.dot(Ns9, xi)*N9

    # print('M3: ', pd.DataFrame(M3))
    # print('M2: ', pd.DataFrame(M2))
    # print('invM2: ', pd.DataFrame(np.linalg.inv(M2)))

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


    alpha, beta = naturalToCartesian(xi,yi)

    Dx=np.sqrt((np.dot(Ns,xi))**2+(np.dot(Ns,yi))**2)/(JacDet)
    Dy=np.sqrt((np.dot(Nr,xi))**2+(np.dot(Nr,yi))**2)/(JacDet)
    # print('Alpha: ', alpha)
    # print('beta: ', beta)
    # print('Dx: ', Dx)
    # print('Dy: ', Dy)
    transformMat = np.array([[np.sin(beta)*Dx, -np.sin(alpha)*Dy], [-np.cos(beta)*Dx, np.cos(alpha)*Dy]])

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

if __name__ == "__main__":

    xi = np.array([5. , 0. , 0. , 5. , 2.5, 0. , 2.5, 5. , 2.5])
    yi = np.array([5. , 5. , 0. , 0. , 5. , 2.5, 0. , 2.5, 2.5])

    ri = -0.7745966692414834
    si = -0.7745966692414834

    N, Nr, Ns = shapeFun9(ri, si)
    xr = np.matmul(Nr,xi)
    xs = np.matmul(Ns,xi)
    yr = np.matmul(Nr,yi)
    ys=np.matmul(Ns,yi)
    Jac = np.array([[xr, xs],[yr, ys]])
    JacInv = np.linalg.inv(Jac)
    JacDet = np.linalg.det(Jac)

    N, Bf,Bc, detJ, NrsVal, J = getQuadraticShapeFunctions(9,ri, si, xi, yi)





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

#unit testing
if __name__ == "__main__":
    pass


    # %matplotlib
    # import matplotlib.pyplot as plt

    # # import sample plate and show it
    # import sampleGeometry
    # from displayModel import *
    # sampleModel=sampleGeometry.importModel()
    # plotInputGeometry(sampleModel)

    # # create mesh
    # from generateMesh import *
    # meshInput1=MeshInput(showGmshMesh=False, elementType='QUAD', nEdgeNodes=33)
    # #meshInput1=MeshInput(showGmshMesh=True, elementType='QUAD', meshSize=5e-2)
    # generateMesh(sampleModel, meshInput1)

    # # compute
    # from solveModel import *
    # solveModel(sampleModel)
    # # display results
    # plotResults(sampleModel)

    # coords = sampleModel.results.outPos
    # z = sampleModel.results.wVert
    # x = coords[:,0]
    # y=coords[:,1]
    # fig = plt.figure()

    # ax = fig.gca(projection='3d')
    # t=10
    # ax.plot_trisurf(x*10,y*10,z,cmap=cm.jet)


    # def axisEqual3D(ax):
    #     extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    #     sz = extents[:,1] - extents[:,0]
    #     centers = np.mean(extents, axis=1)
    #     # maxsize = max(abs(sz))
    #     r = sz[2]/2*2.6
    #     for ctr, dim in zip([centers[2]], 'z'):
    #         getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    # # extents = np.array([getattr(ax, 'get_zlim'))
    # # sz = extents[:,1] - extents[:,0]
    # # centers = np.mean(extents, axis=1)
    # # getattr(ax, 'set_zlim')()
    # axisEqual3D(ax)
    # # Hide grid lines
    # ax.grid(False)

    # # Hide axes ticks
    # ax.set_xticks([])
    # ax.set_yticks([])
    # ax.set_zticks([])

    # plt.show()

    # # from mpl_toolkits.mplot3d import Axes3D # 3D plot
    # # x= sampleModel.results.bendingMomentsPositions[:,0]
    # # y = sampleModel.results.bendingMomentsPositions[:,1]
    # # z= sampleModel.results.bendingMoments[:,0]
    # # #z= sampleModel.results.shearForce[:,0]

    # # # print('x, y: ',sampleModel.results.bendingMomentsPositions)
    # # # print(z[0::10])

    # # fig = plt.figure(

    # # ax = fig.gca(projection='3d')
    # # ax.plot_trisurf(x,y,z,cmap=cm.jet)
    # # plt.show()




# # shapefunctions, devectorized:

#     nPoints = ri.size
#     nodeCoordinates = np.zeros((elemType, 2))
#     nodeCoordinates[:,0]=xi
#     nodeCoordinates[:,1]=yi

#     if elemType==2:
#         # Define shape functions
#         N1 = lambda r, s: 0.5*(1-r)
#         N2 = lambda r, s: 0.5*(1+r)

#         # Form the shape function matrix
#         Nfun= lambda r, s: [N1(r,s), N2(r,s)]
#         Nval=np.array(Nfun(ri, si))
#         Nval=np.moveaxis(Nval, -1, 0)

#         N=np.zeros((nPoints,3, 3*elemType))
#         N[0, 0::3]=Nval
#         N[1, 1::3]=Nval
#         N[2, 2::3]=Nval

#         # Define shape function derivatives, derive deformation matrix
#         N1r = lambda r, s: -0.25*(1-s)
#         N2r = lambda r, s: 0.25*(1-s)
#         N3r = lambda r, s: 0.25*(1+s)
#         N4r = lambda r, s: -0.25*(1+s)

#         N1s = lambda r, s: -0.25*(1-r)
#         N2s = lambda r, s: -0.25*(1+r)
#         N3s = lambda r, s: 0.25*(1+r)
#         N4s = lambda r, s: 0.25*(1-r)

#         NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)],[N4r(r, s), N4s(r, s)]])
#         NrsVal=np.array(NrsFun(ri,si))

#         # matmul treat NrsVal as stack of matrixes residing in the LAST 2 indexes

#         J=np.matmul(nodeCoordinates.transpose(), NrsVal)
#         detJ = np.linalg.det(J)
#         invJ = np.linalg.inv(J)

#         NrsVal = np.matmul(NrsVal,invJ)
#         Bf=np.zeros((3,3*elemType))
#         Bf[0,1::3]=NrsVal[:,0]
#         Bf[1,2::3]=NrsVal[:,1]
#         Bf[2,1::3]=NrsVal[:,1]
#         Bf[2,2::3]=NrsVal[:,0]

#         Bc=np.zeros((nPoints,2,3*elemType))
#         Bc[0,0::3]=NrsVal[:,0]
#         Bc[0,1::3]=Nval
#         Bc[1,0::3]=NrsVal[:,1]
#         Bc[1,2::3]=Nval

#     elif elemType==3:
#         # Define shape functions
#         N1 = lambda r, s: r
#         N2 = lambda r, s: s
#         N3 = lambda r, s: 1-r-s

#         # Form the shape function matrix
#         Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s)]
#         Nval=np.array(Nfun(ri, si))

#         N=np.zeros((3, 3*elemType))
#         N[0, 0::3]=Nval
#         N[1, 1::3]=Nval
#         N[2, 2::3]=Nval

#         # Define shape function derivatives, derive deformation matrix
#         N1r = lambda r, s: 1*np.ones(len(r))
#         N2r = lambda r, s: 0*np.ones(len(r))
#         N3r = lambda r, s: -1*np.ones(len(r))

#         N1s = lambda r, s: 0*np.ones(len(r))
#         N2s = lambda r, s: 1*np.ones(len(r))
#         N3s = lambda r, s: -1*np.ones(len(r))

#         NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)]])
#         NrsVal=np.array(NrsFun(ri,si))

#         # Jacobian matrix
#         J=np.matmul(nodeCoordinates.transpose(), NrsVal)
#         detJ = np.linalg.det(J)
#         invJ = np.linalg.inv(J)

#         NrsVal = np.matmul(NrsVal,invJ)
#         Bf=np.zeros((3,3*elemType))
        
#         Bf[0,1::3]=NrsVal[:,0]
#         Bf[1,2::3]=NrsVal[:,1]
#         Bf[2,1::3]=NrsVal[:,1]
#         Bf[2,2::3]=NrsVal[:,0]

#         Bc=np.zeros((2,3*elemType))
#         Bc[0,0::3]=NrsVal[:,0]
#         Bc[0,1::3]=Nval[:]
#         Bc[1,0::3]=NrsVal[:,1]
#         Bc[1,2::3]=Nval[:]

#     elif elemType==4:
#         # Define shape functions
#         N1 = lambda r, s: 0.25*(1-r)*(1-s)
#         N2 = lambda r, s: 0.25*(1+r)*(1-s)
#         N3 = lambda r, s: 0.25*(1+r)*(1+s)
#         N4 = lambda r, s: 0.25*(1-r)*(1+s)

#         # Form the shape function matrix
#         Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s), N4(r,s)]
#         Nval=np.array(Nfun(ri, si))


#         N=np.zeros((nPoints,3, 3*elemType))
#         N[:,0, 0::3]=Nval
#         N[:,1, 1::3]=Nval
#         N[:,2, 2::3]=Nval

#         # Define shape function derivatives, derive deformation matrix
#         N1r = lambda r, s: -0.25*(1-s)
#         N2r = lambda r, s: 0.25*(1-s)
#         N3r = lambda r, s: 0.25*(1+s)
#         N4r = lambda r, s: -0.25*(1+s)

#         N1s = lambda r, s: -0.25*(1-r)
#         N2s = lambda r, s: -0.25*(1+r)
#         N3s = lambda r, s: 0.25*(1+r)
#         N4s = lambda r, s: 0.25*(1-r)

#         NrsFun = lambda r,s: np.array([[N1r(r, s), N1s(r, s)], [N2r(r, s), N2s(r, s)], [N3r(r, s), N3s(r, s)],[N4r(r, s), N4s(r, s)]])
#         NrsVal=np.array(NrsFun(ri,si))

#         # matmul treat NrsVal as stack of matrixes residing in the LAST 2 indexes

#         J=np.matmul(nodeCoordinates.transpose(), NrsVal)
#         detJ = np.linalg.det(J)
#         invJ = np.linalg.inv(J)

#         NrsVal = np.matmul(NrsVal,invJ)
#         Bf=np.zeros((3,3*elemType))
#         Bf[0,1::3]=NrsVal[:,0]
#         Bf[1,2::3]=NrsVal[:,1]
#         Bf[2,1::3]=NrsVal[:,1]
#         Bf[2,2::3]=NrsVal[:,0]

#         Bc=np.zeros((2,3*elemType))
#         Bc[0,0::3]=NrsVal[:,0]
#         Bc[0,1::3]=Nval
#         Bc[1,0::3]=NrsVal[:,1]
#         Bc[1,2::3]=Nval