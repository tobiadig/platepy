'''
-----------------------------------------------------------
Generates and stores the mesh of a plate Model.
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 30.05.2021
'''
#%% Basic modules
import numpy as np
import pandas as pd
import copy
import gmsh # To create CAD model and mesh

def generateMesh(self,showGmshMesh=False,showGmshGeometryBeforeMeshing = False, elementDefinition='MITC-4-N', \
    meshSize=6e-1, nEdgeNodes=0, order='linear', meshDistortion = False, distVal = 100,\
        deactivateRotation=False):
    ''' Generates mesh and stores it in the plateModel object according to the selected options. Gmsh model has to be already 
        initialized and structural elements have to be added to the model.
        ~~~~~~~~~~~~~~~~~~~
        Input:
        ~~~~~~~~~~~~~~~~~~~
        * **self**: plateModel object. 
        * **showGmshMesh = False**: If True, the model is shown through the built-in fltk terminal (after the mesh generation). 
        * **showGmshMesh = True**: If True, the model is shown through the built-in fltk terminal (before the mesh generation). 
        * **elementDefinition = None**: String defining the desired FE-element in the following form: "type-nNodes-integration". 

            - type: DB for displacement-based elements or MITC. 
            - nNodes: number of nodes ( currently 3, 4 or 9). 
            - integration: Desired Gauss-quadrature for the calculation of the stiffness matrixes. R for reduced or N for normal. 

        * **meshSize = 8e-1**: target mesh size around the point entities. If nEdgeNodes > 0, meshSize is ignored. 
        * **nEdgeNodes = 0**: Prescribes the number of nodes on each edge. Plate must be rectangular.
        * **order = "linear"**: Prescribes the order of the elements. "linear for first order elements (default) and "quadratic" for second order elements. 
        * **meshDistortion = False**: Boolean, if True a mesh distortion function is applied. 
        * **distVal = 100**: Severeness of the mesh distortion function. 
        ~~~~~~~~~~~~~~~~~~~
        Return:
        ~~~~~~~~~~~~~~~~~~~
    '''

    elementType, elementShape, elementIntegration = _getElementDefinition(elementDefinition)
    gmsh.model.mesh.clear()

    # Set recombination rules
    if elementShape == 4 or elementShape == 9 :
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3) #0: simple, 1: blossom (default), 2: simple full-quad, 3: blossom full-quad
        for i in range(0, len(self.plates)):
            gmsh.model.geo.mesh.setRecombine(2, self.plates[i].tag)
    elif elementShape != 3:
        raise Exception
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    gmsh.model.geo.synchronize()

    # Manual assignment of edge nodes
    try:
        pointTags=gmsh.model.getEntities(0)
        if nEdgeNodes>0:
            gmsh.model.geo.mesh.setTransfiniteCurve(1, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(2, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(3, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(4, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
            gmsh.model.geo.synchronize()
        # automatically distort mesh with Gmsh, the results are pretty bad
        # elif nEdgeNodes>0 and meshDistortion:
            # gmsh.model.geo.mesh.setTransfiniteCurve(1, nEdgeNodes, "Progression", progVal)
        #     gmsh.model.geo.mesh.setTransfiniteCurve(2, nEdgeNodes, "Progression", progVal)
        #     gmsh.model.geo.mesh.setTransfiniteCurve(3, nEdgeNodes, "Progression", progVal)
        #     gmsh.model.geo.mesh.setTransfiniteCurve(4, nEdgeNodes, "Progression", progVal)
        #     gmsh.model.geo.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
        #     gmsh.model.geo.synchronize()
        else:
            gmsh.model.mesh.setSize(pointTags,meshSize)
    except:
        print('manual assignment of edge nodes failed')
        raise

    if showGmshGeometryBeforeMeshing:
        gmsh.fltk.run()

    # Mesh generation
    try:
        gmsh.model.mesh.generate()
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.geo.synchronize()
    except:
        print('gmsh mesh generation failed')
        raise

    # select the right order of the elements
    if order == 'quadratic':
        gmsh.model.mesh.setOrder(2)
        gmsh.model.geo.synchronize()
    elif order !='linear':
        raise Exception('order not recognised')
    
    # eventually open fltk UI
    if showGmshMesh:
        gmsh.fltk.run()

    # Generates nodes array
    nodeTagsModel , nodeCoords, _ = gmsh.model.mesh.getNodes()

    nodesArray = np.array(nodeCoords).reshape((-1,3))
    nodesArrayPd = pd.DataFrame(nodesArray, index = nodeTagsModel) # Why a dataframe? Nodes are not always continuous since sometimes are removed with removeDuplicateNodes. Using a pandas dataframe allows to uniquely call a node with his tag
    nodesRotationsPd = pd.DataFrame(np.zeros(nodeTagsModel.size), index =nodeTagsModel)
    gmshToCoherentNodesNumeration = pd.DataFrame(range(0,len(nodeTagsModel)), index = nodeTagsModel)
    # gmshNumeration is the one with node Tags, CoherentNumeration is from 0 to nNodes-1

    # distort the mesh if required
    if meshDistortion:
        nodesArrayPd = _distortMesh(nodesArrayPd, distVal)

    gmshModel = gmsh.model
    platesList = self.plates
    elementsList, getElementByTagDictionary = _getElementsList(gmshModel,platesList, elementType, elementShape,elementIntegration,gmshToCoherentNodesNumeration,nodesArrayPd)
    plateElementsList = copy.deepcopy(elementsList)  # useful if there is the need to distinguish plate elements from the downstand bem elements

    #assemble 2D elements for line load
    loadsList = self.loads
    loadsList = _getLineLoadForces(gmshModel, loadsList,gmshToCoherentNodesNumeration,nodesArrayPd)

    #generate BCs and nodes directions by iterating over wall segments
    wallList = self.walls
    columnsList = self.columns
    BCs, nodesRotationsPd,elasticallySupportedNodes = _getBCsArray(gmshModel,wallList,columnsList,nodesRotationsPd,deactivateRotation)
    elasticallySupportedNodes = gmshToCoherentNodesNumeration.loc[elasticallySupportedNodes].to_numpy()

    nextNode = int(np.max(nodeTagsModel)+1)
    downStandBeamsList = self.downStandBeams
    myPlate = self.plates[0]
    downStandBeamsList, gmshToCoherentNodesNumeration,nodesArray,nodesArrayPd,nodesRotationsPd =\
    _createNewNodesForDownStandBeams(gmshModel,nodesArray,nodesArrayPd,nodesRotationsPd,downStandBeamsList,gmshToCoherentNodesNumeration,nextNode)
    downStandBeamsList, AmatList, elementsList,nodesRotationsPd = \
    _getDownStandBeamsElements(gmshModel,nodesRotationsPd,nodesArray,elementsList,\
    downStandBeamsList,elementType,myPlate,gmshToCoherentNodesNumeration,nodesArrayPd)

    # Store everything into the mesh object
    self.mesh = _Mesh(nodesArrayPd,nodesRotationsPd, elementsList, BCs, AmatList, plateElementsList, getElementByTagDictionary,elasticallySupportedNodes)

def setMesh(self, nodesArray, elements, BCs,elementDefinition = None, load = None):
    ''' Allows to manually define node positions, elements connectivity, boundary conditions and loads. \n
        Input: \n
            * self: plateModel object.\n
            * nodesArray: nNodes x 3 numpy array. Columns are the x-y-z coordinates of each node.\n
            * elements: (nElements x nElementNodes) connectivity matrix. Each row is an element, the columns contain the node tags which build the element. \n
            * BCs: n x 4 numpy array. n is the number of restrained nodes, the first column is the node tag, the other column represent the condition of the three DOFs (1 is blocked, 0 is free). \n
            * load = None: n x 4 numpy array. n is the number of loaded nodes, the first column is the node tag, the other column represent the magnitude of the load for the relative DOF. \n
        Return: \n
            * - \n
    '''
    elementType, elementShape, elementIntegration = _getElementDefinition(elementDefinition)
    elementsList = []
    nNodes = nodesArray.shape[0]
    platesList = self.plates
    k=1
    for element in elements:
        newElement = _Element()
        newElement.tag = k
        newElement.nNodes  = len(element)
        newElement.connectivity  = element
        newElement.shape=elementShape
        newElement.Db = platesList[0].Df 
        newElement.Ds = platesList[0].Dc 
        newElement.type = elementType
        newElement.integration = elementIntegration
        newElement.coherentConnectivity = pd.DataFrame(element-1)
        newElement.coordinates = np.zeros((len(element), 3))
        newElement.coordinates[:,0:2] = nodesArray[element-1, :]
        newElement.whichPlate  = 0
        elementsList.append(newElement)
        k+=1

    nodesRotationsPd = pd.DataFrame(np.zeros((nNodes, 1)), index=range(1, nNodes+1))
    nodesArrayPd =pd.DataFrame(nodesArray, index=range(1, nNodes+1) )
    self.mesh = _Mesh(nodesArrayPd,nodesRotationsPd, elementsList, BCs, [], elementsList, [],[])
    self.mesh.load = load

def _getElementDefinition(elementDefinition):
    '''
        Extracts information from the elementDefinition String. \n
        Input: \n
        * elementDefinition : String defining the desired FE-element in form: "type-nNodes-integration". \n
            \t * type: DB for displacement-based elements or MITC. \n
            \t * nNodes: number of nodes ( currently 3, 4 or 9). \n
            \t * integration: Desired Gauss-quadrature for the calculation of the stiffness matrixes. R for reduced or N for normal. \n

        Return: \n
        * elementType \n
        * elementShape \n
        * elementIntegration
    '''
    temp = elementDefinition.split('-')
    elementType = temp[0]
    elementShape = int(temp[1])
    elementIntegration = temp[2]

    return elementType, elementShape, elementIntegration

def _distortMesh(nodesArray, alpha):
    '''
    Displace the nodes in nodesArray according to a distortion function. \n
    Input: \n
        * nodesArray: Pandas dataframe of shape n x 3 with the x-y-z coordinates of each node  \n
        * alpha: Variable in the distortion algorithm, controlling the severeness of the distortion.\n
    Return: \n
        *   newNodesArray: Pandas dataframe with the distorted nodes
    '''
    myIndex = nodesArray.index.to_numpy()
    nodesArrayNumpy = nodesArray.to_numpy()
    v1=np.ones(nodesArrayNumpy.shape[0])
    x0 = nodesArrayNumpy[:,0]
    y0 = nodesArrayNumpy[:,1]
    a=np.max(x0)
    newNodes = np.zeros(nodesArray.shape)
    newNodes[:,0] = x0+(v1-np.abs(2*x0/a-1))*(2*y0/a-v1)*alpha
    newNodes[:,1] = y0+2*(v1-np.abs(2*y0/a-1))*(2*x0/a-v1)*alpha
    newNodesArray = pd.DataFrame(newNodes, index = myIndex)
    return newNodesArray

def _getElementsList(gmshModel,platesList, elementType, elementShape,elementIntegration,gmshToCoherentNodesNumeration,nodesArrayPd):
    ''' Creates a list of all elements objects of the model.\n
        Input: \n
            * gmshModel: model object of the gmsh library. \n
            * platesList: List of plate objects. \n
            * elementType, elementShape, elementIntegration: Components of the elementDefinition. \n
            * gmshToCoherentNodesNumeration: Pandas dataframe where the indexes are the node tags assigned by gmsh, the values are an array from 0 to nNodes-1. \n
            * nodesArrayPd: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates. \n
        Return: \n
            * elementsList: list of all element objects of the model. \n
            * getElementByTagDictionary: dictionary where the keys are the tags assigned by gmsh, the values are the respective element objects.
    '''
    elementsList = []
    getElementByTagDictionary ={}

    # Assign all the information to each element object and save it in the elementsList
    for i in range(0, len(platesList)):
        _, elemTags, _ = gmshModel.mesh.getElements(2,platesList[i].tag)
        for elemTag in elemTags[0]:
            _ , nodeTags = gmshModel.mesh.getElement(elemTag)
            newElement = _Element()
            newElement.tag = elemTag
            newElement.whichPlate = i
            newElement.Db = platesList[i].Df 
            newElement.Ds = platesList[i].Dc 
            newElement.type = elementType
            newElement.shape = elementShape
            newElement.integration = elementIntegration
            newElement.nNodes  = len(nodeTags)
            newElement.connectivity  = nodeTags
            newElement.coherentConnectivity = gmshToCoherentNodesNumeration.loc[nodeTags]
            newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()

            elementsList.append(newElement)
            getElementByTagDictionary[elemTag] = newElement
    return elementsList, getElementByTagDictionary

def _getLineLoadForces(gmshModel, loadsList, gmshToCoherentNodesNumeration, nodesArrayPd):
    ''' Creates a list of all elements objects of the model.\n
        Input: \n
            * gmshModel: model object of the gmsh library. \n
            * loadsList: List of load objects. \n
            * gmshToCoherentNodesNumeration: Pandas dataframe where the indexes are the node tags assigned by gmsh, the values are an array from 0 to nNodes-1. \n
            * nodesArrayPd: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates. \n
        Return: \n
            * loadsList: list of all load objects of the model.
    '''
    for p in loadsList:
        if p.case == 'line':
            tags = gmshModel.getEntitiesForPhysicalGroup(p.physicalGroup[0],p.physicalGroup[1])
            _, elementTags, nodeTags = gmshModel.mesh.getElements(1,tags[0])   
            elements1DList = []
            for elemTag in elementTags[0]:
                elementType, nodeTags = gmshModel.mesh.getElement(elemTag)
                newElement = _Element()
                newElement.tag = elemTag
                newElement.nNodes  = len(nodeTags)
                newElement.connectivity  = nodeTags
                newElement.coherentConnectivity = gmshToCoherentNodesNumeration.loc[nodeTags]
                newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()
                newElement.whichPlate  = 1  
                elements1DList.append(newElement)
            p.elements1DList = elements1DList
        elif p.case == 'point':
            tags = gmshModel.getEntitiesForPhysicalGroup(p.physicalGroup[0],p.physicalGroup[1])
            _, elementTags, nodeTags = gmshModel.mesh.getElements(0,tags[0])   
            p.pointLoadNode = nodeTags[0]

    return loadsList

def _getBCsArray(gmshModel,wallList,columnsList,nodesRotationsPd,deactivateRotation):
    ''' Creates array where the rows are the constrained nodes, the first columns is the node tag and the 
        other 3 columns define if the correspondend DOF is restrained or free (see supportCondition). 
        Additionally also the nodes rotations are stored.\n
        Input: \n
            * gmshModel: model object of the gmsh library. \n
            * wallList: List of wall objects. \n
            * columnsList: List of column objects. \n
            * nodesRotationsPd: pandas dataframe where the indexes are the node tags, the values are 0. \n 

        Return: \n
            * BCs: n x 4 numpy array. n is the number of restrained nodes, the first column is the node tag, the other column represent the condition of the three DOFs (1 is blocked, 0 is free). \n
            * nodesRotationsPd: pandas dataframe where the indexes are the node tags, the values are the node rotation in radians. 
            * elasticallySupportedNodes: 
    '''
    BCsDic = {}
    elasticallySupportedNodes = []
    for wall in wallList:
        dim = wall.physicalGroup[0]
        nodeTags, _ = gmshModel.mesh.getNodesForPhysicalGroup(dim,wall.physicalGroup[1])
        enitiesTags = gmshModel.getEntitiesForPhysicalGroup(dim,wall.physicalGroup[1])
        for wallLine in enitiesTags:
            nodeTags, coord, _ = gmshModel.mesh.getNodes(dim, wallLine, includeBoundary=True)
            coord = coord.reshape(-1,3)
            # start and end nodes are blocked (why?)
            p1Tag = nodeTags[-2]
            BCsDic[p1Tag] = np.array([1,1,1])
            p2Tag = nodeTags[-1]
            BCsDic[p2Tag] = np.array([1,1,1])
            lineDirection = coord[-1,:] -coord[-2,:]
            if deactivateRotation:
                rotationKiller = 0
            else:
                rotationKiller = 1
            if lineDirection[0]==0 and lineDirection[1]>0:
                lineRot = np.pi/2*rotationKiller
            elif lineDirection[0]==0 and lineDirection[1]<0:
                lineRot = np.pi/2*3*rotationKiller
            else:
                lineRot = np.arctan(lineDirection[1]/lineDirection[0])*rotationKiller
            nodesRotationsPd.loc[nodeTags[:-2]] = lineRot
            for node in nodeTags[:-2]:
                BCsDic[node] = wall.support
                if wall.isElasticallySupported:
                    elasticallySupportedNodes.append(node)
    for col in columnsList:
        dim = col.physicalGroup[0]
        nodeTags, _ = gmshModel.mesh.getNodesForPhysicalGroup(dim,col.physicalGroup[1])
        for node in nodeTags:
            BCsDic[node] = col.support
    BCs = np.zeros((len(BCsDic), 4))
    for count, a in enumerate(BCsDic):
        BCs[count, 0] = a
        BCs[count,1:] = BCsDic[a]

    elasticallySupportedNodes = np.array(elasticallySupportedNodes)
    return BCs, nodesRotationsPd,elasticallySupportedNodes

def _createNewNodesForDownStandBeams(gmshModel,nodesArray,nodesArrayPd,nodesRotationsPd,downStandBeamsList,gmshToCoherentNodesNumeration,nextNode):
    ''' Iterates of the downstand beams and creates the new nodes. The array of nodes and rotations are expanded and returned accordingly. \n
        Input: \n
            * gmshModel: model object of the gmsh library.\n
            * nodesArray: nNodes x 3 numpy array. Columns are the x-y-z coordinates of each node.\n
            * nodesArrayPd: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates.\n
            * nodesRotationsPd: pandas dataframe where the indexes are the node tags, the values are the node rotation in radians\n
            * downStandBeamsList: List of downStandBeam objects. \n
            * gmshToCoherentNodesNumeration: Pandas dataframe where the indexes are the node tags assigned by gmsh, the values are an array from 0 to nNodes-1.\n
            * nextNode: nNodes + 1. This will be the tag of the first beam node.\n
        Return: \n
            * downStandBeamsList: Updated input. \n
            * gmshToCoherentNodesNumeration: Updated input. \n
            * nodesArray: Updated input. \n
            * nodesArrayPd: Updated input. \n
            * nodesRotationsPd: Updated input.
    '''
    nextCoherentNode = int(nodesArray.shape[0])
    for uz in downStandBeamsList:
        dim = uz.physicalGroup[0]
        nodesUZ, coord = gmshModel.mesh.getNodesForPhysicalGroup(dim,uz.physicalGroup[1])
        nNodesUZ = len(nodesUZ)
        newNodesUZ = np.array(range(nextNode, nextNode+nNodesUZ), dtype=int)
        uz.newNodesUZ = newNodesUZ
        coherentNodesUZ = np.array(range(nextCoherentNode, nextCoherentNode+nNodesUZ), dtype=int)
        uz.coherentNodesUZ = coherentNodesUZ
        tempDF = pd.DataFrame(coherentNodesUZ, index=newNodesUZ)
        gmshToCoherentNodesNumeration = gmshToCoherentNodesNumeration.append(tempDF)
        gmshToCoherentNodesNumeration.index = gmshToCoherentNodesNumeration.index.astype(int)
        nextCoherentNode += nNodesUZ
        nextNode = nextNode+nNodesUZ
        tempNodesArray = nodesArrayPd.loc[nodesUZ].to_numpy()
        nodesArray = np.append(nodesArray, tempNodesArray, axis=0)
        nodesArrayPd = nodesArrayPd.append(pd.DataFrame(tempNodesArray, index=newNodesUZ))
        nodesArrayPd.index = nodesArrayPd.index.astype(int)
        uzNodesToNodesNumeration = pd.DataFrame(newNodesUZ, index = nodesUZ)
        uzNodesToNodesNumeration.index = uzNodesToNodesNumeration.index.astype(int)
        uz.uzNodesToNodesNumeration = uzNodesToNodesNumeration
        coherentNodesPlate = gmshToCoherentNodesNumeration.loc[nodesUZ].to_numpy()[:,0]
        uz.coherentNodesPlate = coherentNodesPlate
        nodesRotationsPd = nodesRotationsPd.append(pd.DataFrame(np.zeros(nNodesUZ), index =newNodesUZ))
        nodesRotationsPd.index = nodesRotationsPd.index.astype(int)
    return downStandBeamsList, gmshToCoherentNodesNumeration,nodesArray,nodesArrayPd,nodesRotationsPd

def _getDownStandBeamsElements(gmshModel,nodesRotationsPd,nodesArray,elementsList,\
    downStandBeamsList,elementType,myPlate,gmshToCoherentNodesNumeration,nodesArrayPd):
    ''' Add the downstand beam elements to the model elements list. 
    Additionally, creates the multi-constraint matrix which connects plate and beam DOFs. \n
        Input: \n
            * gmshModel: model object of the gmsh library.\n
            * nodesRotationsPd: pandas dataframe where the indexes are the node tags, the values are the node rotation in radians\n
            * nodesArray: nNodes x 3 numpy array. Columns are the x-y-z coordinates of each node.\n
            * elementsList: List of containing the elements objects of the plate. \n 
            * downStandBeamsList: List of downStandBeam objects. \n
            * elementType: String, can be "DB" or "MITC". \n
            * myPlate: plate object the downstand beam is connected to. \n
            * gmshToCoherentNodesNumeration: Pandas dataframe where the indexes are the node tags assigned by gmsh, the values are an array from 0 to nNodes-1.\n
            * nodesArrayPd: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates.\n
        Return: \n
            * downStandBeamsList: Updated input. \n
            * AmatList: List of MultiFreedom contraint matrix of shape nConstraint x nDOFs. nConstraint is 3 x number of downstandBeam nodes. \n
            * elementsList: Updated input. \n
            * nodesRotationsPd: Updated input.
    '''
    AmatList = []
    uzElementsList = []
    for uz in downStandBeamsList:
        newNodesUZ = uz.newNodesUZ
        uzNodesToNodesNumeration = uz.uzNodesToNodesNumeration 
        coherentNodesPlate = uz.coherentNodesPlate
        coherentNodesUZ = uz.coherentNodesUZ
        dim = uz.physicalGroup[0]

        enitiesTags = gmshModel.getEntitiesForPhysicalGroup(dim,uz.physicalGroup[1])
        for uzLine in enitiesTags:
            # nodeTags, coord, _ = gmsh.model.mesh.getNodes(dim, uzLine, includeBoundary=True)
            # coord = coord.reshape(-1,3)
            elementTypes, elementTags, nodeTags =gmshModel.mesh.getElements(dim,uzLine)
            for elemTag in elementTags[0]:
                _, nodeTags = gmshModel.mesh.getElement(elemTag)
                newElement = _Element()
                newElement.tag = elemTag
                newElement.correspondingPlateElements = nodeTags
                newElement.nNodes  = len(nodeTags)
                realNodeTags = uzNodesToNodesNumeration.loc[nodeTags].to_numpy()[:,0]
                newElement.connectivity  = realNodeTags
                newElement.coherentConnectivity = gmshToCoherentNodesNumeration.loc[realNodeTags]
                realCoherentNodeTags=gmshToCoherentNodesNumeration.loc[realNodeTags].to_numpy()[:,0]
                newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()
                newElement.whichPlate  = 1
                newElement.shape =2
                newElement.type = 'timo'
                newElement.integration = 'R'
                elementsList.append(newElement)
                uzElementsList.append(newElement)

                coord = nodesArray[realCoherentNodeTags,0:2]
                lineDirection = coord[1,:] -coord[0,:]
                lineRot = _getVectorDirection(lineDirection)
                plateRots = nodesRotationsPd.loc[nodeTags].to_numpy()
                # print('before: ',nodesRotationsPd)
                
                nodesRotationsPd.loc[realNodeTags] = lineRot
                # nodesRotationsPd = assignNumpyArrayToDataFrame(nodesRotationsPd, realNodeTags, plateRots, lineRot)
                # print('after: ',nodesRotationsPd)

        nDofs = nodesArray.shape[0]*3
        nConstraints = len(newNodesUZ)*3
        hp = myPlate.thickness
        hb = uz.crossSection.thickness
        A = _getMFCMatrix(coherentNodesPlate,coherentNodesUZ,nodesRotationsPd,nDofs,nConstraints,elementType, hb, hp)
        AmatList.append(A)

        uz.elementsList = uzElementsList
    return downStandBeamsList, AmatList, elementsList,nodesRotationsPd

def _getVectorDirection(lineDirection):
    ''' Return direction of the given vector in radians.\n
        Input: \n
            * lineDirection: Numpy array of length 2 with x and y coordinates. \n
        Return: \n
            * lineRot: Direction of the input vector in radians.
    '''
    if lineDirection[0]==0 and lineDirection[1]>0:
        lineRot = np.pi/2
    elif lineDirection[0]==0 and lineDirection[1]<0:
        lineRot = np.pi/2*3
    elif lineDirection[1]==0 and lineDirection[0]>0:
        lineRot = 0
    elif lineDirection[1]==0 and lineDirection[0]<0:
        lineRot = np.pi
    else:
        lineRot = np.arctan(lineDirection[1]/lineDirection[0])
    return lineRot

def _getMFCMatrix(coherentNodesPlate,coherentNodesUZ,nodesRotationsPd,nDofs,nConstraints,elementType, hb, hp):
    ''' Computes the MultiFreedom matrix A. \n
        Input: \n
            * coherentNodesPlate: Nodes of the plate (coherently numerated from 0 to nNodes-1).\n
            * coherentNodesUZ: Nodes of the downStand Beam (coherently numerated from 0 to nNodes-1)\n
            * nodesArrayPd: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates.\n
            * nDOFs: Number of degrees of freedom, including the downstand beam nodes. \n
            * nConstraints: Number of constrained degrees of freedom (= number of downstand beam nodes * 3). \n
            * elementType: String, can be "DB" or "MITC". \n
            * hb: thickness of the downstand beam (without thickness of the plate). \n
            * hp: thickness of the plate. \n
        Return: \n
            * A: (nConstraints x nDOFS) Multifreedom numpy matrix. \n
    '''
    A=np.zeros((nConstraints, nDofs))
    for i, plateNode in enumerate(coherentNodesPlate):
        uzNode = coherentNodesUZ[i]
        nodeRotation = nodesRotationsPd.loc[uzNode+1].to_numpy()[0]


        # print(nodesRotationsPd)
        a1 = np.zeros(nDofs)
        a2 = np.zeros(nDofs)
        a3 = np.zeros(nDofs)

        if elementType == 'DB':
            correspondingRotationDOF = 1
            mult = -1
        elif elementType == 'MITC' and (nodeRotation ==0):
            correspondingRotationDOF = 2
            mult = -1
        elif elementType == 'MITC' and (nodeRotation > 3.0 and nodeRotation<3.3):
            correspondingRotationDOF = 2
            mult = 1
        elif elementType == 'MITC' and (nodeRotation > 1.5 and nodeRotation<1.7):
            correspondingRotationDOF = 1
            mult = -1
        elif elementType == 'MITC' and (nodeRotation > 4.6 and nodeRotation<4.9):
            correspondingRotationDOF = 1
            mult = 1
        else:
            raise TypeError('elementype not defined?')

        a1[plateNode*3] = -1
        a1[uzNode*3+1] = 1
        a2[plateNode*3+correspondingRotationDOF] = -1*mult 
        a2[uzNode*3+2] = 1
        a3[plateNode*3+correspondingRotationDOF] = -(hb+hp)*0.5*mult
        a3[uzNode*3] = 1
        A[3*i:3*i+3, :] = np.array([a1, a2, a3])
    return A

class _Mesh:
    '''
        Stores all informations regarding the mesh. \n
        Input: \n
            * nodesArray: Pandas Dataframe where the indexes are the node tags assigned by gmsh, the values are a nNodes x 3 array with the x-y-z coordinates.\n
            * nodesRotations: pandas dataframe where the indexes are the node tags, the values are the node rotation in radians.\n
            * elementsList: list of element objects. \n
            * BCs: n x 4 numpy array. n is the number of restrained nodes, the first column is the node tag, the other column represent the condition of the three DOFs (1 is blocked, 0 is free). \n
            * AmatList:  List of MultiFreedom contraint matrix of shape nConstraint x nDOFs. nConstraint is 3 x number of downstandBeam nodes. \n
            * plateElementsList: List of element objects before adding the downstand beam elements. \n
            * getElementByTagDictionary: dictionary where the keys are the tags assigned by gmsh, the values are the respective element objects. \n
        Return: \n
        * -
    '''
    def __init__(self,nodesArray, nodesRotations, elementsList, BCs, AmatList, plateElementsList, getElementByTagDictionary,elasticallySupportedNodes):
        self.nodesArray = nodesArray  #pandas!
        self.nodesRotations = nodesRotations
        self.elementsList=elementsList
        self.BCs = BCs
        self.load = None
        self.AmatList = AmatList
        self.getElementByTagDictionary = getElementByTagDictionary
        self.plateElementsList = plateElementsList
        self.elasticallySupportedNodes = elasticallySupportedNodes

class _Element:
    '''
        Stores all information regarding an element
        Atributes: \n
            * tag: tag automatically attributed by gmsh. \n
            * shape: number of nodes.\n
            * nNodes: number of nodes.\n
            * connectivity: node tags defining the element.\n
            * coordinates: x-y-z coordinates of the nodes. \n
            * whichPlate: to which plate belogs the element. \n
            * BbMat: Bending train matrix.\n
            * rotationMatrix: \n
            * coherentConnectivity: coherent node tags. \n
            * type: "DB" or "MITC". \n
            * integration: "R" for reduced, "N" for normal Gauss quadrature. \n
            * correspondingPlateElements: \n
            * Db: Bending material matrix, coming from the plate object. \n
            * Ds: Shear material matrix, coming from the plate object. \n
    '''
    def __init__(self):
        self.tag  = None
        self.shape = None
        self.nNodes  = None
        self.connectivity  = None
        self.coordinates  = None
        self.whichPlate  = None
        self.BbMat = None
        self.rotationMatrix = None
        self.coherentConnectivity = None #rearranges nodes with a sequential nummeration
        self.type = None
        self.integration = None
        self.correspondingPlateElements = None
        self.Db = None
        self.Ds = None

