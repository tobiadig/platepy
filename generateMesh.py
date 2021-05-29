''' Module Information
-----------------------------------------------------------
Purpose of module: generates the mesh and stores it as attribute of the plateModel class
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''
#%% Basic modules
import numpy as np
import pandas as pd
import copy
import gmsh # To create CAD model and mesh

def generateMesh(self,showGmshMesh=False,showGmshGeometryBeforeMeshing = False, elementDefinition=None, \
    meshSize=5e-2, nEdgeNodes=0, order='linear', meshDistortion = False, distVal = 100,\
        deactivateRotation=False):
    ''' Generates mesh and stores it in the plateModel object according to the selected options.\n
            Input: \n
            * self: plateModel object. \n
            * showGmshMesh = False: If True, the model is shown through the built-in fltk terminal (after the mesh generation). \n
            * showGmshMesh = True: If True, the model is shown through the built-in fltk terminal (before the mesh generation). \n
            * elementDefinition = None: String defining the desired FE-element in the following form: "type-nNodes-integration". \n
            \t * type: DB for displacement-based elements or MITC. \n
            \t * nNodes: number of nodes ( currently 3, 4 or 9). \n
            \t * integration: Desired Gauss-quadrature for the calculation of the stiffness matrixes. R for reduced or N for normal. \n
            * meshSize = 5e-2: target mesh size around the point entities. If nEdgeNodes > 0, meshSize is ignored. \n
            * nEdgeNodes = 0: Prescribes the number of nodes on each edge. Plate must be rectangular.
            * order = "linear": Prescribes the order of the elements. "linear for first 
            order elements (default) and "quadratic" for second order elements. \n
            * meshDistortion = False: Boolean, if True a mesh distortion function is applied. \n
            * distVal = Severeness of the mesh distortion function. \n
            Return: \n
            *   -
    '''

    elementType, elementShape, elementIntegration = _getElementDefinition(elementDefinition)
    gmsh.model.mesh.clear()

    # Set recombination rules
    if elementShape == 4 or elementShape == 9 :
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2) #0: simple, 1: blossom (default), 2: simple full-quad, 3: blossom full-quad
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
        #     gmsh.model.geo.mesh.setTransfiniteCurve(1, nEdgeNodes, "Progression", progVal)
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
    #gmshNumeration is the one with node Tags, CoherentNumeration is from 0 to nNodes-1

    # distort the mesh if required
    if meshDistortion:
        nodesArrayPd = distortMesh(nodesArrayPd, distVal)
    
    elementsList = []
    getElementByTagDictionary ={}

    # Assign all the information to each element object and save it in the elementsList
    for i in range(0, len(self.plates)):
        _, elemTags, _ = gmsh.model.mesh.getElements(2,self.plates[i].tag)
        for elemTag in elemTags[0]:
            _ , nodeTags = gmsh.model.mesh.getElement(elemTag)
            newElement = Element()
            newElement.tag = elemTag
            newElement.whichPlate = i
            newElement.Db = self.plates[i].Df 
            newElement.Ds = self.plates[i].Dc 
            newElement.type = elementType
            newElement.shape = elementShape
            newElement.integration = elementIntegration
            newElement.nNodes  = len(nodeTags)
            newElement.connectivity  = nodeTags
            newElement.coherentConnectivity = gmshToCoherentNodesNumeration.loc[nodeTags]
            newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()

            elementsList.append(newElement)
            getElementByTagDictionary[elemTag] = newElement

    plateElementsList = copy.deepcopy(elementsList)  # usefull if there is the need to distinguish plate elements from the downstand bem elements

    #assemble 2D elements for line load
    for p in self.loads:
        if p.case == 'line':
            tags = gmsh.model.getEntitiesForPhysicalGroup(p.physicalGroup[0],p.physicalGroup[1])
            _, elementTags, nodeTags = gmsh.model.mesh.getElements(1,tags[0])   
            elements1DList = []
            for elemTag in elementTags[0]:
                elementType, nodeTags = gmsh.model.mesh.getElement(elemTag)
                newElement = Element()
                newElement.tag = elemTag
                newElement.nNodes  = len(nodeTags)
                newElement.connectivity  = nodeTags
                newElement.coherentConnectivity = gmshToCoherentNodesNumeration.loc[nodeTags]
                newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()
                newElement.whichPlate  = 1  
                elements1DList.append(newElement)
            p.elements1DList = elements1DList

    #generate BCs and nodes directions by iterating over wall segments
    wallList = self.walls
    gmshModel = gmsh.model
    columnsList = self.columns
    BCs, nodesRotationsPd = getBCsArray(gmshModel,wallList,columnsList,nodesRotationsPd,deactivateRotation)



    nextNode = int(np.max(nodeTagsModel)+1)
    nextCoherentNode = int(nodesArray.shape[0])

    AmatList = []
    for uz in self.downStandBeams:
        uzElementsList = []
        dim = uz.physicalGroup[0]
        nodesUZ, coord = gmsh.model.mesh.getNodesForPhysicalGroup(dim,uz.physicalGroup[1])
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

    for uz in self.downStandBeams:
        newNodesUZ = uz.newNodesUZ
        uzNodesToNodesNumeration = uz.uzNodesToNodesNumeration 
        coherentNodesPlate = uz.coherentNodesPlate
        coherentNodesUZ = uz.coherentNodesUZ
        dim = uz.physicalGroup[0]
        nDofs = nodesArray.shape[0]*3
        nConstraints = len(newNodesUZ)*3
        A=np.zeros((nConstraints, nDofs))

        hp = self.plates[0].thickness
        hb = uz.thickness
        enitiesTags = gmsh.model.getEntitiesForPhysicalGroup(dim,uz.physicalGroup[1])
        for uzLine in enitiesTags:
            # nodeTags, coord, _ = gmsh.model.mesh.getNodes(dim, uzLine, includeBoundary=True)
            # coord = coord.reshape(-1,3)
            elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim,uzLine)
            for elemTag in elementTags[0]:
                _, nodeTags = gmsh.model.mesh.getElement(elemTag)
                newElement = Element()
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

                p1Tag = realCoherentNodeTags[0]

                p2Tag = realCoherentNodeTags[1]
                coord = nodesArray[realCoherentNodeTags,0:2]

                lineDirection = coord[1,:] -coord[0,:]

                #particular directions (np.arctan not defined)
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


                plateRots = nodesRotationsPd.loc[nodeTags].to_numpy()
                # print('before: ',nodesRotationsPd)
                
                nodesRotationsPd.loc[realNodeTags] = lineRot
                # nodesRotationsPd = assignNumpyArrayToDataFrame(nodesRotationsPd, realNodeTags, plateRots, lineRot)
                # print('after: ',nodesRotationsPd)

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
        
        AmatList.append(A)

        uz.elementsList = uzElementsList

    self.mesh = Mesh(nodesArrayPd,nodesRotationsPd, elementsList, BCs)
    self.mesh.AmatList = AmatList
    self.mesh.plateElementsList = plateElementsList
    self.mesh.getElementByTagDictionary = getElementByTagDictionary












































def getBCsArray(gmshModel,wallList,columnsList,nodesRotationsPd,deactivateRotation):
    BCsDic = {}
    for wall in wallList:
        dim = wall.physicalGroup[0]
        nodeTags, _ = gmshModel.mesh.getNodesForPhysicalGroup(dim,wall.physicalGroup[1])
        enitiesTags = gmshModel.getEntitiesForPhysicalGroup(dim,wall.physicalGroup[1])

        for wallLine in enitiesTags:
            nodeTags, coord, _ = gmshModel.mesh.getNodes(dim, wallLine, includeBoundary=True)
            coord = coord.reshape(-1,3)
            # start and end nodes are blocked
            p1Tag = nodeTags[-2]
            BCsDic[p1Tag] = np.array([1,1,1])
            p2Tag = nodeTags[-1]
            BCsDic[p2Tag] = np.array([1,1,1])
            lineDirection = coord[-1,:] -coord[-2,:]

            #particular directions (np.arctan not defined)
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
            
            nodesRotationsPd.loc[nodeTags[:]] = lineRot

            for node in nodeTags[:]:
                BCsDic[node] = wall.support

    for col in columnsList:
        dim = col.physicalGroup[0]
        nodeTags, _ = gmshModel.mesh.getNodesForPhysicalGroup(dim,col.physicalGroup[1])

        for node in nodeTags:
            BCsDic[node] = col.support

    BCs = np.zeros((len(BCsDic), 4))
    for count, a in enumerate(BCsDic):
        BCs[count, 0] = a
        BCs[count,1:] = BCsDic[a]

    return BCs, nodesRotationsPd


def setMesh(self, nodesArray, elements, BCs, load = None):
    elementsList = []
    nNodes = nodesArray.shape[0]
    k=1
    for element in elements:
        newElement = Element()
        newElement.tag = k

        newElement.nNodes  = len(element)
        newElement.connectivity  = element
        newElement.shape=4
        newElement.coherentConnectivity = pd.DataFrame(element-1)

        newElement.coordinates = np.zeros((len(element), 3))
        newElement.coordinates[:,0:2] = nodesArray[element-1, :]

        newElement.whichPlate  = 1  
        elementsList.append(newElement)
        k+=1

    nodesRotationsPd = pd.DataFrame(np.zeros((nNodes, 1)), index=range(1, nNodes+1))
    nodesArrayPd =pd.DataFrame(nodesArray, index=range(1, nNodes+1) )
    self.mesh = Mesh(nodesArrayPd,nodesRotationsPd, elementsList, BCs)
    self.mesh.load = load


class Mesh:
    '''
        stores informations about nodes coordinates and rotation, element connectivity and boundary conditions
    '''
    def __init__(self,nodesArray, nodesRotations, elementsList, BCs):
        self.nodesArray = nodesArray  #pandas!
        self.nodesRotations = nodesRotations
        self.elementsList=elementsList
        self.BCs = BCs
        self.load = None
        self.AmatList = None
        self.getElementByTagDictionary = None

class Element:
    '''
        stores all information regarding a particular element
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

def distortMesh(nodesArray, alpha):
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

def assignNumpyArrayToDataFrame(xDF, indexesToModify, plateRotations, uzRotation):
    index = xDF.index
    values = xDF.to_numpy()
    nIndexes = len(indexesToModify)
    values[indexesToModify-1] = uzRotation * np.ones((nIndexes,1)) - plateRotations

    xDF = pd.DataFrame(values, index=index)
    return xDF

def _getElementDefinition(elementDefinition):
    '''
        Extracts information from the elementDefinition String. \n
        Input: \n
        * elementDefinition : String defining the desired FE-element in the following form: "type-nNodes-integration". \n
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