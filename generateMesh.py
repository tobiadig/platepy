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
def distortMesh(nodesArray, alpha):
    myIndex = nodesArray.index.to_numpy()
    
    nodesArrayNumpy = nodesArray.to_numpy()
    v1=np.ones(nodesArrayNumpy.shape[0])
    x0 = nodesArrayNumpy[:,0]
    y0 = nodesArrayNumpy[:,1]
    a=np.max(x0)
    
    newNodes = np.zeros(nodesArray.shape)
    newNodes[:,0] = x0+(v1-np.abs(2*x0/a-1))*(2*y0/a-v1)*alpha
    newNodes[:,1] = y0+2*(v1-np.abs(2*y0/a-1))*(2*x0/a-v1)*alpha
    
    # xMask = np.logical_or(x0==0, x0==a)
    # yMask = np.logical_or(y0==0, y0==a)
    # newNodes[xMask,0] = x0[xMask]
    # newNodes[yMask,1] = y0[yMask]
    newNodesArray = pd.DataFrame(newNodes, index = myIndex)
    return newNodesArray

def generateMesh(self,showGmshMesh=False,showGmshGeometryBeforeMeshing = False, elementDefinition=None, meshSize=5e-2, nEdgeNodes = 0, order='linear', meshDistortion = False, distVal = 100, deactivateRotation=False):
    ''' Input/Output descriptions
        self: PlateModel class, where the geometry is initialized
        meshInput: options for the creation of the mesh
        return
        the method generates a mesh using gmsh library. The mesh is stored in the PlateModel class.
        following information are also generated and stored:
        nodes coordinates and orientation, element connectivity, boundary conditions
    '''
    temp = elementDefinition.split('-')

    elementType = temp[0]
    elementShape = int(temp[1])
    elementIntegration = temp[2]


    gmsh.model.mesh.clear()
    if elementShape == 4 or elementShape == 9 :
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2) #0: simple, 1: blossom (default), 2: simple full-quad, 3: blossom full-quad
        for i in range(0, len(self.plates)):
            gmsh.model.geo.mesh.setRecombine(2, self.plates[i].tag)
    elif elementShape != 3:
        raise Exception
    
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    gmsh.model.geo.synchronize()

    # manual assignment of edge nodes
    try:
        pointTags=gmsh.model.getEntities(0)
        if nEdgeNodes>0:
            gmsh.model.geo.mesh.setTransfiniteCurve(1, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(2, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(3, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(4, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
            gmsh.model.geo.synchronize()
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

    # mesh generation
    if showGmshGeometryBeforeMeshing:
        gmsh.fltk.run()
    try:
        gmsh.model.mesh.generate()
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.geo.synchronize()
    except:
        print('gmsh mesh generation failed')
        raise


    if order == 'quadratic':
        gmsh.model.mesh.setOrder(2)
        gmsh.model.geo.synchronize()
    elif order !='linear':
        raise Exception('order not recognised')
    
    # eventually open fltk UI
    if showGmshMesh:
        gmsh.fltk.run()

    #generate nodes array
    nodeTagsModel , nodeCoords, _ = gmsh.model.mesh.getNodes()

    nodesArray = np.array(nodeCoords).reshape((-1,3))
    # print('old: ', nodesArray)tanjadahlia
    nodesArrayPd = pd.DataFrame(nodesArray, index = nodeTagsModel)
    nodesRotationsPd = pd.DataFrame(np.zeros(nodeTagsModel.size), index =nodeTagsModel)
    gmshToCoherentNodesNumeration = pd.DataFrame(range(0,len(nodeTagsModel)), index = nodeTagsModel)

    # distort the mesh a bit
    if meshDistortion:
        nodesArrayPd = distortMesh(nodesArrayPd, distVal)
    
    elementsList = []
    # _, elemTags, _ = gmsh.model.mesh.getElements(2)
    k=1
    getElementByTagDictionary ={}

    for i in range(0, len(self.plates)):
        _, elemTags, _ = gmsh.model.mesh.getElements(2,self.plates[i].tag)
        for elemTag in elemTags[0]:
            elemType, nodeTags = gmsh.model.mesh.getElement(elemTag)
            newElement = Element()
            newElement.tag = elemTag
            newElement.whichPlate = i
            # print('element ',elemTag,' belongs to plate ',i)
            newElement.type = elementType
            newElement.shape = elementShape
            newElement.integration = elementIntegration

            newElement.nNodes  = len(nodeTags)
            newElement.connectivity  = nodeTags

            newElement.coherentConnectivity = gmshToCoherentNodesNumeration.loc[nodeTags]

            # old = nodesArray[nodeTags-1,:]

            newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()
            elementsList.append(newElement)
            getElementByTagDictionary[elemTag] = newElement
            k+=1

    plateElementsList = copy.deepcopy(elementsList)

    #assemble 2D elements for line load

    k=1
    for p in self.loads:
        if p.case == 'line':
            tags = gmsh.model.getEntitiesForPhysicalGroup(p.physicalGroup[0],p.physicalGroup[1])
            elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(1,tags[0])   
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
                k+=1
            p.elements1DList = elements1DList

    #generate BCs and nodes directions by iterating over wall segments
    
    BCsDic = {}
    for wall in self.walls:
        dim = wall.physicalGroup[0]
        nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim,wall.physicalGroup[1])
        enitiesTags = gmsh.model.getEntitiesForPhysicalGroup(dim,wall.physicalGroup[1])

        for wallLine in enitiesTags:
            nodeTags, coord, _ = gmsh.model.mesh.getNodes(dim, wallLine, includeBoundary=True)
            coord = coord.reshape(-1,3)
            # start and en nodes are blocked
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
                BCsDic[node] = wall.support.supportCondition
        # print('nodesRotation in mesh gen: ', nodesRotationsPd)

    for col in self.columns:
        dim = col.physicalGroup[0]
        nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim,col.physicalGroup[1])

        for node in nodeTags:
            BCsDic[node] = col.support.supportCondition

    BCs = np.zeros((len(BCsDic), 4))
    for count, a in enumerate(BCsDic):
        BCs[count, 0] = a
        BCs[count,1:] = BCsDic[a]

    nextNode = int(np.max(nodeTagsModel)+1)
    nextCoherentNode = int(nodesArray.shape[0])

    # print('gmshToCoherentNodesNumeration: ',gmshToCoherentNodesNumeration)
    # print('nodesArray: ', nodesArray)
    # print('nodesArrayPd: ', nodesArrayPd)
    # print('nodesRotationsPd: ', nodesRotationsPd)
    AmatList = []
    for uz in self.downStandBeams:
        uzElementsList = []
        dim = uz.physicalGroup[0]
        nodesUZ, coord = gmsh.model.mesh.getNodesForPhysicalGroup(dim,uz.physicalGroup[1])
        # print('nodesUZ: ', nodesUZ)
        nNodesUZ = len(nodesUZ)
        newNodesUZ = np.array(range(nextNode, nextNode+nNodesUZ), dtype=int)
        uz.newNodesUZ = newNodesUZ
        # print('newNodesUZ: ', newNodesUZ)
        coherentNodesUZ = np.array(range(nextCoherentNode, nextCoherentNode+nNodesUZ), dtype=int)
        uz.coherentNodesUZ = coherentNodesUZ
        tempDF = pd.DataFrame(coherentNodesUZ, index=newNodesUZ)

        gmshToCoherentNodesNumeration = gmshToCoherentNodesNumeration.append(tempDF)
        gmshToCoherentNodesNumeration.index = gmshToCoherentNodesNumeration.index.astype(int)
        nextCoherentNode += nNodesUZ
        nextNode = nextNode+nNodesUZ
        tempNodesArray = nodesArrayPd.loc[nodesUZ].to_numpy()
        # print('tempNodesArray ',tempNodesArray)
        nodesArray = np.append(nodesArray, tempNodesArray, axis=0)
        # print('nodesArray ',nodesArray)
        nodesArrayPd = nodesArrayPd.append(pd.DataFrame(tempNodesArray, index=newNodesUZ))
        # print('nodesArrayPd ',nodesArrayPd)
        nodesArrayPd.index = nodesArrayPd.index.astype(int)
        uzNodesToNodesNumeration = pd.DataFrame(newNodesUZ, index = nodesUZ)

        uzNodesToNodesNumeration.index = uzNodesToNodesNumeration.index.astype(int)
        # print('uzNodesToNodesNumeration',uzNodesToNodesNumeration)
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
        #NEED TO CREATE THE CONSTRAINT MATRIX BEAM - PLATE!
        nDofs = nodesArray.shape[0]*3
        nConstraints = len(newNodesUZ)*3
        A=np.zeros((nConstraints, nDofs))

        hp = self.plates[0].thickness
        hb = uz.thickness

        for i, plateNode in enumerate(coherentNodesPlate):
            uzNode = coherentNodesUZ[i]
            a1 = np.zeros(nDofs)
            a2 = np.zeros(nDofs)
            a3 = np.zeros(nDofs)

            if elementType == 'DB':
                correspondingRotationDOF = 1
                mult = -1
            elif elementType == 'MITC':
                correspondingRotationDOF = 2
                mult = -1

            a1[plateNode*3] = -1
            a1[uzNode*3+1] = 1
            a2[plateNode*3+correspondingRotationDOF] = -1*mult 
            a2[uzNode*3+2] = 1
            a3[plateNode*3+correspondingRotationDOF] = -(hb+hp)*0.5*mult
            a3[uzNode*3] = 1
            A[3*i:3*i+3, :] = np.array([a1, a2, a3])
        
        AmatList.append(A)

        enitiesTags = gmsh.model.getEntitiesForPhysicalGroup(dim,uz.physicalGroup[1])

        for uzLine in enitiesTags:
            # nodeTags, coord, _ = gmsh.model.mesh.getNodes(dim, uzLine, includeBoundary=True)
            # coord = coord.reshape(-1,3)
            elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim,uzLine)
            for elemTag in elementTags[0]:
                elementType, nodeTags = gmsh.model.mesh.getElement(elemTag)
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
                    lineRot = np.pi/2*0
                elif lineDirection[0]==0 and lineDirection[1]<0:
                    lineRot = np.pi/2*3*0
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

        uz.elementsList = uzElementsList

        # print('nodesRotation in mesh gen: ', nodesRotationsPd)
    # print('gmshToCoherentNodesNumeration: ',gmshToCoherentNodesNumeration)
    # print('nodesArray: ', nodesArray)
    # print('nodesArrayPd: ', nodesArrayPd)
    # print('nodesRotationsPd: ', nodesRotationsPd)
    self.mesh = Mesh(nodesArrayPd,nodesRotationsPd, elementsList, BCs)
    self.mesh.AmatList = AmatList
    self.mesh.plateElementsList = plateElementsList
    self.mesh.getElementByTagDictionary = getElementByTagDictionary

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


def assignNumpyArrayToDataFrame(xDF, indexesToModify, plateRotations, uzRotation):
    index = xDF.index
    values = xDF.to_numpy()
    nIndexes = len(indexesToModify)
    values[indexesToModify-1] = uzRotation * np.ones((nIndexes,1)) - plateRotations

    xDF = pd.DataFrame(values, index=index)
    return xDF