''' Module Information
-----------------------------------------------------------
Purpose of module: generates the mesh and stores it as attribute of the plateModel class
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''
#%% Basic modules
import numpy as np
import pandas as pd

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

    xMask = np.logical_or(x0==0, x0==a)
    yMask = np.logical_or(y0==0, y0==a)
    newNodes[xMask,0] = x0[xMask]
    newNodes[yMask,1] = y0[yMask]
    newNodesArray = pd.DataFrame(newNodes, index = myIndex)
    return newNodesArray

def generateMesh(self,showGmshMesh=False, elementType = 'QUAD', meshSize=5e-2, nEdgeNodes = 0, order='linear', meshDistortion = False, progVal=1.2):

    ''' Input/Output descriptions
        self: PlateModel class, where the geometry is initialized
        meshInput: options for the creation of the mesh

        return
        the method generates a mesh using gmsh library. The mesh is stored in the PlateModel class.
        following information are also generated and stored:
        nodes coordinates and orientation, element connectivity, boundary conditions
    '''
    
    gmsh.model.mesh.clear()
    if elementType == 'QUAD':
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1) #0: simple, 1: blossom (default), 2: simple full-quad, 3: blossom full-quad
        gmsh.model.geo.mesh.setRecombine(2, 1)
    elif elementType != 'TRI':
        raise TypeError(meshInput.elementType, "not recognised")
    
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    gmsh.model.geo.synchronize()

    # manual assignment of edge nodes
    try:
        pointTags=gmsh.model.getEntities(0)
        if nEdgeNodes>0 and not(meshDistortion):
            gmsh.model.geo.mesh.setTransfiniteCurve(1, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(2, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(3, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(4, nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
            gmsh.model.geo.synchronize()
        elif nEdgeNodes>0 and meshDistortion:
            gmsh.model.geo.mesh.setTransfiniteCurve(1, nEdgeNodes, "Progression", progVal)
            gmsh.model.geo.mesh.setTransfiniteCurve(2, nEdgeNodes, "Progression", progVal)
            gmsh.model.geo.mesh.setTransfiniteCurve(3, nEdgeNodes, "Progression", progVal)
            gmsh.model.geo.mesh.setTransfiniteCurve(4, nEdgeNodes, "Progression", progVal)
            gmsh.model.geo.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
            gmsh.model.geo.synchronize()
        else:
            gmsh.model.mesh.setSize(pointTags,meshSize)
    except:
        print('manual assignment of edge nodes failed')
        raise


    # mesh generation
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
    # print('old: ', nodesArray)
    nodesArrayPd = pd.DataFrame(nodesArray, index = nodeTagsModel)
    nodesRotationsPd = pd.DataFrame(np.zeros(nodeTagsModel.size), index =nodeTagsModel)
    gmshToCoherentNodesNumeration = pd.DataFrame(range(0,len(nodeTagsModel)), index = nodeTagsModel)

    # distort the mesh a bit
    alpha = 120
    nodesArrayPd = distortMesh(nodesArrayPd, alpha)

    elementsList = []
    _, elemTags, _ = gmsh.model.mesh.getElements(2)  

    k=1
    for elemTag in elemTags[0]:
        elemType, nodeTags = gmsh.model.mesh.getElement(elemTag)




        newElement = Element()
        newElement.tag = elemTag

        if elementType == 'QUAD':
            newElement.shape = 4
        else:
            newElement.shape =3       

        newElement.nNodes  = len(nodeTags)
        newElement.connectivity  = nodeTags

        newElement.coherentConnectivity = gmshToCoherentNodesNumeration.loc[nodeTags]

        # old = nodesArray[nodeTags-1,:]

        newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()

        newElement.whichPlate  = 1 
        elementsList.append(newElement)
        k+=1

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
            if lineDirection[0]==0 and lineDirection[1]>0:
                lineRot = np.pi/2
            elif lineDirection[0]==0 and lineDirection[1]<0:
                lineRot = np.pi/2*3
            else:
                lineRot = np.arctan(lineDirection[1]/lineDirection[0])
            
            nodesRotationsPd.loc[nodeTags[:-2]] = lineRot

            for node in nodeTags[:-2]:
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

        #generate BCs
    # iterate over walls and columns
    # isBCSinitiated = False
    # for wall in self.walls:
    #     nodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(wall.physicalGroup[0],wall.physicalGroup[1])
    #     BCsTemp = np.zeros((len(nodeTags), 4))
    #     BCsTemp[:,0]=nodeTags
    #     BCsTemp[:,1]=np.ones(len(nodeTags))*wall.support.supportCondition[0]
    #     BCsTemp[:,2]=np.ones(len(nodeTags))*wall.support.supportCondition[1]
    #     BCsTemp[:,3]=np.ones(len(nodeTags))*wall.support.supportCondition[2]  

    #     if not(isBCSinitiated):
    #         isBCSinitiated=True
    #         BCs = BCsTemp
    #     else:
    #         BCs = np.append(BCs,BCsTemp, axis=0)

    # store result in a mesh class and then in plateModel
    # print('nodesRotationsPd: ', nodesRotationsPd)

    
    self.mesh = Mesh(nodesArrayPd,nodesRotationsPd, elementsList, BCs)

def setMesh(self, nodesArray, elements, BCs, load = None):
    elementsList = []
    nNodes = nodesArray.shape[0]
    k=1
    for element in elements:
        newElement = Element()
   
        newElement.tag = k

        newElement.nNodes  = len(element)
        newElement.connectivity  = element
        newElement.shape=len(element)
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
