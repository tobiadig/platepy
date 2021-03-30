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

def generateMesh(self,meshInput):
    ''' Input/Output descriptions
        self: PlateModel class, where the geometry is initialized
        meshInput: options for the creation of the mesh

        return
        the method generates a mesh using gmsh library. The mesh is stored in the PlateModel class.
        following information are also generated and stored:
        nodes coordinates and orientation, element connectivity, boundary conditions
    '''
    if meshInput.elementType == 'QUAD':
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1) #0: simple, 1: blossom (default), 2: simple full-quad, 3: blossom full-quad
        gmsh.model.geo.mesh.setRecombine(2, 1)
    elif meshInput.elementType != 'TRI':
        raise TypeError(meshInput.elementType, "not recognised")
    
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    gmsh.model.geo.synchronize()

    # manual assignment of edge nodes
    try:
        pointTags=gmsh.model.getEntities(0)
        if meshInput.nEdgeNodes>0:
            gmsh.model.geo.mesh.setTransfiniteCurve(1, meshInput.nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(2, meshInput.nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(3, meshInput.nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteCurve(4, meshInput.nEdgeNodes)
            gmsh.model.geo.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
            gmsh.model.geo.synchronize()
        else:
            gmsh.model.mesh.setSize(pointTags,meshInput.meshSize)
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
    
    # eventually open fltk UI
    if meshInput.showGmshMesh:
        gmsh.fltk.run()

    #generate nodes array
    nodeTagsModel , nodeCoords, _ = gmsh.model.mesh.getNodes()

    nodesArray = np.array(nodeCoords).reshape((-1,3))
    # print('old: ', nodesArray)
    nodesArrayPd = pd.DataFrame(nodesArray, index = nodeTagsModel)
    nodesRotationsPd = pd.DataFrame(np.zeros(nodeTagsModel.size), index =nodeTagsModel)



    # generate elements list of class Element
    # iterate over 2D surfaces and get nodeTags of each mesh-element
    elementsList = []
    _, elemTags, _ = gmsh.model.mesh.getElements(2)  #TODO: implement more plates

    k=1
    for elemTag in elemTags[0]:
        elementType, nodeTags = gmsh.model.mesh.getElement(elemTag)
        newElement = Element()
        newElement.tag = elemTag

        newElement.nNodes  = len(nodeTags)
        newElement.connectivity  = nodeTags

        # old = nodesArray[nodeTags-1,:]

        newElement.coordinates = nodesArrayPd.loc[nodeTags].to_numpy()

        newElement.whichPlate  = 1  #TODO: implement more plates
        elementsList.append(newElement)
        k+=1

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
    self.mesh = Mesh(nodesArray,nodesRotationsPd, elementsList, BCs)

class MeshInput:
    '''
        options for the creation of the model mesh
    '''
    def __init__(self, showGmshMesh=False, elementType = 'QUAD', meshSize=5e-2, nEdgeNodes = 0):
        self.showGmshMesh = showGmshMesh
        self.elementType = elementType
        self.meshSize = meshSize
        self.nEdgeNodes = nEdgeNodes

class Mesh:
    '''
        stores informations about nodes coordinates and rotation, element connectivity and boundary conditions
    '''
    def __init__(self,nodesArray, nodesRotations, elementsList, BCs):
        self.nodesArray = nodesArray
        self.nodesRotations = nodesRotations
        self.elementsList=elementsList
        self.BCs = BCs

class Element:
    '''
        stores all information regarding a particular element
    '''
    def __init__(self):
        self.tag  = None
        self.nNodes  = None
        self.connectivity  = None
        self.coordinates  = None
        self.whichPlate  = None
        self.BbMat = None
        self.rotationMatrix = None
