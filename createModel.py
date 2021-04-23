''' Module Information
-----------------------------------------------------------
Purpose of module: define the classes and methods needed to initialize the model adn define the geometry
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''

# basic modules
import numpy as np
import matplotlib.pyplot as plt

import gmsh # To create CAD model and mesh

class PlateModel:   
    '''
        contains all the information about the current model (geometry, physical properties, mesh, results)
        contains following methods to add structural elements and loads to the model:
        addPlate, addWall, addColumn, addLoad
    '''
    def __init__(self, name):
        self.name = name
        self.loads = []
        self.plates = []
        self.walls = []
        self.columns = []
        self.joints = []
        self.downStandBeams = []
        self.isInitialized = False
        self.mesh = None
        self.results = None
        self.axes = {}

    def addPlate(self, newPlate):  
        '''
        add a plate element of class "Plate" to the model
        '''
        if not(self.isInitialized):
            gmsh.initialize()
            gmsh.option.setNumber("General.Verbosity",0)
            gmsh.model.add(self.name)
            self.isInitialized =True
        
        self.plates.append(newPlate)    # list of the plates contained in the model
        
        #create points according to the given coordinates
        nPoints = len(newPlate.outlineCoords)
        pointTags = [0]*nPoints
        i=0
        for newPoint in newPlate.outlineCoords:
            pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)
            i=i+1

        #create lines
        nLines = nPoints*1-1
        linesTags = [0]*nLines
        for i in range(0, nLines):
            linesTags[i] = gmsh.model.geo.addLine(pointTags[i], pointTags[i+1])
        curveLoopTag = gmsh.model.geo.addCurveLoop(linesTags)
        planeSurfaceTag= gmsh.model.geo.addPlaneSurface([curveLoopTag])
        newPlate.elementComposition.append(planeSurfaceTag)

        # create physical group
        physicalTag = gmsh.model.addPhysicalGroup(2, [planeSurfaceTag])
        newPlate.physicalGroup = (2, physicalTag) # assign a tuple to the plate with: (dim, physTag)

        gmsh.model.geo.synchronize()

    def addWall(self, newWall):
        '''
            add a plate element of class "Plate" to the model
        '''
        if not(self.isInitialized):
            gmsh.initialize()
            gmsh.model.add(self.name)
            self.isInitialized =True

        self.walls.append(newWall)
        
        #create points according to the given coordinates
        nPoints = len(newWall.outlineCoords)
        pointTags = [0]*nPoints
        i=0
        for newPoint in newWall.outlineCoords:
            pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)
            i=i+1

        #create lines
        nLines = nPoints*1-1
        linesTags = [0]*nLines
        for i in range(0, nLines):
            linesTags[i] = gmsh.model.geo.addLine(pointTags[i], pointTags[i+1])
        
        #create physical group
        physicalTag = gmsh.model.addPhysicalGroup(1, linesTags)
        newWall.physicalGroup = (1, physicalTag)

        gmsh.model.geo.synchronize()
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        entities = np.array(gmsh.model.getEntities(1))
        for line in linesTags:
            if line in entities[:,1]:
                gmsh.model.mesh.embed(1, [line], 2, 1) #TODO: the wall should recognize to which plate it belongs. temporarly just plate 1
        
        gmsh.model.geo.synchronize()

    def addColumn(self, newColumn):
        '''
            add column to the plateMode
        '''
        if not(self.isInitialized):
            gmsh.initialize()
            gmsh.model.add(self.name)
            self.isInitialized =True

        self.columns.append(newColumn)

        #create points according to the given coordinates
        nPoints = len(newColumn.outlineCoords)
        pointTags = [0]*nPoints
        i=0
        for newPoint in newColumn.outlineCoords:
            pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)
            i=i+1

        #creates physical group
        physicalTag = gmsh.model.addPhysicalGroup(0, pointTags)
        newColumn.physicalGroup = (0, physicalTag)

        gmsh.model.geo.synchronize()
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        if newColumn.isInPlate:
            gmsh.model.mesh.embed(0, pointTags, 2, 1)#TODO: column should recognize to which plate it belongs. temporarly just plate 1
            gmsh.model.geo.synchronize()
        else:
            gmsh.model.mesh.embed(0, pointTags, 1, 3)#TODO: column should recognize to which line it belongs
            gmsh.model.geo.synchronize()

    def addLoad(self, newLoad):
        '''
            add a load to the model
        '''
        self.loads.append(newLoad)
        
        if newLoad.case =='area':
            pass
        elif newLoad.case == 'line':
            nPoints = len(newLoad.outlineCoords)
            pointTags = [0]*nPoints
            i=0
            for newPoint in newLoad.outlineCoords:
                pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)
                i=i+1

            #create lines
            nLines = nPoints*1-1
            linesTags = [0]*nLines
            for i in range(0, nLines):
                linesTags[i] = gmsh.model.geo.addLine(pointTags[i], pointTags[i+1])

            #create physical group
            physicalTag = gmsh.model.addPhysicalGroup(1, linesTags)
            newLoad.physicalGroup = (1, physicalTag)

            gmsh.model.geo.synchronize()
            gmsh.model.geo.removeAllDuplicates()
            gmsh.model.geo.synchronize()
        elif newLoad.case == 'nodes':
            pass
    def clearMesh(self):
        self.mesh = None

class Plate:
    def __init__(self, inputDict):
        self.outlineCoords = inputDict["outlineCoords"]
        self.thickness = inputDict["thickness"]
        self.surfaceLevel = inputDict["surfaceLevel"]
        self.body = inputDict["body"]
        self.stiffnessFactor = inputDict["stiffnessFactor"]
        self.physicalGroup = None
        self.elementComposition = []
        self.nodeComposition = []

        E=self.body.eModule
        G=self.body.gModule
        nu=self.body.nu
        h=self.thickness
        self.Df=h**3/12*E/(1-nu**2)*np.array([[1, nu, 0],
                                                [nu, 1, 0],
                                                [0, 0, (1-nu)/2]])

        self.Dc = 5/6*h*np.array([[G,0],[0,G]]) #alpha = 5/6 according to ferreira p. 162

    def plot(self, axGeometry):
        coords = self.outlineCoords
        axGeometry.plot(coords[:,0],coords[:,1], color='gray')

class Wall:
    def __init__(self, inputDict):
        self.outlineCoords = inputDict["outlineCoords"]
        self.high = inputDict["high"]
        self.body = inputDict["body"]
        self.support = inputDict["support"]
        self.thickness = inputDict["thickness"]
        self.physicalGroup = None
        self.elementComposition = []
        self.nodeComposition = None

    def plot(self, axGeometry):
        coords = self.outlineCoords
        axGeometry.plot(coords[:,0],coords[:,1], color='g')
        for i in range(0, coords.shape[0]-1):
            p1 = coords[i]
            p2 = coords[i+1]
            d = self.thickness/2

            wallDir = p2 - p1
            length = np.dot(wallDir, wallDir)**0.5
            normWallDir = -np.array([wallDir[1], wallDir[0]])/length
            l1 = np.zeros((2,2))
            l2 = np.zeros((2,2))

            l1[0,:] = p1 + d*normWallDir
            l1[1,:] = p2 + d*normWallDir

            l2[0,:] = p1 - d*normWallDir
            l2[1,:] = p2 - d*normWallDir
            linWidth = 0.7
            axGeometry.plot(l1[:,0],l1[:,1], color='g', linewidth = linWidth)
            axGeometry.plot(l2[:,0],l2[:,1], color='g', linewidth = linWidth)

class Column:
    def __init__(self, inputDict, isInPlate = False):
        self.outlineCoords = inputDict["outlineCoords"]
        self.high = inputDict["high"]
        self.body = inputDict["body"]
        self.support = inputDict["support"]
        self.crossSection = inputDict["crossSection"]
        self.width = inputDict["width"]
        self.physicalGroup = None
        self.elementComposition = None
        self.nodeComposition = None
        self.isInPlate = isInPlate

    def plot(self, ax):
        x = self.outlineCoords[0,0]
        y= self.outlineCoords[0,1]
        d = self.width/2
        x1=x-d
        x2=x+d
        y1 = y-d
        y2=y+d

        xi = np.array([x1, x2, x2, x1, x1])
        yi= np.array([y1, y1, y2, y2, y1])

        ax.plot(xi,yi, color='b')
        ax.scatter(x,y,marker="D")

class Concrete:
    def __init__(self, inputDict):
        self.eModule = inputDict["eModule"]
        self.gModule = inputDict["gModule"]
        self.density = inputDict["density"]
        self.nu = inputDict["nu"]

class Support:
    def __init__(self,supportCondition):
        self.supportCondition = supportCondition

class CrossSection:
    def __init__(self, A, Iy, Iz):
        self.A = A
        self.Iy = Iy
        self.Iz = Iz

class Load:
    def __init__(self,case, magnitude):
        self.magnitude = magnitude
        self.case=case
        self.outlineCoords = np.array([])
        self.physicalGroup = None
        self.elements1DList = None
        self.nodePattern = None


