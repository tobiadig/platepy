''' Defines the classes and methods needed to initialize the model and define the geometry.


Copywrite Tobia Diggelmann (ETH Zurich) 29.05.2021
'''

# basic modules
import numpy as np
import gmsh # To create CAD model and mesh

class PlateModel:   
    '''
    A PlateModel object is used as basis for the calculation of FE-model. Is used to store structural elements, mesh and results.
    '''
    def __init__(self):
        gmsh.initialize()
        gmsh.model.remove()
        gmsh.option.setNumber("General.Verbosity",0)
        self.loads = []
        self.plates = []
        self.walls = []
        self.columns = []
        self.joints = []
        self.downStandBeams = []
        self.isInitialized = False
        self.mesh = None
        self.resultsInformation = None
        self.axes = {}
    def addPlate(self, newPlate):  
        '''
        Add a plate object to the plateModel.

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        **newPlate**: `Plate` object to be added.

        '''
        self.plates.append(newPlate)

        # Create points in the Gmsh model
        nPoints = len(newPlate.outlineCoords)
        pointTags = [0]*nPoints
        for i, newPoint in enumerate(newPlate.outlineCoords):
            pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)

        # Create lines in the Gmsh model
        nLines = nPoints*1-1
        linesTags = [0]*nLines
        for i in range(0, nLines):
            linesTags[i] = gmsh.model.geo.addLine(pointTags[i], pointTags[i+1])
        
        # Create surface in the Gmsh model
        curveLoopTag = gmsh.model.geo.addCurveLoop(linesTags)
        planeSurfaceTag= gmsh.model.geo.addPlaneSurface([curveLoopTag])
        newPlate.tag=planeSurfaceTag
        
        # Assign surface to a physical group
        physicalTag = gmsh.model.addPhysicalGroup(2, [planeSurfaceTag])
        newPlate.physicalGroup = (2, physicalTag) # assign a tuple to the plate with: (dim, physTag)
        
        gmsh.model.geo.synchronize()
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

    def addWall(self, newWall):
        '''
        Add a wall to the plateModel. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        **newWall**: `Wall` object to be added. 

        '''
        self.walls.append(newWall)
        
        # Create points in the Gmsh model
        nPoints = len(newWall.outlineCoords)
        pointTags = [0]*nPoints
        i=0
        for newPoint in newWall.outlineCoords:
            pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)
            i=i+1

        # Create lines in the Gmsh model
        nLines = nPoints*1-1
        linesTags = [0]*nLines
        for i in range(0, nLines):
            linesTags[i] = gmsh.model.geo.addLine(pointTags[i], pointTags[i+1])
        
        # Assign line to a physical group
        physicalTag = gmsh.model.addPhysicalGroup(1, linesTags)
        newWall.physicalGroup = (1, physicalTag)

        gmsh.model.geo.synchronize()
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        # Embed lines into the surface (walls embedded into the plate)
        entities = np.array(gmsh.model.getEntities(1))
        for line in linesTags:
            if line in entities[:,1]:
                gmsh.model.mesh.embed(1, [line], 2, 1) 
        gmsh.model.geo.synchronize()

    def addDownStandBeam(self, newDSB):
        '''
        Add a downstand beam to the plateModel. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        **newDSB**: `DownStandBeam` object to be added. 

        '''
        self.downStandBeams.append(newDSB)
        
        # Create points in the Gmsh model
        nPoints = len(newDSB.outlineCoords)
        pointTags = [0]*nPoints
        i=0
        for newPoint in newDSB.outlineCoords:
            pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)
            i=i+1

        # Create lines in the Gmsh model
        nLines = nPoints*1-1
        linesTags = [0]*nLines
        for i in range(0, nLines):
            linesTags[i] = gmsh.model.geo.addLine(pointTags[i], pointTags[i+1])
        
        # Assign line to a physical group
        physicalTag = gmsh.model.addPhysicalGroup(1, linesTags)
        newDSB.physicalGroup = (1, physicalTag) #(dimension, physicalTag)

        gmsh.model.geo.synchronize()
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        entities = np.array(gmsh.model.getEntities(1))
        for line in linesTags:
            if line in entities[:,1]:
                gmsh.model.mesh.embed(1, [line], 2, 1) 

    def addColumn(self, newColumn):
        '''
        Add a column to the plateModel. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        **newColumn**: `Column` object to be added.

        '''
        self.columns.append(newColumn)

        # Create points in the Gmsh model
        nPoints = len(newColumn.outlineCoords)
        pointTags = [0]*nPoints
        i=0
        for newPoint in newColumn.outlineCoords:
            pointTags[i] = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)
            i=i+1

        # Assign point to a physical group
        physicalTag = gmsh.model.addPhysicalGroup(0, pointTags)
        newColumn.physicalGroup = (0, physicalTag)

        gmsh.model.geo.synchronize()
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        if newColumn.isInPlate:
            gmsh.model.mesh.embed(0, pointTags, 2, 1)
            gmsh.model.geo.synchronize()
        else:
            pass # I didn't manage to embed a point into a line,
                # user has to define a plate corner where the column will be...
            # gmsh.model.mesh.embed(0, pointTags, 1, 3)
            # gmsh.model.geo.synchronize()

    def addLoad(self, newLoad):
        '''
        Add a load to the plateModel.

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        **newLoad**: `Load` object to be added.

        '''
        self.loads.append(newLoad)
        
        if newLoad.case =='area':
            pass # in this case the force is by default applied in the solve module
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
        elif newLoad.case == 'point':

            newPoint = newLoad.outlineCoords
            pointTags = gmsh.model.geo.addPoint(newPoint[0], newPoint[1], 0)

            #create physical group
            physicalTag = gmsh.model.addPhysicalGroup(0, [pointTags])
            newLoad.physicalGroup = (0, physicalTag)
            

            gmsh.model.geo.synchronize()
            gmsh.model.mesh.embed(0,[pointTags], 2, 1)
            gmsh.model.geo.synchronize()
            gmsh.model.geo.removeAllDuplicates()
            gmsh.model.geo.synchronize()

    def clearMesh(self):
        '''Method to delete the current mesh of the plateModel.'''
        self.mesh = None

class Plate:
    def __init__(self, inputDict, isUnterZug=False, t=0):
        '''
        A plate object contains all characteristics regarding geometry and material.

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **inputDict**: dictionary with following entries:
            - "outlineCoords": n x 2 numpy array with n x-y couples, representing the boundaries of the plate  
            - "thickness": thickness of the plate 
            - "body": `Concrete` object 
        * **isUnterZug = False**: Boolean, True if the plate aims to modell a downstand beam (default is False)
        * **t = 0**: If the plate aims to modell a downstand beam, thickness of the sorrounding plate (default is 0) 

        '''
        self.outlineCoords = inputDict["outlineCoords"]
        self.thickness = inputDict["thickness"]
        self.body = inputDict["body"]
        self.physicalGroup = None
        self.elementComposition = []
        self.nodeComposition = []
        self.tag=None
        self.isUnterZug = isUnterZug

        E=self.body.eModule
        G=self.body.gModule
        nu=self.body.nu
        h=self.thickness
        if not(isUnterZug):
            self.Df=h**3/12*E/(1-nu**2)*np.array([[1, nu, 0],
                                                    [nu, 1, 0],
                                                    [0, 0, (1-nu)/2]])

            self.Dc = 5/6*h*np.array([[G,0],[0,G]]) #alpha = 5/6 according to Ferreira p. 162
        elif isUnterZug:
            self.Dc = 5/6*h*np.array([[G,0],[0,G]])
            e_uz = (h-t)/2
            hModified = (h**3+ 6*h*e_uz**2)**(1/3)
            self.Df =  1/12*E*np.array([[(hModified)**3, 0,    0],
                                        [0,        (t)**3, 0],
                                        [0,         0,   t**3*1/2]])

    def _plot(self, axGeometry):
        '''
        Plots a plate object in the axGeometry axis. \n
        Input: \n
        *   axGeometry: Target axis for the plot \n
        Return: \n
        *   -
        '''
        coords = self.outlineCoords
        axGeometry.plot(coords[:,0],coords[:,1], color='gray')
        if self.isUnterZug:
            axGeometry.fill_between(coords[0:2,0], coords[0:2,1], y2=coords[2:4,1],color='yellow')

class Wall:
    def __init__(self, inputDict, verticalDisplacement = False):
        '''
        A wall object contains all characteristics regarding geometry and support conditions. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **inputDict**: dictionary with following entries:
            - "outlineCoords": n x 2 numpy array with n x-y couples, representing the points defining the wall's outline.  
            - "thickness": thickness of the wall 
            - "support": Numpy array of length 3, each element is either 1 or 0 to block or leave free the relative degree of freedom. 

        '''
        self.outlineCoords = inputDict["outlineCoords"]
        self.support = inputDict["support"]
        self.thickness = inputDict["thickness"]
        self.physicalGroup = None
        self.elementComposition = []
        self.nodeComposition = None
        self.isElasticallySupported = verticalDisplacement
        if verticalDisplacement:
            try:
                self.high = inputDict["effective high"]
                self.body = inputDict["body"]
            except:
                Exception('Please define wall body and high')

    def _plot(self, axGeometry):
        '''
        Plots a wall object in the axGeometry axis. \n
        Input: \n
        *   axGeometry: Target axis for the plot \n
        Return: \n
        *   -
        '''
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
            axGeometry.plot(np.array([l1[0,0],l2[0,0]]),np.array([l1[0,1],l2[0,1]]), color='g', linewidth = linWidth)
            axGeometry.plot(np.array([l1[1,0],l2[1,0]]),np.array([l1[1,1],l2[1,1]]), color='g', linewidth = linWidth)

class Column:
    def __init__(self, inputDict, isInPlate = False):
        '''
        A column object contains all characteristics regarding geometry and support conditions. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **inputDict**: dictionary with following entries:
            - "outlineCoords": 1 x 2 numpy array with the x and y coordinates of the column. 
            - "width": Width of the column (square-shaped). 
            - "support": Numpy array of length 3, each element is either 1 or 0 to block or leave free the relative degree of freedom. 
        * **isInPlate = False**: Boolean, True if the columns is positioned inside a plate (default is False) .

        '''
        self.outlineCoords = inputDict["outlineCoords"]
        self.support = inputDict["support"]
        self.width = inputDict["width"]
        self.physicalGroup = None
        self.elementComposition = None
        self.nodeComposition = None
        self.isInPlate = isInPlate

    def _plot(self, ax):
        '''
        Plots a column object in the axGeometry axis. \n
        Input: \n
        *   axGeometry: Target axis for the plot \n
        Return: \n
        *   -
        '''
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
        ax.scatter(x,y,marker="D", color='b')

class DownStandBeam:
    def __init__(self, inputDict):
        '''
        A downStandBeam object contains all characteristics regarding geometry and material. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **inputDict**: dictionary with following entries: 
            - "outlineCoords": n x 2 numpy array with n x-y couples, representing the points defining the downstand beam's outline.  
            - "body": `Concrete` object.
            - "crossSection": `CrossSection` object. 

        '''
        self.outlineCoords = inputDict["outlineCoords"]
        self.body = inputDict["body"]
        self.crossSection = inputDict["crossSection"]
        self.physicalGroup = None
        self.elementsList = []
        self.nodeComposition = None
        self.elementsList = None
        self.Amat = None
        self.uzNodesToNodesNumeration = None
        self.coherentNodesPlate = None
        self.coherentNodesUZ = None
        self.newNodesUZ = None

    def _plot(self, axGeometry):
        '''
        Plots a downStandBeam object in the axGeometry axis. \n
        Input: \n
        *   axGeometry: Target axis for the plot \n
        Return: \n
        *   -
        '''
        coords = self.outlineCoords
        axGeometry.plot(coords[:,0],coords[:,1], color='grey', linestyle='dashed')
        for i in range(0, coords.shape[0]-1):
            p1 = coords[i]
            p2 = coords[i+1]
            d = self.crossSection.width/2
            
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
            axGeometry.plot(l1[:,0],l1[:,1], color='grey', linewidth = linWidth)
            axGeometry.plot(l2[:,0],l2[:,1], color='grey', linewidth = linWidth)
            if l1[0,0] != l1[1,0]:
                axGeometry.fill_between(l1[:,0], l2[:,1], l1[:,1],color='yellow')
            else:
                axGeometry.fill_between(np.array([l1[0,0], l2[0,0]]), np.array([l1[1,1], l2[1,1]]), np.array([l1[0,1], l2[0,1]]),color='yellow')

class Concrete:
    def __init__(self, inputDict):
        '''
        A concrete object contains all characteristics regarding the material used for plates and downstand beams.

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **inputDict**: dictionary with following entries:
            - "eModule": Elastic modulus E.  
            - "gModule": Shear modulus G. 
            - "nu": Poisson's ratio. 

        '''
        self.eModule = inputDict["eModule"]
        self.gModule = inputDict["gModule"]
        self.nu = inputDict["nu"]

class StandardConcrete(Concrete):
    def __init__(self, concreteType):
        '''
        Defines a concrete object from a repository of standard concretes.

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **concreteType**: string defining the standard concrete type. Possibilities are: 
            - "C25_30" 
            - "C30_37" 
            - "unit" 

        '''
        dic = {}
        if concreteType == 'C25_30':
            self.eModule = 32.1*1e6 #kN/m2
            self.gModule = 14.36*1e6 #kN/m2
            self.nu = 0.17
        elif concreteType == 'C30_37':
            self.eModule = 33.6*1e6 #kN/m2
            self.gModule = 14.36*1e6 #kN/m2
            self.nu = 0.17
        elif concreteType == 'unit':
            self.eModule = 10.92*1e6 #kN/m2
            self.gModule = 10.92/2.6*1e6 #kN/m2
            self.nu = 0.3
        else:
            raise ValueError('Concrete type not recognised')
        

class CrossSection:
    def __init__(self, A, Iy, Iz, b,h):
        '''
        A crossSection object contains all the geometrical information used in downstand beams. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **A**: Area of the cross section. 
        * **Iy**: Second moment of area in respect to the y-axis.
        * **Iz**: Second moment of area in respect to the z-axis.
        * **b**: Width of the structural element.
        * **h** : High of the structural element.

        '''        
        self.A = A
        self.Iy = Iy
        self.Iz = Iz
        self.width = b
        self.thickness = h

class Load:
    def __init__(self,case, magnitude):
        '''
        A load object contains all information which define a load, including magnitude, type and position. 

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **case**: String defining the type of load. Acceptable values are: 
            - "line": Line load, outline is to be additionally defined. 
            - "area": Constant load distributed on the entire structure.
            - "point": Concentrated load, position is to be additionally defined. 
        * **magnitude**: Numpy array of length 3, each element define the magnitude of the applied load for the relative degree of freedom.

        '''        
        self.magnitude = magnitude
        self.case=case
        self.outlineCoords = np.array([])
        self.physicalGroup = None
        self.elements1DList = None
        self.nodePattern = None
        self.pointLoadNode = None

