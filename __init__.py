'''
*Made with ❤️ by Tobia Diggelmann (ETH Zurich), June 2021*

Provides

* Methods to easily model structural components for plates via the Gmsh API (see https://gmsh.info/).
* Mesh generation via Gmsh API.
* FEM calculation of displacements, rotations, bending moments, shear forces and required bending moment resistance according to SIA262.
* Results visualization with isolines and 3D plots.
* Analytical solutions according to Kirchhoff plate theory.



The present documentation provides the user with an overview of the package and its functioning. For a more *hands-on* approach to getting started, refer to the tutorials in appendix [A. Tutorials](#A. Tutorials). The tutorial are also already available in the package distribution in the *tutorials* folder.


1. Introduction
====================

The purpose of the present package, developed during a master thesis at the ETH Zurich, is to provide a simple set of tools to model and compute plate models using Finite Elements. The package is based on the `Gmsh API`, used to create geometrical models and to discretize the system into elements. The package has been developed for users with knowledge about civil engineering, while every aspect regarding finite elements has been automatized as much as possible. This allows for a basic FE-modelling without direct knowledge of FEM. However, all kind of parameters regarding meshing and computation can be customized. 

2. How to get platepy
========================

Platepy is available on https://pypi.org/project/platepy/ , and can therefore be obtained through the terminal with the command (without the $ sign):

`$ pip install platepy`

Alternatively, if problems with accessing the package through the global path are encountered, the packages can be accessed on the GitHub repository: https://github.com/tobiadig/platepy . All downloaded modules should be in a folder named "platepy" and located in the working directory.

When the platepy package has been downloaded (via pip, or via GitHub and located in the working directory) the functions in the package can be accessed in a Python script with the command `from platepy import *`.

Following packages are required for the correct functioning of platepy:

* numpy
* scipy
* matplotlib
* pandas
* Gmsh
* tqdm

All these packages can be installed through pip with the command (without the $ sign)

`$ pip install numpy, scipy, matplotlib, pandas, gmsh, tqdm`

NOTE: On windows, in order to access the packages, the system has to be restarted after installation (to update the global path to the packages).


3. How to use platepy
=============================

Once the functions of platepy are accessible, the workflow is as follows:

1. Create a model with structural components (slabs, walls, columns etc.) and loads.
2. Generate a mesh.
3. Solve the system using FEM.
4. Display the results.

Each steps is presented in the following sections, examples are provided in [A. Tutorials](#A. Tutorials).


3.1. Create a model
----------------------------

** Create model with the Gmsh API **


The package can be used to model a plate static system. The plate model is represented by a `PlateModel` class. A model can be inizialized as follows:


    exampleModel = PlateModel()

The user has then to define all the components which compose the static system together with their properties, and the loads acting on the system. The components have then to be added to the `PlateModel` object. The available structural components are:

* plates
* walls
* columns
* downstand beams

The procedure to add a structural component is:

1. Create a dictionary. Create the required entries, which are documented in the submodule [platepy.createModel](#create).
2. Instantiate the structural component's class with the created dictionary.
3. Add the object to the `PlateModel` with the specific method of the component.

The general structure is therefore as follows: ::

    inputDic = {}
    inputDic['entry1'] = value1
    inputDic['entry2'] = value2
    <component>Object = <Component>Class(inputDic)
    exampleModel.add<Component>(<component>Object)

The loads are defined similarly with the `Load` class.

NOTE: It is advisable to check if the geometry has been correctly defined before moving to the next steps of the analysis. To inspect the geometry call the `plotInputGeometry function` defined in [displayModel](#displayModel).

NOTE: Plate outline has to be a closed polygon, points defining walls and downstand beams **have** to lie inside the plate or to coincide with corners of the plate. ** Walls and downstand beams starting or ending on a plate's edge will raise an error**. In this cases, add points to the plate's outline coordinates accordingly.

** Create model manually **

It is also possible for the user to manually define the nodes positions, the element connectivity, blocked DOFs and nodal forces. See the tutorial [Manually generate mesh](#Manually generate mesh)


3.2. Generate mesh
----------------------------

** Generate mesh automatically with Gmsh**

The mesh of the created model can be automatically generated with the Gmsh API. Call the `generateMesh` function defined in [generateMesh](#generateMesh). The most critical parameter is the mesh size. By default the meshSize is 0.6m, which works fine for plates with a size in the range of 8-12m. It is therefore important to inspect the mesh before the next steps of the analysis and to check if the mesh size is too coarse (risk of inaccurate results) or too fine (risk of the calculation taking too long). On an normal laptop, a model with up to 4'000 elements is generally sufficiently fast (ca. 400 elements/s on standard laptops).

The mesh can be inspected by setting the parameter `showGmshMesh` to True.

NOTE: The automatic generation of the mesh is generally a critical step, since not always the Gmsh library is able to generate the elements and the errors raised are quite cryptic. If a mesh cannot be generated, try following steps:

1. Check if the geometry has been correctly defined. The geometry can also be inspected with by setting the `showGmshGeometryBeforeMeshing` option of the generateMesh function to True. 

    a. Plate outline has to be a closed polygon.

    b. Points defining walls and downstand beams **have** to lie inside the plate or to coincide with corners of the plate. ** Walls and downstand beams starting or ending on a plate's edge will raise an error**. In this cases, add points to the plate's outline coordinates accordingly.

    c. Columns have either to lie inside the plate (with the `isInPlate` option set to True), or to lie on a plate's corner. ** Columns lying on a plate's edge will raise an error**.

2. By very sharp angles or very close structural elements, it is possible that Gmsh fails to generate a mesh. Try a finer mesh or simplify the model.
3. If the meshSize is too coarse compared to the model's size, an error will rise. 


<a name="Generate mesh manually"></a>**Generate mesh manually**

It is also possible for the user to manually define the nodes positions, the element connectivity, blocked DOFs and nodal forces. This is particularly useful if a simple, yet precise mesh is required (for example for a patch test). An example is shown in [tutorial 3](#t3).


3.3. Solve
----------------------------
Once the mesh has been successfully generated, the system can be solved by calling the `solveModel` function of the submodule [solveModel](#solveModel). The results have to be scaled according to the dimensions used in the input value. It is advised to use dimensions consistently in order to get the right results (for example meters and kN). To get displacements in mm, moments in kNm and shear forces in kN, the result scale tuple will be (0.001,1,1). The user can also choose the position where the internal forces should be evaluated with the `internalForcePosition` parameter. By default internal forces are evaluated at the centre of each element.

After the calculation, the user can also compute the values along an arbitrary line using the `computeBeamComponents` function. Single or multiple specific points can be evaluated with the `evaluateAtPoints` function.

3.4. Display results
----------------------------

The submodule [displayModel](#displayModel) contains all functions which allow to display the results. Platepy allows to:

* Display the input geometry with the `plotInputGeometry` function.
* Display the mesh with `plotMesh`. This is an alternative to the build-in Gmsh GUI, activated by setting the `showGmshMesh` parameter of the `generateMesh` function equal to True. plotMesh also allow to save the image.
* Display displacements and internal forces with the `plotResults` function. Refer to the submodule [displayModel](#displayModel) for details.
* Display results along an arbitrary line, computed previously with the `computeBeamComponents` function.

The created images can be manually saved through the plt.show() interface, or automatically saved in SVG format by setting the `saveToSVG` parameter equal to True.



<a name="A. Tutorials"></a> A. Tutorials
=============================================


<a name="t1"></a>A.1. Creating a simple slab geometry and computing the FEM solution 
------------------------------------------------------------------------

    # ------------------------------------------------------------------------------
    #
    # platepy tutorial 1
    #
    # Creating a simple slab geometry and computing the FEM solution
    #
    # ------------------------------------------------------------------------------

    # The program is entirely defined in the `platepy.py' module (which contains the
    # full documentation of all the function). The easier way is to directly import
    # all function with the * command:
    from platepy import *
    import numpy as np

    # The model will now be initialized
    tutorialModel = PlateModel()

    # To define the structural components of the model a dictionary has to be created.
    # The required entries of the dictionary are described in the documentation of the 
    # relative class.

    # let's define the concrete by defining the dictionary and the required entries:
    # (attention to the units! by default all lengths are considered to be in meters,
    # forces in kN).
    concreteDict = {}
    concreteDict["eModule"] = 32.1*1e6 #kN/m2
    concreteDict["gModule"] =  14.36*1e6 #kN/m2
    concreteDict["nu"] = 0.17
    C25_30 = Concrete(concreteDict)

    # alternatively, a standard concrete type can be used:
    C30_37 = StandardConcrete("C30_37")


    # Let's now define the structural components with the same procedure.
    # First a square plate with a 10m side:
    plateDict = {}
    plateDict["outlineCoords"]=np.array([[0,0], [10,0],[10,10], [0,10], [0, 0]])
    plateDict["thickness"] = 0.3
    plateDict["body"]=C30_37
    plate = Plate(plateDict)

    # and let's add the plate to the model:
    tutorialModel.addPlate(plate)

    # now we define two walls on the left and right sides of the plate:
    # the wall dictionary can be re-used by only changing the outline.
    wallDict = {}
    wallDict["outlineCoords"] = np.array([[0,0],[0,10]])
    wallDict["thickness"] = 0.5 # m
    # The entries in the support array define is the relative DOF is free (0)
    # or blocked (1). The first entry is the vertical displacement, the second
    # is the rotation on the axis in the direction of the wall, the third entry
    # is the rotation on the axis perpendicular to the wall. For simply supported
    # hard boundary condition:
    wallDict["support"] = np.array([1, 0, 1])
    wall1 = Wall(wallDict)

    # the second wall:
    wallDict["outlineCoords"] = np.array([[10,0],[10,10]])
    wall2 = Wall(wallDict)

    # and add to the model:
    tutorialModel.addWall(wall1)
    tutorialModel.addWall(wall2)

    # add a column in the centre:
    columnDict = {}
    columnDict["outlineCoords"] = np.array([[5,5]])
    columnDict["support"] = np.array([1, 0, 0])
    columnDict["width"] = 0.5
    col1 = Column(columnDict, isInPlate=True) # if the plate is inside the plate, set "isInPlate = True"
    tutorialModel.addColumn(col1)

    # Lastly, a constant load distributed over the plate:
    load = Load('area', np.array([-10,0,0]))
    tutorialModel.addLoad(load)

    # The model is ready! Let's check the geometry and then 
    # generate the mesh and see if the fineness is
    # appropriate or has the be adjusted
    plotInputGeometry(tutorialModel)
    plt.show()
    generateMesh(tutorialModel, showGmshMesh=True)

    # The mesh has been generated correctly, let's solve the model.
    # The scale for the displacements has been changed to 1e3, in
    # order to transform m to mm
    solveModel(tutorialModel, resultsScales=(1e3,1,1))
    #
    # Now the results can be plotted, for example
    # Vertical displacements, bending moments and shear in x-direction
    plotResults(tutorialModel, ['vDisp', 'Mx', 'Vx'])
    plt.show()



<a name="t2"></a>A.2. Downstand beams and plotting results over a line ('Schnitte')
----------------------------------------------------------------------

    # ------------------------------------------------------------------------------
    #
    # platepy tutorial 2
    #
    # Downstand beams and plotting results over a line ('Schnitte')
    #
    # ------------------------------------------------------------------------------

    # Let's create a model with a square plate supported by two walls on opposite
    # Edges, but this time with a downstand beam spanning between the walls.
    # Important is to add coordinates in the outline of walls and plates for
    #  where the downstand beam will come!

    from platepy import *
    import numpy as np
    tutorialModel = PlateModel()

    C30_37 = StandardConcrete("C25_30")

    plateDict = {}
    # points [10,5] and [0,5] have been added
    plateDict["outlineCoords"]=np.array([[0,0], [10,0],[10,5],[10,10], [0,10],[0,5], [0, 0]])
    plateDict["thickness"] = 0.3
    plateDict["body"]=C30_37
    plate = Plate(plateDict)
    tutorialModel.addPlate(plate)

    wallDict = {}
    wallDict["outlineCoords"] = np.array([[0,0],[0,5],[0,10]])
    wallDict["thickness"] = 0.5 # m
    wallDict["support"] = np.array([1, 0, 1])
    wall1 = Wall(wallDict)
    wallDict["outlineCoords"] = np.array([[10,0],[10,5],[10,10]])
    wall2 = Wall(wallDict)
    tutorialModel.addWall(wall1)
    tutorialModel.addWall(wall2)

    # Let's now define a downstand beam spanning between the two walls
    dsbDic = {}
    dsbDic['outlineCoords'] = np.array([[0,5],[10,5]])
    dsbDic['body'] = C30_37
    # The cross section object stores information about the dsb:
    # Area, Iy, Iz, width and high (m and mm^4)
    dsbDic['crossSection'] = CrossSection(0.3*0.3, 0.3*4/12, 0, 0.3, 0.3)
    dsb = DownStandBeam(dsbDic)
    tutorialModel.addDownStandBeam(dsb)

    # Add the load:
    load = Load('area', np.array([-10,0,0]))
    tutorialModel.addLoad(load)

    # As in tutorial 1, let's check the geometry, generate the mesh and solve
    plotInputGeometry(tutorialModel)
    plt.show()
    generateMesh(tutorialModel)
    solveModel(tutorialModel, resultsScales=(1e3,1,1))

    # Instead of plotting the results with isolines, let's inspect a vertical cut
    # for Vertical displacements, bending moments and shear in y-direction

    # The cut is defined by starting and end coordinates, the number of points to be evaluated and 
    # the scale of the result
    computeBeamComponents(tutorialModel, (5,0), (5,10), 150, resultsScales=(1e3,1,1))

    # the results can be now plotted with the function plotBeamComponent
    # the suffix '_line' at the end of the desired plot has to be added
    plotBeamComponent(tutorialModel, ['vDisp_line', 'My_line', 'Vy_line'])
    plt.show()

<a name="t3"></a>A.3. Manually generate mesh
------------------------------------------------------------------------

    # ------------------------------------------------------------------------------
    #
    # platepy tutorial 3
    #
    # Manually create a simple mesh
    #
    # ------------------------------------------------------------------------------


    # The present tutorial shows how to manually define nodes, element, boundary conditions
    # and nodal loads. This can be useful when a simple, yet precise defined mesh is required
    # For example for a patch test.

    # Let's prepare the model and the material
    import numpy as np
    from platepy import *

    tutorialModel = PlateModel()
    C25_30 = StandardConcrete("C25_30")

    # the plate has still to be assigned to the model in order to specify the material properties
    # of the elements
    plateDict = {}
    plateDict["outlineCoords"]=np.array([[0,0], [10,0]]) # the outline is irrelevant, since 
                                                            # the shape will be manually defined later
    plateDict["thickness"] = 0.1
    plateDict["body"]=C25_30
    plate = Plate(plateDict)
    tutorialModel.addPlate(plate)

    # now, instead of defining loads, structural components and generate the mesh, all steps are
    # performed manually.

    # Firstly, the array of nodes coordinates has to be defined.
    nodesArray = np.array([[0,0],    # Node 1
                            [10, 0], # Node 2
                            [10,10], # Node 3
                            [0,10],  # Node 4
                            [2, 2],  # Node 5
                            [8,3],   # Node 6
                            [8,7],   # Node 7
                            [4,7]])  # Node 8

    #Then the elements connectivity
    elements = np.array([[1,2,6,5], # Element 1
                        [2,3,7,6],  # Element 2
                        [7,3,4,8],  # Element 3
                        [1,5,8,4],  # Element 4
                        [5,6,7,8]]) # Element 5

    #Then, the boundary condition. The first column is the node of reference, the other columns
    # define if their relative DOF is blocked or free (vDisp/xRot/yRot). 
    # Let's clamp the left-hand side and block all rotations.
    BCs = np.array([[1, 1, 1, 1],
                    [4, 1, 1, 1],
                    [2, 0, 1, 1],
                    [3, 0, 1, 1],
                    [5, 0, 1, 1],
                    [6, 0, 1, 1],
                    [7, 0, 1, 1],
                    [8, 0, 1, 1]])

    # The nodal forces have to be defined. Let's put two vertical loads on the right-hand side.
    forces = np.array([[2, -5, 0, 0],
                        [3, -5, 0, 0]])

    # A load object is created, then the "nodePattern" attribute is assigned with 
    # the above defined array of nodal forces
    forcePattern = Load('nodes', np.array([0,0,0])) # the magnitude is not used
    forcePattern.nodePattern=forces
    tutorialModel.addLoad(forcePattern)

    # nodes, elements and boundary conditions are assigned to the model with the setMesh function
    # differently from the generateMesh function, there is no default element type
    setMesh(tutorialModel, nodesArray, elements, BCs,  elementDefinition='MITC-4-N')

    solveModel(tutorialModel, resultsScales=(1e6,1,1))

    # Let's now look at the displacement in the individual nodes with the "text+mesh" plot type
    plotResults(tutorialModel,plotType='text+mesh',valuesToPlotList=['vDisp'])
    plt.show()

<a name="t4"></a>A.4. Compute analytical solutions
------------------------------------------------------------------------

    # ------------------------------------------------------------------------------
    #
    # platepy tutorial 4
    #
    # Compute analytical solution according to Timoshenko's formulas
    #
    # ------------------------------------------------------------------------------

    # This tutorial shows how to use the analyticPlateSolutions module
    # We will access the displacement at the location (x=-2, y=1) (relative 
    # to the plate's center) of a simply supported rectangular plate with shape
    # a=10, b=12

    from platepy import *

    # First, let's define the input options as attributes of the POpts class
    # The possible entries are described in the documentation
    pOpts = POpts()
    pOpts.shape="rectangular"
    pOpts.depth = "thin"
    pOpts.support = "simplySupported"
    pOpts.geometry = (10,12)  # width a and high b
    pOpts.material = Material(10920, 0.3, 0.1) #E, nu and thickness

    # Here the Load options are defined, as attributes of the LOpts class
    # The possible entries are described in the documentation
    lOpts = LOpts()
    lOpts.type = "distributed"
    lOpts.magnitude = 1

    # The number of term for the Taylor expantion has to be defined.
    # 40 terms is usually more than enough
    sOpts = SOpts()
    sOpts.nTerms = 40

    #  inPos is a nPointsx2 numpy array containing the x-y coordinates of
    # the points we want to evaluate. In our case we only need on point.
    xPos = -2
    yPos = 1
    inPos=np.array([[xPos,yPos]])  # x-y coordinates

    quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)

    # values is a (nPoints, nOutputs) array of the quantities values evaluated at outPos.
    # the output values correspond to the True entries of the array quantities.
    # the entries of quantities correspond to: (vDisp, Mx, My, Mxy, Vx, Vy).
    # In this case, all values are available apart from the twisting moment Mxy, therefore
    # values has shape (1,5).

    vDisp_mm = values[0,0]*1000

    print("The displacement at (", xPos,",",yPos,") corresponds to ", vDisp_mm, ' mm')


'''

from .analyticPlateSolutions import *
from .createModel import *
from .generateMesh import *
from .solveModel import *
from .displayModel import *

