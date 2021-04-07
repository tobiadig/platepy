''' Module Information
-----------------------------------------------------------
Purpose of module: example of geometry definition of a simply supported rectangular slab
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''

# basic modules
import numpy as np
from createModel import *

def importModel():
    ConcreteDict = {}
    ConcreteDict["eModule"] = 10920 #kN/m2
    ConcreteDict["gModule"] =  10920/(2*(1+0.3))#kN/m2
    ConcreteDict["density"] = 2.5 # t/m3
    ConcreteDict["nu"] = 0.3
    C25_30 = Concrete(ConcreteDict)

    distributedLoad = Load(np.array([-1, 0, 0]))

    a=1
    b=1
    h=0.1
    plateDict = {}
    plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])

    # plateDict["outlineCoords"]=np.array([[0,0], [1,0], [1,1], [2,1], [2,2], [0, 2], [0,0]])

    plateDict["thickness"] = h
    plateDict["surfaceLevel"] = 0
    plateDict["body"]=C25_30
    plateDict["stiffnessFactor"] = 1
    plate1 = Plate(plateDict)

    wallDict = {}
    # wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
    wallDict["outlineCoords"] = np.array([[0.25,0.25], [0.75,0.25], [0.75,0.75], [0.25,0.75], [0.25,0.25]])
    # wallDict["outlineCoords"] = np.array([[0,0], [a,0]])
    # wallDict["outlineCoords"] = np.array([[0,0], [a,0]])
    
    wallDict["high"] = 3 # m
    wallDict["body"] = C25_30
    wallDict["support"] = Support(np.array([1, 1, 1]))
    wallDict["thickness"] = 0.05 # m
    wall1 = Wall(wallDict)

    wallDict["outlineCoords"] = np.array([[0,b], [a,b]])
    wall2 = Wall(wallDict)

    # wallDict["outlineCoords"] = np.array([[0,b], [a,b], [a,2*b]])
    # wall1 = Wall(wallDict)


    columnDict = {}
    columnDict["outlineCoords"] = np.array([[0.5*a,b*0.55]])
    columnDict["high"] = 3
    columnDict["body"] = C25_30
    columnDict["support"] = Support(np.array([1, 0, 0]))
    columnDict["crossSection"] = None
    columnDict["width"] = 0.05

    # columnDict["outlineCoords"] = np.array([[0.15*a,b*0.35]])
    col1 = Column(columnDict,isInPlate = True)

    columnDict["outlineCoords"] = np.array([[a*0.5,b]])
    col2 = Column(columnDict)

    columnDict["outlineCoords"] = np.array([[a,b*0.5]])
    col3 = Column(columnDict)

    columnDict["outlineCoords"] = np.array([[0.5*a,0]])
    col4 = Column(columnDict)
    columnDict["outlineCoords"] = np.array([[0.5*a,0.5]])
    col5 = Column(columnDict, isInPlate=True)

    

    firstModel = PlateModel("plateModel1")
    
    firstModel.addPlate(plate1)
    firstModel.addWall(wall1)
    # firstModel.addWall(wall2)
    # firstModel.addColumn(col1)
    # firstModel.addColumn(col2)
    # firstModel.addColumn(col3)
    # firstModel.addColumn(col4)
    # firstModel.addColumn(col5)
    firstModel.addLoad(distributedLoad)
    

    return firstModel


