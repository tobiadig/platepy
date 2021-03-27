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
    supportCondition = Support(np.array([1, 1, 0]))

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
    wallDict["outlineCoords"] = np.array([[0,0], [a,0]])
    # wallDict["outlineCoords"] = np.array([[0,0], [a,0]])
    
    wallDict["high"] = 3 # m
    wallDict["body"] = C25_30
    wallDict["support"] = supportCondition
    wallDict["thickness"] = 0.05 # m
    wall1 = Wall(wallDict)

    wallDict["outlineCoords"] = np.array([[0,b], [a,b]])
    wall2 = Wall(wallDict)

    # wallDict["outlineCoords"] = np.array([[0,b], [a,b], [a,2*b]])
    # wall1 = Wall(wallDict)


    columnDict = {}
    columnDict["outlineCoords"] = np.array([[0.0*a,b*0.5]])
    columnDict["high"] = 3
    columnDict["body"] = C25_30
    columnDict["support"] = Support(np.array([1, 0, 0]))
    columnDict["crossSection"] = None
    columnDict["width"] = 0.05
    col1 = Column(columnDict)

    firstModel = PlateModel("plateModel1")
    firstModel.addPlate(plate1)
    firstModel.addWall(wall1)
    firstModel.addWall(wall2)

    firstModel.addLoad(distributedLoad)
    firstModel.addColumn(col1)

    return firstModel


