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
    ConcreteDict["eModule"] = 32.1*1e6 #kN/m2
    ConcreteDict["gModule"] =  14.36*1e6 #kN/m2
    ConcreteDict["density"] = 2.5 # t/m3
    ConcreteDict["nu"] = 0.17
    C25_30 = Concrete(ConcreteDict)

    distributedLoad = Load('area',np.array([-1, 0, 0]))

    a=10
    b=10
    h=0.1
    plateDict = {}
    plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])
    plateDict["thickness"] = h
    plateDict["surfaceLevel"] = 0
    plateDict["body"]=C25_30
    plateDict["stiffnessFactor"] = 1
    plate1 = Plate(plateDict)

    wallDict = {}
    wallDict["outlineCoords"] = np.array([[0,0], [a,0],[a,b],[0,b],[0,0]])
    wallDict["high"] = 3 # m
    wallDict["body"] = C25_30
    wallDict["support"] = Support(np.array([1, 1, 1]))
    wallDict["thickness"] = 0.5 # m
    wall1 = Wall(wallDict)

    # columnDict = {}
    # columnDict["outlineCoords"] = np.array([[0,b*2]])
    # columnDict["high"] = 3
    # columnDict["body"] = C25_30
    # columnDict["support"] = Support(np.array([1, 0, 0]))
    # columnDict["crossSection"] = None
    # columnDict["width"] = 0.5

    # col1 = Column(columnDict,isInPlate = False)

    firstModel = PlateModel("plateModel1")
    firstModel.addPlate(plate1)
    firstModel.addWall(wall1)

    firstModel.addLoad(distributedLoad)
    
    return firstModel