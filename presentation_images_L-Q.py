#%%

import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress
import os

ConcreteDict = {}

ConcreteDict["eModule"] = 10920 #kN/m2
ConcreteDict["gModule"] =  10920/(2*1.3) #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.0
C25_30 = Concrete(ConcreteDict)

a=4
b=1
h=0.001
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('line',np.array([-1*0,1,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
# wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 0, 1]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

# patchTestModel.addLoad(lineLoad)

import matplotlib.pyplot as plt
import pandas as pd
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *

# 1 d kragarm
nodesArray = np.array([[0,0],
                        [1, 0], 
                        [2,0],
                        [3,0], 
                        [4, 0], 
                        [4,1],
                        [3, 1], 
                        [2,1],
                        [1,1], 
                        [0, 1]])


elements = np.array([[1,2,9,10],
                    [2,3,8,9],
                    [3,4,7,8],
                    [4,5,6,7]])

BCs = np.array([[1,1,1,1],
                [10,1,1,1]])

forces = np.array([[5, -0.5*0, 0.5, 0],
                    [6, -0.5*0, 0.5, 0]])

forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
#% img settings
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\L -distorted\\"

figSize=2.0
figRatio = 4.2/1.2
#%  regular mesh with numbers
plotMesh(patchTestModel, plotStrucElements=False)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'regular_mesh_LR.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#%%
import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress
import os

ConcreteDict = {}

ConcreteDict["eModule"] = 10920 #kN/m2
ConcreteDict["gModule"] =  10920/(2*1.3) #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.0
C25_30 = Concrete(ConcreteDict)

a=4
b=1
h=0.001
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('line',np.array([-1*0,1,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
# wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 0, 1]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

# patchTestModel.addLoad(lineLoad)

import matplotlib.pyplot as plt
import pandas as pd
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *

# 1 d kragarm
nodesArray = np.array([[0,0],
                        [1.5, 0], 
                        [2,0],
                        [3.3,0], 
                        [4, 0], 
                        [4,1],
                        [3, 1], 
                        [2.5,1],
                        [0.5,1], 
                        [0, 1]])


elements = np.array([[1,2,9,10],
                    [2,3,8,9],
                    [3,4,7,8],
                    [4,5,6,7]])

BCs = np.array([[1,1,1,1],
                [10,1,1,1]])

forces = np.array([[5, -0.5*0, 0.5, 0],
                    [6, -0.5*0, 0.5, 0]])

forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
#% img settings
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\L -distorted\\"

figSize=2.0
figRatio = 4.2/1.2
#%  distorted mesh with numbers
plotMesh(patchTestModel, plotStrucElements=False)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'distorted_mesh_LR.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% distorted mesh w/o numbers
plotMesh(patchTestModel, plotStrucElements=False, plotNodes=False, plotPoints=True)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'distorted_mesh_LR_No-Nodes.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - Mx
forces = np.array([[5, 0, 0.5, 0.0],
                    [6, 0, 0.5, 0.0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces

patchTestModel.loads[0] = forcePattern

elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[0], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'LR_constantMoment_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'LR_constantMoment_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[0], internalForcePosition = 'center', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'LR_constantMoment_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - Mx
forces = np.array([[5, -0.5, 0, 0.0],
                    [6, -0.5, 0, 0.0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces

patchTestModel.loads[0] = forcePattern
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[0], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'LR_constantShear_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'LR_constantShear_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[0], internalForcePosition = 'center', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'LR_constantShear_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)


#%%#############################################################
# Quadratic element

figSize=2.5
figRatio = 4.2/1.2
#%
import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress

ConcreteDict = {}
# ConcreteDict["eModule"] = 30e6 #kN/m2
# ConcreteDict["gModule"] =  15e6 #kN/m2
# ConcreteDict["density"] = 2.5 # t/m3
# ConcreteDict["nu"] = 0.0

ConcreteDict["eModule"] = 10920 #kN/m2
ConcreteDict["gModule"] =  10920/(2*1.3) #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.0
C25_30 = Concrete(ConcreteDict)
a=4
b=1
h=0.001
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('line',np.array([-1,-1*0,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
# wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 1]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

# patchTestModel.addLoad(lineLoad)

import matplotlib.pyplot as plt
import pandas as pd
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *

# 1 d kragarm
nodesArray = np.array([[0,0],  #1
                        [1, 0], #2
                        [2,0],#3
                        [3,0], #4
                        [4, 0], #5
                        [4,1],#6
                        [3, 1], #7
                        [2,1],#8
                        [1,1], #9
                        [0, 1], #10
                        [0.5,0], #11
                        [1.5,0],#12
                        [2.5,0],#13
                        [3.5,0],#14
                        [4,0.5],#15
                        [3.5,1],#16
                        [2.5,1],#17
                        [1.5,1],#18
                        [0.5,1],#19
                        [0,0.5],#20
                        [0.5,0.5],#21
                        [1,0.5],#22
                        [1.5,0.5],#23
                        [2,0.5],#24
                        [2.5,0.5],#25
                        [3,0.5],#26
                        [3.5,0.5]])#27

elements = np.array([[1,2,9,10, 11, 22, 19, 20, 21],
                    [2,3,8,9, 12, 24, 18, 22, 23],
                    [3,4,7,8, 13, 26, 17, 24, 25],
                    [4,5,6,7, 14, 15, 16, 26, 27]])

BCs = np.array([[1,1,1,1],
                [10,1,1,1],
                [20,1,1,1]])

forces = np.array([[5, -0.0, 1/6, 0],
                    [6, -0.00, 1/6, 0],
                    [15,-0.0,2/3,0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
#%  regular mesh with numbers
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\Q -regular\\"
plotMesh(patchTestModel, plotStrucElements=False)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'regular_mesh_QR.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% regular mesh w/o numbers
plotMesh(patchTestModel, plotStrucElements=False, plotNodes=False, plotPoints=True)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'regular_mesh_QR_No-Nodes.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - Mx

elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantMoment_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantMoment_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'intPoints', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantMoment_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - Mx
forces = np.array([[5, -1/6,0, 0],
                    [6, -1/6,0, 0],
                    [15,-2/3,0,0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces

patchTestModel.loads[0] = forcePattern
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantShear_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantShear_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'intPoints', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantShear_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)


#%%#############################################################
# Quadratic element DISTORTED MESH
figSize=2.2
figRatio = 4.2/1.2
#%
import numpy as np
from createModel import *
from solveModel import *
from generateMesh import *
from scipy.stats import linregress

ConcreteDict = {}
# ConcreteDict["eModule"] = 30e6 #kN/m2
# ConcreteDict["gModule"] =  15e6 #kN/m2
# ConcreteDict["density"] = 2.5 # t/m3
# ConcreteDict["nu"] = 0.0

ConcreteDict["eModule"] = 10920 #kN/m2
ConcreteDict["gModule"] =  10920/(2*1.3) #kN/m2
ConcreteDict["density"] = 2.5 # t/m3
ConcreteDict["nu"] = 0.0
C25_30 = Concrete(ConcreteDict)
a=4
b=1
h=0.001
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('line',np.array([-1,-1*0,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
# wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 1]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

# patchTestModel.addLoad(lineLoad)

import matplotlib.pyplot as plt
import pandas as pd
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
p1x = 0
p2x = 1.5
p3x = 2
p4x= 3.3
p5x = 4
p6x = 4
p7x = 3
p8x = 2.5
p9x = 0.5
p10x = 0

p20 = 0
p22 = 0.5*(p2x+p9x)
p24 = 0.5*(p3x+p8x)
p26 = 0.5*(p4x+p7x)
p15 = 4

# 1 d kragarm
nodesArray = np.array([[p1x,0],  #1
                        [p2x, 0], #2
                        [p3x,0],#3
                        [p4x,0], #4
                        [p5x, 0], #5
                        [p6x,1],#6
                        [p7x, 1], #7
                        [p8x,1],#8
                        [p9x,1], #9
                        [p10x, 1], #10
                        [0.5*(p1x+p2x),0], #11
                        [0.5*(p2x+p3x),0],#12
                        [0.5*(p3x+p4x),0],#13
                        [0.5*(p4x+p5x),0],#14
                        [0.5*(p5x+p6x),0.5],#15
                        [0.5*(p7x+p6x),1],#16
                        [0.5*(p8x+p7x),1],#17
                        [0.5*(p9x+p8x),1],#18
                        [0.5*(p10x+p9x),1],#19
                        [p20,0.5],#20
                        [0.5*(p20+p22),0.5],#21
                        [p22,0.5],#22
                        [0.5*(p24+p22),0.5],#23
                        [p24,0.5],#24
                        [0.5*(p24+p26),0.5],#25
                        [p26,0.5],#26
                        [0.5*(p26+p15),0.5]])#27

elements = np.array([[1,2,9,10, 11, 22, 19, 20, 21],
                    [2,3,8,9, 12, 24, 18, 22, 23],
                    [3,4,7,8, 13, 26, 17, 24, 25],
                    [4,5,6,7, 14, 15, 16, 26, 27]])

BCs = np.array([[1,1,1,1],
                [10,1,1,1],
                [20,1,1,1]])

forces = np.array([[5, -0.0, 1/6, 0],
                    [6, -0.00, 1/6, 0],
                    [15,-0.0,2/3,0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
#%  regular mesh with numbers
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\Q -distorted\\"
plotMesh(patchTestModel, plotStrucElements=False)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'distorted_mesh_QR.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% regular mesh w/o numbers
plotMesh(patchTestModel, plotStrucElements=False, plotNodes=False, plotPoints=True)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'distorted_mesh_QR_No-Nodes.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - Mx

elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantMoment_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantMoment_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'intPoints', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantMoment_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - Mx
forces = np.array([[5, -1/6,0, 0],
                    [6, -1/6,0, 0],
                    [15,-2/3,0,0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces

patchTestModel.loads[0] = forcePattern
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantShear_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantShear_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[2], internalForcePosition = 'intPoints', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(-0.1, 4.1)
currentAx.set_ylim(-0.1, 1.1)
plt.axis('off')
imgName = 'QR_constantShear_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)


#%% plot 1-d plot

# x = patchTestModel.results.internalForcesPositions[:,0]
# y = patchTestModel.results.internalForcesPositions[:,1]
# z = patchTestModel.results.shearForces[:,0]
# z = np.abs(z)

# ax = plt.subplot()
# ax.plot(x,z, color = 'k')
# ax.grid()
# plt.yscale('log')
# ax.set_xlim(0.49, 0.51)
# # ax.set_ylim(0.6,1.4)
# plt.xlabel('x [-]')
# plt.ylabel('Vx [-]')

# plt.show()
