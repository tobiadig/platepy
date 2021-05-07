#%%  MITC 4 constant moment
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


nodesArray = np.array([[0,0],
                        [10, 0], 
                        [10,10],
                        [0,10], 
                        [2, 2], 
                        [8,3], 
                        [8,7],
                        [4,7]])

elements = np.array([[1,2,6,5],
                    [2,3,7,6],
                    [7,3,4,8],
                    [1,5,8,4],
                    [5,6,7,8]])



BCs = np.array([[1, 1, 0, 1],
                [4, 1, 0, 1],
                [2, 0, 0, 0],
                [3, 0, 0, 0]])

forces = np.array([[2, 0, 0, 5],
                    [3, 0, 0, 5]])

forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
#% img settings
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC4\\"

figSize=4.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1
#%  distorted mesh with numbers
plotMesh(patchTestModel, plotStrucElements=False)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'mesh_numbered.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% distorted mesh w/o numbers
plotMesh(patchTestModel, plotStrucElements=False, plotNodes=False, plotPoints=True)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'mesh.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - Mx
# forces = np.array([[5, 0, 0.5, 0.0],
#                     [6, 0, 0.5, 0.0]])
# forcePattern = Load('nodes', np.array([0,0,0]))
# forcePattern.nodePattern=forces

# patchTestModel.loads[0] = forcePattern

elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[1], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantMoment_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant moment - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantMoment_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[1], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantMoment_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#%% MITC4 constant shear########################################
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC4\\"
figSize=4.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1

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


nodesArray = np.array([[0,0],
                        [10, 0], 
                        [10,10],
                        [0,10], 
                        [2, 2], 
                        [8,3], 
                        [8,7],
                        [4,7]])

elements = np.array([[1,2,6,5],
                    [2,3,7,6],
                    [7,3,4,8],
                    [1,5,8,4],
                    [5,6,7,8]])



BCs = np.array([[1, 1, 1, 1],
                [4, 1, 1, 1],
                [2, 0, 1, 1],
                [3, 0, 1, 1],
                [5, 0, 1, 1],
                [6, 0, 1, 1],
                [7, 0, 1, 1],
                [8, 0, 1, 1]])

forces = np.array([[2, -5, 0, 0],
                    [3, -5, 0, 0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e0, elementDefinition=elemTypes[1], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantShear_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantShear_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[1], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])

currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantShear_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#%% MITC4 constant torsion########################################
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC4\\"
figSize=4.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1

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


nodesArray = np.array([[0,0],
                        [10, 0], 
                        [10,10],
                        [0,10], 
                        [2, 2], 
                        [8,3], 
                        [8,7],
                        [4,7]])

elements = np.array([[1,2,6,5],
                    [2,3,7,6],
                    [7,3,4,8],
                    [1,5,8,4],
                    [5,6,7,8]])



BCs = np.array([[1, 1, 0, 0],
                [4, 1, 0, 0],
                [2, 1, 0, 0]])

forces = np.array([[3, -1, 0, 0]])
forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[1], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantT_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

# torsion
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['xy'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantT_Mxy.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantT_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[1], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantT_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)



#%% MITC9 constant shear########################################
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC9\\"
figSize=6.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1

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


nodesArray = np.array([[0,0], #1
                        [10, 0], #2
                        [10,10], #3
                        [0,10], #4
                        [2, 2], #5
                        [8,3], #6
                        [8,7], #7
                        [4,7], #8
                        [5,0], #9
                        [1,1], #10
                        [5,1.25], #11
                        [9,1.5], #12
                        [5,2.5], #13
                        [0,5], #14
                        [1.5,4.75], #15
                        [3,4.5], #16
                        [5.5,4.75], #17
                        [8,5], #18
                        [9,5], #19
                        [10,5], #20
                        [2,8.5], #21
                        [6,7], #22
                        [9,8.5], #23
                        [5.5,8.5], #24
                        [5,10]]) #25

elements = np.array([[1,2,6,5,9,12,13,10,11],
                    [6,2,3,7,12,20,23,18,19],
                    [8,7,3,4,22,23,25,21,24],
                    [1,5,8,4,10,16,21,14,15],
                    [5,6,7,8,13,16,22,18,17]])


BCs = np.array([[1, 1, 1, 1],
                [4, 1, 1, 1],
                [14, 1, 1, 1],
                [2,0,1,1],
                [3,0,1,1],
                [5,0,1,1],
                [6,0,1,1],
                [7,0,1,1],
                [8,0,1,1]])

fTot = 10

forces = np.array([[2, -fTot/6, 0, 0],
                    [3, -fTot/6, 0, 0],
                    [20,-fTot*2/3, 0, 0]])

forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)
plotMesh(patchTestModel, plotStrucElements=False)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'mesh_numbered.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% distorted mesh w/o numbers
plotMesh(patchTestModel, plotStrucElements=False, plotNodes=False, plotPoints=True)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'mesh.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantV_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantV_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantV_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#%% MITC9 constant moment
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC9\\"
figSize=6.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1

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

nodesArray = np.array([[0,0], #1
                        [10, 0], #2
                        [10,10], #3
                        [0,10], #4
                        [2, 2], #5
                        [8,3], #6
                        [8,7], #7
                        [4,7], #8
                        [5,0], #9
                        [1,1], #10
                        [5,1.25], #11
                        [9,1.5], #12
                        [5,2.5], #13
                        [0,5], #14
                        [1.5,4.75], #15
                        [3,4.5], #16
                        [5.5,4.75], #17
                        [8,5], #18
                        [9,5], #19
                        [10,5], #20
                        [2,8.5], #21
                        [6,7], #22
                        [9,8.5], #23
                        [5.5,8.5], #24
                        [5,10]]) #25

elements = np.array([[1,2,6,5,9,12,13,10,11],
                    [6,2,3,7,12,20,23,18,19],
                    [8,7,3,4,22,23,25,21,24],
                    [1,5,8,4,10,16,21,14,15],
                    [5,6,7,8,13,16,22,18,17]])


BCs = np.array([[1, 1, 1, 0],
                [4, 1, 1, 0],
                [14, 1, 1, 0]])

fTot = 10

forces = np.array([[2,0,0, -fTot/6],
                    [3, 0,0,-fTot/6],
                    [20,0,0,-fTot*2/3]])

forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)

# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-14, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantM_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantM_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantM_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)


#%% MITC9 constant torsion
myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC9\\"
figSize=6.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1

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


nodesArray = np.array([[0,0], #1
                        [10, 0], #2
                        [10,10], #3
                        [0,10], #4
                        [2, 2], #5
                        [8,3], #6
                        [8,7], #7
                        [4,7], #8
                        [5,0], #9
                        [1,1], #10
                        [5,1.25], #11
                        [9,1.5], #12
                        [5,2.5], #13
                        [0,5], #14
                        [1.5,4.75], #15
                        [3,4.5], #16
                        [5.5,4.75], #17
                        [8,5], #18
                        [9,5], #19
                        [10,5], #20
                        [2,8.5], #21
                        [6,7], #22
                        [9,8.5], #23
                        [5.5,8.5], #24
                        [5,10]]) #25

elements = np.array([[1,2,6,5,9,12,13,10,11],
                    [6,2,3,7,12,20,23,18,19],
                    [8,7,3,4,22,23,25,21,24],
                    [1,5,8,4,10,16,21,14,15],
                    [5,6,7,8,13,16,22,18,17]])


BCs = np.array([[1, 1, 0, 0],
                [4, 1, 0, 0],
                [2, 1, 0, 0],
                [14, 1, 0, 0],
                [9, 1, 0, 0]])

fTot = 10

forces = np.array([[3, -10,0,0]])

forcePattern = Load('nodes', np.array([0,0,0]))
forcePattern.nodePattern=forces
patchTestModel.addLoad(forcePattern)

setMesh(patchTestModel, nodesArray, elements, BCs)

# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)

plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['xy'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantT_Mxy.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#% constant shear - displacement
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=True, bendingMomentsToPlot=[] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantT_w.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% constant moment - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantT_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#MITC9 regular mesh constant shear

myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC9 -regular\\"
figSize=6.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1

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

a=10
b=10
h=0.001
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('line',np.array([-1,0,0]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
# wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 0]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

patchTestModel.addLoad(lineLoad)

import matplotlib.pyplot as plt
import pandas as pd
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
generateMesh(patchTestModel, showGmshMesh=False, elementType='QUAD', nEdgeNodes=3, order='quadratic')

plotMesh(patchTestModel, plotStrucElements=False)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'mesh_numbered.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
#% distorted mesh w/o numbers
plotMesh(patchTestModel, plotStrucElements=False, plotNodes=False, plotPoints=True)
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'mesh.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)
# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']

#% constant shear - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=[] ,shearForcesToPlot=['x'])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantV_Vx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#MITC9 regular mesh

myPath = "C:\\Users\\Diggelmann\\Desktop\\MSc thesis\\presentazioni\\02_midterm\\final_images\\MITC9 -regular\\"
figSize=6.0
figRatio = 1
xlim1 = -0.1
xlim2 = 10.1
ylim1 = -0.1
ylim2 = 10.1

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

a=10
b=10
h=0.001
plateDict = {}
plateDict["outlineCoords"]=np.array([[0,0], [a,0], [a,b], [0,b],[0,0]])

plateDict["thickness"] = h
plateDict["surfaceLevel"] = 0
plateDict["body"]=C25_30
plateDict["stiffnessFactor"] = 1
plate1 = Plate(plateDict)

lineLoad = Load('line',np.array([0,0,-1]))
lineLoad.outlineCoords=np.array([[a,0], [a,b]])

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0], [0,b]])
# wallDict["outlineCoords"] = np.array([[0,0], [a,0], [a,b], [0,b], [0,0]])
wallDict["high"] = 3 # m
wallDict["body"] = C25_30
wallDict["support"] = Support(np.array([1, 1, 0]))
wallDict["thickness"] = 0.05 # m
wall1 = Wall(wallDict)

patchTestModel = PlateModel("plateModel1")
patchTestModel.addPlate(plate1)
patchTestModel.addWall(wall1)

patchTestModel.addLoad(lineLoad)

import matplotlib.pyplot as plt
import pandas as pd
# pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
from displayModel import *

# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
generateMesh(patchTestModel, showGmshMesh=False, elementType='QUAD', nEdgeNodes=3, order='quadratic')


# compute
# ElemType: Quadrangluar or Triangular + Linear or Quadratic or MITC + Reduced or Normal Integration
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']

#% constant shear - Vx
solveModel(patchTestModel, resultsScaleIntForces = (1, 1e0), resultsScaleVertDisp = 1e-7, elementDefinition=elemTypes[3], internalForcePosition = 'nodes', smoothedValues =False)
plotResults(patchTestModel,displacementPlot='text+mesh', verticalDisplacement=False, bendingMomentsToPlot=['x'] ,shearForcesToPlot=[])
currentFig = plt.gcf()
currentAx = plt.gca()
currentFig.set_size_inches(figRatio*figSize,1*figSize)
currentAx.set_xlim(xlim1, xlim2)
currentAx.set_ylim(ylim1, ylim2)
plt.axis('off')
imgName = 'constantM_Mx.svg'
savingPath = os.path.join(myPath, imgName)
plt.savefig(savingPath)