# ------------------------------------------------------------------------------
#
# platepy demo 2
#
# T-shaped plate with downstand beam
#
# ------------------------------------------------------------------------------
from platepy import *
import numpy as np
import matplotlib.pyplot as plt

C25_30Dic ={}
C25_30Dic['eModule'] = 32
C25_30Dic['gModule'] = 14
C25_30Dic['nu'] = 0.3
C25_30 = Concrete(C25_30Dic)

plateDic = {}
plateDic['thickness'] = 0.1
plateDic['body'] = C25_30
plateDic['outlineCoords'] = np.array([[0,0],[10,0],[10,10],[20,10],[20,20],[-10,20],[-10,10],[0,10],[0,0]])
plate = Plate(plateDic)

wallDic = {}
wallDic['outlineCoords'] = np.array([[0,0],[10,0]])
wallDic['thickness'] = 0.5
wallDic['support'] = np.array([1,0,1])
wall1 = Wall(wallDic)

wallDic['outlineCoords'] = np.array([[-10,20],[20,20]])
wall2 = Wall(wallDic)

colDic = {}
colDic['outlineCoords'] = np.array([[-10,10]])
colDic['width'] = 0.5
colDic['support'] = np.array([1,0,0])
col1 = Column(colDic)

colDic['outlineCoords'] = np.array([[0,10]])
col2 = Column(colDic)

colDic['outlineCoords'] = np.array([[10,10]])
col3 = Column(colDic)

colDic['outlineCoords'] = np.array([[20,10]])
col4 = Column(colDic)

dsbDic = {}
dsbDic['outlineCoords'] = np.array([[0,10], [10,10]])
dsbDic['body'] = C25_30
hBeam = 0.2
bBeam = 0.2
dsbDic['crossSection'] =CrossSection(hBeam*bBeam, hBeam**3*bBeam/12, 0, bBeam, hBeam)
dsb = DownStandBeam(dsbDic)

load = Load('area', np.array([-1,0,0]))

demoModel = PlateModel()
demoModel = PlateModel()
demoModel.addPlate(plate)
demoModel.addWall(wall1)
demoModel.addWall(wall2)
demoModel.addColumn(col1)
demoModel.addColumn(col2)
demoModel.addColumn(col3)
demoModel.addColumn(col4)
demoModel.addLoad(load)
demoModel.addDownStandBeam(dsb)

generateMesh(demoModel)
solveModel(demoModel)
plotResults(demoModel, valuesToPlotList=['vDisp'])
plt.show()