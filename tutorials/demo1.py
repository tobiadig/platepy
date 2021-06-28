# ------------------------------------------------------------------------------
#
# platepy demo 1
#
# Complicated geometry
#
# ------------------------------------------------------------------------------

from platepy import *


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
plateDict["outlineCoords"]=np.array([[0,0], [a,0],[a,b], [2*a,b],[2*a,1.5*b], [2*a, 2*b], [0,2*b], [0,0]])
plateDict["thickness"] = h
plateDict["body"]=C25_30
plate1 = Plate(plateDict)

wallDict = {}
wallDict["outlineCoords"] = np.array([[0,0],[a,0]])
wallDict["support"] = np.array([1, 0, 1])
wallDict["thickness"] = 0.5 # m
wall1 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[0.5*a,0.5*b],[0.5*a,b]])
wall2 = Wall(wallDict)

wallDict["outlineCoords"] = np.array([[2*a,b],[2*a,1.5*b],[1.5*a, 1.5*b]])
wall3 = Wall(wallDict)

columnDict = {}
columnDict["outlineCoords"] = np.array([[0.5*a,1.5*b]])
columnDict["body"] = C25_30
columnDict["support"] = np.array([1, 0, 0])
columnDict["crossSection"] = None
columnDict["width"] = 1
col1 = Column(columnDict, isInPlate=True)

columnDict["outlineCoords"] = np.array([[2*a,2*b]])
col2 = Column(columnDict)

columnDict["outlineCoords"] = np.array([[0,2*b]])
col3 = Column(columnDict)

firstModel = PlateModel()
firstModel.addPlate(plate1)

firstModel.addWall(wall1)
firstModel.addWall(wall2)
firstModel.addWall(wall3)

firstModel.addColumn(col1)
firstModel.addColumn(col2)
firstModel.addColumn(col3)

firstModel.addLoad(distributedLoad)

# import sample plate and show it

generateMesh(firstModel, showGmshMesh=False, elementDefinition='MITC-4-N', meshSize=6e-1)

# plotMesh(firstModel, plotNodes = False, plotStrucElements = False, plotPoints = True)
# compute

solveModel(firstModel, resultsScales = (1e3,1,1), internalForcePosition = 'center')
# plotInputGeometry(firstModel)
# plotMesh(firstModel, plotStrucElements=False, plotNodes=False, plotPoints=True)
plt.show()
#%%
# startCoord = (1.5*a,b)
# endCoord = (0.5*a,2*b)

# nEvaluationPoints = 1000
# computeBeamComponents(firstModel, startCoord, endCoord, nEvaluationPoints,resultsScales = (1,1, 1),integrationWidth = 0, nIntegrationPoints =0)

# # possiblePlots vDisp_line, Mx_line,My_line,Mxy_line,Vx_line,Vy_line
# valuesToPlotList = ['vDisp_line', 'Mx_line', 'Vx_line']
# plotBeamComponent(firstModel, valuesToPlotList, plotOnMesh=False)


# possiblePlots: vDisp, Mx, My, Mxy, Vx, Vy
#                 Mx_Rd_+, My_Rd_+, Mx_Rd_-, My_Rd_-, Mx_Rd_max, My_Rd_max
#                 N_DSB, V_DSB, M_DSB
# plotInputGeometry(firstModel)
valuesToPlotList=['vDisp']
plotResults(firstModel, valuesToPlotList, plotType = '3d',saveToSVG=True, saveImage=False)
plt.show()

