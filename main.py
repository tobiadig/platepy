#%%
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
import sampleGeometry
from displayModel import *
sampleModel=sampleGeometry.importModel()
# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
meshInput1=MeshInput(showGmshMesh=False, elementType='QUAD', nEdgeNodes=21)
# meshInput1=MeshInput(showGmshMesh=True, elementType='TRI', meshSize=8e-2)
generateMesh(sampleModel, meshInput1)

# compute
from solveModel import *
solveModel(sampleModel, resultsScaleIntForces = 1e-3, resultsScaleVertDisp = 1e-13)
# wVert = values[:,0]

# for i, coord in enumerate(outPos):
#     sampleModel.geometryInterface.text(coord[0], coord[1], '{:.2f}'.format(wVert[i]*1000))

# display results
# plotResults(sampleModel)
# plt.show()

#%% max moment

plotResults(sampleModel, verticalDisplacement=True, bendingMomentsToPlot=['x', 'y', 'xy'],shearForcesToPlot=['x', 'y'])

M = sampleModel.results.bendingMoments[:,0]

print('Moments: ',np.min(M)*1e-5)
plt.show()
from AnalyticPlateSolutions import *

# analythical vetical displacement, rectangular, simply supported distributed
# all inputs in kN and m
    # pOpts : plate options (e.g. dataclass or named tuple)
    #     pOpts.shape    = "rectangular" | "circular"
    #     pOpts.depth    = "thick" | "thin"
    #     pOpts.support  = "simplySupported" | "clamped"
    #     pOpts.geometry = list of pertinent parameters sufficient to describe geometry
    #     pOpts.material = list of pertinent parameters sufficient to describe material
    # lOpts : load options (e.g. dataclass or named tuple)
    #     lOpts.type      = "concentrated" | "distributed"
    #     lOpts.position  = list of pertinent parameters sufficient to describe position
    #     lOpts.magnitude = magnitude of vertical force
    # sOpts : solution options (e.g. dataclass or named tuple)
    #     sOpts.nTerms = list (e.g. describing the amount of series expansion terms requested)
    # inPos : positions at which output quantities are requested (i.e. array of 2D points)
###################
#ANALYTHICAL SOLUTION
pOpts = POpts()
pOpts.shape="rectangular"
pOpts.depth = "thin"
pOpts.support = "simplySupported"
pOpts.geometry = (1,1)
pOpts.material = Material(10920, 0.3, 0.1) #E, nu and h

lOpts = LOpts()
lOpts.type = "distributed"
lOpts.magnitude = 1

sOpts = SOpts()
sOpts.nTerms = 40

inPos=np.array([[0,0]])

quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
# print('Maximum moment, analythical: {:.2f}'.format(values[0,2]*1000)) 