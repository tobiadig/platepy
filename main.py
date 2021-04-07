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
# meshInput1=MeshInput(showGmshMesh=False, elementType='QUAD', nEdgeNodes=21)
meshInput1=MeshInput(showGmshMesh=True, elementType='TRI', meshSize=8e-2)
generateMesh(sampleModel, meshInput1)

# compute
from solveModel import *
solveModel(sampleModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3)

# display results
plotResults(sampleModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=[],shearForcesToPlot=[])
plt.show()