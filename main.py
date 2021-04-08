#%%
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
import cubusGeometry
from displayModel import *
sampleModel=cubusGeometry.importModel()
# plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
# meshInput1=MeshInput(showGmshMesh=False, elementType='QUAD', nEdgeNodes=31)
meshInput1=MeshInput(showGmshMesh=True, elementType='QUAD', meshSize=5e-1)
generateMesh(sampleModel, meshInput1)

# compute
from solveModel import *
solveModel(sampleModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3)

# display results
plotResults(sampleModel,displacementPlot='isolines', verticalDisplacement=True, bendingMomentsToPlot=['x', 'y', 'xy'],shearForcesToPlot=['x'])
plt.show()