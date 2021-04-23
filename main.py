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
# generateMesh(sampleModel, showGmshMesh=False, elementType='QUAD', meshSize=5e-1)
generateMesh(sampleModel, showGmshMesh=False, elementType='QUAD', nEdgeNodes=3, order ='quadratic')

# compute
from solveModel import *
elemTypes = ['L-R', 'MITC4-N', 'Q-R', 'MITC9-N']
solveModel(sampleModel, resultsScaleIntForces = (1, 1), resultsScaleVertDisp = 1e3,elemType=elemTypes[2], internalForcePosition = 'center')

# display results
plotResults(sampleModel,displacementPlot='text', verticalDisplacement=False, bendingMomentsToPlot=[],shearForcesToPlot=[])
print(pd.DataFrame(sampleModel.results.wVert))
plt.show()