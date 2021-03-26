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
# meshInput1=MeshInput(showGmshMesh=True, elementType='TRI', nEdgeNodes=11)
meshInput1=MeshInput(showGmshMesh=True, elementType='QUAD', meshSize=5e-2)
generateMesh(sampleModel, meshInput1)

# compute
from solveModel import *
solveModel(sampleModel)
# wVert = values[:,0]

# for i, coord in enumerate(outPos):
#     sampleModel.geometryInterface.text(coord[0], coord[1], '{:.2f}'.format(wVert[i]*1000))

# display results
# plotResults(sampleModel)
# plt.show()


#%%
