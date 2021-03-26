#%%
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option("display.max_rows", None, "display.max_columns", None)

# import sample plate and show it
import sampleGeometry
from displayModel import *
sampleModel=sampleGeometry.importModel()
plotInputGeometry(sampleModel)

# create mesh
from generateMesh import *
meshInput1=MeshInput(showGmshMesh=False, elementType='QUAD', nEdgeNodes=11)

# meshInput1=MeshInput(showGmshMesh=True, elementType='QUAD', meshSize=5e-2)
generateMesh(sampleModel, meshInput1)

# compute
from solveModel import *
outPos, values = solveModel(sampleModel)
wVert = values[:,0]
# print(wVert)
# df = pd.DataFrame({'x': outPos[:,0], 'y': outPos[:,1], 'w':wVert})
# print(df)

for i, coord in enumerate(outPos):
    sampleModel.geometryInterface.text(coord[0], coord[1], '{:.2f}'.format(wVert[i]*1000))

# display results
plotResults(sampleModel)
plt.show()


#%%
# coords = sampleModel.results.bendingMomentsPositions
# z = sampleModel.results.bendingMoments[:,0]
# nEdge = int(np.sqrt(z.size))
# e = 1/nEdge/2
# numC = z.size-(nEdge-1)*4
# # print('nedge ',nEdge)
# # print('zsize ', z.size)
# # print('numc: ', numC)


# Llow = e*1.2
# Lup = (1-e)*0.98

# # print('lLow: ', Llow)
# # print('Lup: ', Lup)
# newCoords = np.zeros((numC, 2))
# newZ = np.zeros(numC)
# k=0

# # print('e ', e)
# for i, c in enumerate(coords):
#     # print(1,c[0], c[1])
#     if (c[0]>Llow and c[0]<Lup) and (c[1]>Llow and c[1]<Lup):
#         # print(2,c[0], c[1])
#         newCoords[k, :]=c
#         newZ[k]=z[i]
#         k+=1
# # print('k ',k)
# # print(newZ)

# # %matplotlib
# fig=plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)