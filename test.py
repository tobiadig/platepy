#%%import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
a=np.random.randn(2500).reshape((50,50))

# define the colormap
cmap = plt.get_cmap('PuOr')

# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

# define the bins and normalize and forcing 0 to be part of the colorbar!
bounds = np.arange(np.min(a),np.max(a),.5)
idx=np.searchsorted(bounds,0)
bounds=np.insert(bounds,idx,0)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

plt.imshow(a,interpolation='none',norm=norm,cmap=cmap)
plt.colorbar()
plt.show()
#%%
test = {}

test['triangular']={}
test['triangular'][1]={}
test['triangular'][1]['points'] = [0,1,2]
print(test['triangular'][1]['points'])
#%%
import numpy as np
A=np.array([[1,2], [3,4]])
B =  np.array([[1,2], [3,4]])
C = np.array([[1,2], [3,4]])



#%%
import numpy as np
N1 = lambda r, s: 1*np.ones(len(r))
N3 = lambda r, s: (1-r)*(1-s)
ri=np.array([1,2,3])
si = np.array([1,2,3])
A=np.array([N1(ri,si),N3(ri,si)])
print(A[1,2])

#%%
import numpy as np
a=np.array([1,2,3])
a=a.reshape((3,1))
print(a.shape)

#%%
import numpy as np
a = np.array([[2,2],[2,2]])
b = np.array([1,2])

print(a/b)

#%%

N9 = lambda r, s: (1-r**2)*(1-s**2)
N8 = lambda r,s: 0.5*(1-s**2)*(1+r)-0.5*N9(r,s)
N7 = lambda r, s: 0.5*(1-r**2)*(1-s)-0.5*N9(r,s)
N6 = lambda r, s: 0.5*(1-s**2)*(1-r)-0.5*N9(r,s)
N5 = lambda r, s: 0.5*(1-r**2)*(1+s)-0.5*N9(r,s)
N4 = lambda r, s: 0.25*(1+r)*(1-s) - 0.5*N7(r,s) - 0.5*N8(r,s) - 0.25*N9(r,s)
N3 = lambda r, s: 0.25*(1-r)*(1-s) - 0.5*N7(r,s) - 0.5*N6(r,s) - 0.25*N9(r,s)
N2 = lambda r, s: 0.25*(1-r)*(1+s) - 0.5*N5(r,s) - 0.5*N6(r,s) - 0.25*N9(r,s)
N1 = lambda r, s: 0.25*(1+r)*(1+s) - 0.5*N5(r,s) - 0.5*N8(r,s) - 0.25*N9(r,s)

print(N1(0.7, 0.7))

N1 = lambda r, s: 0.25*(r**2-r)*(s**2-s)
N2 = lambda r, s: 0.25*(r**2+r)*(s**2-s)
N3 = lambda r, s: 0.25*(r**2+r)*(s**2+s)
N4 = lambda r, s: 0.25*(r**2-r)*(s**2+s)
N5 = lambda r, s: 0.5*(s**2-s)*(1-r**2)
N6 = lambda r, s: 0.5*(r**2+r)*(1-s**2)
N7 = lambda r, s: 0.5*(s**2+s)*(1-r**2)
N8 = lambda r, s: 0.5*(r**2-r)*(1-s**2)
N9 = lambda r, s: (1-r**2)*(1-s**2)

print(N1(-0.7, -0.7))
#%%
import numpy as np
a=np.array([[1,2,3],[1,2,3],[1,2,3]])
b=np.array([[1,2,3],[1,2,3],[1,2,3]])

c=np.matmul(a, b)
c2=np.matmul(c, c)
d=a@b@c

print(c2)
print(d)

#%%
from solveModel import *
import numpy as np
import pandas as pd
elemType = 4
ri = np.array([-1])
si = np.array([0])
xi = np.array([0,5,5,0])
yi = np.array([0,0,5, 5])

N, Bf,Bc, detJ= getLinearVectorizedShapeFunctions(elemType,ri, si, xi, yi)

Bc=Bc[0,:,:]
N = N[0,:,:]
# print('N: ', pd.DataFrame(N))
print(pd.DataFrame(Bc))
ri = -1
si = 0
N, Bb,Bs, detJ = getMITCShapefunctions(ri, si, xi, yi)
print(pd.DataFrame(Bs))

#%%
from scipy import sparse    # To work with sparse matrix
import numpy as np
A= sparse.csr_matrix(np.array([[1,2,3], [1,2,3], [1,2,3]]))
B= sparse.csr_matrix(np.array([[2,2,2]]).transpose())

C=A*B

print(C)

#%%
import numpy as np
ri = np.random.rand(4,)
si = np.random.rand(4,)
print(ri.shape)

N1 = lambda r, s: 0.25*(1-r)*(1-s)
N2 = lambda r, s: 0.25*(1+r)*(1-s)
N3 = lambda r, s: 0.25*(1+r)*(1+s)
N4 = lambda r, s: 0.25*(1-r)*(1+s)

# Form the shape function matrix
Nfun= lambda r, s: [N1(r,s), N2(r,s), N3(r,s), N4(r,s)]


Nval=np.array(Nfun(ri, si))
print(Nval)
Nval=np.moveaxis(Nval, -1, 0)
print(Nval)
#%%
GPr = -0.7745966692414834
GPs = -0.7745966692414834
N4r = 0.25*(1-GPs) + 0.5*GPr*(1-GPs) -0.25*(1-GPs**2) + 0.5*GPr*(1-GPs**2)

N4 = 0.25*(1+GPr)*(1-GPs) - 0.25*(1-GPr**2)*(1-GPs) - 0.25*(1-GPs**2)*(1+GPr) - 0.25*(1-GPr**2)*(1-GPs**2)
print('N4: ', N4)
print('N4r: ',N4r)