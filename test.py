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
a='Q-MITC-N'

b=a.split('-')
print(b[2]=='R')