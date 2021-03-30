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