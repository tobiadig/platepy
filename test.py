#%%import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

nodes = np.array([[1,0,0,0],
                    [2,0,1,0],
                    [4, 1,1,0],
                    [5,0,1,0]])
# nodes={}
# nodes[1]= np.array([0,0,0])
# nodes[2]=np.array([0,1,0])
# nodes[4]=np.array([1,1,0])
# nodes[5]=np.array([0,1,0])
nodesPd = pd.DataFrame(nodes[:,1:], index=nodes[:,0])
myNodes = np.array([1,2])

print(nodesPd.loc[myNodes].to_numpy())
