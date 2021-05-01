import numpy as np
from scipy.linalg import block_diag # to create the rotation matrix
def rotMatrix(theta):
    '''
        creates rotation matrix of angle theta for a single node
    '''
    A = np.array([[1., 0, 0],
                  [0, np.cos(theta), -np.sin(theta)],
                  [0, np.sin(theta), np.cos(theta)]], dtype=float)

    return A

def getRotationMatrix(elementType, elemNodesRotations):
    Ri = []
    RiInit = False
    for dofRot in elemNodesRotations:
        if not(RiInit):
            R=rotMatrix(dofRot)
            RiInit=True
        else:
            R = block_diag(R, rotMatrix(dofRot))
    if elementType in ['L', 'Q', 'MITC4']:
        return R
    elif elementType in ['MITC9']:
        indixes = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26]) #ohne 24!
        newR = R[:,indixes]
        newR = newR[indixes,:]
        return newR
    else:
        raise TypeError('Element type in getRotationMatrix not recognised')