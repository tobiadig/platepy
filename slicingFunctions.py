import numpy as np

def getKCoeff(elementType, coherentElemNodes):
    nNodes = len(coherentElemNodes)
    if elementType in ['L', 'Q', 'MITC4']:
        kCoeff = np.zeros((3*nNodes),dtype=int)
        for i in range(0,3):
            kCoeff[0+i::3]=coherentElemNodes*3+i
        discartedDOF = None
    elif elementType in ['MITC9']:
        kCoeffTemp = np.zeros((3*nNodes),dtype=int)
        for i in range(0,3):
            kCoeffTemp[0+i::3]=coherentElemNodes*3+i

        kCoeff = np.zeros((3*nNodes-1),dtype=int)
        kCoeff[:-2] = kCoeffTemp[:-3]
        kCoeff[-2:] = kCoeffTemp[-2:]
        discartedDOF = kCoeffTemp[-3]
    else:
        raise TypeError('Element type in getCoeff not recognised')
    return kCoeff, discartedDOF


def getRowsColumns(kCoeff, nMatrixDofs):
    rows = np.zeros(nMatrixDofs,dtype=int)
    columns = np.zeros(nMatrixDofs,dtype=int)
    c=0
    for j in kCoeff:
        for i in kCoeff:
            rows[c] = i
            columns[c] = j
            c+=1

    return rows, columns