''' Module Information
-----------------------------------------------------------
Purpose of module: display given geometry and results of a plateModel class using matplotlib
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''
# Basic modules
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D # 3D plot
from matplotlib import cm   # contour plot

def plotInputGeometry(self, figaspect = 1):
    w, h = plt.figaspect(1.)
    mult = 1.5
    fig1, axGeometry = plt.subplots(figsize=(w*mult, h*mult))

    for plate in self.plates:
        plate.plot(axGeometry)  

    for wall in self.walls:
        wall.plot(axGeometry)

    for column in self.columns:
        column.plot(axGeometry)

    self.axes['InputGeometry'] = axGeometry
    return axGeometry

def plotResults(self, verticalDisplacement = True,displacementPlot = 'isolines', bendingMomentsToPlot = [], shearForcesToPlot = []):
    outAxis = []
    valPoss = ['x', 'y', 'xy']
    
    outVal = np.zeros(6, dtype=bool)
    outVal[0]=verticalDisplacement
    for i, a in enumerate(valPoss):
        if a in bendingMomentsToPlot:
            outVal[i+1] = True
        if a in shearForcesToPlot:
            outVal[i+4] = True

    if outVal[0]: #plot vertical displacements
        if displacementPlot == 'isolines':
            axVertDisp = plotInputGeometry(self)
        
            x= self.results.outPos[:,0]
            y = self.results.outPos[:,1]
            z= self.results.wVert*1000
            cs = plt.tricontour(x,y,z,colors='r')

            axVertDisp.clabel(cs)
            xLim = np.array([np.min(x), np.max(x)])
            yLim = np.array([np.min(y), np.max(y)])
            a=xLim[1]-xLim[0]
            b= yLim[1]-yLim[0]

            # fig.set_size_inches(10*a,10*b)
            marginWidth = 0.1
            axVertDisp.set_xlim(xLim[0]-marginWidth*a, xLim[1]+marginWidth*a)
            axVertDisp.set_ylim(yLim[0]-marginWidth*b, yLim[1]+marginWidth*b)
            zMaxString = '{:.2f}'.format(self.results.wMax[2])
            axVertDisp.text(self.results.wMax[0],self.results.wMax[1], zMaxString,color='r', bbox=dict(facecolor='none', edgecolor='red'))
        elif displacementPlot == '3d':
            axVertDisp = fig.gca(projection='3d')
            axVertDisp.plot_trisurf(self.results.outPos[:,0],self.results.outPos[:,1],self.results.wVert,cmap=cm.jet)
        else:
            raise TypeError('type of plot does not exist')

        self.axes['VerticalDisplacement'] = axVertDisp
        outAxis.append(axVertDisp)

