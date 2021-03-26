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

def plotInputGeometry(self):
    w, h = plt.figaspect(1.)
    mult = 1.5
    fig1, axGeometry = plt.subplots(figsize=(w*mult, h*mult))
    fig2, axResults = plt.subplots(figsize=(w*mult, h*mult))

    
    for plate in self.plates:
        plate.plot(axGeometry)
        plate.plot(axResults)      

    for wall in self.walls:
        wall.plot(axGeometry)
        wall.plot(axResults)

    for column in self.columns:
        column.plot(axGeometry)
        column.plot(axResults)

    self.geometryInterface = axGeometry
    self.resultsInterface = axResults
    return axGeometry

def plotResults(self):
    axResults = self.resultsInterface

    x= self.results.outPos[:,0]
    y = self.results.outPos[:,1]
    z= self.results.wVert*1000

    # ax = fig.gca()
    # ax.plot_trisurf(self.results.outPos[:,0],self.results.outPos[:,1],self.results.wVert,cmap=cm.jet)
    cs = plt.tricontour(x,y,z,colors='r')

    axResults.clabel(cs)
    # for plate in self.plates:
    #     plate.plot()        

    # for wall in self.walls:
    #     wall.plot()

    # for column in self.columns:
    #     column.plot()

    # fig = plt.gcf()
    # ax = plt.gca()
    xLim = np.array([np.min(x), np.max(x)])
    yLim = np.array([np.min(y), np.max(y)])
    a=xLim[1]-xLim[0]
    b= yLim[1]-yLim[0]

    # fig.set_size_inches(10*a,10*b)
    marginWidth = 0.1
    axResults.set_xlim(xLim[0]-marginWidth*a, xLim[1]+marginWidth*a)
    axResults.set_ylim(yLim[0]-marginWidth*b, yLim[1]+marginWidth*b)
    zMaxString = '{:.2f}'.format(self.results.wMax[2])
    axResults.text(self.results.wMax[0],self.results.wMax[1], zMaxString,color='r', bbox=dict(facecolor='none', edgecolor='red'))
    self.resultsInterface = axResults
    return axResults
