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
            fig = plt.figure()
            axVertDisp = fig.gca(projection='3d')
            axVertDisp.plot_trisurf(self.results.outPos[:,0],self.results.outPos[:,1],self.results.wVert,cmap=cm.jet)
        else:
            raise TypeError('type of plot does not exist')

        self.axes['VerticalDisplacement'] = axVertDisp
        outAxis.append(axVertDisp)

    if outVal[1]:
        newCoords, newZ = momentsPlotPreparation(self, 0)

        fig=plt.figure()
        axMx = fig.gca(projection='3d')
        axMx.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)
        self.axes['Mx'] = axMx
        outAxis.append(axMx)
    if outVal[2]:
        newCoords, newZ = momentsPlotPreparation(self, 1)

        fig=plt.figure()
        axMy = fig.gca(projection='3d')
        axMy.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)
        self.axes['My'] = axMy
        outAxis.append(axMy)


    if outVal[3]:
        newCoords, newZ = momentsPlotPreparation(self, 2)

        fig=plt.figure()
        axMxy = fig.gca(projection='3d')
        axMxy.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)
        self.axes['Mxy'] = axMxy
        outAxis.append(axMxy)

    

def momentsPlotPreparation(self, nMoment):
    coords = self.results.bendingMomentsPositions
    z = self.results.bendingMoments[:,nMoment]
    nEdge = int(np.sqrt(z.size))
    e = 1/nEdge/2
    numC = z.size-(nEdge-1)*4
    # print('nedge ',nEdge)
    # print('zsize ', z.size)
    # print('numc: ', numC)
    Llow = e*1.2
    Lup = (1-e)*0.98

    # print('lLow: ', Llow)
    # print('Lup: ', Lup)
    newCoords = np.zeros((numC, 2))
    newZ = np.zeros(numC)
    k=0

    # print('e ', e)
    for i, c in enumerate(coords):
        # print(1,c[0], c[1])
        if (c[0]>Llow and c[0]<Lup) and (c[1]>Llow and c[1]<Lup):
            # print(2,c[0], c[1])
            newCoords[k, :]=c
            newZ[k]=z[i]
            k+=1
    # print('k ',k)
    # print(newZ)

    return newCoords, newZ
