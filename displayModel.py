''' Module Information
-----------------------------------------------------------
Purpose of module: display given geometry and results of a plateModel class using matplotlib
-----------------------------------------------------------
- Copywrite Tobia Diggelmann (ETH Zurich) 24.03.2021
'''
# Basic modules
import numpy as np
import matplotlib
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
        x= self.results.outPos[:,0]
        y = self.results.outPos[:,1]
        z= self.results.wVert

        if displacementPlot == 'isolines':
            axVertDisp = myIsoPlot(self,x,y,z, theTitle='Vertical Displacements')


        elif displacementPlot == '3d':
            fig = plt.figure()
            axVertDisp = fig.gca(projection='3d')
            axVertDisp.plot_trisurf(self.results.outPos[:,0],self.results.outPos[:,1],self.results.wVert,cmap=cm.jet)

        # elif displacementPlot == 'table':
            # for i, coord in enumerate(outPos):
                # sampleModel.geometryInterface.text(coord[0], coord[1], '{:.2f}'.format(wVert[i]*1000))
        else:
            raise TypeError('type of plot does not exist')


        self.axes['Mx'] = axVertDisp
        outAxis.append(axVertDisp)

    if outVal[1]: #plot Mx
        x = self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z = self.results.bendingMoments[:,0]
        if displacementPlot == 'isolines':
            axMx = myIsoPlot(self,x,y,z, theTitle='M_x')

        elif displacementPlot == '3d':
            fig=plt.figure()
            axMx = fig.gca(projection='3d')
            axMx.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)

        else:
            raise TypeError('type of plot does not exist')

        self.axes['Mx'] = axMx
        outAxis.append(axMx)        

    if outVal[2]: # plot My
        x= self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z= self.results.bendingMoments[:,1]
        if displacementPlot == 'isolines':
            axMy = myIsoPlot(self,x,y,z,theTitle='M_y')

        elif displacementPlot == '3d':
            fig=plt.figure()
            axMy = fig.gca(projection='3d')
            axMy.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)
        else:
            raise TypeError('type of plot does not exist')
        self.axes['My'] = axMy
        outAxis.append(axMy)

    if outVal[3]:  # plot Mxy
        x= self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z= self.results.bendingMoments[:,2]
        if displacementPlot == 'isolines':
            axMxy = myIsoPlot(self,x,y,z,theTitle='M_{xy}')
           
        elif displacementPlot == '3d':
            fig=plt.figure()
            axMxy = fig.gca(projection='3d')
            axMxy.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)
        else:
            raise TypeError('type of plot does not exist')
        self.axes['Mxy'] = axMxy
        outAxis.append(axMxy)

    if outVal[4]: # plot Vx
        x = self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z = self.results.shearForces[:,0]
        if displacementPlot == 'isolines':
            axVx = myIsoPlot(self,x,y,z,theTitle='V_x')

        elif displacementPlot == '3d':
            fig=plt.figure()
            axVx = fig.gca(projection='3d')
            axVx.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)

        else:
            raise TypeError('type of plot does not exist')

        self.axes['Vx'] = axVx
        outAxis.append(axVx)

    if outVal[5]: # plot Vy
        x = self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z = self.results.shearForces[:,1]
        if displacementPlot == 'isolines':
            axVy = myIsoPlot(self,x,y,z,theTitle='V_y')

        elif displacementPlot == '3d':
            fig=plt.figure()
            axVy = fig.gca(projection='3d')
            axVy.plot_trisurf(newCoords[:,0],newCoords[:,1],newZ,cmap=cm.jet)

        else:
            raise TypeError('type of plot does not exist')

        self.axes['Vy'] = axVy
        outAxis.append(axVy)


def myIsoPlot(self,x,y,z, theTitle = ''):
    outAx = plotInputGeometry(self)

    red = np.array([255/256, 45/256, 0/256, 1])
    grey = np.array([128/256, 128/256, 128/256, 1])
    blue = np.array([0/256, 55/256, 255/256, 1])
    newcolors = [red,grey,blue]
    mycmp = matplotlib.colors.ListedColormap(newcolors)
    bounds = np.array([-1e70, -0.00001,0.000001, 1e70])

    norm = matplotlib.colors.BoundaryNorm(bounds, 3)

    # myColorMap = matplotlib.colors.Colormap('redAndBlue')
    # myColorMap.set_over(self, color='r')
    # myColorMap.set_under(self, color='b')

    # cs = plt.tricontour(x,y,z,colors='r')
    cs = plt.tricontour(x,y,z,cmap=mycmp, norm=norm)

    outAx.clabel(cs, fmt='%1.1f')
    xLim = np.array([np.min(x), np.max(x)])
    yLim = np.array([np.min(y), np.max(y)])
    a=xLim[1]-xLim[0]
    b= yLim[1]-yLim[0]

    # fig.set_size_inches(10*a,10*b)
    marginWidth = 0.1
    outAx.set_xlim(xLim[0]-marginWidth*a, xLim[1]+marginWidth*a)
    outAx.set_ylim(yLim[0]-marginWidth*b, yLim[1]+marginWidth*b)

    iMMax = np.argmax(z)
    iMMin = np.argmin(z)

    zMinString = '{:.2f}'.format(z[iMMin])
    if np.abs(z[iMMin])>0.1:
        outAx.text(x[iMMin],y[iMMin], zMinString,color='r', bbox=dict(facecolor='w', edgecolor='red'), zorder=1000)

    zMaxString = '{:.2f}'.format(z[iMMax])
    if np.abs(z[iMMax])>0.1:
        outAx.text(x[iMMax],y[iMMax], zMaxString,color='b', bbox=dict(facecolor='w', edgecolor='blue'), zorder=1000)

    outAx.set_title(theTitle)
    return outAx



























# def momentsPlotPreparation(self, nMoment):
#     coords = self.results.bendingMomentsPositions
#     z = self.results.bendingMoments[:,nMoment]
#     nEdge = int(np.sqrt(z.size))
#     e = 1/nEdge/2
#     numC = z.size-(nEdge-1)*4
#     # print('nedge ',nEdge)
#     # print('zsize ', z.size)
#     # print('numc: ', numC)
#     a=10000
#     Llow = e*1.2*a
#     Lup = (1-e)*0.98*a

#     # print('lLow: ', Llow)
#     # print('Lup: ', Lup)
#     newCoords = np.zeros((numC, 2))
#     newZ = np.zeros(numC)
#     k=0

#     # print('e ', e)
#     for i, c in enumerate(coords):
#         # print(1,c[0], c[1])
#         if (c[0]>Llow and c[0]<Lup) and (c[1]>Llow and c[1]<Lup):
#             # print(2,c[0], c[1])
#             newCoords[k, :]=c
#             newZ[k]=z[i]
#             k+=1
#     # print('k ',k)
#     # print(newZ)

#     return newCoords, newZ
