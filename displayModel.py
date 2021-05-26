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
import base64
from io import BytesIO
from tqdm import tqdm

def plotInputGeometry(self, figaspect = 1):
    w, h = plt.figaspect(1.)
    mult = 1.5
    fig, axGeometry = plt.subplots(figsize=(w*mult, h*mult))

    for plate in self.plates:
        plate.plot(axGeometry)  

    for wall in self.walls:
        wall.plot(axGeometry)

    for column in self.columns:
        column.plot(axGeometry)
        
    for uz in self.downStandBeams:
        uz.plot(axGeometry)

    self.axes['InputGeometry'] = axGeometry
    return fig,axGeometry

def plotResults(self, verticalDisplacement = True,displacementPlot = 'isolines', bendingMomentsToPlot = [], shearForcesToPlot = [], saveImage=False):
    outAxis = []
    outFig = []
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

        theTitle='w'
        fig, axOut = plotInternalForce(self,displacementPlot,theTitle, x,y,z,saveImage)

        self.axes[theTitle] = axOut
        outAxis.append(axOut)

    if outVal[1]: #plot Mx
        x = self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z = self.results.bendingMoments[:,0]
        theTitle='Mx'
        fig, axOut = plotInternalForce(self,displacementPlot,theTitle, x,y,z,saveImage)

        self.axes[theTitle] = axOut
        outAxis.append(axOut)

    if outVal[2]: # plot My
        x= self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z= self.results.bendingMoments[:,1]
        theTitle='My'
        fig, axOut = plotInternalForce(self,displacementPlot,theTitle, x,y,z,saveImage)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    if outVal[3]:  # plot Mxy
        x= self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z= self.results.bendingMoments[:,2]
        theTitle='Mxy'
        fig, axOut = plotInternalForce(self,displacementPlot,theTitle, x,y,z,saveImage)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    if outVal[4]: # plot Vx
        x = self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z = self.results.shearForces[:,0]
        theTitle='Vx'
        fig, axOut = plotInternalForce(self,displacementPlot,theTitle, x,y,z,saveImage)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    if outVal[5]: # plot Vy
        x = self.results.internalForcesPositions[:,0]
        y = self.results.internalForcesPositions[:,1]
        z = self.results.shearForces[:,1]
        theTitle='Vy'
        fig, axOut = plotInternalForce(self,displacementPlot,theTitle, x,y,z,saveImage)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    return outFig, outAxis

def plotInternalForce(self,plotType,theTitle, x,y,z,saveImage):
    if plotType == 'isolines':
        fig,axOut = myIsoPlot(self,x,y,z,theTitle=theTitle)
        if saveImage:
            buf = BytesIO()
            fig.savefig(buf, format="png")
            data5 = base64.b64encode(buf.getbuffer()).decode("ascii")
            # plt.savefig(r'C:\Users\Diggelmann\Desktop\FEMFlask\static\images\new_plot.png')
            # outFig.append(data5)
    elif plotType == '3d':
        fig=plt.figure()
        axOut = fig.gca(projection='3d')
        axOut.plot_trisurf(x,y,z,cmap=cm.jet)
    elif plotType == 'text':
        fig, axOut = myTextPlot(self, x,y,z,theTitle=theTitle)
    elif plotType == 'text+mesh':
        fig, axOut = myTextOnMeshPlot(self,x,y,z, theTitle = theTitle)
    else:
        raise TypeError('type of plot does not exist')

    return fig, axOut

def myTextPlot(self,x,y,z, theTitle = ''):
    fig,outAx = plotInputGeometry(self)
    for i,a in enumerate(z):
        outAx.text(x[i], y[i], '{:.4f}'.format(a)) # number of Nachkomastellen to be displayed
        outAx.scatter(x[i], y[i], marker = ".", facecolor = "k")

    outAx.set_title(theTitle)
    return fig, outAx

def myTextOnMeshPlot(self,x,y,z, theTitle = ''):
    fig, outAx = plotMesh(self, plotNodes=False, plotStrucElements=False)
    for i,a in enumerate(z):
        outAx.text(x[i]+0.03, y[i]+0.03, '{:.2f}'.format(a)) # number of Nachkomastellen to be displayed
        outAx.scatter(x[i], y[i], marker = ".", facecolor = "r",zorder=0)

    outAx.set_title(theTitle)
    return fig, outAx


def myIsoPlot(self,x,y,z, theTitle = ''):
    fig,outAx = plotInputGeometry(self)

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

    zMinString = '{:.3f}'.format(z[iMMin])
    if np.abs(z[iMMin])>0.1:
        outAx.text(x[iMMin],y[iMMin], zMinString,color='r', bbox=dict(facecolor='w', edgecolor='red'), zorder=1000)

    zMaxString = '{:.3f}'.format(z[iMMax])
    if np.abs(z[iMMax])>0.1:
        outAx.text(x[iMMax],y[iMMax], zMaxString,color='b', bbox=dict(facecolor='w', edgecolor='blue'), zorder=1000)

    outAx.set_title(theTitle)
    return fig,outAx

def plotMesh(self, plotNodes = True, plotStrucElements = True, plotPoints = False):
    if plotStrucElements:
        fig,outAx = plotInputGeometry(self)
    else:
        fig, outAx = plt.subplots()
    nodes = self.mesh.nodesArray.to_numpy()

    #plot lines:
    elementsList = self.mesh.elementsList
    print('Plotting mesh')
    for element in tqdm(elementsList):
        elemNodes = element.coherentConnectivity.to_numpy()[:,0]
        # nEdges = element.shape
        nEdges = 4
        

        elemCoords = element.coordinates
        if len(elemCoords)<nEdges:
            continue
        xValues = np.zeros((nEdges+1,))
        yValues = np.zeros((nEdges+1,))

        xValues[0:-1] = elemCoords[0:nEdges,0]
        xValues[-1] = elemCoords[0,0]

        yValues[0:-1] = elemCoords[0:nEdges,1]
        yValues[-1] = elemCoords[0,1]
        outAx.plot(xValues, yValues, color='k', zorder=-1)
    #plot nodes:
    k=1
    if plotNodes:
        for node in nodes:
            outAx.scatter(node[0], node[1], facecolor='r', marker='.')
            outAx.text(node[0]+0.03, node[1]+0.03, k)
            k+=1
    if plotPoints:
        k=1
        for node in nodes:
            outAx.scatter(node[0], node[1], facecolor='r', marker='.')
            # outAx.text(node[0], node[1], k)
            k+=1
    return fig, outAx 

def plotBeamComponent(self,lineName, verticalDisplacement = True, bendingMomentsToPlot = [], shearForcesToPlot = [], plotOnMesh=False):
    outAxis = []
    outFig = []
    valPoss = ['x', 'y', 'xy']
    
    outVal = np.zeros(6, dtype=bool)
    outVal[0]=verticalDisplacement
    for i, a in enumerate(valPoss):
        if a in bendingMomentsToPlot:
            outVal[i+1] = True
        if a in shearForcesToPlot:
            outVal[i+4] = True

    schnitt = self.results.schnittList[lineName]
    x=schnitt.arrayEvaluationPoints[:,0]
    y=schnitt.arrayEvaluationPoints[:,1]

    if outVal[0]: #plot vertical displacements
        z= schnitt.verticalDisplacements

        theTitle='w'
        fig, axOut = plotSchnittValues(self,theTitle, x,y,z,plotOnMesh)

        self.axes[theTitle] = axOut
        outAxis.append(axOut)

    if outVal[1]: #plot Mx
        z = schnitt.bendingMoments[:,0]

        theTitle='Mx'
        fig, axOut = plotSchnittValues(self,theTitle, x,y,z,plotOnMesh)

        self.axes[theTitle] = axOut
        outAxis.append(axOut)

    if outVal[2]: # plot My

        z= schnitt.bendingMoments[:,1]
        theTitle='My'
        fig, axOut = plotSchnittValues(self,theTitle, x,y,z,plotOnMesh)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    if outVal[3]:  # plot Mxy

        z= schnitt.bendingMoments[:,2]
        theTitle='Mxy'
        fig, axOut = plotSchnittValues(self,theTitle, x,y,z,plotOnMesh)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    if outVal[4]: # plot Vx

        z = schnitt.shearForces[:,0]
        theTitle='Vx'
        fig, axOut = plotSchnittValues(self,theTitle, x,y,z,plotOnMesh)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    if outVal[5]: # plot Vy

        z = schnitt.shearForces[:,1]
        theTitle='Vy'
        fig, axOut = plotSchnittValues(self,theTitle, x,y,z,plotOnMesh)

        self.axes[theTitle] = axOut
        outAxis.append(axOut) 

    return outFig, outAxis


def plotSchnittValues(self,theTitle, x,y,z,plotOnMesh):
    if not plotOnMesh:
        fig,outAx = plotInputGeometry(self)
    else:
        fig, outAx = plotMesh(self, plotNodes=False, plotStrucElements=False, plotPoints =True)

    iMMax = np.argmax(np.abs(z))
    # iMMin = np.argmin(z)


    magValue = 2/np.max(np.abs(z))
    zNorm = z*magValue

    lineDir = np.array([(x[-1]-x[0]),(y[-1]-y[0])])
    lineDir = lineDir/np.sqrt(lineDir[0]**2+lineDir[1]**2)
    normLineDir = np.array([lineDir[1], -lineDir[0]])
    zPoints = np.zeros((x.shape[0],2))
    zPoints[:,0] = x+zNorm*normLineDir[0]
    zPoints[:,1] = y+zNorm*normLineDir[1]

    for i in range(0,x.shape[0]):
        outAx.plot(np.array([x[i], zPoints[i,0]]),np.array([y[i], zPoints[i,1]]), color = 'grey')

    outAx.plot(x,y,color='k')
    outAx.plot(zPoints[:,0], zPoints[:,1], color='k')





    # zMinString = '{:.3f}'.format(z[iMMin])
    # if np.abs(z[iMMin])>0.1:
    #     outAx.text(zPoints[iMMin,0],zPoints[iMMin,1], zMinString,color='r', bbox=dict(facecolor='w', edgecolor='red'), zorder=1000)

    zMaxString = '{:.3f}'.format(z[iMMax])
    if np.abs(z[iMMax])>0.1:
        outAx.text(zPoints[iMMax,0],zPoints[iMMax,1], zMaxString,color='k', bbox=dict(facecolor='w', edgecolor='grey'), zorder=1000)




    # xLim = np.array([np.min(x), np.max(x)])
    # yLim = np.array([np.min(y), np.max(y)])
    # a=xLim[1]-xLim[0]
    # b= yLim[1]-yLim[0]

    # # fig.set_size_inches(10*a,10*b)
    # marginWidth = 0.1
    # outAx.set_xlim(xLim[0]-marginWidth*a, xLim[1]+marginWidth*a)
    # outAx.set_ylim(yLim[0]-marginWidth*b, yLim[1]+marginWidth*b)


    outAx.set_title(theTitle)
    return fig,outAx




