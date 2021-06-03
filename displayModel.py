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
import matplotlib.tri as tri

from mpl_toolkits.mplot3d import Axes3D # 3D plot
from matplotlib import cm   # contour plot
import base64
from io import BytesIO
from tqdm import tqdm

def plotInputGeometry(self, figaspect = 1):
    w, h = plt.figaspect(figaspect)
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
    for p in self.loads:
        if p.case == 'point':
            mySize = 100
            axGeometry.scatter(p.outlineCoords[0], p.outlineCoords[1],marker="x", color='r', s=mySize )
            froceString = 'F = '+str(p.magnitude[0])
            axGeometry.text(p.outlineCoords[0]*1.1, p.outlineCoords[1]*1.1,froceString, fontsize = 15)

    self.axes['InputGeometry'] = axGeometry
    return fig,axGeometry

def plotResults(self, valuesToPlotList, plotType = 'isolines',saveToSVG=False, saveImage=False):
    outAxis = []
    outFig = []
    resultsDictionary = self.resultsInformation.resultsDictionary
    for valueToPlot in valuesToPlotList:
        theTitle=valueToPlot
        resultToPlot = resultsDictionary[valueToPlot]
        myZ = resultToPlot.z*resultToPlot.resultScale
        
        fig, axOut = plotInternalForce(self,plotType,theTitle, resultToPlot.x,resultToPlot.y,myZ,saveImage)
        self.axes[theTitle] = valueToPlot
        outAxis.append(axOut)
        if saveToSVG:
            fig.savefig(theTitle+'.svg')
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
        percentageTrianglesToMantain = 98
        triang = removeTrianglesOutsidePlate(x,y,percentageTrianglesToMantain)
        axOut = fig.gca(projection='3d')
        axOut.plot_trisurf(triang,z,cmap=plt.get_cmap('winter'))
        axOut.grid(False)
        # Hide axes ticks
        axOut.set_xticks([])
        axOut.set_yticks([])
        axOut.set_zticks([])
    elif plotType == 'text':
        fig, axOut = myTextPlot(self, x,y,z,theTitle=theTitle)
    elif plotType == 'text+mesh':
        fig, axOut = myTextOnMeshPlot(self,x,y,z, theTitle = theTitle)
    else:
        raise TypeError('type of plot does not exist')

    return fig, axOut

def removeTrianglesOutsidePlate(x,y,percentageTrianglesToMantain):
    triang = tri.Triangulation(x, y)
    triangles = triang.triangles
    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
    alpha = np.percentile(maxi, percentageTrianglesToMantain)

    triang.set_mask(maxi > alpha)
    return triang

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
    percentageTrianglesToMantain = 98    
    triang = removeTrianglesOutsidePlate(x,y,percentageTrianglesToMantain)
    cs = plt.tricontour(triang,z,cmap=mycmp, norm=norm)

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

def plotBeamComponent(self, valuesToPlotList, plotOnMesh=False):
    outAxis = []
    outFig = []
    resultsDictionary = self.resultsInformation.resultsDictionary
    for valueToPlot in valuesToPlotList:
        theTitle=valueToPlot
        resultToPlot = resultsDictionary[valueToPlot]
        myZ = resultToPlot.z*resultToPlot.resultScale
        fig, axOut = plotSchnittValues(self,theTitle, resultToPlot.x,resultToPlot.y,myZ,plotOnMesh)
        self.axes[valueToPlot] = valueToPlot
        outAxis.append(axOut)
    return outFig, outAxis

def plotSchnittValues(self,theTitle, x,y,z,plotOnMesh):
    maxVal = np.max(self.plates[0].outlineCoords)
    if not plotOnMesh:
        fig,outAx = plotInputGeometry(self)
    else:
        fig, outAx = plotMesh(self, plotNodes=False, plotStrucElements=False, plotPoints =True)

    iMMax = np.argmax(np.abs(z))
    # iMMin = np.argmin(z)
    magValue = 0.2*maxVal/np.max(np.abs(z))
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

    zMaxString = '{:.3f}'.format(z[iMMax])
    if np.abs(z[iMMax])>0.1:
        outAx.text(zPoints[iMMax,0],zPoints[iMMax,1], zMaxString,color='k', bbox=dict(facecolor='w', edgecolor='grey'), zorder=1000)

    outAx.set_title(theTitle)
    return fig,outAx




