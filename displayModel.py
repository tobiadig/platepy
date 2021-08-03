'''
<a name="displayModel"></a> Displays geometry and results of a plateModel class using `matplotlib`

Copywrite Tobia Diggelmann (ETH Zurich) 09.06.2021
'''
# Basic modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import base64
from io import BytesIO
from tqdm import tqdm

def plotInputGeometry(self, figaspect = 1):
    '''
    Display the structural elements contained in a `PlateModel` object.

    ~~~~~~~~~~~~~~~~~~~
    INPUT
    ~~~~~~~~~~~~~~~~~~~

    *   **self**: `PlateModel` object with initialized structural components.
    *   **figaspect = 1**: aspect ratio of the output figure (width:high).

    ~~~~~~~~~~~~~~~~~~~
    RETURN
    ~~~~~~~~~~~~~~~~~~~
    '''
    w, h = plt.figaspect(figaspect)
    mult = 1.5
    # fig = plt.figure()
    # axGeometry = plt.axes()
    fig, axGeometry = plt.subplots(figsize=(w*mult, h*mult))

    for plate in self.plates:
        plate._plot(axGeometry)

    for wall in self.walls:
        wall._plot(axGeometry)

    for uz in self.downStandBeams:
        uz._plot(axGeometry)

    for column in self.columns:
        column._plot(axGeometry)

    self.axes['InputGeometry'] = axGeometry
    return fig,axGeometry

def plotResults(self, valuesToPlotList, plotType = 'isolines',saveToSVG=False, saveBase64=False):
    '''
        Displays the requested result over the model's structural components.

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **self**: Solved `PlateModel` object with results stored.
        * **valuesToPlotList**: List of strings representing the requested output plots. Possible values are:
            - "vDisp": Vertical displacements.
            - "Mx": Bending moments in x-direction.
            - "My": Bending moments in y-direction.
            - "Mxy": Twisting moments.
            - "Vx": Shear forces in x-direction.
            - "Vy": Shear forces in y-direction.
            - "Mx_Rd_+": Required positive bending resistance in x-direction.
            - "My_Rd_+": Required positive bending resistance in y-direction.
            - "Mx_Rd_-": Required negative bending resistance in x-direction.
            - "My_Rd_-": Required negative bending resistance in y-direction.
            - "Mx_Rd_max": Required maximum bending resistance in x-direction.
            - "My_Rd_max": Required maximum bending resistance in y-direction.

        * **plotType = "isolines"**: String representing the type of plot. Possible values are:

            - "isolines"(default): Lines with the same output-value are extrapolated and plotted.
            - "3d": 3-dimensional plot.
            - "text": The values at the data-points are printed as text near to the data-points.
            - "text+mesh": The values at the data-points are printed as text near to the data-points on the model's mesh.

        * **saveToSVG = False**: if True, saves the plot as SVG in the current folder.
        * **saveBase64 = False**: Experimental.

        ~~~~~~~~~~~~~~~~~~~
        RETURN
        ~~~~~~~~~~~~~~~~~~~

        *   **outFig**, **outAxis**
    '''
    outAxis = []
    outFig = []
    resultsDictionary = self.resultsInformation.resultsDictionary
    for valueToPlot in valuesToPlotList:
        theTitle=valueToPlot
        resultToPlot = resultsDictionary[valueToPlot]
        myZ = resultToPlot.z*resultToPlot.resultScale

        fig, axOut = _plotInternalForce(self,plotType,theTitle, resultToPlot.x,resultToPlot.y,myZ,saveBase64)
        self.axes[theTitle] = valueToPlot
        outAxis.append(axOut)
        outFig.append(fig)
        if saveToSVG:
            fig.savefig(theTitle+'.svg')
    return outFig, outAxis

def plotMesh(self, plotNodes = True, plotStrucElements = True, plotPoints = False, saveBase64=False):
    '''
    Plot nodes and elements given an initialized mesh.

    ~~~~~~~~~~~~~~~~~~~
    INPUT
    ~~~~~~~~~~~~~~~~~~~

    * **self**: `PlateModel` object with initialized mesh.
    * **plotNodes = True**: If True, plots the nodes with numeration.
    * **plotStrucElements = True**: If True, also plots the underlying structural components.
    * **plotPoints = False**: If True, plots the nodes without numeration.

    ~~~~~~~~~~~~~~~~~~~
    RETURN
    ~~~~~~~~~~~~~~~~~~~

    * **fig**, **outAx**
    '''
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
    outAx.set_xticks([])
    outAx.set_yticks([])

    if saveBase64:
        buf = BytesIO()
        fig.savefig(buf, format="png")
        data5 = base64.b64encode(buf.getbuffer()).decode("ascii")
        fig = data5

    return fig, outAx

def plotBeamComponent(self, valuesToPlotList, plotOnMesh=False, saveToSVG = False):
    '''
    Display the requested line-results over the model structural components.

    ~~~~~~~~~~~~~~~~~~~
    INPUT
    ~~~~~~~~~~~~~~~~~~~

    * **self**: Solved `PlateModel` object with line-results stored.
    * **valuesToPlotList**: List of strings representing the requested output plots. Possible values are:

        - "vDisp_line": Vertical displacements on a line.
        - "Mx_line": Bending moments in x-direction on a line.
        - "My_line": Bending moments in y-direction on a line.
        - "Mxy_line": Twisting moment on a lines.
        - "Vx_line": Shear forces in x-direction on a line.
        - "Vy_line": Shear forces in y-direction on a line.

    * **plotOnMesh = False**: If True, the line results are plotted over the underlying mesh.
    * **saveToSVG = False**: if True, saves the plot as SVG in the current folder.
    ~~~~~~~~~~~~~~~~~~~
    RETURN
    ~~~~~~~~~~~~~~~~~~~

    *   **outFig**, **outAxis**
    '''
    outAxis = []
    outFig = []
    resultsDictionary = self.resultsInformation.resultsDictionary
    for valueToPlot in valuesToPlotList:
        if _valuesAreNotBeamValues(valueToPlot):
            raise Exception('add "_line" at the end of the plot string!')
        theTitle=valueToPlot
        resultToPlot = resultsDictionary[valueToPlot]
        myZ = resultToPlot.z*resultToPlot.resultScale
        fig, axOut = _plotSchnittValues(self,theTitle, resultToPlot.x,resultToPlot.y,myZ,plotOnMesh)
        self.axes[valueToPlot] = valueToPlot
        outAxis.append(axOut)
        if saveToSVG:
            fig.savefig(theTitle+'.svg')

    # Hide axes ticks
    axOut.set_xticks([])
    axOut.set_yticks([])

    return outFig, outAxis

def _plotInternalForce(self,plotType,theTitle, x,y,z,saveBase64):
    '''
        Plot the given z values at the x-y coordinates. \n
        Input: \n
        * self: plateModel object with stored geometry.\n
        * plotType: String representing the type of plot. \n
        * theTitle: String to be displayed as title, by default is the same string as in plotType. \n
        * x, y: Coordinates of the data points. \n
        * z: Values at the data points. \n
        * saveBase64: experimental.\n
        Return: \n
        * fig, axOut: Figure and axis with the requested plot.
    '''
    if plotType == 'isolines':
        fig,axOut = _myIsoPlot(self,x,y,z,theTitle=theTitle)
        if saveBase64:
            buf = BytesIO()
            fig.savefig(buf, format="png")
            data5 = base64.b64encode(buf.getbuffer()).decode("ascii")
            fig = data5
            # plt.savefig(r'C:\Users\Diggelmann\Desktop\FEMFlask\static\images\new_plot.png')
            # outFig.append(data5)
    elif plotType == '3d':
        fig=plt.figure()
        percentageTrianglesToMantain = 99
        triang = _removeTrianglesOutsidePlate(x,y,percentageTrianglesToMantain)
        axOut = fig.gca(projection='3d')
        axOut.plot_trisurf(triang,z,cmap=plt.get_cmap('jet'))
        # axOut.plot_trisurf(triang,z,cmap=plt.get_cmap('winter_r'))
        # axOut.plot_trisurf(triang,z,cmap=plt.get_cmap('viridis'))
        axOut.pbaspect = [1, 1, 0.5]
        axOut.grid(False)
        # Hide axes ticks
        axOut.set_xticks([])
        axOut.set_yticks([])
        axOut.set_zticks([])
    elif plotType == 'text':
        fig, axOut = _myTextPlot(self, x,y,z,theTitle=theTitle)
    elif plotType == 'text+mesh':
        fig, axOut = _myTextOnMeshPlot(self,x,y,z, theTitle = theTitle)
    else:
        raise TypeError('type of plot does not exist')

    return fig, axOut

def _myTextPlot(self,x,y,z, theTitle = ''):
    '''
        The values at the data-points are printed as text near to the data-points. \n
        Input: \n
        * self: plateModel object with stored geometry.\n
        * x, y: Coordinates of the data points. \n
        * z: Values at the data points. \n
        * theTitle: String to be displayed as title.\n
        Return: \n
        * fig, outAx
    '''
    fig,outAx = plotInputGeometry(self)
    for i,a in enumerate(z):
        outAx.text(x[i], y[i], '{:.4f}'.format(a)) # number of Nachkomastellen to be displayed
        outAx.scatter(x[i], y[i], marker = ".", facecolor = "k")

    outAx.set_title(theTitle)
    return fig, outAx

def _myTextOnMeshPlot(self,x,y,z, theTitle = ''):
    '''
        The values at the data-points are printed as text near to the data-points on the overlying mesh. \n
        Input: \n
        * self: plateModel object with stored geometry.\n
        * x, y: Coordinates of the data points. \n
        * z: Values at the data points. \n
        * theTitle: String to be displayed as title.\n
        Return: \n
        * fig, outAx
    '''
    fig, outAx = plotMesh(self, plotNodes=False, plotStrucElements=False)
    for i,a in enumerate(z):
        outAx.text(x[i]+0.03, y[i]+0.03, '{:.2f}'.format(a)) # number of Nachkomastellen to be displayed
        outAx.scatter(x[i], y[i], marker = ".", facecolor = "r",zorder=0)

    outAx.set_title(theTitle)
    return fig, outAx

def _myIsoPlot(self,x,y,z, theTitle = ''):
    '''
        Input: \n
        * self: plateModel object with stored geometry.\n
        * x, y: Coordinates of the data points. \n
        * z: Values at the data points. \n
        * theTitle: String to be displayed as title.\n
        Return: \n
        * fig, outAx
    '''
    fig,outAx = plotInputGeometry(self)

    red = np.array([255/256, 45/256, 0/256, 1])
    grey = np.array([128/256, 128/256, 128/256, 1])
    blue = np.array([0/256, 55/256, 255/256, 1])
    newcolors = [red,grey,blue]
    mycmp = matplotlib.colors.ListedColormap(newcolors)
    bounds = np.array([-1e70, -0.00001,0.000001, 1e70])

    norm = matplotlib.colors.BoundaryNorm(bounds, 3)
    percentageTrianglesToMantain = 98
    triang = _removeTrianglesOutsidePlate(x,y,percentageTrianglesToMantain)
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
    # Hide axes ticks
    outAx.set_xticks([])
    outAx.set_yticks([])

    return fig,outAx

def _plotSchnittValues(self,theTitle, x,y,z,plotOnMesh):
    '''
        Plot line-result. \n
        Input: \n
        * self: plate Object. \n
        * theTitle: String to be displayed as the title of the plot.
        * x, y: Coordinates of the data points. \n
        * z: Values at the data points. \n
        * plotOnMesh = False: If True, the line results are plotted over the underlying mesh.\n
        Return: \n
        * fig,outAx
    '''
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
    if np.abs(z[iMMax])>0.0:
        outAx.text(zPoints[iMMax,0],zPoints[iMMax,1], zMaxString,color='k', bbox=dict(facecolor='w', edgecolor='grey'), zorder=1000)

    outAx.set_title(theTitle)
    return fig,outAx

def _removeTrianglesOutsidePlate(x,y,percentageTrianglesToMantain):
    '''
        Remove the larger alpha percentage of triangles.
        Input: \n
        * x, y: Coordinates of the data points. \n
        * percentageTrianglesToMantain: 100 - alpha. \n
        Return: \n
        * triang: triangles object.
    '''
    triang = tri.Triangulation(x, y)
    triangles = triang.triangles
    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
    alpha = np.percentile(maxi, percentageTrianglesToMantain)

    triang.set_mask(maxi > alpha)
    return triang

def _valuesAreNotBeamValues(valueToPlot):
    if '_line' in valueToPlot:
        return False
    else:
        return True
