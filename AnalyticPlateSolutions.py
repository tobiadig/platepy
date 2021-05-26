''' Module Information
-----------------------------------------------------------
Purpose of module: Provide analytic plate bending solutions
-----------------------------------------------------------
- Solutions stem from Timoshenko's book: "Plates and Shells"
- Copywrite Tobia Diggelmann (ETH Zurich) and Adrian Egger (Cubus AG), 13.03.2021
'''

#%% Import only necessary modules
import numpy as np
import warnings                                                                #for unimplemented cases
from collections import namedtuple                                             #for generating test data if __name__ = __main__

#%% Driver dispatching to the correct subroutine
def AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos):
    ''' Input/Output descriptions
    Parameters
    ----------
    pOpts : plate options (e.g. dataclass or named tuple)
        pOpts.shape    = "rectangular" | "circular"
        pOpts.depth    = "thick" | "thin"
        pOpts.support  = "simplySupported" | "clamped"
        pOpts.geometry = list of pertinent parameters sufficient to describe geometry
        pOpts.material = list of pertinent parameters sufficient to describe material
    lOpts : load options (e.g. dataclass or named tuple)
        lOpts.type      = "concentrated" | "distributed"
        lOpts.position  = list of pertinent parameters sufficient to describe position
        lOpts.magnitude = magnitude of vertical force
    sOpts : solution options (e.g. dataclass or named tuple)
        sOpts.nTerms = list (e.g. describing the amount of series expansion terms requested)
    inPos : positions at which output quantities are requested (i.e. array of 2D points)

    Returns
    -------
    quantities = boolean list (of size 8) of calculated and therefore returned outputs (possible: Wz, Rx, Ry, Mx, Mx, Mxy, Vx, Vy)
    values     = array (of size nReturnedQuantities x length(inPos) ) containing quanties calculated at inPos
    outPos     = positions where output quantities are calculated (i.e. array of 2D points)
    '''
    
    #compute D : The plate bending constant
    # E = pOpts.material[0] ;nu = pOpts.material[1]; h = pOpts.material[2]
    # pOpts.material[3] = E*h**3/(12*(1-nu**2))                                  
    
    if   pOpts.shape == "rectangular" and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "concentrated":
        quantities, values, outPos = Rect_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "rectangular" and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "distributed":
        quantities, values, outPos = Rect_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos)
    
    elif pOpts.shape == "rectangular" and pOpts.depth == "thin"  and pOpts.support == "clamped"         and lOpts.type == "concentrated":
        quantities, values, outPos = Rect_thin_cSupport_Concntr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "rectangular" and pOpts.depth == "thin"  and pOpts.support == "clamped"         and lOpts.type == "distributed":
        quantities, values, outPos = Rect_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos)
    
    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "concentrated":
        quantities, values, outPos = Circ_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "distributed":
        quantities, values, outPos = Circ_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos)

    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "clamped"         and lOpts.type == "concentrated":
        quantities, values, outPos = Circ_thin_cSupport_Concntr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "clamped"         and lOpts.type == "distributed":
        quantities, values, outPos = Circ_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos)

    elif pOpts.shape == "rectangular"  and pOpts.depth == "thin"  and pOpts.support == "sSupportwBeams"  and lOpts.type == "distributed":
        quantities, values, outPos = Rect_thin_sSupportwBeams_Distr(pOpts, lOpts, sOpts, inPos)
    
    else:
        raise ValueError('Requested combination for analytic solution either not feasible or implemented yet')    
    
    return quantities, values, outPos

#%% individual analytic solution

class Material:
    def __init__(self, E, nu, h):
        self.E = E
        self.nu = nu
        self.h = h
        self.D = E*h**3/(12*(1-nu**2))
        self.myLambda = None

class POpts:
    def __init__(self):
        self.shape = None
        self.depth = None
        self.support = None
        self.geometry = None
        self.material = None

class LOpts:
    def __init__(self):
        self.type = None
        self.position = None
        self.magnitude = None

class SOpts:
    def __init__(self, nTerms = 20):
        self.nTerms = nTerms


def Rect_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    print('simplySupported')
    quantities=[True, False, False, True, True, False, False, False]
    #Actually: Wz, Rx, Ry, Mx, Mx, Mxy, Vx, Vy)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    nOutputs = 5
    x = inPos[:,0]
    y = inPos[:,1]
    q0 = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    nTerms = sOpts.nTerms

    nPoints = len(x)
    v1=np.ones(nPoints)
    x = x+ v1*a/2

    # x += v1*a/2
    # print('x',x)
    # print('y',y)
    # print('q0',q0)
    # print('D',D)
    # print('nu',nu)
    # print('a',a)
    # print('b',b)
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,2, nOutputs))
    w = np.zeros(nPoints)
    Mxy0Factor = np.zeros(nPoints)
    Myy0Factor = np.zeros(nPoints)
    QxFactor = np.zeros(nPoints)
    QyFactor = np.zeros(nPoints)

    for m in range(1,nTerms*2+1,2):
        alfa = m*np.pi*b/(2*a)
        Am = -2*(alfa*np.tanh(alfa)+2)/((np.pi*m)**5*np.cosh(alfa))
        Bm = 2/((np.pi*m)**5*np.cosh(alfa))

        Mxy0Factor += m**2*(2*nu*Bm-(1-nu)*Am)*np.sin(m*np.pi*x/a)
        Myy0Factor += m**2*(2*Bm+(1-nu)*Am)*np.sin(m*np.pi*x/a) 
        QxFactor += m**3*Bm*np.cosh(m*np.pi*y/a)*np.cos(m*np.pi*x/a) 
        QyFactor += m**3*Bm*np.sinh(m*np.pi*y/a)*np.sin(m*np.pi*x/a)
        w+= (v1*4/(np.pi**5*m**5)+Am*np.cosh(m*np.pi*y/a)+Bm*np.pi*m*y/a*np.sinh(m*np.pi*y/a))*np.sin(m*np.pi*x/a)

    Mxy0 = q0*x*(a*v1-x)/2 - q0*a**2*np.pi**2*Mxy0Factor
    Myy0 = nu*q0*x*(a*v1-x)/2 - q0*a**2*np.pi**2*Myy0Factor
    Qx = q0*(a*v1-2*x)/2 - 2*np.pi**3*q0*a*QxFactor
    Qy = -2*np.pi**3*q0*a*QyFactor
    w=w*q0*a**4/D

    values[:,0] = w
    outPos[:,:,0] = inPos
    values[:,1] = Mxy0
    outPos[:,0,1] = inPos[:,0]
    values[:,2] = Myy0
    outPos[:,0,2] = inPos[:,0]
    values[:,3] = Qx
    outPos[:,:,3] = inPos
    values[:,4] = Qy  
    outPos[:,:,4] = inPos

    return quantities, values, outPos


def Rect_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos):
# Rectangular_SimplySupported_Concentrated(x,y,r,s,P,D,nu,a,b,outQuantities,xTerms=10, yTerms=10): #2)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    #Actually: Wz, Rx, Ry, Mx, Mx, Mxy, Vx, Vy)
    quantities=[True, False, False, True, True, True, False, False]
    nOutputs = 4
    x = inPos[:,0]
    y = inPos[:,1]
    r = lOpts.position[0]
    P = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    nTerms = sOpts.nTerms
    nPoints = len(x)
    v1=np.ones(nPoints)
    values = np.zeros((nPoints, nOutputs))
    outPos=np.zeros((nPoints,2,nOutputs))
    w = np.zeros(nPoints)
    Mx = np.zeros(nPoints)
    My = np.zeros(nPoints)
    Mxy = np.zeros(nPoints)
    y1=b*v1-y

    for m in range(1,nTerms*2+1,2):
        beta = m*np.pi*b/(a) # io
        temp=(v1+beta/np.tanh(beta)*v1-beta*y1/b/np.tanh(beta*y1/b)-beta*s/b/np.tanh(beta*s/b)*v1)*np.sinh(beta*s/b)*np.sinh(beta*y1/b)*np.sin(m*np.pi*r/a)*np.sin(m*np.pi*x/a)/(m**3*np.sinh(beta))
        w+=temp
        # temp =[v1+nu*v1+(1-nu)*m*np.pi*y/a]
        # print(temp.shape)
        Mx += m**(-1)*np.sin(m*np.pi*r/a)*np.sin(m*np.pi*x/a)*(v1+nu*v1+(1-nu)*m*np.pi*y/a)*np.exp(-m*np.pi*y/a)
        My += 1/m*np.sin(m*np.pi*r/a)*np.sin(m*np.pi*x/a)*(v1+nu*v1-(1-nu)*m*np.pi*y/a)*np.exp(-m*np.pi*y/a)
        Mxy += np.sin(m*np.pi*r/a)*np.cos(m*np.pi*x/a)*np.exp(-m*np.pi*y/a)

    values[:,0]=w*P*a**2/(np.pi**3*D)
    outPos[:,:,0] = inPos
    values[:,1]=Mx*P/2/a
    outPos[:,:,1] = inPos
    values[:,2]=My*P/2/a
    outPos[:,:,2] = inPos
    values[:,3]=Mxy*-P/(2*a)*y*(1-nu)
    outPos[:,:,3] = inPos

    return quantities, values, outPos



# ----------
# pOpts : plate options (e.g. dataclass or named tuple)
#     pOpts.shape    = "rectangular" | "circular"
#     pOpts.depth    = "thick" | "thin"
#     pOpts.support  = "simplySupported" | "clamped"
#     pOpts.geometry = list of pertinent parameters sufficient to describe geometry


#     pOpts.material = list of pertinent parameters sufficient to describe material
#                         E = material[0]
# nu = material[1]
# h = material[2]
# D = E*h**3/(12*(1-nu**2))

# lOpts : load options (e.g. dataclass or named tuple)
#     lOpts.type      = "concentrated" | "distributed"
#     lOpts.position  = list of pertinent parameters sufficient to describe position

#     lOpts.magnitude = magnitude of vertical force
# sOpts : solution options (e.g. dataclass or named tuple)
#     sOpts.nTerms = list (e.g. describing the amount of series expansion terms requested)
# inPos : positions at which output quantities are requested (i.e. array of 2D points)

# Returns
# -------
# quantities = boolean list (of size 8) of calculated and therefore returned outputs (possible: Wz, Rx, Ry, Mx, Mx, Mxy, Vx, Vy)
# values     = array (of size nReturnedQuantities x length(inPos) ) containing quanties calculated at inPos
# outPos     = positions where output quantities are calculated (i.e. array of 2D points)
# '''

def Rect_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    print('clamped')

    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    #Actually: Wz, Rx, Ry, Mx, Mx, Mxy, Vx, Vy)
    quantities=[True, False, False, False, False, False, False, False]
    nOutputs = 1
    x = inPos[:,0]
    y = inPos[:,1]
    P = lOpts.magnitude
    D = pOpts.material.D
    print('D: ', D)
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    print('a: ', a, ' and b: ', b)
    nTerms = sOpts.nTerms
    nPoints = len(x)
    v1=np.ones(nPoints)
    values = np.zeros((nPoints, nOutputs))
    outPos=np.zeros((nPoints,2,nOutputs))
    wClamped = np.zeros(nPoints)
    Mx = np.zeros(nPoints)
    My = np.zeros(nPoints)
    Mxy = np.zeros(nPoints)
    y1=b*v1-y

    for m in range(1,nTerms*2+1,2):
        alpha = m*np.pi*b/(2*a)
        wClamped = wClamped + (-1)**((m-1)/2)/m**5*np.cos(m*np.pi*x/a)*(v1-(alpha*np.tanh(alpha)+2)/(2*np.cosh(alpha))*np.cosh(m*np.pi*y/a)+v1*(m*np.pi*y)/(2*np.cosh(alpha)*a)*np.sinh(m*np.pi*y/a))
    wClamped=wClamped*4*P*a**4/(np.pi**5*D)
    print('w: ', wClamped)

    values[:,0]=wClamped
    outPos[:,:,0] = inPos

    return quantities, values, outPos


def Rect_thin_cSupport_Concntr(pOpts, lOpts, sOpts, inPos):
    pass
def Circ_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    pass
def Circ_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos):
    pass
def Circ_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    pass
def Circ_thin_cSupport_Concntr(pOpts, lOpts, sOpts, inPos):
    pass
def Rect_thin_sSupportwBeams_Distr(pOpts, lOpts, sOpts, inPos):

    quantities=[True, False, False, False, False, False, False, False]
    #Actually: Wz, Rx, Ry, Mx, Mx, Mxy, Vx, Vy)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    nOutputs = 1
    x = inPos[:,0]
    y = inPos[:,1]
    q0 = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    myLambda = pOpts.material.myLambda
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    nTerms = sOpts.nTerms
    nPoints = len(x)

    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,2, nOutputs))
    w = np.zeros(nPoints)

    v1=np.ones(nPoints)
    for m in range(1,nTerms*2+1,2):
        alpham = m*np.pi*b/(2*a)
        Am = 4/(m**5*np.pi**5)*\
            (nu*(1+nu)*np.sinh(alpham)-nu*(1-nu)*alpham*np.cosh(alpham)-m*np.pi*myLambda*(2*np.cosh(alpham) +alpham*np.sinh(alpham)))/ \
                ((3+nu)*(1-nu)*np.sinh(alpham)*np.cosh(alpham)-(1-nu)**2*alpham+2*m*np.pi*myLambda*np.cosh(alpham)**2)

        Bm = 4/(m**5*np.pi**5)*nu*(1-nu)*np.sinh(alpham)+m*np.pi*myLambda*np.cosh(alpham)/((3+nu)*(1-nu)*np.sinh(alpham)*np.cosh(alpham)-(1-nu)**2*alpham+2*m*np.pi*myLambda*np.cosh(alpham)**2)

        w +=q0*a**4/D* ((4/(m**5*np.pi**5)*v1 +Am*np.cosh(m*np.pi*y/a) + Bm*m*np.pi*y/a*np.sinh(m*np.pi*y/a))*np.sin(m*np.pi*x/a))

    values[:,0] = w
    outPos[:,:,0] = inPos
    
    return quantities, values, outPos

#%% for unit testing! This code is executed if this python file is run!
if __name__ == "__main__":
    #provide some demo input and verify output
    PlateOptions = namedtuple('PlateOptions', 'shape depth support geometry material')
    pOpts = PlateOptions("rectangular", "thin", "simplySupported", [10, 10], [30e9, 0.3, 0.2, 22e6])
    
    LoadOptions = namedtuple('LoadOptions', 'type position magnitude')
    lOpts = LoadOptions("distributed", [5, 5], 1e3)

    SolverOptions = namedtuple('SolverOptions', 'nTerms')
    sOpts = SolverOptions([10, 10])
    
    inPos = np.random.rand(5,2)

    quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
    print(quantities); print(values); print(outPos)
    
#%% Old code base below   
#------------------------------------------------------------------------------               

def Rectangular_SimplySupported_Distributed(x,y,q0,D,nu,a,b,outQuantities,xTerms=10, yTerms=10): #1) 
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    # check the inputs
    if (outQuantities[1] or outQuantities[2]):
        warnings.warn('The values of Mx and My are only available at the center of the plate in y-direction ( y = b/2). Passed y-values are ignored')

    if outQuantities[3]:
        raise TypeError('Mxy value not available')

    nOutputs = np.sum(outQuantities)
    nPoints = len(x)
    v1=np.ones(nPoints)
    x = x+ v1*a/2

    # x += v1*a/2
    # print('x',x)
    # print('y',y)
    # print('q0',q0)
    # print('D',D)
    # print('nu',nu)
    # print('a',a)
    # print('b',b)
    result = np.zeros((nPoints, outQuantities.size))
    w = np.zeros(nPoints)
    Mxy0Factor = np.zeros(nPoints)
    Myy0Factor = np.zeros(nPoints)
    QxFactor = np.zeros(nPoints)
    QyFactor = np.zeros(nPoints)

    for m in range(1,xTerms*2+1,2):
        
        Am = -2*(alfa*np.tanh(alfa)+2)/((np.pi*m)**5*np.cosh(alfa))
        Bm = 2/((np.pi*m)**5*np.cosh(alfa))

        Mxy0Factor += m**2*(2*nu*Bm-(1-nu)*Am)*np.sin(m*np.pi*x/a)
        Myy0Factor += m**2*(2*Bm+(1-nu)*Am)*np.sin(m*np.pi*x/a) 
        QxFactor += m**3*Bm*np.cosh(m*np.pi*y/a)*np.cos(m*np.pi*x/a) 
        QyFactor += m**3*Bm*np.sinh(m*np.pi*y/a)*np.sin(m*np.pi*x/a)
        w+= (v1*4/(np.pi**5*m**5)+Am*np.cosh(m*np.pi*y/a)+Bm*np.pi*m*y/a*np.sinh(m*np.pi*y/a))*np.sin(m*np.pi*x/a)

    Mxy0 = q0*x*(a*v1-x)/2 - q0*a**2*np.pi**2*Mxy0Factor
    Myy0 = nu*q0*x*(a*v1-x)/2 - q0*a**2*np.pi**2*Myy0Factor
    Qx = q0*(a*v1-2*x)/2 - 2*np.pi**3*q0*a*QxFactor
    Qy = -2*np.pi**3*q0*a*QyFactor
    w=w*q0*a**4/D

    result[:,0] = w
    result[:,1] = Mxy0
    result[:,2] = Myy0
    result[:,4] = Qx
    result[:,5] = Qy  

    return result[:,outQuantities]

def Rectangular_SimplySupported_Concentrated(x,y,r,s,P,D,nu,a,b,outQuantities,xTerms=10, yTerms=10): #2)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    # check the inputs
    if outQuantities[4] or outQuantities[5]:
        raise TypeError('Q values are not available')
    if b/a<10:
        warnings.warn('This theory is appliable only for very long inputs, check if this applies to your case')
    if outQuantities[1] or outQuantities[2] or outQuantities[3]:
        warnings.warn('The moments can be calculated only by \eta=0, eta values will be ignored')
        warnings.warn('Moments in points too close to the concentrated load could be imprecise, consider increasing the number of order')
    nOutputs = np.sum(outQuantities)
    nPoints = len(x)
    v1=np.ones(nPoints)
    result = np.zeros((nPoints, outQuantities.size))
    w = np.zeros(nPoints)
    Mx = np.zeros(nPoints)
    My = np.zeros(nPoints)
    Mxy = np.zeros(nPoints)
    y1=b*v1-y

    for m in range(1,xTerms*2+1,2):
        beta = m*np.pi*b/(a) # io
        temp=(v1+beta/np.tanh(beta)*v1-beta*y1/b/np.tanh(beta*y1/b)-beta*s/b/np.tanh(beta*s/b)*v1)*np.sinh(beta*s/b)*np.sinh(beta*y1/b)*np.sin(m*np.pi*r/a)*np.sin(m*np.pi*x/a)/(m**3*np.sinh(beta))
        w+=temp
        # temp =[v1+nu*v1+(1-nu)*m*np.pi*y/a]
        # print(temp.shape)
        Mx += m**(-1)*np.sin(m*np.pi*r/a)*np.sin(m*np.pi*x/a)*(v1+nu*v1+(1-nu)*m*np.pi*y/a)*np.exp(-m*np.pi*y/a)
        My += 1/m*np.sin(m*np.pi*r/a)*np.sin(m*np.pi*x/a)*(v1+nu*v1-(1-nu)*m*np.pi*y/a)*np.exp(-m*np.pi*y/a)
        Mxy += np.sin(m*np.pi*r/a)*np.cos(m*np.pi*x/a)*np.exp(-m*np.pi*y/a)

    result[:,0]=w*P*a**2/(np.pi**3*D)
    result[:,1]=Mx*P/2/a
    result[:,2]=My*P/2/a
    result[:,3]=Mxy*-P/(2*a)*y*(1-nu)

    return result[:,outQuantities]

#rectangular plates with all edges built in, distributed load p .197
def Rectangular_Clamped_Distributed(x,y,q,D,nu,a,b,outQuantities,xTerms=10): #3)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    wTemp = 0
    for m in range(1,xTerms*2+1,2):
        alfa = m*np.pi*b/(2*a)
        wTemp = wTemp + (-1)**((m-1)/2)/m**5*np.cos(m*np.pi*x/a)*(1-(alpha*np.tanh(alpha)+2)/(2*np.cosh(alpha)*np.cosh(m*np.pi.y/a)+(m*np.pi*y)/(2*np.cosh(alpha)*a)*np.sinh(m*np.pi*y/a)))

    w=wTemp*4*q*a**4/(np.pi**5*D)


    

def Rectangular_Clamped_Concentrated(x,y,r,s,P,D,nu,a,b,outQuantities,xTerms=100, yTerms=10):  #4)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    # check the inputs
    if outQuantities[1] or  outQuantities[2] or outQuantities[3] or outQuantities[4] or outQuantities[5]:
        raise TypeError('for this plate only w is available')

    nOutputs = np.sum(outQuantities)
    nPoints = len(x)
    v1=np.ones(nPoints)
    result = np.zeros((nPoints, outQuantities.size))
    w=np.zeros(nPoints)
    for m in range(1,xTerms*2+1,2):
        alfa = m*np.pi*b/(2*a)
        w+=m**(-3)*np.cos(m*np.pi*x/a)*((np.tanh(alfa)-alfa/(np.cosh(alfa)**2))* \
            np.cosh(m*np.pi*y)-np.sinh(m*np.pi*y/a)- \
                m*np.pi*y/a*np.tanh(alfa)*np.sinh(m*np.pi*y/a)+m*np.pi*y/a*np.cosh(m*np.pi*y/a))

    result[:,0]=w*P*a**2/(2*np.pi**3*D)
    return result[:,outQuantities]

def Rectangular_sAtEdges_Distributed(x,y,q,D,nu,a,b,outQuantities, xTerms =10): #5)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    return TypeError('td: to be implemented')

def Circular_SimplySupported_Distributed(r,q,D,nu,a,outQuantities): #6)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    # check the inputs
    if outQuantities[5] or  outQuantities[3]:
        raise TypeError('Mxy and Qy values are not available')

    if outQuantities[1] or  outQuantities[2]:
        warnings.warn('Mx and My correspond to the moment in radial and tangential directions.')

    nOutputs = np.sum(outQuantities)
    nPoints = len(r)
    v1=np.ones(nPoints)
    result = np.zeros((nPoints, outQuantities.size))

    result[:,4]=q*r/2
    result[:,0]=q*(v1*a**2-r**2)/(64*D)*((5+nu)/(1+nu)*(v1*a**2-r**2))
    result[:,1] = q/16*(3+nu)*(v1*a**2-r**2)
    result[:,2] = q/16*(v1*a**2*(3+nu)-r**2*(1+nu*3))
    return result[:,outQuantities]

def Circular_SimplySupported_Concentrated(r,xi,P,D,nu,a,outQuantities): #7)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    # check the inputs
    if outQuantities[5] or  outQuantities[3]:
        raise TypeError('Mxy and Qy values are not available')

    if not(xi==0):
        warnings.warn('The concentrated load is applied in the center of the plate, xi value will be ignored')
    nOutputs = np.sum(outQuantities)
    nPoints = len(r)
    v1=np.ones(nPoints)
    result = np.zeros((nPoints, outQuantities.size))

    result[:,4]=P*r/2
    result[:,0]= P/(16*np.pi*D)*((3+nu)/(1+nu)*(v1*a**2-r**2)+2*r**2*np.log10(r/a))
    result[:,1] = P/(4*np.pi)*(1+nu)*np.log10(a/r)
    result[:,2] = P/(4*np.pi)*((1+nu)*np.log10(a/r)+1-nu)
    return result[:,outQuantities]

def Circular_Clamped_Distributed(r,q,D,nu,a,outQuantities): #8)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    # check the inputs
    if outQuantities[5] or  outQuantities[3]:
        raise TypeError('Mxy and Qy values are not available')
    if outQuantities[1] or  outQuantities[2]:
        warnings.warn('Mx and My correspond to the moment in radial and tangential directions.')

    nOutputs = np.sum(outQuantities)
    nPoints = len(r)
    v1=np.ones(nPoints)
    result = np.zeros((nPoints, outQuantities.size))

    result[:,4]=q*r/2
    result[:,0]=q/(64*D)*(v1*a**2-r**2)**2
    result[:,1] = q/16*(a**2*(1-nu)-r**2*(3+nu))
    result[:,2] = q/16*(a**2*(1-nu)-r**2*(1+nu*3))
    return result[:,outQuantities]

def Circular_Clamped_Concentrated(r,theta,xi,P,D,nu,a,outQuantities): #9)
    # outQuantities: [w, Mx, My, Mxy, Qx, Qy]
    notJustW = (outQuantities[1] or  outQuantities[2] or outQuantities[3] or outQuantities[4] or outQuantities[5])
    nOutputs = np.sum(outQuantities)
    nPoints = len(r)
    v1=np.ones(nPoints)
    result = np.zeros((nPoints, outQuantities.size))
    x=r/a

    if not(notJustW):
        result[:,0] = P*a**2/(16*np.pi*D)*((v1-x**2)*(v1-xi**2)+(x**2+xi**2-2*x*xi*np.cos(theta))*np.log10((x**2+xi**2-2*x*xi*np.cos(theta))/(v1+x**2*xi**2-2*x*xi*np.cos(theta))))
    elif notJustW and xi:
        result[:,4]=q*r/2
        result[:,0]= P*r**2/(8*np.pi*D)*np.log10(r/a)+P/(16*np.pi*D)*(v1*a**2-r**2)
        result[:,1] = P/(4*np.pi)*((1+nu)*np.log10(a/r)-v1)
        result[:,2] = P/(4*np.pi)*((1+nu)*np.log10(a/r)-nu*v1)
    else:
        raise TypeError('By placing the concentrated load outside the center only the vertical displacement can be calculated. Put xi = 0 or just call for w.')

    return result[:,outQuantities]

def TimshenkoPlateSolution(points,plateType,geometry, material,loadCond,loadVal, outQuantities, \
                            supportCond = 'sSupported', thin = True, loadCoordinates = True ):
    #define common parameters for all calculation
    E = material[0]
    nu = material[1]
    h = material[2]
    D = E*h**3/(12*(1-nu**2))

    if plateType == 'rectangular':
        #check inputs
        if  len(points)<2:
            raise IndexError('For rectangular plate both x and y coords required')
        if len(geometry)<2:
            raise IndexError('For rectangular plate both a and b dimensions required')
        
        x = points[0]
        y = points[1]
        a = geometry[0]
        b = geometry[1]

        if supportCond == 'sSupported':
            if loadCond == 'distributed':
                return Rectangular_SimplySupported_Distributed(x,y,loadVal,D,nu,a,b,outQuantities)
            
            elif loadCond == 'concentrated':
                #check inputs
                if type(loadCoordinates)==bool:
                    raise TypeError('coordinates of the concentrated load missing')

                r=loadCoordinates[0]
                s=loadCoordinates[1]

                return Rectangular_SimplySupported_Concentrated(x,y,r,s,loadVal,D,nu,a,b,outQuantities)

            else:
                raise ValueError('Load condition must be either "distributed" or "concentrated"')
        
        elif supportCond == 'clamped':
            
            if loadCond == 'distributed':
                return Rectangular_Clamped_Distributed(x,y,loadVal,D,nu,a,b,outQuantities)
        
            elif loadCond == 'concentrated':
                #check inputs
                if type(loadCoordinates)==bool:
                    raise TypeError('coordinates of the concentrated load missing')

                r=loadCoordinates[0]
                s=loadCoordinates[1]

                return Rectangular_Clamped_Concentrated(x,y,r,s,loadVal,D,nu,a,b,outQuantities)
            else:
                raise ValueError('Load condition must be either "distributed" or "concentrated"')
        
        elif supportCond == 'sAtCorners':

            return Rectangular_sAtEdges_Distributed(x,y,loadVal,D,nu,a,b,outQuantities)

        else:
            raise ValueError('Support condition must be either "clamped", "sAtCorners" or "sSupported" (default)')


    elif plateType == 'circular':
        #define inputs common for circular plates
        #check inputs
        r=points
        a=geometry

        if supportCond == 'sSupported':
            if loadCond == 'distributed':
                return Circular_SimplySupported_Distributed(r,loadVal,D,nu,a,outQuantities)
        
            elif loadCond == 'concentrated':
                xi=loadCoordinates
                return Circular_SimplySupported_Concentrated(r,xi,loadVal,D,nu,a,outQuantities)

            else:
                raise ValueError('Load condition must be either "distributed" or "concentrated"')
        
        elif supportCond == 'clamped':
            if loadCond == 'distributed':
                return Circular_Clamped_Distributed(r,loadVal,D,nu,a,outQuantities)
            
            elif loadCond == 'concentrated':
                theta=loadCoordinates[1]
                xi=loadCoordinates[0]
                return Circular_Clamped_Concentrated(r,theta,xi,loadVal,D,nu,a,outQuantities)

            else:
                ValueError('Load condition must be either "distributed" or "concentrated"')

        else:
            raise ValueError('Support condition must be either "clamped" or "sSupported" (default)')

    else:
        raise ValueError('Plate type must be either "rectangular" or "circular"')
        
        