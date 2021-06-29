''' Provide analytic plate bending solutions. Solutions stem from Timoshenko's book: "Plates and Shells" (1959). 

Copywrite Tobia Diggelmann (ETH Zurich) and Adrian Egger (Cubus AG), 11.06.2021
'''
import numpy as np
from scipy.interpolate import interp1d
#test
def AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos):
    '''
    Compute the analytical solution of a plate problem according to Timoshenko's book "Theory of plates and shells" (1959).

    ~~~~~~~~~~~~~~~~~~~
    INPUT
    ~~~~~~~~~~~~~~~~~~~

    * **pOpts** : object with following attributes 
        - pOpts.shape    = "rectangular" or "circular".
        - pOpts.depth    = "thick" or "thin".
        - pOpts.support  = "simplySupported" or "clamped".
        - pOpts.geometry = for rectangles: tuple with (a,b) dimensions, for circles: radius r.
        - pOpts.material = material object. 
    * **lOpts** : object with following attributes
        - lOpts.type      = "concentrated" or "distributed"
        - lOpts.position  = if the load is concentrated, tuple with (x,y) coordinates.
        - lOpts.magnitude = magnitude of vertical force.
    * **sOpts** : solution options 
        - sOpts.nTerms = Integer, number of terms in the Taylor expansion.
    * **inPos** : Numpy array with shape (n, 2) containing the positions at which output quantities are requested.

    ~~~~~~~~~~~~~~~~~~~
    RETURN
    ~~~~~~~~~~~~~~~~~~~

    * **quantities** = boolean list (of size 8) of calculated and therefore returned outputs (possible: Wz, Rx, Ry, Mx, Mx, Mxy, Vx, Vy). 
    * **values**     = Numpy array of size (length(inPos),nReturnedQuantities ) containing quantities calculated at inPos. 
    * **outPos**     = Numpy array of size (length(inPos), 2, nReturnedQuantities) containing positions where output quantities are calculated.
    '''

    if   pOpts.shape == "rectangular" and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "concentrated":
        quantities, values, outPos = _Rect_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "rectangular" and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "distributed":
        quantities, values, outPos = _Rect_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos)

    elif pOpts.shape == "rectangular" and pOpts.depth == "thin"  and pOpts.support == "clamped"         and lOpts.type == "distributed":
        quantities, values, outPos = _Rect_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos)
    
    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "concentrated":
        quantities, values, outPos = _Circ_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "simplySupported" and lOpts.type == "distributed":
        quantities, values, outPos = _Circ_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos)

    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "clamped"         and lOpts.type == "concentrated":
        quantities, values, outPos = _Circ_thin_cSupport_Concntr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "circular"    and pOpts.depth == "thin"  and pOpts.support == "clamped"         and lOpts.type == "distributed":
        quantities, values, outPos = _Circ_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos)

    elif pOpts.shape == "rectangular"  and pOpts.depth == "thin"  and pOpts.support == "sSupportwBeams"  and lOpts.type == "distributed":
        quantities, values, outPos = _Rect_thin_sSupportwBeams_Distr(pOpts, lOpts, sOpts, inPos)
    elif pOpts.shape == "rectangular"  and pOpts.depth == "thin"  and pOpts.support == "columnswBeams"  and lOpts.type == "distributed":
        w = _Rect_thin_columnswBeams_Distr(pOpts, lOpts, sOpts, inPos)
        return w
    
    else:
        raise ValueError('Requested combination for analytic solution either not feasible or implemented yet')
    return quantities, values, outPos

class Material:
    def __init__(self, E, nu, h):
        '''
        Stores material properties for the pOpts class.

        ~~~~~~~~~~~~~~~~~~~
        INPUT
        ~~~~~~~~~~~~~~~~~~~

        * **E**: Young's modulus.
        * **nu**: Poisson's ratio.
        * **h**: Plate's thickness.

        ~~~~~~~~~~~~~~~~~~~
        ATTRIBUTES
        ~~~~~~~~~~~~~~~~~~~

        * **myLambda**: For analytical models containing beam stiffeners, ratio between beam and plate stiffness (EI/aD).

        '''
        self.E = E
        self.nu = nu
        self.h = h
        self.D = E*h**3/(12*(1-nu**2))
        self.myLambda = None

class POpts:
    def __init__(self):
        '''
        Plate options.

        ~~~~~~~~~~~~~~~~~~~
        ATTRIBUTES
        ~~~~~~~~~~~~~~~~~~~

        * **pOpts.shape**    = "rectangular" or "circular".
        * **pOpts.depth**    = "thick" or "thin".
        * **pOpts.support**  = "simplySupported" or "clamped".
        * **pOpts.geometry** = for rectangles: tuple with (a,b) dimensions, for circles: radius r.
        * **pOpts.material** = material object. 

        '''
        self.shape = None
        self.depth = None
        self.support = None
        self.geometry = None
        self.material = None

class LOpts:
    def __init__(self):
        '''
        Load options.

        ~~~~~~~~~~~~~~~~~~~
        ATTRIBUTES
        ~~~~~~~~~~~~~~~~~~~

        * **lOpts.type**      = "concentrated" or "distributed"
        * **lOpts.position**  = if the load is concentrated, tuple with (x,y) coordinates.
        * **lOpts.magnitude** = magnitude of vertical force.
        
        '''
        self.type = None
        self.position = None
        self.magnitude = None

class SOpts:
    def __init__(self, nTerms = 20):
        '''
        Solution options.

        ~~~~~~~~~~~~~~~~~~~
        ATTRIBUTES
        ~~~~~~~~~~~~~~~~~~~

        **sOpts.nTerms = 20** = number of terms of the Taylor's expansions, by default 20 (usually is enough)."
        
        '''
        self.nTerms = nTerms

def _Rect_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    quantities=[True, True, True, False, True, True]
    #Actually: Wz, Mx, My, Mxy, Vx, Vy)
    nOutputs = 5
    x = inPos[:,0]
    y = inPos[:,1]
    nPoints = len(x)
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,2, nOutputs))
    q0 = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    nTerms = sOpts.nTerms
    v1=np.ones(nPoints)
    x = x+ v1*a/2
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

def _Rect_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos):
    ##Actually: Wz, Mx, My, Mxy, Vx, Vy)
    nPoints = len(inPos[:,0])
    v1=np.ones(nPoints)
    quantities=[True, True, True, True, False, False]
    nOutputs = 4
    P = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    x = inPos[:,0]+v1*0.5*a
    y = inPos[:,1]+v1*0.5*b
    r = lOpts.position[0]+v1*0.5*a
    s = lOpts.position[1]+v1*0.5*b
    nTerms = sOpts.nTerms
    

    values = np.zeros((nPoints, nOutputs))
    outPos=np.zeros((nPoints,2,nOutputs))
    w = np.zeros(nPoints)
    Mx = np.zeros(nPoints)
    My = np.zeros(nPoints)
    Mxy = np.zeros(nPoints)
    y1=b*v1-y
    for m in range(1,nTerms*2+1,2):
        beta = m*np.pi*b/(a) # io
        temp=(v1+beta/np.tanh(beta)*v1-beta*y1/b/np.tanh(beta*y1/b)-beta*s/b/np.tanh(beta*s/b)*v1)\
            *np.sinh(beta*s/b)*np.sinh(beta*y1/b)*np.sin(m*np.pi*r/a)*np.sin(m*np.pi*x/a)/(m**3*np.sinh(beta))
        w+=temp
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

def _Rect_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    quantities=[True, True, True, False, False, False]
    #Actually: Wz, Mx, My, Mxy, Vx, Vy)
    nOutputs = 3
    nPoints = 2
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,2, nOutputs))
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    timoTable = np.array([[1, 0.00126, -0.0513, -0.513, 0.0231, 0.231],
                            [1.1, 0.00150, -0.0581, -0.0538, 0.0264, 0.0231],
                            [1.2, 0.00172, -0.0639, -0.0554, 0.0299, 0.0228],
                            [1.3, 0.00191, -0.0687, -0.0563, 0.0327, 0.0222],
                            [1.4, 0.00207, -0.0726, -0.0568, 0.0349, 0.0212],
                            [1.5, 0.00220, -0.0757, -0.0570, 0.0368, 0.0203],
                            [1.6, 0.00230, -0.0780, -0.0571, 0.0381, 0.0193],
                            [1.7, 0.00238, -0.0799, -0.0571, 0.0392, 0.0182],
                            [1.8, 0.00245, -0.0812, -0.0571, 0.0401, 0.0174],
                            [1.9, 0.00249, -0.0822, -0.0571, 0.0407, 0.0165],
                            [2.0, 0.00254, -0.0829, -0.0571, 0.0412, 0.0158],
                            [1000, 0.00260, -0.0833, -0.571, 0.0417, 0.0125]])
    outPos[1,0,1] = a/2
    outPos[1, 1, 2] = b / 2
    f = interp1d(timoTable[:,0], timoTable[:,1:], axis=0)
    vals = f(b/a)
    values[0, 0] = vals[0]
    values[1,1] = vals[1]
    values[1,2] = vals[3]
    values[0,1] = vals[4]
    values[0,2] = vals[5]
    return quantities, values, outPos

def _Rect_thin_sSupportwBeams_Distr(pOpts, lOpts, sOpts, inPos):
    quantities=[True, False, False, False, False, False]
    #Actually: Wz, Mx, Mx, Mxy, Vx, Vy)
    nOutputs = 1
    x = inPos[:,0]
    y = inPos[:,1]
    nPoints =1
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,2, nOutputs))
    q0 = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    myLambda = pOpts.material.myLambda
    a = pOpts.geometry[0]
    b = pOpts.geometry[1]
    nTerms = sOpts.nTerms
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

def _Rect_thin_columnswBeams_Distr(pOpts, lOpts, sOpts, inPos):
    quantities=[True, True, True, False, False, False]
    #Actually: Wz,  Mx, Mx, Mxy, Vx, Vy)
    nOutputs = 3
    nPoints =1
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,2, nOutputs))
    myLambda = pOpts.material.myLambda
    gammas = np.array([0,0.5, 1, 2, 3, 4, 5, 10, 25, 50, 100, 1000])
    wTimo = np.array([0.0257, 0.01174, 0.00873, 0.00668, 0.00588, 0.00546, 0.00519, 0.00464, 0.00429, 0.00418, 0.00412, 0.00406])
    mTimo = np.array([0.1109, 0.0691, 0.0601, 0.0539, 0.0515, 0.0502, 0.0494, 0.0477, 0.0467, 0.0463, 0.0462, 0.0460])
    f1 = interp1d(gammas, wTimo, kind='linear')
    f2= interp1d(gammas, mTimo, kind='linear')

    values[0,0] = f1(myLambda)
    values[0,1] = f2(myLambda)
    values[0,2] = f2(myLambda)

    return quantities, values, outPos

def _Circ_thin_sSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    quantities=[True, True, True, False, True, False]
    #Actually: Wz, Mx, My, Mxy, Vx, Vy)
    nOutputs = 4
    r = inPos
    nPoints = len(r)
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,1, nOutputs))
    q0 = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    v1=np.ones(nPoints)
    values[:,0] = q0*(a**2*v1-r**2)/(64*D)*((5+nu)/(1+nu)*a**2*v1-r**2)
    outPos[:,:,0] = inPos
    values[:,1] = q0/16*(3+nu)*(a**2*v1-r**2)
    outPos[:,0,1] = inPos
    values[:,2] = q0/16*(a**2*(3+nu)*v1-r**2*(1+3*nu))
    outPos[:,0,2] = inPos
    values[:,3] = q0*r/2
    outPos[:,:,3] = inPos
    return quantities, values, outPos

def _Circ_thin_sSupport_Concntr(pOpts, lOpts, sOpts, inPos):
    quantities=[True, True, True, False, False, False]
    #Actually: Wz, Mx, My, Mxy, Vx, Vy)
    nOutputs = 3
    r = inPos
    nPoints = len(r)
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,1, nOutputs))
    P = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    v1=np.ones(nPoints)
    values[:,0] = P/(16*np.pi*D)*((3+nu)/(1+nu)*(a**2*v1-r**2) + 2*r**2*np.log10(r/a))
    outPos[:,:,0] = inPos
    values[:,1] = P/(4*np.pi)*(1+nu)*np.log10(a/r)
    outPos[:,0,1] = inPos
    values[:,2] = P/(4*np.pi)*((1+nu)*np.log10(a/r)+(1+nu)*v1)
    outPos[:,0,2] = inPos
    return quantities, values, outPos

def _Circ_thin_cSupport_Distrib(pOpts, lOpts, sOpts, inPos):
    quantities=[True, True, True, False, True, False]
    #Actually: Wz, Mx, My, Mxy, Vx, Vy)
    nOutputs = 4
    r = inPos
    nPoints = len(r)
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,1, nOutputs))
    q0 = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    v1=np.ones(nPoints)
    values[:,0] = q0/(64*D)*(a**2*v1-r**2)**2
    outPos[:,:,0] = inPos
    values[:,1] = q0/16*(a**2*(1+nu)*v1-r**2*(3+nu))
    outPos[:,0,1] = inPos
    values[:,2] = q0/16*(a**2*(1+nu)*v1-r**2*(1+3*nu))
    outPos[:,0,2] = inPos
    values[:,3] = q0*r/2
    outPos[:,:,3] = inPos
    return quantities, values, outPos

def _Circ_thin_cSupport_Concntr(pOpts, lOpts, sOpts, inPos):
    quantities=[True, False, False, False, False, False]
    #Actually: Wz, Mx, My, Mxy, Vx, Vy)
    nOutputs = 1
    r = inPos
    nPoints = len(r)
    values = np.zeros((nPoints, nOutputs))
    outPos = np.zeros((nPoints,1, nOutputs))
    P = lOpts.magnitude
    D = pOpts.material.D
    nu = pOpts.material.nu
    a = pOpts.geometry[0]
    b= lOpts.position[0]
    theta = lOpts.position[1]
    x=r/a
    xi = b/a
    v1=np.ones(nPoints)
    values[:,0] = P*a**2/(16*np.pi*D)((v1-x**2)*(1-xi**2)+(x**2+xi**2-2*x*xi*np.cos(theta)*np.log10((x**2+xi**2*v1-2*x*xi*np.cos(theta))/(v1+x**2*xi**2-2*x*xi*np.cos(theta)))))
    outPos[:,:,0] = inPos
    return quantities, values, outPos

#%% for unit testing. This code is executed if this python file is run
if __name__ == "__main__":
    pOpts = POpts()
    pOpts.shape="rectangular"
    pOpts.depth = "thin"
    pOpts.support = "clamped" 
    pOpts.geometry = (1,1)
    pOpts.material = Material(10920, 0.3, 0.1) #E, nu and h

    lOpts = LOpts()
    lOpts.type = "distributed"
    lOpts.magnitude = 1

    sOpts = SOpts()
    sOpts.nTerms = 4

    inPos=np.array([[0,0]])

    quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
    print('outpos: ', outPos)
    print('values: ', values)