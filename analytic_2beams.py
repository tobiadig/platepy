import numpy as np
def TwoBeamsAnalythical(pOpts, lOpts, sOpts, inPos):
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

    nPoints = len(x)
    v1=np.ones(nPoints)
    for m in range(1,nTerms*2+1,2):
        alpham = m*np.pi*b/(2*a)
        Am = 4/(m**5*np.pi**5)*nu*(1+nu)*np.sinh(alpham)-nu*(1-nu)*alpham*np.cosh(alpham)-m*np.pi*(2*np.cosh(alpham) +alpham*np.sinh(alpham))/...
        ((3+nu)*(1-nu)*np.sinh(alpham)*np.cosh(alpham)-(1-nu)**2*alpham+2*m*np.pi*myLambda*np.cosh(alpham)**2)
        Bm = 4/(m**5*np.pi**5)*nu*(1-nu)*np.sinh(alpham)+m*np.pi*myLambda*np.cosh(alpham)/...
        ((3+nu)*(1-nu)*np.sinh(alpham)*np.cosh(alpham)-(1-nu)**2*alpham+2*m*np.pi*myLambda*np.cosh(alpham)**2)

        w = w +q0*a**4/D* (4/(m**5*np.pi**5)*v1 +Am*np.cosh(m*np.pi*y/a) + Bm*m*np.pi*y/a*np.sinh(m*np.pi*y/a)*np.sin(m*np.pi*x/a))

    values[:,0] = w
    outPos[:,:,0] = inPos
    
    return quantities, values, outPos
