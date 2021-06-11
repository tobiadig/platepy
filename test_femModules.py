from .analyticPlateSolutions import *
import numpy as np

def test_Rect_ss_Distrib():
    #Table 8
    pOpts = POpts()
    pOpts.shape="rectangular"
    pOpts.depth = "thin"
    pOpts.support = "simplySupported" 
    a = 1
    b = 1.5
    pOpts.geometry = (a,b)
    pOpts.material = Material(10920, 0.3, 0.1) #E, nu and h
    lOpts = LOpts()
    lOpts.type = "distributed"
    lOpts.magnitude = 1
    sOpts = SOpts()
    sOpts.nTerms = 40
    inPos=np.array([[0,0],[-a/2, 0], [0, -b/2]])
    quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
    assert round(values[0,0],5) == 0.00772
    assert round(values[0,1],4) == 0.0812
    assert round(values[0,2],4) == 0.0498
    assert round(values[1,3],3) == 0.424
    assert round(values[2,4],3) == 0.364

def test_rect_ss_conc():
    # Table 23
    a=1
    b=1.0
    pOpts = POpts()
    pOpts.shape="rectangular"
    pOpts.depth = "thin"
    pOpts.support = "simplySupported" 
    pOpts.geometry = (a,b)
    pOpts.material = Material(10920, 0.3, 0.1) #E, nu and h
    lOpts = LOpts()
    lOpts.type = "concentrated"
    lOpts.position=(0, 0)
    lOpts.magnitude = 1
    sOpts = SOpts()
    sOpts.nTerms = 20
    inPos=np.array([[0,0]])
    quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
    assert round(values[0,0],5) == 0.01160

    b = 1.4
    pOpts.geometry = (a,b)
    quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)
    assert round(values[0,0],5) == 0.01486
