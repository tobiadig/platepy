# ------------------------------------------------------------------------------
#
# platepy tutorial 4
#
# Compute analytical solution according to Timoshenko's formulas
#
# ------------------------------------------------------------------------------

# This tutorial shows how to use the analyticPlateSolutions module
# We will access the displacement at the location (x=-2, y=1) (relative 
# to the plate's center) of a simply supported rectangular plate with shape
# a=10, b=12

from platepy import *

# First, let's define the input options as attributes of the POpts class
# The possible entries are described in the documentation
pOpts = POpts()
pOpts.shape="rectangular"
pOpts.depth = "thin"
pOpts.support = "simplySupported"
pOpts.geometry = (10,12)  # width a and high b
pOpts.material = Material(10920, 0.3, 0.1) #E, nu and thickness

# Here the Load options are defined, as attributes of the LOpts class
# The possible entries are described in the documentation
lOpts = LOpts()
lOpts.type = "distributed"
lOpts.magnitude = 1

# The number of term for the Taylor expantion has to be defined.
# 40 terms is usually more than enough
sOpts = SOpts()
sOpts.nTerms = 40

#  inPos is a nPointsx2 numpy array containing the x-y coordinates of
# the points we want to evaluate. In our case we only need on point.
xPos = -2
yPos = 1
inPos=np.array([[xPos,yPos]])  # x-y coordinates

quantities, values, outPos = AnalyticPlateSolutions(pOpts, lOpts, sOpts, inPos)

# values is a (nPoints, nOutputs) array of the quantities values evaluated at outPos.
# the output values correspond to the True entries of the array quantities.
# the entries of quantities correspond to: (vDisp, Mx, My, Mxy, Vx, Vy).
# In this case, all values are available apart from the twisting moment Mxy, therefore
# values has shape (1,5).

vDisp_mm = values[0,0]*1000

print("The displacement at (", xPos,",",yPos,") corresponds to ", vDisp_mm, ' mm')
