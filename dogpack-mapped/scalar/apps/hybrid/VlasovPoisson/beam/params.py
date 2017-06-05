# List of the basic parameters for this problem:

import math

# Parameters for the model:
a   = 1.
vth = 1.
eta = 0.25

# Note: eta = omega / omega0

# Parameters derived from the model:
omega0 = 2.*vth/(a*eta)   # == 8.0
omega  = eta * omega0     # definition of eta

# solve omega0 = sqrt{ n0 / (2*(1-eta^2) ) } for n0:
n0 = 2*(1-eta**2)*omega0**2

# Time stepping parameters:
# Tfinal = 1.68 / omega0

T  = (2.0 * math.pi) / omega0
dt = T / 250. # time step used by Besse and Sonnendrucker (2005) paper

Tfinal = 1.68 / omega0
Tfinal = 2.0 * T

# TODO - this doesn't agree with what's in the paper!!!
#
# They say, "The computation runs on 128 processors of a Compaq Alphaserver
# ES45 1GHz and takes around 5 hours when the final time step Tfinal = 2 T."
#
# They define T = 2 pi / omega0, which means that with our definition of
# omega0, we would have a different final time!?
