def SampleBasis1(spts, space_order ):
    """Sample basis functions on mesh.

    This routine samples the (orthogonal) basis functions on a mesh.  It
    returns a single vector that can be used to convert modal to nodal values.

    Parameters
    ----------

        space_order : int
        Order of accuracy in space.  We use a total of kmax basis functions.
            In 1D, space_order == kmax.
            In 2D,        kmax == space_order*(space_order+1)/2.
        Note: We do not use a full Tensor product of Legendre polynomials in DoGPack.

        spts   : ndarray
        location of points (on canonical element) to be sampled.
        spts[:] \in [-1,1].

    Returns
    -------

        phi   :  ndarray
        Basis functions evaluated at each point.  phi.shape = (numpts, kmax)

    Notes
    -----

    This routine uses the orthonormal Legendre polynomials.  That is, when
    integrated with a weighted average,

        0.5 * \int_{-1}^1 \phi^k \phi^k = \delta_{k,l},

    where \delta_{k,l} is the Kronecker delta function.
    """

    import numpy as np
    sqrt = np.sqrt

    # Quick error check:
    if( kmax < 0 or kmax > 6 ):
        raise( ValueError("Error: kmax = %d." % kmax ) )

    sq3  = sqrt(3.0)
    sq5  = sqrt(5.0)
    sq7  = sqrt(7.0)
    sq11 = sqrt(11.0)

    phi = np.zeros( (len(spts), kmax), dtype=np.float64 )
    for (n1,xi) in enumerate( spts ):

        xi2 = xi*xi
        xi3 = xi2*xi
        xi4 = xi3*xi
        xi5 = xi4*xi

        phi[n1,0] = 1.0

        if(kmax>1):
            phi[n1,1] = sq3*xi

        if(kmax>2):
            phi[n1,2] = sq5*(1.5*xi2 - 0.5)

        if(kmax>3):
            phi[n1,3] = sq7*(2.5*xi3 - 1.5*xi)

        if(kmax>4):
            phi[n1,4] = 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0

        if(kmax>5):
            phi[n1,5] = (63.0/8.0)*sq11 * ( xi5 - (10.0/9.0)*xi3 + (5.0/21.0)*xi )

    return phi
    #-------------------------------------------------------------------------#

def SampleNonOrthonormalLegendrePolys(spts, kmax ):
    """Sample Non-orthonormal Legendre polynomials on an element.

    This routine samples the (non-orthonormal) basis functions on a mesh.  It
    returns a single vector that can be used to convert modal to nodal values.

    Parameters
    ----------

        kmax   : int
        Maximum number of basis functions to be used.

        spts   : ndarray
        location of points (on canonical element) to be sampled.
        spts[:] \in [-1,1].

    Returns
    -------

        phi   :  ndarray
        Basis functions evaluated at each point.  phi.shape = (numpts, kmax)

    Notes
    -----

    The Legendre polynomials satisfy the recursive relationship

        (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x),

    where P_0 = 1, and P_1 = x.  They satisfy the orthogonality relationship,

        \int_{-1}^{1} P_m(x) P_n(x)\,dx = \\frac{2}{2n+1} \delta_{mn}.
    """

    import numpy as np
    sqrt = np.sqrt

    # Quick error check:
    if( kmax < 0 ):
        raise( ValueError("Error: kmax = %d." % kmax ) )

    phi = np.zeros( (len(spts), kmax), dtype=np.float64 )

    # First polynomial is the constant function
    phi[:,0] = 1.0 + 0*spts

    # Second polynomial (if needed)
    if( kmax > 1 ):
        phi[:,1] = spts.copy()

    # Every other polynomial follows a recursive definition
    for n in range(1, kmax-1 ):
        phi[:,n+1] = ( (2.0*n+1.)*spts*phi[:,n] - n * phi[:,n-1] ) / (n+1.0)

    return phi

def LegQuad(N, (a,b)=(-1.,1.), tol=1e-15, maxiter=5 ):
    """Legendre quadrature.

xi,wgts = LegQuad( N [, (a,b)=(-1.,1.)] [, tol=1e-15] [, maxiter=5] )

Compute the Legendre-Gauss quadrature points of order N with
        wgts(x) = 1 ,    x in (-1,1).
The resulting quadrature points are then mapped to xi in [a,b].

With only one input value (N), the [a,b] default to the interval [-1,1].

Parameters
----------

   N       : int
   Number of qudrature points.  Should be a positive.

   [a,b]   : tuple
   interval of integration.  (Default [a,b] = [-1,1] )

   tol     : float
   tolerance used for Newton iteration.  (Default: 1e-15 )

   maxiter : int
   maximum number of iterations used in Newton iteration.  ( Default: 5 )
   
Returns
-------

   xi   : ndarray
   Quadrature points, in interval [-1,1]

   wgts : ndarray
   Quadrature weights.

Notes
-----

Integration is then accomplished by \sum_i f( xi_i ) w_i.

Examples
--------

>>> import numpy as np
>>> import LegQuad
>>> f = lambda x: np.exp(-x)
>>> xi, wgts = LegQuad.LegQuad(5)
>>> print( xi )
[-0.90617985, -0.53846931,  0.        ,  0.53846931, 0.90617985]
>>> print( wgts )
[ 0.23692689,  0.47862867,  0.56888889, 0.47862867,  0.23692689]
>>> print( sum( wgts ) )
2.0
>>> # Now, approximate \int_{-1}^1 f(x) dx
>>> print( sum( f(xi)*wgts ) )
2.35040238646

Description of the algorithm
----------------------------

The zeros of the Nth Legendre Polynomial P_{N}(x) are computed using:

    P_{n+1} = [(2n+1)/(n+1)]x P_{n} - [n/(n+1)] P_{n-1}.


Using a very accurate initial guess, the zeros are found by standard
application of Newton's method:

                  xnew = xold - P_{N}(xold) / P_{N}'(xold)

The derivative P_{N}'(x) is also found using recursion:

  P_{N}' = N/(1-x^2) ( -x P_{N} + P_{N-1} )

Weights are then obtained by evaluating

              w_i = c_i / ( (1-x^2)*(P_{N}'(x))^2 )

at x = x_i.  Instead of using c_i, the weights are normalized according to

    sum{w_i} =  \int_{-1}^{1} dx = 2.

Original Author: Mathew F. Causley, Michigan State University Fall 2012.
Translated into Python by David C. Seal, Michigan State University Fall 2014.
"""

    import numpy as np

    #-------------------Initial guess is O(N^-4) accurate.--------------------%
    #
    # THIS IS THE PART THAT IS SPECIAL ABOUT THIS CODE.
    #
    #-------------------------------------------------------------------------%
    th  = np.arange(4*N-1, 2, -4 ) * np.pi/(4.*N+2.)             # Iniitial guess for theta
    x   = np.cos( th + (2.0*np.cos(th)/np.sin(th))/(4.*N+2.)**2 )  # Initial guess for the zeros/nodes

    #------------Use recursion to evaluate P_{N}(x) and P_{N}'(x)-------------%
    for k in range( maxiter ):

        Pnm1 = np.ones(N)    #P_0 = 1
        Pn   = x             #P_1 = x

        for n in range(1, N): #compute P_{n+1}(x) at the approximate roots x=x_i.
            Pnp1 = ( (2.*n+1.)*x*Pn- n*Pnm1 )/(n+1.)
            Pnm1 = Pn
            Pn   = Pnp1

        #--------Construct Q = P_{N}' and update x with Newton's method-------%
        Q  = N*(-x*Pn + Pnm1)/( (1.-x**2) )
        x  = x - Pn / Q

        if( np.linalg.norm(Pn/Q, np.inf )<tol ):
            break


    #--------------Map from [-1,1] to [a,b], and compute weights--------------%
    xi   = 0.5*(1.0-x)*a + 0.5*(1.+x)*b
    wgts = 1./( (1.-x**2)*Q**2 )
    wgts = (wgts/sum(wgts))*(b-a)

    return xi, wgts
    #-------------------------------------------------------------------------#

def SampleBasisCart2(spts, space_order ):
    """Sample basis functions on mesh.

    This routine samples the (orthogonal) basis functions on a mesh.  It
    returns a single vector that can be used to convert modal to nodal values.

    Parameters
    ----------

        spts   : ndarray
        location of points (on canonical element) to be sampled.  The
        canonical variables are
            spts[:,0] =  xi[:] \in [-1,1].
            spts[:,1] = eta[:] \in [-1,1].

        space_order : int
        Order of accuracy in space.  We use a total of kmax basis functions.
            In 1D, space_order == kmax.
            In 2D,        kmax == space_order*(space_order+1)/2.
        Note: We do not use a full Tensor product of Legendre polynomials in DoGPack.

    Returns
    -------

        phi   :  ndarray
        Basis functions evaluated at each point.  phi.shape = (numpts, kmax)

    Notes
    -----

    This routine uses the orthonormal Legendre polynomials.  That is, when
    integrated with a weighted average,

        0.5 * \int_{-1}^1 \phi^k \phi^k = \delta_{k,l},

    where \delta_{k,l} is the Kronecker delta function.
    """

    # Quick error check:
    if( space_order < 0 or space_order > 5 ):
        raise( ValueError("Error: space_order = %d." % space_order ) )

    import numpy as np
    sqrt = np.sqrt

    sq3 = sqrt(3.0)
    sq5 = sqrt(5.0)
    sq7 = sqrt(7.0)

    LegVals = np.zeros( spts.shape, dtype=np.float64 )

    numpts, ndims = spts.shape
    for m in range(0, mpts ):

      xi  = spts[m,0]
      eta = spts[m,1]

      xi2 = xi*xi
      xi3 = xi2*xi
      xi4 = xi3*xi

      eta2 = eta*eta
      eta3 = eta2*eta
      eta4 = eta3*eta

      if (space_order==5):
        LegVals[ 0,m]  = 1.0      
        LegVals[ 1,m]  = sq3*xi
        LegVals[ 2,m]  = sq3*eta      
        LegVals[ 3,m]  = 3.0*xi*eta 
        LegVals[ 4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[ 5,m]  = sq5*(1.5*eta2 - 0.5)      
        LegVals[ 6,m]  = sq3*sq5*eta*(1.5*xi2 - 0.5)      
        LegVals[ 7,m]  = sq3*sq5*xi*(1.5*eta2 - 0.5)      
        LegVals[ 8,m]  = sq7*(2.5*xi3 - 1.5*xi)
        LegVals[ 9,m]  = sq7*(2.5*eta3 - 1.5*eta)
        LegVals[10,m] = sq3*sq7*(2.5*xi3 - 1.5*xi)*eta
        LegVals[11,m] = sq3*sq7*(2.5*eta3 - 1.5*eta)*xi
        LegVals[12,m] = 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0)      
        LegVals[13,m] = 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0      
        LegVals[14,m] = 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0
        continue

      elif (space_order==4):
        LegVals[0,m]  = 1.0      
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta      
        LegVals[3,m]  = 3.0*xi*eta 
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)      
        LegVals[6,m]  = sq3*sq5*eta*(1.5*xi2 - 0.5)      
        LegVals[7,m]  = sq3*sq5*xi*(1.5*eta2 - 0.5)      
        LegVals[8,m]  = sq7*(2.5*xi3 - 1.5*xi)
        LegVals[9,m]  = sq7*(2.5*eta3 - 1.5*eta)
        continue

      elif (space_order==3):
        LegVals[0,m]  = 1.0
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta
        LegVals[3,m]  = 3.0*xi*eta
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)
        continue

      elif (space_order==2):
        LegVals[0,m]  = 1.0  
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta
        continue

      elif (space_order==1):
        LegVals[0,m]  = 1.0

    return LegVals
    #-------------------------------------------------------------------------#

if __name__== '__main__':
    """Test basic operations on Legendre polynomials.
    """

    import numpy as np

    Nbig,kmax,tol = (20,5, 1e-14 )
    spts, wgts = LegQuad( Nbig )

    phiDoGPack  = SampleBasis1(spts, kmax )
    phiLegendre = SampleNonOrthonormalLegendrePolys(spts, kmax )

    # Check for orthonormality of the basis functions
    for k1 in range(kmax):
        for k2 in range(kmax):
            val = 0.5*sum( phiDoGPack[:,k1]*phiDoGPack[:,k2]*wgts )
            if( k1 == k2 ):
                assert( abs( val - 1.0 ) < tol )
            else:
                assert( abs( val ) < tol )

    # Check for orthogonality of Legendre polynomials
    for k1 in range(kmax):
        for k2 in range(kmax):
            val = ( 2.0/(2.0*k1+1.0) )**(-1)*sum( phiLegendre[:,k1]*phiLegendre[:,k2]*wgts )
            if( k1 == k2 ):
                assert( abs(val-1.0) < tol )
            else:
                assert( abs( val ) < tol )

    # Check for conversion between Legendre and our basis functions
    for k1 in range(kmax):
        phitest = np.sqrt( (2.0*k1+1.0) ) * phiLegendre[:,k1] - phiDoGPack[:,k1]
        assert( max( abs( phitest ) ) < tol )

    

