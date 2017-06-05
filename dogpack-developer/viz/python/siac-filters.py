def q_init_test( x ):
    """Initial conditions.

This routine can be used to test Gaussian quadrature rules.
"""

    import numpy as np
    return np.sin( 2.0*np.pi*x )

def construct_pre_post_process_errors( q, time, speed, q_init, xbounds, N=6 ):
    """
CONSTRUCT_PRE_POST_PRCOESS_ERRORS.  Evaluate the solution and compute errors.

    Parameters
    ----------

        q       : ndarray
        vector containing (DoGPack) basis coefficnets.
        q.shape = (mx, kmax )

        t       : float
        time of final solution

        q_init  : function handle
        Function describing the initial conditions.

        xbounds : tuple
        tuple with two numbers describing lower and upper bounds of domain.

        N       : int
        Optional.  Integer describing how many quadrature points to use to
        evaluate integral describing the errors.

    Returns
    -------

        Error   : ndarray
        Vector containing errors evaluated at each point.


        erinf   : float
        Vector describing the L2 error of the solution, computed by evaluated
        exact-approx at each quadrature point.
        
"""

    import numpy as no
    from LegendrePolys import SampleBasis1, LegQuad

    # Size of current unkown
    mx, kmax = q.shape

    # Grid information
    xleft  = xbounds[0]
    xright = xbounds[1]
    dx     = ( xright - xleft ) / mx


    # Vector of cell centers for each element.  This adds a single ghost cell to
    # the left and right of the domain.
    x = np.linspace( xleft+0.5*dx, xright-0.5*dx, mx )
    
    Approx = []          # Pointwise approximation to solution
    Exact  = []          # Exact solution
    Y      = []          # Grid points
    er2, erinf = (0.,0.) # errors

    # Gauss-legendre polynomial at each quadrature point
    wg, xg = LegQuad(N)
    pl     = SampleBasis1(xg, kmax )

    # main loop over entire domain
    for i in range( mx ):

        # Integrate, and add the error over grid cell x_i to the solution
        for j in range( N ):

            xr  = x[i] + xg[j]*dx*0.5
            xrt = xr - speed*time

            # Approximation
            uapprox = 0
            for kk in range( kmax ):
                uapprox = uapprox + q[i, kk]*pl[j, kk]
            Approx.append( uapprox )

            # Exact
            exact = q_init( xrt )
            Exact.append( exact )
            er2 = er2 + (exact-uapprox)**2*dx*wg[j]

    # Pointwise error, and L-inf error
    Error = abs( np.array(Exact) - np.array(Approx) );
    erinf = max( np.array( Error) )

    return Error, erinf

def SymmetricPostprocess( q, time, kmax, xbounds ):
    """Postprocess the solution.

Calculate post-processed solution, and evaluate it.  This version only uses the symmetric post-processor

Parameters
----------

    q   : ndarray
    Vector of unkowns.  q.shape = (mx, kmax).  
    
    time  : float
    Time for solution.  This is used to construct exact solution.

    xbounds : tuple
    Size of computational domain.

Returns
-------

  er2pp      : TODO
  l2 error

  erinfpp    : TODO l-inf error
"""

    import numpy as no
    from LegendrePolys import SampleBasis1, LegQuad

    # Size of current unkown
    mx, kmax = q.shape
    mpp      = kmax-1

    # Grid information
    xleft  = xbounds[0]
    xright = xbounds[1]
    dx     = ( xright - xleft ) / mx


    # Vector of cell centers for each element.  This adds a single ghost cell to
    # the left and right of the domain.
    x = linspace( xleft-0.5*dx, xright+0.5*dx, Elmts+2 )
    
    Approx = []          # Pointwise approximation to solution
    Exact  = []          # Exact solution
    Y      = []          # Grid points
    er2, erinf = (0.,0.) # errors

    # in the symmetric part
    rspline = 4*mpp
    morder  = mpp+1
    numbs1  = rspline+1

    # shift function
    kwide = int((3*mpp+1)/2)
    rx=0.
    ccsym = np.zeros( (2*kwide+1,mpp+1,6) )

    # Extract the conversion coeffient matrix
    # TODO - rewrite these!
    if (mpp==1):
        C1to6k1
    elif(mpp==2):
        C1to6k2
    elif(mpp==3):
        C1to6k3
    elif(mpp==4):
        C1to6k4
    else:
        fprintf('error, mpp>6\n')

    # DEFINE THE COEFFEICIENTS ccsym used for constructing post-processed solution.
    #
    # This section converts the `monomial' basis to Legendre basis.
    #
    ccmod = np.zeros( (2*kwide+1,mpp+1,6) )
    for kk in range(2*kwide+1):

        ccmod[kk,0,:] =    ccsym[kk,0,:]
        ccmod[kk,1,:] = 2.*ccsym[kk,1,:]
        if (mpp>=2):
           ccmod[kk,2,:] = 6.*ccsym[kk,2,:]-0.5*ccsym[kk,0,:]
        if (mpp>=3):
           ccmod[kk,3,:] = 20.*ccsym[kk,3,:]-3.0*ccsym[kk,1,:]
        if (mpp>=4):
           ccmod[kk,4,:] = 70.*ccsym[kk,4,:]-15.*ccsym[kk,2,:]+0.375*ccsym[kk,0,:]
        if (mpp>=5):
           ccmod[kk,5,:] = 252.*ccsym[kk,5,:]-70.*cc[kk,3,:]+3.75*ccsym[kk,1,:] 
        if (mpp>=6):
           ccmod[kk,6,:] = 924.*ccsym[kk,6,:]-315.*ccsym[kk,4,:]+26.25*ccsym[kk,2,:] - 0.3125*ccsym[kk,0,mmp,:]

    # ----------------------------------------------------------------------------
    # Compute errors when compared with exact solution (sampled at GL-points)
    # ----------------------------------------------------------------------------
    er2pp, errinfpp = (0.,0.)

    # Each of these will eventually have size ( (Elmts+1)*6, 1 )
    Ypp      = []      # Mesh - sampled at each quadrature point
    Exactpp  = []      # Exact solution - sampled at each quadrature point
    ppApprox = []      # Symmetric, post-processed solution sampled at each quadrature point

    # Gaussian quadrature weights and points, (using exactly 6 points)
    wg, xg = LegQuad(N)

    for i in range(1, mx ): # modified to loop over interior points only (-DS)

        # Loop over each Gaussian quadrature point
        for j in range(6):

            xr  = x[i] + xg[j]*dx*0.5
            xrt = xr - speed*time

            # Symmetric Post-processed Approximation
            ppuapprox=0;
            for k in range(morder):
                #for m =-kwide:kwide
                for m in range( -kwide, kwide ):
                    if( (i+m)<=0 ):
                        ppuapprox = ppuapprox + ccmod(m+kwide+1,k,j)*u(k,Elmts+(i+m))
                    elif ((i+m)>=(Elmts+2)):
                        ppuapprox = ppuapprox + ccmod(m+kwide+1,k,j)* u(k,2+i+m-(Elmts+2))
                    else:
                        ppuapprox = ppuapprox + ccmod(m+kwide+1,k,j)*u(k,i+m)

            #ppApprox=[ppApprox;ppuapprox]
            ppApprox.append( ppuapprox )

            # Exact Solution
            exact   = q_init( xrt );

            Exactpp.append( exact )
            Ypp.append( xr )

            # Running total for l2-error
            er2pp=er2pp+(exact-ppuapprox)^2*dx*wg(j);

    # Pointwise error, and L-inf error
    Error = abs( np.array(Exactpp) - np.array(ppApprox) );
    erinf = max( np.array( ppError) )

    # L2 error
    er2pp   = np.sqrt(0.5*er2pp)

if __name__=='__main__':
    
    import numpy as np
    from LegendrePolys import SampleBasis1, LegQuad

    kmax   = 5          # Number of basis functions to use
    N      = 27         # Quadrature order
    niter  = 5          # Number of iterations to run

    # Grid information and polynomial order
    xleft, xright = (0., 1.)

    # Basis functions evaluated on the canonial element
    xg, wg = LegQuad(N)
    phi    = SampleBasis1(xg, kmax )

    for n in range(niter):

        mx     = int(20*2**(n-1))
        dx     = ( xright - xleft ) / mx
        x      = np.linspace( xleft+0.5*dx, xright-0.5*dx, mx )

        # Construct a projection onto the initial conditions
        q = np.zeros( (mx,kmax) )
        for i in range(mx):

            for k in range(kmax):

                # Integrate over grid cell x_i to find the solution
                qval = 0.
                for mp in range( N ):

                    # Initial conditions and quadrature point
                    xr   = x[i] + xg[mp]*dx*0.5
                    qval = qval + q_init_test( xr )*phi[mp,k]*0.5*wg[mp]

                q[i,k] = qval

        if( n == 0 ):
            Error, erinf = construct_pre_post_process_errors( q, 1.0, 1.0, q_init_test, (xleft,xright) )
            print('erinf = %2.5e; log2( ratio ) = ---' % erinf )
        else:
            Error_old = Error
            erinf_old = erinf
            Error, erinf = construct_pre_post_process_errors( q, 1.0, 1.0, q_init_test, (xleft,xright) )
            print('erinf = %2.5e; log2( ratio ) = %2.3f' % ( erinf, np.log2( erinf_old / erinf ) ) )



