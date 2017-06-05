// Legendre2d_Unst.cpp.
//
// This file describes common operations performed one the 2D (Unstructured) basis
// functions.
//
// See also: L2Project_Unst, L2ProjectGrad_Unst and $DOGPACK/lib/2d/cart/Legendre2d.

#include "dogdefs.h"                // for tensors and unsupported_value_error
#include "MonomialsToLegendre.h"    // Converts monomials to Legendre

// -------------------------------------------------------------------------- //
// Routine for evaluating DG basis function at a list of points.
//
// Input:
// ------
//
//      spts( 1:numpts, 1:2 )  - points (in canonical coordinates)
//
//          xi  = spts( :, 1 ) - the "x"-coordinate, after mapping
//          eta = spts( :, 2 ) - the "y"-coordinate, after mapping
//
// Returns:
// --------
//
//      phi( 1:numpts, 1:kmax ) - basis functions evaluated at each point.  It
//                                is expected that kmax is consisten with 
//                                dogParamsCart2.get_kmax()
//
// See also: SetLegendreGrad.
// -------------------------------------------------------------------------- //
void SetLegendreAtPoints_Unst(const dTensor2& spts, dTensor2& phi)
{

    // Number of points and number of basis functions
    const int mpoints = phi.getsize(1);
    const int    kmax = phi.getsize(2);

    // monomial basis (non-orthogonal)
    dTensor2   mu(mpoints, kmax); 

    // Loop over each quadrature point and construct monomial polys
    for (int m=1; m<=mpoints; m++)
    {
        // coordinates
        const double xi   = spts.get(m,1);      
        const double xi2  = xi*xi;
        const double xi3  = xi2*xi;
        const double xi4  = xi3*xi;
        const double eta  = spts.get(m,2);
        const double eta2 = eta*eta;
        const double eta3 = eta2*eta;
        const double eta4 = eta3*eta;      

        // monomials basis (non-orthogonal)
        switch( kmax )
        {
            case 15:  // fifth order		    		    
                mu.set(m, 15, eta4     );
                mu.set(m, 14, xi4      );
                mu.set(m, 13, xi2*eta2 );
                mu.set(m, 12, eta3*xi  );
                mu.set(m, 11, xi3*eta  );

            case 10:  // fourth order
                mu.set(m, 10, eta3     );
                mu.set(m, 9,  xi3      );
                mu.set(m, 8,  xi*eta2  );
                mu.set(m, 7,  eta*xi2  );

            case 6:  // third order
                mu.set(m, 6,  eta2     );
                mu.set(m, 5,  xi2      );
                mu.set(m, 4,  xi*eta   );		    

            case 3:  // second order		    
                mu.set(m, 3, eta       );
                mu.set(m, 2, xi        );

            case 1:  // first order
                mu.set(m, 1, 1.0       );

                break;		    
        }
    }

    // Loop over each quadrature point and construct Legendre polys.
    // 
    // That is, convert the monomial basis functions (mu) to orthonormal basis
    // functions (phi) through the change of basis matrix, Mmat.
    //
    // See: MonomialsToLegendre.h.
    for (int m=1; m<=mpoints; m++)    
    for (int i=1; i<=kmax; i++)
    {
        double tmp = 0.0;
        for (int j=1; j<=i; j++)
        {  tmp += Mmat[i-1][j-1]*mu.get(m,j);  }
        phi.set(m,i, tmp );      
    }

}

// -------------------------------------------------------------------------- //
// Routine for evaluating the derivatives of the DG basis function at a list 
// of points.
//
// Input:
// ------
//
//   spts( 1:numpts, 1:2 )  - points (in canonical coordinates)
//
//          xi  = xpts( :, 1 ) 
//          eta = xpts( :, 2 )
//
// Returns:
// --------
//
// Derivatives of the basis functions:
//
//      phi_xi ( 1:numpts, 1:kmax ) - \pd{ phi^k }{ xi }
//      phi_eta( 1:numpts, 1:kmax ) - \pd{ phi^k }{ eta }
//
// In order to convert to a "physical" derivative, one needs to include the
// Jacobian for the change of basis from (x,y) <-> (xi, eta).
//
// See also: SetLegendreAtPoints_Unst, L2ProjectGrad_Unst.
// -------------------------------------------------------------------------- //
void SetLegendreGrad_Unst( const dTensor2& spts, dTensor2& phi_xi, dTensor2& phi_eta )
{

//  const int space_order = dogParams.get_space_order();
//  const int kmax        = dogParams.get_kmax();

    const int mpoints     = phi_xi.getsize(1);
    const int kmax_fout   = phi_xi.getsize(2);

    dTensor2   mu_xi(mpoints, kmax_fout);    //  xi-derivative of monomial basis (non-orthogonal)
    dTensor2  mu_eta(mpoints, kmax_fout);    // eta-derivative of monomial basis (non-orthogonal)

    // Construct monomial polys
    for (int m=1; m<=mpoints; m++)
    {
        // coordinates
        const double xi   = spts.get(m,1);      
        const double xi2  = xi*xi;
        const double xi3  = xi2*xi;
        const double xi4  = xi3*xi;
        const double eta  = spts.get(m,2);
        const double eta2 = eta*eta;
        const double eta3 = eta2*eta;
        const double eta4 = eta3*eta;      

        // Gradient of monomial basis functions at each gaussian quadrature point
        switch( kmax_fout )
        {

            case 15:  // fifth order
                mu_xi.set( m,15,  0.0         );
                mu_xi.set( m,14,  4.0*xi3     );
                mu_xi.set( m,13,  2.0*xi*eta2 );
                mu_xi.set( m,12,  eta3        );
                mu_xi.set( m,11,  3.0*xi2*eta );

                mu_eta.set( m,15, 4.0*eta3    );
                mu_eta.set( m,14, 0.0         );
                mu_eta.set( m,13, 2.0*xi2*eta );
                mu_eta.set( m,12, 3.0*eta2*xi );
                mu_eta.set( m,11, xi3 );

            case 10:  // fourth order
                mu_xi.set( m,10,  0.0        );
                mu_xi.set( m,9,   3.0*xi2    );			
                mu_xi.set( m,8,   eta2       );
                mu_xi.set( m,7,   2.0*eta*xi );

                mu_eta.set( m,10, 3.0*eta2   );
                mu_eta.set( m,9,  0.0        );
                mu_eta.set( m,8,  2.0*eta*xi );
                mu_eta.set( m,7,  xi2        );

            case 6:  // third order
                mu_xi.set( m,6,  0.0      );
                mu_xi.set( m,5,  2.0*xi   );			
                mu_xi.set( m,4,  eta      );

                mu_eta.set( m,6,  2.0*eta );			
                mu_eta.set( m,5,  0.0     );
                mu_eta.set( m,4,  xi      );

            case 3:  // second order
                mu_xi.set( m,3,  0.0 );
                mu_xi.set( m,2,  1.0 );

                mu_eta.set( m,3, 1.0 );
                mu_eta.set( m,2, 0.0 );

            case 1:  // first order
                mu_xi.set( m,1,  0.0 );

                mu_eta.set( m,1, 0.0 );
                break;
        }

        // Loop over each quadrature point and construct Legendre polys
        for (int i=1; i<=kmax_fout; i++)
        {
            double tmp1 = 0.0;
            double tmp2 = 0.0;
            for (int j=1; j<=i; j++)
            {  
                tmp1 = tmp1 + Mmat[i-1][j-1]*mu_xi.get(m,j);  
                tmp2 = tmp2 + Mmat[i-1][j-1]*mu_eta.get(m,j);
            }

            phi_xi.set(m,i,  tmp1 );
            phi_eta.set(m,i, tmp2 );
        }

    }


}

// -------------------------------------------------------------------------- //
// Routine for evaluating the second derivative of the DG basis function.
//
// Input:
// ------
//
//   spts( 1:numpts, 1:2 )  - points (in canonical coordinates)
//
//          xi  = xpts( :, 1 ) 
//          eta = xpts( :, 2 )
//
// Returns:
// --------
//
// Derivatives of the basis functions:
//
//      phi_xi2  ( 1:numpts, 1:kmax ) - \pd2{ phi^k }{ xi^2   }
//      phi_xieta( 1:numpts, 1:kmay ) - \pd2{ phi^k }{ xi eta }
//      phi_eta2 ( 1:numpts, 1:kmay ) - \pd2{ phi^k }{ eta^2  }
//
// These elements still need to be mapped to their physical derivatives through
// the use of a Jacobian.
//
// See also: SetLegendreAtPoints.
// -------------------------------------------------------------------------- //
void LegendreDiff2_Unst(const dTensor2& spts, 
    dTensor2* phi_xi2, dTensor2* phi_xieta, dTensor2* phi_eta2 )
{

//  const int space_order = dogParams.get_space_order();
//  const int kmax        = dogParams.get_kmax();
    const int mpoints     = phi_xi2->getsize(1);
    const int kmax_fout   = phi_xi2->getsize(2);

    // Derivatives for each monomial (non-orthogonal)
    dTensor2   mu_xi2   (mpoints, kmax_fout);
    dTensor2   mu_xieta (mpoints, kmax_fout);
    dTensor2  mu_eta2   (mpoints, kmax_fout);

    // Construct monomial polys
    for (int m=1; m<=mpoints; m++)
    {

        // grid point (x,y) <-> (xi, eta)
        const double xi   = spts.get(m,1);
        const double eta  = spts.get(m,2);
        const double xi2  = xi*xi;
        const double xi3  = xi*xi2;
        const double xi4  = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     

        // Gradient of monomial basis functions at each gaussian quadrature point
        switch( kmax_fout )
        {

            case 15:  // fifth order
                mu_xi2.set( m,15,  0.0         );
                mu_xi2.set( m,14, 12.0*xi2     );
                mu_xi2.set( m,13,  2.0*eta2    );
                mu_xi2.set( m,12,  0.0         );
                mu_xi2.set( m,11,  6.0*xi *eta );

                mu_xieta.set( m,15, 0.0         );
                mu_xieta.set( m,14, 0.0         );
                mu_xieta.set( m,13, 4.0*xi *eta );
                mu_xieta.set( m,12, 3.0*eta2    );
                mu_xieta.set( m,11, 3.0*xi2     );


                mu_eta2.set( m,15,12.0*eta2    );
                mu_eta2.set( m,14, 0.0         );
                mu_eta2.set( m,13, 2.0*xi2     );
                mu_eta2.set( m,12, 6.0*eta *xi );
                mu_eta2.set( m,11, 0.0         );


            case 10:  // fourth order

                mu_xi2.set( m,10,  0.0        );
                mu_xi2.set( m,9,   6.0*xi     );			
                mu_xi2.set( m,8,   0.0        );
                mu_xi2.set( m,7,   2.0*eta    );

                mu_xieta.set( m,10,  0.0     );
                mu_xieta.set( m,9,   0.0     );			
                mu_xieta.set( m,8,   2.0*eta );
                mu_xieta.set( m,7,   2.0*xi  );


                mu_eta2.set( m,10, 6.0*eta  );
                mu_eta2.set( m,9,  0.0      );
                mu_eta2.set( m,8,  2.0*xi   );
                mu_eta2.set( m,7,  0.0      );


            case 6:  // third order
                mu_xi2.set( m,6,  0.0       );
                mu_xi2.set( m,5,  2.0       );			
                mu_xi2.set( m,4,  0.0       );

                mu_xieta.set    ( m,6,  0.0 );			
                mu_xieta.set    ( m,5,  0.0 );
                mu_xieta.set    ( m,4,  1.0 );

                mu_eta2.set     ( m,6,  2.0 );			
                mu_eta2.set     ( m,5,  0.0 );
                mu_eta2.set     ( m,4,  0.0 );


            case 3:  // second order
                mu_xi2.set      ( m,3,  0.0 );
                mu_xi2.set      ( m,2,  0.0 );

                mu_xieta.set    ( m,3, 0.0 );
                mu_xieta.set    ( m,2, 0.0 );

                mu_eta2.set     ( m,3, 0.0 );
                mu_eta2.set     ( m,2, 0.0 );

            case 1:  // first order
                mu_xi2.  set( m,1, 0.0 );
                mu_xieta.set( m,1, 0.0 );
                mu_eta2. set( m,1, 0.0 );
                break;
        }

        // Loop over each quadrature point and construct Legendre polys
        for (int i=1; i<=kmax_fout; i++)
        {
            double tmp[] = {0.0,0.0,0.0};
            for (int j=1; j<=i; j++)
            {  
                tmp[0] += Mmat[i-1][j-1]*mu_xi2.  get(m,j);  
                tmp[1] += Mmat[i-1][j-1]*mu_xieta.get(m,j);  
                tmp[2] += Mmat[i-1][j-1]*mu_eta2. get(m,j);  
            }

            phi_xi2  ->set (m,i, tmp[0] );
            phi_xieta->set (m,i, tmp[1] );
            phi_eta2 ->set (m,i, tmp[2] );
        }

    }
}

// -------------------------------------------------------------------------- //
// Set (2D) quadrature weights and points.  These are the points used for any
// L2Project call.  
//
// NOTE: The quadrature points used for projecting onto the gradient are 
//       different!
//
// Quadrature weights and points.  Our Canonical triangle has
// coordinates:
//
//      (-1/3,-1/3), ( 2/3, -1/3 ), and ( -1/3, 2/3 ).
//
// The vertices are oriented in counterclockwise order 
// (in order to follow a right-hand rule for cross products).
//
// See "High degree efficient symmetrical Gaussian quadrature rules for the
// triangle", D. A. Dunavant Internat. J. Numer. Methods Engrg. 21 (1985),
// no. 6, 1129–1148.  DOI: 10.1002/nme.1620210612
//
// For futher description, see:
//      http://people.sc.fsu.edu/~jburkardt/f_src/dunavant/dunavant.html
// 
// See also: $DOGPACK/lib/Quadrature.h, setQuadPointsGrad_Unst.
// ------------------------------------------------------------------------- //
void setQuadPoints_Unst(int QuadOrder, dTensor1& wgts, dTensor2& spts)
{

    switch( QuadOrder )
    {

        case 1:
            spts.set(1,1, 0.0 );
            spts.set(1,2, 0.0 );

            wgts.set(1, 0.5 );
            break;

        case 2:
            spts.set(1,1,  1.0/3.0 );
            spts.set(1,2, -1.0/6.0 );

            spts.set(2,1, -1.0/6.0 );
            spts.set(2,2, -1.0/6.0 );

            spts.set(3,1, -1.0/6.0 );
            spts.set(3,2,  1.0/3.0 );

            wgts.set(1, 1.0/6.0 );
            wgts.set(2, 1.0/6.0 );
            wgts.set(3, 1.0/6.0 );
            break;

        case 3:
            spts.set(1,1,  0.112615157582632 );
            spts.set(1,2,  0.112615157582632 );

            spts.set(2,1, -0.225230315165263 );
            spts.set(2,2,  0.112615157582632 );

            spts.set(3,1,  0.112615157582632 );
            spts.set(3,2, -0.225230315165263 );

            spts.set(4,1, -0.241757119823562 );
            spts.set(4,2, -0.241757119823562 );

            spts.set(5,1,  0.483514239647126 );
            spts.set(5,2, -0.241757119823562 );

            spts.set(6,1, -0.241757119823562 );
            spts.set(6,2,  0.483514239647126 );

            wgts.set(1, 0.1116907948390055 );
            wgts.set(2, 0.1116907948390055 );
            wgts.set(3, 0.1116907948390055 );
            wgts.set(4, 0.0549758718276610 );
            wgts.set(5, 0.0549758718276610 );
            wgts.set(6, 0.0549758718276610 );
            break;

        case 4:
            spts.set(1,1,  -0.084046588162423 );
            spts.set(1,2,  -0.084046588162423 );

            spts.set(2,1,   0.168093176324846 );
            spts.set(2,2,  -0.084046588162423 );

            spts.set(3,1,  -0.084046588162423 );
            spts.set(3,2,   0.168093176324846 );

            spts.set(4,1,  -0.270244318841831 );
            spts.set(4,2,  -0.270244318841831 );

            spts.set(5,1,   0.540488637683663 );
            spts.set(5,2,  -0.270244318841831 );

            spts.set(6,1,  -0.270244318841831 );
            spts.set(6,2,   0.540488637683663 );

            spts.set(7,1,  -0.280188283488516 );
            spts.set(7,2,  -0.022980882299549 );

            spts.set(8,1,  -0.280188283488516 );
            spts.set(8,2,   0.303169165788066 );

            spts.set(9,1,  -0.022980882299549 );
            spts.set(9,2,   0.303169165788067 );

            spts.set(10,1, -0.022980882299549 );
            spts.set(10,2, -0.280188283488516 );

            spts.set(11,1,  0.303169165788066 );
            spts.set(11,2, -0.022980882299549 );

            spts.set(12,1,  0.303169165788066 );
            spts.set(12,2, -0.280188283488516 );

            wgts.set(1,  0.0583931378631895 );
            wgts.set(2,  0.0583931378631895 );
            wgts.set(3,  0.0583931378631895 );
            wgts.set(4,  0.0254224531851035 );
            wgts.set(5,  0.0254224531851035 );
            wgts.set(6,  0.0254224531851035 );
            wgts.set(7,  0.0414255378091870 );
            wgts.set(8,  0.0414255378091870 );
            wgts.set(9,  0.0414255378091870 );
            wgts.set(10, 0.0414255378091870 );
            wgts.set(11, 0.0414255378091870 );
            wgts.set(12, 0.0414255378091870 );
            break;

        case 5:
            spts.set(1,1,   0.000000000000000 );
            spts.set(1,2,   0.000000000000000 );

            spts.set(2,1,   0.125959254959390 );
            spts.set(2,2,   0.125959254959390 );

            spts.set(3,1,  -0.251918509918779 );
            spts.set(3,2,   0.125959254959390 );

            spts.set(4,1,   0.125959254959390 );
            spts.set(4,2,  -0.251918509918779 );

            spts.set(5,1,  -0.162764025581573 );
            spts.set(5,2,  -0.162764025581573 );

            spts.set(6,1,   0.325528051163147 );
            spts.set(6,2,  -0.162764025581573 );

            spts.set(7,1,  -0.162764025581573 );
            spts.set(7,2,   0.325528051163147 );

            spts.set(8,1,  -0.282786105016302 );
            spts.set(8,2,  -0.282786105016302 );

            spts.set(9,1,   0.565572210032605 );
            spts.set(9,2,  -0.282786105016302 );

            spts.set(10,1, -0.282786105016302 );
            spts.set(10,2,  0.565572210032605 );

            spts.set(11,1, -0.324938555923375 );
            spts.set(11,2, -0.070220503698695 );

            spts.set(12,1, -0.324938555923375 );
            spts.set(12,2,  0.395159059622071 );

            spts.set(13,1, -0.070220503698695 );
            spts.set(13,2, -0.324938555923375 );

            spts.set(14,1, -0.070220503698695 );
            spts.set(14,2,  0.395159059622071 );

            spts.set(15,1,  0.395159059622071 );
            spts.set(15,2, -0.324938555923375 );

            spts.set(16,1,  0.395159059622071 );
            spts.set(16,2, -0.070220503698695 );

            wgts.set(1,  0.0721578038388935 );
            wgts.set(2,  0.0475458171336425 );
            wgts.set(3,  0.0475458171336425 );
            wgts.set(4,  0.0475458171336425 );
            wgts.set(5,  0.0516086852673590 );
            wgts.set(6,  0.0516086852673590 );
            wgts.set(7,  0.0516086852673590 );
            wgts.set(8,  0.0162292488115990 );
            wgts.set(9,  0.0162292488115990 );
            wgts.set(10, 0.0162292488115990 );
            wgts.set(11, 0.0136151570872175 );
            wgts.set(12, 0.0136151570872175 );
            wgts.set(13, 0.0136151570872175 );
            wgts.set(14, 0.0136151570872175 );
            wgts.set(15, 0.0136151570872175 );
            wgts.set(16, 0.0136151570872175 );
            break;
    }

}

// -------------------------------------------------------------------------- //
// Set (2D) quadrature weights and points.  These are the points used for calls
// to L2ProjectGrad_Unst.
//
// NOTE: The quadrature points used for regular projections are 
//       different!
//
// Quadrature weights and points.  Our Canonical triangle has
// coordinates:
//
//      (-1/3,-1/3), ( 2/3, -1/3 ), and ( -1/3, 2/3 ).
//
// The vertices are oriented in counterclockwise order 
// (in order to follow a right-hand rule for cross products).
//
// See "High degree efficient symmetrical Gaussian quadrature rules for the
// triangle", D. A. Dunavant Internat. J. Numer. Methods Engrg. 21 (1985),
// no. 6, 1129–1148.  DOI: 10.1002/nme.1620210612
//
// For futher description, see:
//      http://people.sc.fsu.edu/~jburkardt/f_src/dunavant/dunavant.html
// 
// See also: $DOGPACK/lib/Quadrature.h, setQuadPoints_Unst.
// ------------------------------------------------------------------------- //
void setQuadPointsGrad_Unst(int QuadOrder, dTensor1& wgts, dTensor2& spts)
{

    switch ( QuadOrder )
    {
        case 2:
            spts.set(1,1, 0.0 );
            spts.set(1,2, 0.0 );

            wgts.set(1, 0.5 );
            break;

        case 3:
            spts.set(1,1,  0.112615157582632 );
            spts.set(1,2,  0.112615157582632 );

            spts.set(2,1, -0.225230315165263 );
            spts.set(2,2,  0.112615157582632 );

            spts.set(3,1,  0.112615157582632 );
            spts.set(3,2, -0.225230315165263 );

            spts.set(4,1, -0.241757119823562 );
            spts.set(4,2, -0.241757119823562 );

            spts.set(5,1,  0.483514239647126 );
            spts.set(5,2, -0.241757119823562 );

            spts.set(6,1, -0.241757119823562 );
            spts.set(6,2,  0.483514239647126 );

            wgts.set(1, 0.1116907948390055 );
            wgts.set(2, 0.1116907948390055 );
            wgts.set(3, 0.1116907948390055 );
            wgts.set(4, 0.0549758718276610 );
            wgts.set(5, 0.0549758718276610 );
            wgts.set(6, 0.0549758718276610 );
            break;

        case 4:
            spts.set(1,1,   0.000000000000000 );
            spts.set(1,2,   0.000000000000000 );

            spts.set(2,1,   0.136808730771782 );
            spts.set(2,2,   0.136808730771782 );

            spts.set(3,1,  -0.273617461543563 );
            spts.set(3,2,   0.136808730771782 );

            spts.set(4,1,   0.136808730771782 );
            spts.set(4,2,  -0.273617461543563 );

            spts.set(5,1,  -0.232046826009877 );
            spts.set(5,2,  -0.232046826009877 );

            spts.set(6,1,   0.464093652019754 );
            spts.set(6,2,  -0.232046826009877 );

            spts.set(7,1,  -0.232046826009877 );
            spts.set(7,2,   0.464093652019754 );	 

            wgts.set(1,  0.1125000000000000 );
            wgts.set(2,  0.0661970763942530 );
            wgts.set(3,  0.0661970763942530 );
            wgts.set(4,  0.0661970763942530 );
            wgts.set(5,  0.0629695902724135 );
            wgts.set(6,  0.0629695902724135 );
            wgts.set(7,  0.0629695902724135 );
            break;

        case 5:
            spts.set(1,1,   0.000000000000000 );
            spts.set(1,2,   0.000000000000000 );

            spts.set(2,1,   0.125959254959390 );
            spts.set(2,2,   0.125959254959390 );

            spts.set(3,1,  -0.251918509918779 );
            spts.set(3,2,   0.125959254959390 );

            spts.set(4,1,   0.125959254959390 );
            spts.set(4,2,  -0.251918509918779 );

            spts.set(5,1,  -0.162764025581573 );
            spts.set(5,2,  -0.162764025581573 );

            spts.set(6,1,   0.325528051163147 );
            spts.set(6,2,  -0.162764025581573 );

            spts.set(7,1,  -0.162764025581573 );
            spts.set(7,2,   0.325528051163147 );

            spts.set(8,1,  -0.282786105016302 );
            spts.set(8,2,  -0.282786105016302 );

            spts.set(9,1,   0.565572210032605 );
            spts.set(9,2,  -0.282786105016302 );

            spts.set(10,1, -0.282786105016302 );
            spts.set(10,2,  0.565572210032605 );

            spts.set(11,1, -0.324938555923375 );
            spts.set(11,2, -0.070220503698695 );

            spts.set(12,1, -0.324938555923375 );
            spts.set(12,2,  0.395159059622071 );

            spts.set(13,1, -0.070220503698695 );
            spts.set(13,2, -0.324938555923375 );

            spts.set(14,1, -0.070220503698695 );
            spts.set(14,2,  0.395159059622071 );

            spts.set(15,1,  0.395159059622071 );
            spts.set(15,2, -0.324938555923375 );

            spts.set(16,1,  0.395159059622071 );
            spts.set(16,2, -0.070220503698695 );

            wgts.set(1,  0.0721578038388935 );
            wgts.set(2,  0.0475458171336425 );
            wgts.set(3,  0.0475458171336425 );
            wgts.set(4,  0.0475458171336425 );
            wgts.set(5,  0.0516086852673590 );
            wgts.set(6,  0.0516086852673590 );
            wgts.set(7,  0.0516086852673590 );
            wgts.set(8,  0.0162292488115990 );
            wgts.set(9,  0.0162292488115990 );
            wgts.set(10, 0.0162292488115990 );
            wgts.set(11, 0.0136151570872175 );
            wgts.set(12, 0.0136151570872175 );
            wgts.set(13, 0.0136151570872175 );
            wgts.set(14, 0.0136151570872175 );
            wgts.set(15, 0.0136151570872175 );
            wgts.set(16, 0.0136151570872175 );
            break;
    }


}
