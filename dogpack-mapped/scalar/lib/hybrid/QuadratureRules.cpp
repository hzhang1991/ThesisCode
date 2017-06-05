#include "QuadratureRules.h"

// Desctructor (remove any allocated memory)
QuadratureRules::~QuadratureRules()
{

    delete wgt_unst;
    delete wgt_cart;
    delete wgt;

    delete spts_unst;
    delete spts_cart;
    delete spts;

    delete phi_unst;
    delete phi_cart;
    delete phi;

}

void QuadratureRules::init( int sorder )
{

    // Order of quadrature rules
    const int QuadOrder = sorder;
    assert_ge(QuadOrder,1);
    assert_le(QuadOrder,5);

    this->sorder = sorder;

    int num_quad_cart[] = {1,4,9,16,25};  // == M^2
    int num_quad_unst[] = {1,3,6,12,16};  // slightly smaller than M^2
    int kmax2d_vec[]    = {1,3,6,10,15};  // number of polynomials for each 2D problem

    // number of polynomials for each 2D problem
    this->kmax2d = kmax2d_vec[ sorder-1 ];

    // --------------------------------------------------------------------- //
    //
    // Set up the 2D, Cartesian, quadrature points, weights and polynomials
    //
    // --------------------------------------------------------------------- //

    const int s2 = sorder*sorder;
    const int mpoints1d = sorder;

    // Set up the Cartesian points and weights:
    this->spts_cart = new dTensor2 (s2, 2 );
    this->wgt_cart  = new dTensor1 (s2 );

    const int mpoints_cart = s2;
    const int md2 = sorder/2;

    // 1D Quadrature points (Gauss-Legendre)
    dTensor1 w1d(mpoints1d),x1d(mpoints1d);
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);
    setGaussPoints1d( w1d, x1d );

    // Tensor product Gaussian Quadrature
    int k=0;
    for (int m1=1; m1<=(mpoints1d); m1++)
    for (int m2=1; m2<=(mpoints1d); m2++)
    {
            k = k+1;
            wgt_cart->set(k,  w1d.get(m1)*w1d.get(m2) );

            spts_cart->set(k,1, x1d.get(m1) );
            spts_cart->set(k,2, x1d.get(m2) );
    }

    // Evaluate the cartesian basis functions at each quadrature point
    phi_cart = new dTensor2( mpoints_cart, kmax2d );
    for (int m=1; m<=mpoints_cart; m++)
    {
        // coordinates
        const double xi   = spts_cart->get(m,1);      
        const double xi2  = xi*xi;
        const double xi3  = xi2*xi;
        const double xi4  = xi3*xi;
        const double eta  = spts_cart->get(m,2);
        const double eta2 = eta*eta;
        const double eta3 = eta2*eta;
        const double eta4 = eta3*eta;      

        // Basis Functions evaluated at each quadrature point:
        switch( this->kmax2d )
        {
            case 15:  // fifth order                                 
                phi_cart->set( m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi_cart->set( m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
                phi_cart->set( m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
                phi_cart->set( m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                phi_cart->set( m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

            case 10:  // fourth order
                phi_cart->set( m,10, sq7*(2.5*eta3 - 1.5*eta) );
                phi_cart->set( m,9,  sq7*(2.5*xi3 - 1.5*xi) );
                phi_cart->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi_cart->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

            case 6:  // third order
                phi_cart->set( m,6,  sq5*(1.5*eta2 - 0.5) );
                phi_cart->set( m,5,  sq5*(1.5*xi2 - 0.5) );
                phi_cart->set( m,4,  3.0*xi*eta );                  

            case 3:  // second order                
                phi_cart->set( m,3, sq3*eta );
                phi_cart->set( m,2, sq3*xi  );

            case 1:  // first order
                phi_cart->set( m,1, 1.0 );
                break;

            default:
                unsupported_value_error(this->kmax2d);
        }
        
    }

    // --------------------------------------------------------------------- //
    //
    // Set up the 2D, unstructured, quadrature points, weights and polynoms
    //
    // --------------------------------------------------------------------- //

    // Set up the unstructured points and weights
    int mpoints_unst;
    switch ( QuadOrder )
    {
        case 1:
            mpoints_unst = 1;
            break;

        case 2:
            mpoints_unst = 3;
            break;

        case 3:
            mpoints_unst = 6;
            break;

        case 4:
            mpoints_unst = 12;
            break;

        case 5:	     
            mpoints_unst = 16;
            break;
    }

    // monomial basis (non-orthogonal)
    dTensor2      mu(mpoints_unst, this->kmax2d ); 

    spts_unst = new dTensor2( mpoints_unst, 2 );
    wgt_unst  = new dTensor1( mpoints_unst );
    phi_unst  = new dTensor2( mpoints_unst, this->kmax2d );

    // ----------------------------------------------
    // Set unstructured quadrature weights and points
    // ----------------------------------------------
    switch ( QuadOrder )
    {

        case 1:
            spts_unst->set(1,1, 0.0 );
            spts_unst->set(1,2, 0.0 );

            wgt_unst->set(1,    0.5 );
            break;

        case 2:
            spts_unst->set(1,1,  1.0/3.0 );
            spts_unst->set(1,2, -1.0/6.0 );

            spts_unst->set(2,1, -1.0/6.0 );
            spts_unst->set(2,2, -1.0/6.0 );

            spts_unst->set(3,1, -1.0/6.0 );
            spts_unst->set(3,2,  1.0/3.0 );

            wgt_unst->set(1, 1.0/6.0 );
            wgt_unst->set(2, 1.0/6.0 );
            wgt_unst->set(3, 1.0/6.0 );
            break;

        case 3:
            spts_unst->set(1,1,  0.112615157582632 );
            spts_unst->set(1,2,  0.112615157582632 );

            spts_unst->set(2,1, -0.225230315165263 );
            spts_unst->set(2,2,  0.112615157582632 );

            spts_unst->set(3,1,  0.112615157582632 );
            spts_unst->set(3,2, -0.225230315165263 );

            spts_unst->set(4,1, -0.241757119823562 );
            spts_unst->set(4,2, -0.241757119823562 );

            spts_unst->set(5,1,  0.483514239647126 );
            spts_unst->set(5,2, -0.241757119823562 );

            spts_unst->set(6,1, -0.241757119823562 );
            spts_unst->set(6,2,  0.483514239647126 );

            wgt_unst->set(1, 0.1116907948390055 );
            wgt_unst->set(2, 0.1116907948390055 );
            wgt_unst->set(3, 0.1116907948390055 );
            wgt_unst->set(4, 0.0549758718276610 );
            wgt_unst->set(5, 0.0549758718276610 );
            wgt_unst->set(6, 0.0549758718276610 );
            break;

        case 4:
            spts_unst->set(1,1,  -0.084046588162423 );
            spts_unst->set(1,2,  -0.084046588162423 );

            spts_unst->set(2,1,   0.168093176324846 );
            spts_unst->set(2,2,  -0.084046588162423 );

            spts_unst->set(3,1,  -0.084046588162423 );
            spts_unst->set(3,2,   0.168093176324846 );

            spts_unst->set(4,1,  -0.270244318841831 );
            spts_unst->set(4,2,  -0.270244318841831 );

            spts_unst->set(5,1,   0.540488637683663 );
            spts_unst->set(5,2,  -0.270244318841831 );

            spts_unst->set(6,1,  -0.270244318841831 );
            spts_unst->set(6,2,   0.540488637683663 );

            spts_unst->set(7,1,  -0.280188283488516 );
            spts_unst->set(7,2,  -0.022980882299549 );

            spts_unst->set(8,1,  -0.280188283488516 );
            spts_unst->set(8,2,   0.303169165788066 );

            spts_unst->set(9,1,  -0.022980882299549 );
            spts_unst->set(9,2,   0.303169165788067 );

            spts_unst->set(10,1, -0.022980882299549 );
            spts_unst->set(10,2, -0.280188283488516 );

            spts_unst->set(11,1,  0.303169165788066 );
            spts_unst->set(11,2, -0.022980882299549 );

            spts_unst->set(12,1,  0.303169165788066 );
            spts_unst->set(12,2, -0.280188283488516 );

            wgt_unst->set(1,  0.0583931378631895 );
            wgt_unst->set(2,  0.0583931378631895 );
            wgt_unst->set(3,  0.0583931378631895 );
            wgt_unst->set(4,  0.0254224531851035 );
            wgt_unst->set(5,  0.0254224531851035 );
            wgt_unst->set(6,  0.0254224531851035 );
            wgt_unst->set(7,  0.0414255378091870 );
            wgt_unst->set(8,  0.0414255378091870 );
            wgt_unst->set(9,  0.0414255378091870 );
            wgt_unst->set(10, 0.0414255378091870 );
            wgt_unst->set(11, 0.0414255378091870 );
            wgt_unst->set(12, 0.0414255378091870 );
            break;

        case 5:
            spts_unst->set(1,1,   0.000000000000000 );
            spts_unst->set(1,2,   0.000000000000000 );

            spts_unst->set(2,1,   0.125959254959390 );
            spts_unst->set(2,2,   0.125959254959390 );

            spts_unst->set(3,1,  -0.251918509918779 );
            spts_unst->set(3,2,   0.125959254959390 );

            spts_unst->set(4,1,   0.125959254959390 );
            spts_unst->set(4,2,  -0.251918509918779 );

            spts_unst->set(5,1,  -0.162764025581573 );
            spts_unst->set(5,2,  -0.162764025581573 );

            spts_unst->set(6,1,   0.325528051163147 );
            spts_unst->set(6,2,  -0.162764025581573 );

            spts_unst->set(7,1,  -0.162764025581573 );
            spts_unst->set(7,2,   0.325528051163147 );

            spts_unst->set(8,1,  -0.282786105016302 );
            spts_unst->set(8,2,  -0.282786105016302 );

            spts_unst->set(9,1,   0.565572210032605 );
            spts_unst->set(9,2,  -0.282786105016302 );

            spts_unst->set(10,1, -0.282786105016302 );
            spts_unst->set(10,2,  0.565572210032605 );

            spts_unst->set(11,1, -0.324938555923375 );
            spts_unst->set(11,2, -0.070220503698695 );

            spts_unst->set(12,1, -0.324938555923375 );
            spts_unst->set(12,2,  0.395159059622071 );

            spts_unst->set(13,1, -0.070220503698695 );
            spts_unst->set(13,2, -0.324938555923375 );

            spts_unst->set(14,1, -0.070220503698695 );
            spts_unst->set(14,2,  0.395159059622071 );

            spts_unst->set(15,1,  0.395159059622071 );
            spts_unst->set(15,2, -0.324938555923375 );

            spts_unst->set(16,1,  0.395159059622071 );
            spts_unst->set(16,2, -0.070220503698695 );

            wgt_unst->set(1,  0.0721578038388935 );
            wgt_unst->set(2,  0.0475458171336425 );
            wgt_unst->set(3,  0.0475458171336425 );
            wgt_unst->set(4,  0.0475458171336425 );
            wgt_unst->set(5,  0.0516086852673590 );
            wgt_unst->set(6,  0.0516086852673590 );
            wgt_unst->set(7,  0.0516086852673590 );
            wgt_unst->set(8,  0.0162292488115990 );
            wgt_unst->set(9,  0.0162292488115990 );
            wgt_unst->set(10, 0.0162292488115990 );
            wgt_unst->set(11, 0.0136151570872175 );
            wgt_unst->set(12, 0.0136151570872175 );
            wgt_unst->set(13, 0.0136151570872175 );
            wgt_unst->set(14, 0.0136151570872175 );
            wgt_unst->set(15, 0.0136151570872175 );
            wgt_unst->set(16, 0.0136151570872175 );
            break;
    }

    // Loop over each quadrature point and construct monomial polys
    for (int m=1; m<=mpoints_unst; m++)
    {
        // coordinates
        const double xi   = spts_unst->get(m,1);      
        const double xi2  = xi*xi;
        const double xi3  = xi2*xi;
        const double xi4  = xi3*xi;
        const double eta  = spts_unst->get(m,2);
        const double eta2 = eta*eta;
        const double eta3 = eta2*eta;
        const double eta4 = eta3*eta;      

        // monomials basis (non-orthogonal)
        switch( this->kmax2d )
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

    // Loop over each quadrature point and construct Legendre polys
    for (int m=1; m<=mpoints_unst; m++)    
    for (int i=1; i<=kmax2d; i++)
    {
        double tmp = 0.0;
        for (int j=1; j<=i; j++)
        {  tmp = tmp + Mmat[i-1][j-1]*mu.get(m,j);  }
        phi_unst->set(m,i, tmp );      
    }

    // --------------------------------------------------------------------- //
    //
    // Set up the 4D Quadrature points, and polynomials
    //
    // --------------------------------------------------------------------- //
    const int mpoints = mpoints_unst*mpoints_cart;

    int kmax_vec[] = {1,5,15};  // number of polynomials for each 4D problem
    kmax = kmax_vec[sorder-1];

    spts = new dTensor3( mpoints_unst, mpoints_cart, 4    );
    wgt  = new dTensor2( mpoints_unst, mpoints_cart       );
    phi  = new dTensor3( mpoints_unst, mpoints_cart, kmax );

    // Tensor product, Gaussian quadrature on weights and points:
    for( int m1=1; m1 <= mpoints_unst; m1++ )
    for( int m2=1; m2 <= mpoints_cart; m2++ )
    {

        wgt->set  (m1, m2, wgt_unst->get(m1)*wgt_cart->get(m2) );

        spts->set(m1, m2, 1, spts_unst->get(m1,1) );
        spts->set(m1, m2, 2, spts_unst->get(m1,2) );

        spts->set(m1, m2, 3, spts_cart->get(m2,1) );
        spts->set(m1, m2, 4, spts_cart->get(m2,2) );

    }

    // Evaluate the basis functions:
    for( int m1=1; m1 <= mpoints_unst; m1++ )
    for( int m2=1; m2 <= mpoints_cart; m2++ )
    {
        // 4D-coordinates
        //
        // The unstructured grid is on (x,y).  Canonical variables are mapped
        // from:
        //          (x,y) <-> (xi,eta);  
        //
        // The structured grid is on (vx,vy). Canonical variables are mapped
        // from:
        //
        //          (vx,vy) <-> (mu,tau);
        //
        const double xi   = spts->get(m1, m2, 1);      
        const double xi2  = xi*xi;

        const double eta  = spts->get(m1, m2, 2);
        const double eta2 = eta*eta;

        const double mu   = spts->get(m1, m2, 3);      
        const double mu2  = mu*mu;

        const double tau  = spts->get(m1, m2, 4);
        const double tau2 = tau*tau;

        switch( sorder )
        {
            case 3:  // third order

                // 'cartesian' terms
                phi->set(m1, m2,  15, phi_cart->get( m2, 6 ) );
                phi->set(m1, m2,  14, phi_cart->get( m2, 5 ) );
                phi->set(m1, m2,  13, phi_cart->get( m2, 4 ) );

                // cross (cartesian + unstructured) terms:
                phi->set(m1, m2,  12, 
                    phi_unst->get( m1, 3 ) * phi_cart->get(m2,3) );
                phi->set(m1, m2,  11, 
                    phi_unst->get( m1, 3 ) * phi_cart->get(m2,2) );
                phi->set(m1, m2,  10, 
                    phi_unst->get( m1, 2 ) * phi_cart->get(m2,3) );
                phi->set(m1, m2,  9, 
                    phi_unst->get( m1, 2 ) * phi_cart->get(m2,2) );

                // 'unstructured' terms
                phi->set(m1, m2, 8, phi_unst->get( m1, 6 ) );
                phi->set(m1, m2, 7, phi_unst->get( m1, 5 ) );
                phi->set(m1, m2, 6, phi_unst->get( m1, 4 ) );

            case 2:  // second order		    

                phi->set(m1, m2, 5, sq3*tau                  ); // phi_cart(3)
                phi->set(m1, m2, 4, sq3*mu                   ); // phi_cart(2)
                phi->set(m1, m2, 3, sq2*sq3*(xi + 2.0*eta)   ); // phi_unst(3)
                phi->set(m1, m2, 2, 3.0*sq2*xi               ); // phi_unst(2)

            case 1:  // first order

                phi->set(m1, m2, 1, 1.0       );
                break;

            default:
                unsupported_value_error(sorder);

        }

    }

}
