#include "Legendre2d.h"

// Set quadrature weights and points
// This method was moved from L2ProjectGrad.cpp
void SetQuadrature(const int mpoints1d, dTensor1& wgt, dTensor2& spts)
{
    {
        const int mpoints = mpoints1d*mpoints1d;
        assert_eq(mpoints, wgt.getsize());
        assert_eq(mpoints, spts.getsize(1));
    }

    // Construct the 1D quadrature weights and points
    dTensor1 w1d(mpoints1d), x1d(mpoints1d);
    setGaussPoints1d(w1d, x1d);

    // Tensor product Gaussian Quadrature
    int k=0;
    for (int m1=1; m1<=(mpoints1d); m1++)
        for (int m2=1; m2<=(mpoints1d); m2++)
        {
            k = k+1;
            wgt.set(k,  w1d.get(m1)*w1d.get(m2) );

            spts.set(k,1, x1d.get(m1) );
            spts.set(k,2, x1d.get(m2) );
        }

}

void Legendre2d::SetLegendreIntervals()
{

}


void SetLegendreAtPoints(const dTensor2& spts, dTensor2& phi)
{
    const int space_order = dogParams.get_space_order();
    const int kmax = phi.getsize(2);
    //assert_eq(kmax,dogParams.get_kmax());
    const int mpoints = spts.getsize(1);
    assert_eq(mpoints,phi.getsize(1));
    for (int m=1; m<=mpoints; m++)
    {
        // grid point (x,y)
        const double xi  = spts.get(m,1);
        const double eta = spts.get(m,2);
        const double xi2 = xi*xi;
        const double xi3 = xi*xi2;
        const double xi4 = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     


        // Legendre basis functions at each gaussian quadrature point in the
        // interval [-1,1]x[-1,1].
        switch( kmax )
        {
            case 25:  // fifth order                                 
                phi.set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
                phi.set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq7*(2.5*xi3 - 1.5*xi) );
                phi.set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq5*(1.5*xi2 - 0.5) );
                phi.set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq5*(1.5*eta2 - 0.5) );
                phi.set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*xi );
                phi.set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*eta );
                phi.set( m,18, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi.set( m,17, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );

            case 16:  // fourth order
                phi.set( m,16, sq7*(2.5*eta3 - 1.5*eta)*sq7*(2.5*xi3 - 1.5*xi) );
                phi.set( m,15, sq7*(2.5*eta3 - 1.5*eta)*sq5*(1.5*xi2 - 0.5) );
                phi.set( m,14, sq7*(2.5*xi3 - 1.5*xi)*sq5*(1.5*eta2 - 0.5) );
                phi.set( m,13, sq7*(2.5*eta3 - 1.5*eta)*sq3*xi );
                phi.set( m,12, sq7*(2.5*xi3 - 1.5*xi)*sq3*eta );
                phi.set( m,11, sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( m,10,  sq7*(2.5*xi3 - 1.5*xi) );

            case 9:  // third order
                phi.set( m,9, 5.0*(1.5*xi2 - 0.5)*(1.5*eta2 - 0.5) );
                phi.set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi.set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );
                phi.set( m,6,  sq5*(1.5*eta2 - 0.5) );
                phi.set( m,5,  sq5*(1.5*xi2 - 0.5) );

            case 4:  // second order                
                phi.set( m,4,  3.0*xi*eta );
                phi.set( m,3, sq3*eta );
                phi.set( m,2, sq3*xi  );

            case 1:  // first order
                phi.set( m,1, 1.0 );

                break;

            default:
                unsupported_value_error(kmax);
        }


    }
}

void set_divfree_polys(int mpoints1d,
        const dTensor2& spts,
        dTensor3& phi_divfree)
{
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double dx2 = dx*dx;
    const double dy2 = dy*dy;
    const double dx3 = dx*dx2;
    const double dy3 = dy*dy2;
    const double dx4 = dx*dx3;
    const double dy4 = dy*dy3;
    const double dx5 = dx*dx4;
    const double dy5 = dy*dy4;
    const double dx6 = dx*dx5;
    const double dy6 = dy*dy5;
    const double length = sqrt(dx*dx + dy*dy);
    const double lx5 = sqrt(5.0*dx*dx + dy*dy);
    const double ly5 = sqrt(dx*dx + 5.0*dy*dy);
    const double l4  = sqrt(3.0*dx4 + 45.0*dx2*dy2 + 35.0*dy4);
    const double l6  = sqrt(7.0*dx6 + 107.0*dx4*dy2 + 107.0*dy4*dx2 + 7.0*dy6);

    assert_eq(mpoints1d,dogParams.get_space_order());
    for (int m=1; m<=phi_divfree.getsize(1); m++)
    {
        // grid point (x,y)
        const double xi  = spts.get(m,1);
        const double eta = spts.get(m,2);
        const double xi2 = xi*xi;
        const double xi3 = xi*xi2;
        const double xi4 = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;

        // Divergence-free Legendre basis functions at each
        // gaussian quadrature point in the interval [-1,1]x[-1,1].
        switch( mpoints1d )
        {
            default:
                unsupported_value_error(mpoints1d);
            case 5:  // fifth order (6 basis vectors)
                eprintf("L2Project_divfree has not been implemented"
                        " for mpoints1d = %d", mpoints1d);
            case 4:  // fourth order (5 basis vectors)
                phi_divfree.set(m,14,1,  0.5*sq7*(5.0*eta3-3.0*eta)  );
                phi_divfree.set(m,14,2,  0.0  );

                phi_divfree.set(m,13,1,  0.0  );
                phi_divfree.set(m,13,2, -0.5*sq7*(5.0*xi3-3.0*xi) );

                phi_divfree.set(m,12,1,  (0.5*sq3*sq5)*(dx/length)*eta*(3.0*xi2-1.0) );
                phi_divfree.set(m,12,2, -(0.5*sq3*sq5)*(dy/length)*xi*(3.0*eta2-1.0) );

                phi_divfree.set(m,11,1,  (0.5*(sq5*sq7*dx*xi)/(l4*l6))
                        *((-21.0+105.0*eta2)*dy4+(-48.0+7.0*xi2+135.0*eta2)*dx2*dy2
                            +(9.0*eta2-3.0)*dx4) );
                phi_divfree.set(m,11,2, -(0.5*(sq5*sq7*dy*eta)/(l4*l6))
                        *((-21.0+35.0*eta2)*dy4+(-48.0+21.0*xi2+45.0*eta2)*dx2*dy2
                            +(-3.0+3.0*eta2)*dx4) );

                phi_divfree.set(m,10,1,  0.5*((sq3*sq7*dx)/(length*l4))*xi
                        *(dy2*(5.0*xi2-5.0)+dx2*(5.0*xi2-3.0)) );
                phi_divfree.set(m,10,2, -0.5*((sq3*sq7*dy)/(length*l4))*eta
                        *(dy2*(15.0*xi2-5.0)+dx2*(15.0*xi2-3.0)) );

            case 3:  // third order (4 basis vectors)
                phi_divfree.set(m,9,1,  0.5*sq5*(3.0*eta2-1.0) );
                phi_divfree.set(m,9,2,  0.0 );

                phi_divfree.set(m,8,1,  0.0 );
                phi_divfree.set(m,8,2, -0.5*sq5*(3.0*xi2-1.0) );

                phi_divfree.set(m,7,1,  3.0*sq5*xi*eta*(dx/lx5) );
                phi_divfree.set(m,7,2, -sq5/2.0*(3.0*eta2-1.0)*(dy/lx5) );

                phi_divfree.set(m,6,1,  sq5/2.0*(3.0*xi2-1.0)*(dx/ly5) );
                phi_divfree.set(m,6,2, -3.0*sq5*xi*eta*(dy/ly5) );

            case 2:  // second order (3 basis vectors)
                phi_divfree.set(m,5,1,  sq3*eta );
                phi_divfree.set(m,5,2,  0.0 );

                phi_divfree.set(m,4,1,  0.0 );
                phi_divfree.set(m,4,2, -sq3*xi );

                phi_divfree.set(m,3,1,  dx/length * sq3*xi  );
                phi_divfree.set(m,3,2, -dy/length * sq3*eta );

            case 1:  // first order (2 basis vectors)
                phi_divfree.set(m,2,1,  1.0 );
                phi_divfree.set(m,2,2,  0.0 );

                phi_divfree.set(m,1,1,  0.0 );
                phi_divfree.set(m,1,2, -1.0 );

                break;
        }
    }
}


// Set values of the Legendre polys on each quadrature point
//
void SetLegendre(const int dummy_space_order,
        const dTensor2& spts,
        dTensor2& phi)
{
    const int space_order = dogParams.get_space_order();
    assert_eq(dummy_space_order,space_order);
    SetLegendreAtPoints(spts, phi);
}

// Set values of the Gradient of the Legendre polys on each quadrature point
void SetLegendreGrad(const double dx, const double dy,
        const dTensor2& spts, dTensor2& phi_x, dTensor2& phi_y)
{
    const int space_order = dogParams.get_space_order();
    const int mpoints = spts.getsize(1);
    const int kmax    = dogParams.get_kmax();

    // quick error check
    {
        assert_eq(phi_x.getsize(1), mpoints);
        assert_eq(kmax, phi_x.getsize(2) );
        assert_eq(kmax, phi_y.getsize(2) );
        const int mpoints1d = int(sqrt(mpoints));//space_order-1;
        assert_eq(mpoints, mpoints1d*mpoints1d);
    }

    //
    const double tmpx = 2.0/dx;
    const double tmpy = 2.0/dy;

    for (int m=1; m<=(mpoints); m++)
    {
        // grid point (x,y)
        const double xi   = spts.get(m,1);
        const double eta  = spts.get(m,2);
        const double xi2  = xi*xi;
        const double xi3  = xi*xi2;
        const double xi4  = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     

        // Gradient of Legendre basis functions at each gaussian quadrature point
      // Gradient of Legendre basis functions at each gaussian quadrature point
      switch( kmax )
        {
	case 25:  // fifth order
          phi_x.set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi) );
          phi_x.set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*tmpx*sq7*(7.5*xi2 - 1.5) );
          phi_x.set( m,23, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi)*sq7*(2.5*eta3 - 1.5*eta) );
          phi_x.set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*tmpx*sq5*3.0*xi );
          phi_x.set( m,21, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi)*sq5*(1.5*eta2 - 0.5) );
          phi_x.set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*tmpx );
          phi_x.set( m,19, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi)*sq3*eta );
          phi_x.set( m,18, 0.0 );
          phi_x.set( m,17, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi) );


          phi_y.set( m,25, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
          phi_y.set( m,24, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*sq7*(2.5*xi3 - 1.5*xi) );
          phi_y.set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*tmpy*sq7*(7.5*eta2 - 1.5) );
          phi_y.set( m,22, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*sq5*(1.5*xi2 - 0.5) );
          phi_y.set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*tmpy*sq5*3.0*eta );
          phi_y.set( m,20, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*sq3*xi );
          phi_y.set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*tmpy );
          phi_y.set( m,18, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta) );
          phi_y.set( m,17, 0.0 );

	case 16:  // fourth order
          phi_x.set( m,16, sq7*(2.5*eta3 - 1.5*eta)*tmpx*sq7*(7.5*xi2 - 1.5)  );
          phi_x.set( m,15, sq7*(2.5*eta3 - 1.5*eta)*tmpx*sq5*3.0*xi );
          phi_x.set( m,14, tmpx*sq7*(7.5*xi2 - 1.5) *sq5*(1.5*eta2 - 0.5) );
          phi_x.set( m,13, sq7*(2.5*eta3 - 1.5*eta)*sq3*tmpx );
          phi_x.set( m,12, tmpx*sq7*(7.5*xi2 - 1.5) *sq3*eta );
          phi_x.set( m,11, 0.0 );
          phi_x.set( m,10,  tmpx*sq7*(7.5*xi2 - 1.5)  );


          phi_y.set( m,16, tmpy*sq7*(7.5*eta2 - 1.5)*sq7*(2.5*xi3 - 1.5*xi) );
          phi_y.set( m,15, tmpy*sq7*(7.5*eta2 - 1.5)*sq5*(1.5*xi2 - 0.5) );
          phi_y.set( m,14, sq7*(2.5*xi3 - 1.5*xi)*tmpy*sq5*3.0*eta );
          phi_y.set( m,13, tmpy*sq7*(7.5*eta2 - 1.5)*sq3*xi );
          phi_y.set( m,12, sq7*(2.5*xi3 - 1.5*xi)*sq3*tmpy );
          phi_y.set( m,11, tmpy*sq7*(7.5*eta2 - 1.5) );
          phi_y.set( m,10,  0.0 );

	case 9:  // third order
          phi_x.set( m,9,  sq5*tmpx*3.0*xi*sq5*(1.5*eta2 - 0.5) );
          phi_x.set( m,8,  sq3*sq5*tmpx*(1.5*eta2 - 0.5) );
          phi_x.set( m,7,  sq3*eta*tmpx*sq5*3.0*xi );
          phi_x.set( m,6,  0.0 );
          phi_x.set( m,5,  tmpx*sq5*3.0*xi );


          phi_y.set( m,9,  sq5*(1.5*xi2 - 0.5)*tmpy*sq5*3.0*eta );
          phi_y.set( m,8,  sq3*xi*tmpy*sq5*3.0*eta );
          phi_y.set( m,7,  sq3*sq5*tmpy*(1.5*xi2 - 0.5) );
          phi_y.set( m,6,  tmpy*sq5*3.0*eta );
          phi_y.set( m,5,  0.0 );

	case 4:  // second order
          phi_x.set( m,4,  3.0*tmpx*eta );
          phi_x.set( m,3, 0.0 );
          phi_x.set( m,2, sq3*tmpx  );

          phi_y.set( m,4,  3.0*xi*tmpy );
          phi_y.set( m,3, sq3*tmpy );
          phi_y.set( m,2, 0.0  );

	case 1:  // first order
	  phi_x.set( m, 1,  0.0 );
	  
	  phi_y.set( m, 1,  0.0 );
	  
	  break;
	  
	default:
	  unsupported_value_error(kmax);
        }
    }//end of constructing legendre polynomials
}

#if 0
// This is no longer used and does not define the appropriate
// set of points for maintaining positivity.
// (It formerly was used in ApplyPosLimiter.cpp)
// moved this from L2Project_divfree.cpp
//
void L2Project_set_gauss_lobatto_points(
        int mpoints, int mpoints1d, dTensor2& spts)
{
    assert_eq(mpoints1d, dogParams.get_space_order());
    assert_eq(mpoints,mpoints1d*mpoints1d);

    dTensor1 x1d(mpoints1d);
    setGaussLobattoPoints1d(x1d);

    // Tensor product Gaussian Quadrature
    int k=0;
    for (int m1=1; m1<=(mpoints1d); m1++)
        for (int m2=1; m2<=(mpoints1d); m2++)
        {
            k++;

            //save gauss quad grid point location on interval [-1,1]^2
            spts.set(k,1, x1d.get(m1) );
            spts.set(k,2, x1d.get(m2) );
        }
}
#endif

// for use after applying limiters
//
// should also create a linear projection
// project_linearly_onto_locally_divergence_free_subspace
// which uses averages instead of minmod.
// This could be used to project the electric field onto the
// affine subspace of functions whose divergence agrees with the
// charge density over permittivity (with highest coefficient
// zeroed).
//
// (moved from ApplyLimiter_divfree.cpp)
void project_onto_locally_divergence_free_subspace(dTensorBC4& q)
{
    const int mx = q.getsize(1);
    const int my = q.getsize(2);
    const int mbc  = q.getmbc();
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double length = sqrt(dx*dx + dy*dy);
    const double lx5    = sqrt(5.0*dx*dx + dy*dy);
    const double ly5    = sqrt(dx*dx + 5.0*dy*dy);
    const int how_many = dogParams.get_how_many_vectors_divfree();

    //#pragma omp parallel for
    for (int k=1; k<=how_many; k++)
    {
        int mcomp = dogParams.get_which_compnt_divfree()[k];

        switch ( dogParams.get_space_order() )
        {
            default:
                unsupported_value_error(dogParams.get_space_order());
            case 1:
                break;
            case 2:  // 2nd order in space

                for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
                    for (int j=(3-mbc); j<=(my+mbc-2); j++)
                    {
                        const double B3 = minmod(q.get(i,j,mcomp,2)/dx,
                                -q.get(i,j,mcomp+1,3)/dy);

                        q.set(i,j,mcomp,  2,  dx*B3 );
                        q.set(i,j,mcomp+1,3, -dy*B3 );
                    }

                break;

            case 3:  // 3rd order in space

                for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
                    for (int j=(3-mbc); j<=(my+mbc-2); j++)
                    {
                        double B3 = minmod(q.get(i,j,mcomp,2)/dx,       -q.get(i,j,mcomp+1,3)/dy);
                        double B6 = minmod(q.get(i,j,mcomp,5)/dx,       -q.get(i,j,mcomp+1,4)/(sq5*dy));
                        double B7 = minmod(q.get(i,j,mcomp,4)/(sq5*dx), -q.get(i,j,mcomp+1,6)/dy);

                        q.set(i,j,mcomp,  2,      dx*B3 );
                        q.set(i,j,mcomp+1,3,     -dy*B3 );

                        q.set(i,j,mcomp,  5,      dx*B6 );
                        q.set(i,j,mcomp+1,4, -sq5*dy*B6 );

                        q.set(i,j,mcomp,  4,  sq5*dx*B7 );
                        q.set(i,j,mcomp+1,6,     -dy*B7 );          
                    }

                break;
        }
    }
}

Legendre2d::~Legendre2d()
{
    delete quadratureWeights;
    delete quadraturePoints;
    delete phi;
    delete wgt_phi_transpose;
    delete phi_min_ptr;
    delete phi_max_ptr;
    delete phi_interval;
}

// To maintain positivity in two dimensions we enforce positivity
// at the tensor product of the gauss-lobatto points with the
// Gaussian quadrature points where the Riemann solver is used.
// It is not necessary to use the gauss-lobatto points in any
// integrals, so the weights are not needed. See Theorem 2.1 in
// [ZhangShu10]:
//
// [ZhangShu10] X. Zhang, C-W. Shu, On positivity preserving high
// order discontinuous Galerkin schemes for compressible Euler
// equations on rectangular meshes, Journal of Computational
// Physics (2010), doi: 10.1016/j.jcp.2010.08.016
void Legendre2d::setPositivityPoints()
{

    const int polynomial_order = dogParams.get_space_order()-1;
    const int mGaussLobattoPoints1d = int(ceil((polynomial_order+3)/2.));
    const int mGaussPoints1d = dogParams.get_space_order();

    // Grab the 1D points
    dTensor1 gaussLobattoPoints1d   ( mGaussLobattoPoints1d  );
    dTensor1 gaussPoints1d          ( mGaussPoints1d         );
    setGaussLobattoPoints1d         ( gaussLobattoPoints1d   );
    setGaussPoints1d                ( gaussPoints1d          );

    const int mpoints = mGaussLobattoPoints1d*mGaussPoints1d;
    // will enforce positivity at two meshes.
    numPositivityPoints = 2*mpoints;
    positivityPoints = new dTensor2(numPositivityPoints,2);
    //
    int k1=0;
    for(int m1=1;m1<=mGaussPoints1d;m1++)
        for(int m2=1;m2<=mGaussLobattoPoints1d;m2++)
        {
            k1++;
            positivityPoints->set(k1,1, gaussPoints1d.get(m1));
            positivityPoints->set(k1,2, gaussLobattoPoints1d.get(m2));
        }
    assert_eq(k1,mpoints);
    int k2=mpoints;
    for(int m1=1;m1<=mGaussPoints1d;m1++)
        for(int m2=1;m2<=mGaussLobattoPoints1d;m2++)
        {
            k2++;
            positivityPoints->set(k2,1, gaussLobattoPoints1d.get(m2));
            positivityPoints->set(k2,2, gaussPoints1d.get(m1));
        }
    assert_eq(k2,2*mpoints);

    phiAtPositivityPoints = new dTensor2(positivityPoints->getsize(1),
            dogParams.get_kmax());
    SetLegendreAtPoints(*positivityPoints,*phiAtPositivityPoints);

    // print out coordinates and phi values at Riemann points
    // (compare with setRiemannPoints() and edge_data::init()).
    if(debug3)
    {
        const int kmax = dogParams.get_kmax();
        int k1=0;
        printf("\n points at y boundaries:");
        for(int m1=1;m1<=mGaussPoints1d;m1++)
            for(int m2=1;m2<=mGaussLobattoPoints1d;m2++)
            {
                k1++;
                if(m2==1||m2==mGaussLobattoPoints1d)
                {
                    printf("\n   phi(%d,*) at (%24.16e, %24.16e ):",k1,
                            positivityPoints->get(k1,1),
                            positivityPoints->get(k1,2));
                    for(int k=1;k<=kmax;k++)
                    {
                        printf("\n    phi(%d,%d) = %24.16e",
                                k1,k,phiAtPositivityPoints->get(k1,k));
                    }
                    printf("\n");
                }
            }
        assert_eq(k1,mpoints);
        int k2=mpoints;
        printf("\n points at x boundaries:");
        for(int m1=1;m1<=mGaussPoints1d;m1++)
            for(int m2=1;m2<=mGaussLobattoPoints1d;m2++)
            {
                k2++;
                if(m2==1||m2==mGaussLobattoPoints1d)
                {
                    printf("\n   phi(%d,*) at (%24.16e, %24.16e ):",k2,
                            positivityPoints->get(k2,1),
                            positivityPoints->get(k2,2));
                    for(int k=1;k<=kmax;k++)
                    {
                        printf("\n    phi(%d,%d) = %24.16e",
                                k2,k,phiAtPositivityPoints->get(k2,k));
                    }
                    printf("\n");
                }
            }
        assert_eq(k2,2*mpoints);
    }
}

void Legendre2d::setRiemannPoints()
{
    const int mGaussPoints1d = dogParams.get_space_order();

    dTensor1 gaussPoints1d(mGaussPoints1d);
    setGaussPoints1d(gaussPoints1d);

    numRiemannPoints = mGaussPoints1d*4;
    riemannPoints = new dTensor2(numRiemannPoints,2);
    //
    int k1=0;
    for(int m1=1;m1<=mGaussPoints1d;m1++) // yr
    {
        k1++;
        riemannPoints->set(k1,1, gaussPoints1d.get(m1));
        riemannPoints->set(k1,2, -1.);
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // yl
    {
        k1++;
        riemannPoints->set(k1,1, gaussPoints1d.get(m1));
        riemannPoints->set(k1,2, 1.);
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // xr
    {
        k1++;
        riemannPoints->set(k1,1, -1.);
        riemannPoints->set(k1,2, gaussPoints1d.get(m1));
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // xl
    {
        k1++;
        riemannPoints->set(k1,1, 1.);
        riemannPoints->set(k1,2, gaussPoints1d.get(m1));
    }
    assert_eq(k1,numRiemannPoints);

    phiAtRiemannPoints = new dTensor2(riemannPoints->getsize(1),
            dogParams.get_kmax());
    SetLegendreAtPoints(*riemannPoints,*phiAtRiemannPoints);

    // print out coordinates and phi values at Riemann points
    // (compare with setPositivityPoints())
    if(debug3)
    {
        const int kmax = dogParams.get_kmax();
        for(int m1=1;m1<=numRiemannPoints;m1++)
        {
            printf("\n   phi(%d,*) at (%24.16e, %24.16e ):",m1,
                    riemannPoints->get(m1,1),
                    riemannPoints->get(m1,2));
            for(int k=1;k<=kmax;k++)
            {
                printf("\n    phi(%d,%d) = %24.16e",
                        m1,k,phiAtRiemannPoints->get(m1,k));
            }
            printf("\n");
        }
    }
}

Legendre2d::Legendre2d():
    quadratureWeights(0),
    quadraturePoints(0),
    phi(0),
    wgt_phi_transpose(0),
    phi_min_ptr(0),
    phi_max_ptr(0),
    phi_interval(0),
    quadratureWeightsCoarse(0), quadraturePointsCoarse(0), phiCoarse(0),
    phi_x(0), phi_y(0),
    numPositivityPoints(0), positivityPoints(0), phiAtPositivityPoints(0),
    numRiemannPoints(0), riemannPoints(0), phiAtRiemannPoints(0),
    edgeData(0)
{
    const int space_order = dogParams.get_space_order();
    const int kmax = dogParams.get_kmax();

    // set quadrature for integration accurate to order space_order
    //
    const int mpoints1d = space_order;
    const int mpoints = mpoints1d*mpoints1d;
    quadratureWeights = new dTensor1(mpoints);
    quadraturePoints = new dTensor2(mpoints,2);
    phi = new dTensor2(mpoints,kmax);
    SetQuadrature(mpoints1d,*quadratureWeights,*quadraturePoints);
    SetLegendre(space_order,*quadraturePoints,*phi);
    wgt_phi_transpose = new dTensor2(kmax,mpoints);
    for(int mp=1;mp<=mpoints;mp++)
        for(int k=1;k<=kmax;k++)
        {
            wgt_phi_transpose->set(k,mp, quadratureWeights->get(mp)*phi->get(mp,k));
        }

    // SHOULD ONLY BE CALLED IF NEEDED.
    // THIS FUNCTIONALITY IS NO LONGER USED ANYWAY.
    /*
       const int kmax_divfree = dogParams.get_kmax_divfree();
       phi_divfree = new dTensor3(mpoints,kmax_divfree,2);
       set_divfree_polys(mpoints1d,*quadraturePoints, *phi_divfree);
     */

    // set quadrature for integration accurate to order space_order-1
    //
    //const int mpoints1dCoarse = iMax(1,space_order-1);
    const int mpoints1dCoarse = space_order-1;
    const int mpointsCoarse = mpoints1dCoarse*mpoints1dCoarse;
    quadratureWeightsCoarse = new dTensor1(mpointsCoarse);
    quadraturePointsCoarse = new dTensor2(mpointsCoarse,2);
    phiCoarse = new dTensor2(mpointsCoarse,kmax);
    SetQuadrature(mpoints1dCoarse,*quadratureWeightsCoarse,*quadraturePointsCoarse);
    SetLegendre(space_order,*quadraturePointsCoarse,*phiCoarse);
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    phi_x = new dTensor2(mpointsCoarse,kmax);
    phi_y = new dTensor2(mpointsCoarse,kmax);
    SetLegendreGrad(dx,dy,*quadraturePointsCoarse,*phi_x,*phi_y);

    edgeData = new edge_data;
    edgeData->init();
    //void SetEdgeData(int morder, int kmax, edge_data& EdgeData);
    //SetEdgeData(space_order,kmax,*edgeData);

    //SetLegendreInfinityNorms(phi_maxabs);
    SetLegendreIntervals();

    setPositivityPoints();

    setRiemannPoints();

    // check that positivity points concur with edgeData
    check_edgeData();
};

void Legendre2d::check_edgeData() const
{
    const int mGaussPoints1d = dogParams.get_space_order();
    const dTensor2& phi = *phiAtRiemannPoints;
    const int kmax = dogParams.get_kmax();
    int k1=0;
    for(int m1=1;m1<=mGaussPoints1d;m1++) // yr
    { k1++;
        for(int k=1;k<=kmax;k++)
            {assert_eq(phi.get(k1,k), edgeData->phi_yr->get(m1,k));}
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // yl
    { k1++;
        for(int k=1;k<=kmax;k++)
            {assert_eq(phi.get(k1,k), edgeData->phi_yl->get(m1,k));}
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // xr
    { k1++;
        for(int k=1;k<=kmax;k++)
            assert_eq(phi.get(k1,k), edgeData->phi_xr->get(m1,k));
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // xl
    { k1++;
        for(int k=1;k<=kmax;k++)
            assert_eq(phi.get(k1,k), edgeData->phi_xl->get(m1,k));
    }
    assert_eq(k1,numRiemannPoints);
}

Legendre2d& Legendre2d::instance() // const
{
    static Legendre2d* legendre2d = new Legendre2d;
    return *legendre2d;
}
