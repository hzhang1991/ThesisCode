#include "Legendre3d.h"

Legendre3d::Legendre3d():
    faceData(0)
{
    faceData = new FaceData;
    faceData->init();
};

Legendre3d::~Legendre3d()
{
}

Legendre3d& Legendre3d::instance()
{
    static Legendre3d* legendre3d = new Legendre3d;
    return *legendre3d;
}



// -----------------------------------------------------------------------
// FUNCTIONS USED TO SET QUADRATURE POINTS AND WEIGHTS, AS WELL AS 
// LEGENDRE FUNCTIONS AND DERIVATIVES AT QUADRATURE POINTS
// -----------------------------------------------------------------------

// Set quadrature weights and points
void SetQuadWgtsPts(const int mpoints1d, 
        dTensor1& wgt, 
        dTensor2& spts)
{
    {
        const int mpoints = mpoints1d*mpoints1d*mpoints1d;
        assert_eq(mpoints, wgt.getsize());
        assert_eq(mpoints, spts.getsize(1));
    }

    // Grab the 1D quadrature points
    dTensor1 w1d(mpoints1d),x1d(mpoints1d);
    setGaussPoints1d( w1d, x1d );

    // Tensor product Gaussian Quadrature
    int k=0;
    for (int m1=1; m1<=(mpoints1d); m1++)
        for (int m2=1; m2<=(mpoints1d); m2++)
            for (int m3=1; m3<=(mpoints1d); m3++)
            {
                k = k+1;
                wgt.set(k,  w1d.get(m1)*w1d.get(m2)*w1d.get(m3) );

                spts.set(k,1, x1d.get(m1) );
                spts.set(k,2, x1d.get(m2) );
                spts.set(k,3, x1d.get(m3) );
            }
}


// Loop over each quadrature point to construct Legendre polys
void SetLegendrePolys(const int mpoints, 
        const int kmax,
        const dTensor2& spts, 
        dTensor2& phi)
{
    assert_eq(mpoints,spts.getsize(1));
    assert_eq(mpoints,phi.getsize(1));
    assert_eq(kmax,phi.getsize(2));

    for (int m=1; m<=mpoints; m++)
    {
        // grid point (x,y,z)
        const double xi   = spts.get(m,1);
        const double eta  = spts.get(m,2);
        const double zeta = spts.get(m,3);
        const double xi2  = xi*xi;
        const double xi3  = xi*xi2;
        const double xi4  = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;
        const double zeta2 = zeta*zeta;
        const double zeta3 = zeta*zeta2;
        const double zeta4 = zeta*zeta3;

        // Legendre basis functions at each gaussian quadrature point in the
        // interval [-1,1]x[-1,1]x[-1,1].
        switch( kmax )
        {
            case 20:  // fourth order
                phi.set( m,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
                phi.set( m,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
                phi.set( m,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
                phi.set( m,17, 3.0*sq3*xi*eta*zeta                 );
                phi.set( m,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
                phi.set( m,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
                phi.set( m,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
                phi.set( m,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
                phi.set( m,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
                phi.set( m,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );

            case 10:  // third order
                phi.set( m,10, onehalf*sq5*(3.0*zeta2-1.0)         );
                phi.set( m,9,  onehalf*sq5*(3.0*eta2-1.0)          );
                phi.set( m,8,  onehalf*sq5*(3.0*xi2-1.0)           );
                phi.set( m,7,  3.0*eta*zeta                        );
                phi.set( m,6,  3.0*xi*zeta                         );
                phi.set( m,5,  3.0*xi*eta                          );

            case 4:  // second order  
                phi.set( m,4,  sq3*zeta                            );
                phi.set( m,3,  sq3*eta                             );
                phi.set( m,2,  sq3*xi                              );

            case 1:  // first order
                phi.set( m,1, 1.0                                  );

                break;

            default:
                unsupported_value_error(kmax);
        }

    }

}


// Loop over each quadrature point to construct Legendre polys
void SetLegendrePolysGrad(const double dx,
        const double dy,
        const double dz,
        const int mpoints, 
        const int kmax, 
        const dTensor2& spts, 
        dTensor2& phi_x,
        dTensor2& phi_y,
        dTensor2& phi_z)
{

    {
        assert_eq(mpoints, spts.getsize(1));
        assert_eq(mpoints, phi_x.getsize(1));
        assert_eq(kmax,    phi_x.getsize(2));
        assert_eq(mpoints, phi_y.getsize(1));
        assert_eq(kmax,    phi_y.getsize(2));
        assert_eq(mpoints, phi_z.getsize(1));
        assert_eq(kmax,    phi_z.getsize(2));
    }

    const double tmpx = 2.0/dx;
    const double tmpy = 2.0/dy;
    const double tmpz = 2.0/dz;

    for (int m=1; m<=(mpoints); m++)
    {
        // grid point (x, y, z)
        const double xi    = spts.get(m, 1);
        const double eta   = spts.get(m, 2);
        const double zeta  = spts.get(m, 3);
        const double xi2   = xi*xi;
        const double xi3   = xi*xi2;
        const double xi4   = xi*xi3;
        const double eta2  = eta*eta;
        const double eta3  = eta*eta2;
        const double eta4  = eta*eta3;
        const double zeta2 = zeta*zeta;
        const double zeta3 = zeta*zeta2;
        const double zeta4 = zeta*zeta3; 

        // Gradient of Legendre basis functions at each gaussian quadrature point
        switch( kmax )
        {
            case 20:  // fourth order
                phi_x.set( m,20, 0.0                                  );
                phi_x.set( m,19, 0.0                                  );
                phi_x.set( m,18, tmpx*3.0*onehalf*sq7*(5.0*xi2-1.0)   );
                phi_x.set( m,17, tmpx*3.0*sq3*eta*zeta                );
                phi_x.set( m,16, 0.0                                  );
                phi_x.set( m,15, tmpx*onehalf*sq3*sq5*(3.0*zeta2-1.0) );
                phi_x.set( m,14, 0.0                                  );
                phi_x.set( m,13, tmpx*onehalf*sq3*sq5*(3.0*eta2-1.0)  );
                phi_x.set( m,12, tmpx*3.0*sq3*sq5*zeta*xi             );
                phi_x.set( m,11, tmpx*3.0*sq3*sq5*eta*xi              );

                phi_y.set( m,20, 0.0                                  );
                phi_y.set( m,19, tmpy*3.0*onehalf*sq7*(5.0*eta2-1.0)  );
                phi_y.set( m,18, 0.0                                  );      
                phi_y.set( m,17, tmpy*3.0*sq3*xi*zeta                 );
                phi_y.set( m,16, tmpy*onehalf*sq3*sq5*(3.0*zeta2-1.0) );
                phi_y.set( m,15, 0.0                                  );
                phi_y.set( m,14, tmpy*3.0*sq3*sq5*eta*zeta            );
                phi_y.set( m,13, tmpy*3.0*sq3*sq5*eta*xi              );
                phi_y.set( m,12, 0.0                                  );          
                phi_y.set( m,11, tmpy*onehalf*sq3*sq5*(3.0*xi2-1.0)   );

                phi_z.set( m,20, tmpz*3.0*onehalf*sq7*(5.0*zeta2-1.0) );
                phi_z.set( m,19, 0.0                                  );
                phi_z.set( m,18, 0.0                                  );      
                phi_z.set( m,17, tmpz*3.0*sq3*xi*eta                  );
                phi_z.set( m,16, tmpz*3.0*sq3*sq5*eta*zeta            );
                phi_z.set( m,15, tmpz*3.0*sq3*sq5*zeta*xi             );
                phi_z.set( m,14, tmpz*onehalf*sq3*sq5*(3.0*eta2-1.0)  );
                phi_z.set( m,13, 0.0                                  );
                phi_z.set( m,12, tmpz*onehalf*sq3*sq5*(3.0*xi2-1.0)   );      
                phi_z.set( m,11, 0.0                                  );

            case 10:  // third order
                phi_x.set( m,10, 0.0                                  );
                phi_x.set( m,9,  0.0                                  );
                phi_x.set( m,8,  tmpx*3.0*sq5*xi                      );
                phi_x.set( m,7,  0.0                                  );
                phi_x.set( m,6,  tmpx*3.0*zeta                        );
                phi_x.set( m,5,  tmpx*3.0*eta                         );

                phi_y.set( m,10, 0.0                                  );
                phi_y.set( m,9,  tmpy*3.0*sq5*eta                     );
                phi_y.set( m,8,  0.0                                  );
                phi_y.set( m,7,  tmpy*3.0*zeta                        );
                phi_y.set( m,6,  0.0                                  );
                phi_y.set( m,5,  tmpy*3.0*xi                          );

                phi_z.set( m,10, tmpz*3.0*sq5*zeta                    );
                phi_z.set( m,9,  0.0                                  );
                phi_z.set( m,8,  0.0                                  );
                phi_z.set( m,7,  tmpz*3.0*eta                         );
                phi_z.set( m,6,  tmpz*3.0*xi                          );
                phi_z.set( m,5,  0.0                                  );

            case 4:  // second order  
                phi_x.set( m,4,  0.0                                  );
                phi_x.set( m,3,  0.0                                  );
                phi_x.set( m,2,  tmpx*sq3                             );

                phi_y.set( m,4,  0.0                                  );
                phi_y.set( m,3,  tmpy*sq3                             );
                phi_y.set( m,2,  0.0                                  );

                phi_z.set( m,4,  tmpz*sq3                             );
                phi_z.set( m,3,  0.0                                  );
                phi_z.set( m,2,  0.0                                  );

            case 1:  // first order
                phi_x.set( m,1,  0.0                                  );

                phi_y.set( m,1,  0.0                                  );

                phi_z.set( m,1,  0.0                                  );

                break;

            default:
                unsupported_value_error(kmax);
        }
    }
}
