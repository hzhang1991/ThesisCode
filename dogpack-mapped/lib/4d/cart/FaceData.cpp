#include "FaceData.h"

FaceData::FaceData()
// Constructor
{
  const int mpoints1d = dogParams.get_space_order();
  const int kmax = dogParams.get_kmax();

  KMAX_MAX = 20;
  
  // 2D Gaussian quadrature weights and points
  wgts2d = new dTensor1(mpoints1d*mpoints1d); wgts2d->setall(0.);
  xpts2d = new dTensor2(mpoints1d*mpoints1d,2); xpts2d->setall(0.);
  
  // Legendre basis functions
  phi_xl = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); phi_xl->setall(0.);
  phi_xr = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); phi_xr->setall(0.);
  phi_yl = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); phi_yl->setall(0.);
  phi_yr = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); phi_yr->setall(0.);
  phi_zl = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); phi_zl->setall(0.);
  phi_zr = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); phi_zr->setall(0.);
  
  // weights times Legendre basis functions
  wght_phi_xl = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); wght_phi_xl->setall(0.);
  wght_phi_xr = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); wght_phi_xr->setall(0.);
  wght_phi_yl = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); wght_phi_yl->setall(0.);
  wght_phi_yr = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); wght_phi_yr->setall(0.);
  wght_phi_zl = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); wght_phi_zl->setall(0.);
  wght_phi_zr = new dTensor2(mpoints1d*mpoints1d,KMAX_MAX); wght_phi_zr->setall(0.);
}

FaceData::~FaceData()
// Destructor
{
  // Gaussian quadrature points and weights
  delete wgts2d;
  delete xpts2d;
  
  // Legendre basis functions
  delete phi_xl;
  delete phi_xr;
  delete phi_yl;
  delete phi_yr;
  delete phi_zl;
  delete phi_zr;
}

// Initial FaceData
void FaceData::init()
{
  FaceData& FaceData = *this;
  const int morder = dogParams.get_space_order();
  const int kmax = dogParams.get_kmax();
  
  double xi,xi2,xi3,xi4;
  double eta,eta2,eta3,eta4;
  double zeta,zeta2,zeta3,zeta4;
  const double dx = dogParamsCart4.get_dx();
  const double dy = dogParamsCart4.get_dy();
  const double dz = dogParamsCart4.get_dz();
  
  // ---------------------------------
  // Quick error check
  // ---------------------------------
  assert_printf(kmax==(morder*(morder+1)*(morder+2))/6, "\n"
                "       morder = %d\n"
                "         kmax = %d\n",
                morder,kmax);
  
  // ---------------------------------
  // Set quadrature weights and points
  // ---------------------------------
  dTensor1 wgts1d(morder);
  dTensor1 xpts1d(morder);

  switch ( morder )
    {
    case 1:
      wgts1d.set(1, 2.0 );
      
      xpts1d.set(1, 0.0 );
      break;
      
    case 2:
      wgts1d.set(1,  1.0 );
      wgts1d.set(2,  1.0 );
      
      xpts1d.set(1, -1.0/sq3 );
      xpts1d.set(2,  1.0/sq3 );
      break;

    case 3:
      wgts1d.set(1,  5.0/9.0 );
      wgts1d.set(2,  8.0/9.0 );
      wgts1d.set(3,  5.0/9.0 );
      
      xpts1d.set(1,  sq3/sq5 );
      xpts1d.set(2,  0.0 );
      xpts1d.set(3, -sq3/sq5 );
      break;
      
    case 4:
      wgts1d.set(1, (18.0 - sq3*sq10)/36.0 );
      wgts1d.set(2, (18.0 + sq3*sq10)/36.0 );
      wgts1d.set(3, wgts1d.get(2) );
      wgts1d.set(4, wgts1d.get(1) );
      
      xpts1d.set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
      xpts1d.set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
      xpts1d.set(3, -xpts1d.get(2) );
      xpts1d.set(4, -xpts1d.get(1) );           
      break;
      
    case 5:      
      wgts1d.set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
      wgts1d.set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
      wgts1d.set(3, 128.0/225.0 );
      wgts1d.set(4, wgts1d.get(2) );
      wgts1d.set(5, wgts1d.get(1) );
      
      xpts1d.set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
      xpts1d.set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
      xpts1d.set(3,  0.0 );
      xpts1d.set(4, -xpts1d.get(2) );
      xpts1d.set(5, -xpts1d.get(1) );
      break;
    }

  assert_eq(wgts2d->getsize(),morder*morder);
  assert_eq(xpts2d->getsize(1),morder*morder);
  assert_eq(xpts2d->getsize(2),2);

  int k=0;
  for (int i=1; i<=morder; i++)
    for (int j=1; j<=morder; j++)
      {
        k=k+1;
        wgts2d->set(k, wgts1d.get(i)*wgts1d.get(j) );
        xpts2d->set(k,1, xpts1d.get(i) );
        xpts2d->set(k,2, xpts1d.get(j) );
      }
  
  // ----------------------------------------------
  // Evaluate Legendre basis functions on the faces
  // ----------------------------------------------
  for (int m=1; m<=(morder*morder); m++)
    {
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Right face (will be Left state in Riemann problem)
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xi   = 1.0;
      eta  = xpts2d->get(m,1);
      zeta = xpts2d->get(m,2);
      
      xi2 = xi*xi;
      xi3 = xi2*xi;
      xi4 = xi3*xi;

      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;

      zeta2 = zeta*zeta;
      zeta3 = zeta2*zeta;
      zeta4 = zeta3*zeta;
      
      phi_xl->set( m,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
      phi_xl->set( m,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
      phi_xl->set( m,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
      phi_xl->set( m,17, 3.0*sq3*xi*eta*zeta                 );
      phi_xl->set( m,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
      phi_xl->set( m,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
      phi_xl->set( m,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
      phi_xl->set( m,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
      phi_xl->set( m,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
      phi_xl->set( m,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );
      phi_xl->set( m,10, onehalf*sq5*(3.0*zeta2-1.0)         );
      phi_xl->set( m,9,  onehalf*sq5*(3.0*eta2-1.0)          );
      phi_xl->set( m,8,  onehalf*sq5*(3.0*xi2-1.0)           );
      phi_xl->set( m,7,  3.0*eta*zeta                        );
      phi_xl->set( m,6,  3.0*xi*zeta                         );
      phi_xl->set( m,5,  3.0*xi*eta                          );      
      phi_xl->set( m,4,  sq3*zeta                            );
      phi_xl->set( m,3,  sq3*eta                             );
      phi_xl->set( m,2,  sq3*xi                              );
      phi_xl->set( m,1,  1.0                                 );
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Left face (will be Right state in Riemann problem)
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xi   = -1.0;
      eta  =  xpts2d->get(m,1);
      zeta =  xpts2d->get(m,2);
      
      xi2 = xi*xi;
      xi3 = xi2*xi;
      xi4 = xi3*xi;

      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;

      zeta2 = zeta*zeta;
      zeta3 = zeta2*zeta;
      zeta4 = zeta3*zeta;

      phi_xr->set( m,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
      phi_xr->set( m,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
      phi_xr->set( m,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
      phi_xr->set( m,17, 3.0*sq3*xi*eta*zeta                 );
      phi_xr->set( m,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
      phi_xr->set( m,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
      phi_xr->set( m,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
      phi_xr->set( m,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
      phi_xr->set( m,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
      phi_xr->set( m,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );
      phi_xr->set( m,10, onehalf*sq5*(3.0*zeta2-1.0)         );
      phi_xr->set( m,9,  onehalf*sq5*(3.0*eta2-1.0)          );
      phi_xr->set( m,8,  onehalf*sq5*(3.0*xi2-1.0)           );
      phi_xr->set( m,7,  3.0*eta*zeta                        );
      phi_xr->set( m,6,  3.0*xi*zeta                         );
      phi_xr->set( m,5,  3.0*xi*eta                          );
      phi_xr->set( m,4,  sq3*zeta                            );
      phi_xr->set( m,3,  sq3*eta                             );
      phi_xr->set( m,2,  sq3*xi                              );
      phi_xr->set( m,1,  1.0                                 );

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Back face (will be Left state in Riemann problem)
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xi   =  xpts2d->get(m,1);
      eta  =  1.0;
      zeta =  xpts2d->get(m,2);
      
      xi2 = xi*xi;
      xi3 = xi2*xi;
      xi4 = xi3*xi;

      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;

      zeta2 = zeta*zeta;
      zeta3 = zeta2*zeta;
      zeta4 = zeta3*zeta;
      
      phi_yl->set( m,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
      phi_yl->set( m,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
      phi_yl->set( m,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
      phi_yl->set( m,17, 3.0*sq3*xi*eta*zeta                 );
      phi_yl->set( m,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
      phi_yl->set( m,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
      phi_yl->set( m,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
      phi_yl->set( m,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
      phi_yl->set( m,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
      phi_yl->set( m,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );
      phi_yl->set( m,10, onehalf*sq5*(3.0*zeta2-1.0)         );
      phi_yl->set( m,9,  onehalf*sq5*(3.0*eta2-1.0)          );
      phi_yl->set( m,8,  onehalf*sq5*(3.0*xi2-1.0)           );
      phi_yl->set( m,7,  3.0*eta*zeta                        );
      phi_yl->set( m,6,  3.0*xi*zeta                         );
      phi_yl->set( m,5,  3.0*xi*eta                          );
      phi_yl->set( m,4,  sq3*zeta                            );
      phi_yl->set( m,3,  sq3*eta                             );
      phi_yl->set( m,2,  sq3*xi                              );
      phi_yl->set( m,1,  1.0                                 );

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Front face (will be Right state in Riemann problem)
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xi   =  xpts2d->get(m,1);
      eta  = -1.0;
      zeta =  xpts2d->get(m,2);
      
      xi2 = xi*xi;
      xi3 = xi2*xi;
      xi4 = xi3*xi;

      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;

      zeta2 = zeta*zeta;
      zeta3 = zeta2*zeta;
      zeta4 = zeta3*zeta;
      
      phi_yr->set( m,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
      phi_yr->set( m,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
      phi_yr->set( m,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
      phi_yr->set( m,17, 3.0*sq3*xi*eta*zeta                 );
      phi_yr->set( m,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
      phi_yr->set( m,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
      phi_yr->set( m,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
      phi_yr->set( m,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
      phi_yr->set( m,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
      phi_yr->set( m,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );
      phi_yr->set( m,10, onehalf*sq5*(3.0*zeta2-1.0)         );
      phi_yr->set( m,9,  onehalf*sq5*(3.0*eta2-1.0)          );
      phi_yr->set( m,8,  onehalf*sq5*(3.0*xi2-1.0)           );
      phi_yr->set( m,7,  3.0*eta*zeta                        );
      phi_yr->set( m,6,  3.0*xi*zeta                         );
      phi_yr->set( m,5,  3.0*xi*eta                          );
      phi_yr->set( m,4,  sq3*zeta                            );
      phi_yr->set( m,3,  sq3*eta                             );
      phi_yr->set( m,2,  sq3*xi                              );
      phi_yr->set( m,1,  1.0                                 );

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Top face (will be Left state in Riemann problem)
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      xi   =  xpts2d->get(m,1);
      eta  =  xpts2d->get(m,2);
      zeta =  1.0;
      
      xi2 = xi*xi;
      xi3 = xi2*xi;
      xi4 = xi3*xi;

      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;

      zeta2 = zeta*zeta;
      zeta3 = zeta2*zeta;
      zeta4 = zeta3*zeta;
      
      phi_zl->set( m,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
      phi_zl->set( m,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
      phi_zl->set( m,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
      phi_zl->set( m,17, 3.0*sq3*xi*eta*zeta                 );
      phi_zl->set( m,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
      phi_zl->set( m,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
      phi_zl->set( m,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
      phi_zl->set( m,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
      phi_zl->set( m,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
      phi_zl->set( m,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );
      phi_zl->set( m,10, onehalf*sq5*(3.0*zeta2-1.0)         );
      phi_zl->set( m,9,  onehalf*sq5*(3.0*eta2-1.0)          );
      phi_zl->set( m,8,  onehalf*sq5*(3.0*xi2-1.0)           );
      phi_zl->set( m,7,  3.0*eta*zeta                        );
      phi_zl->set( m,6,  3.0*xi*zeta                         );
      phi_zl->set( m,5,  3.0*xi*eta                          );
      phi_zl->set( m,4,  sq3*zeta                            );
      phi_zl->set( m,3,  sq3*eta                             );
      phi_zl->set( m,2,  sq3*xi                              );
      phi_zl->set( m,1,  1.0                                 );
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Bottom face (will be Right state in Riemann problem)
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xi   =  xpts2d->get(m,1);
      eta  =  xpts2d->get(m,2);
      zeta = -1.0;
      
      xi2 = xi*xi;
      xi3 = xi2*xi;
      xi4 = xi3*xi;

      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;

      zeta2 = zeta*zeta;
      zeta3 = zeta2*zeta;
      zeta4 = zeta3*zeta;
      
      phi_zr->set( m,20, onehalf*sq7*zeta*(5.0*zeta2-3.0)    );
      phi_zr->set( m,19, onehalf*sq7*eta*(5.0*eta2-3.0)      );
      phi_zr->set( m,18, onehalf*sq7*xi*(5.0*xi2-3.0)        );
      phi_zr->set( m,17, 3.0*sq3*xi*eta*zeta                 );
      phi_zr->set( m,16, onehalf*sq3*sq5*eta*(3.0*zeta2-1.0) );
      phi_zr->set( m,15, onehalf*sq3*sq5*xi*(3.0*zeta2-1.0)  );
      phi_zr->set( m,14, onehalf*sq3*sq5*zeta*(3.0*eta2-1.0) );
      phi_zr->set( m,13, onehalf*sq3*sq5*xi*(3.0*eta2-1.0)   );
      phi_zr->set( m,12, onehalf*sq3*sq5*zeta*(3.0*xi2-1.0)  );
      phi_zr->set( m,11, onehalf*sq3*sq5*eta*(3.0*xi2-1.0)   );
      phi_zr->set( m,10, onehalf*sq5*(3.0*zeta2-1.0)         );
      phi_zr->set( m,9,  onehalf*sq5*(3.0*eta2-1.0)          );
      phi_zr->set( m,8,  onehalf*sq5*(3.0*xi2-1.0)           );
      phi_zr->set( m,7,  3.0*eta*zeta                        );
      phi_zr->set( m,6,  3.0*xi*zeta                         );
      phi_zr->set( m,5,  3.0*xi*eta                          );
      phi_zr->set( m,4,  sq3*zeta                            );
      phi_zr->set( m,3,  sq3*eta                             );
      phi_zr->set( m,2,  sq3*xi                              );
      phi_zr->set( m,1,  1.0                                 );
    }

  for (int m=1; m<=(morder*morder); m++)
    for (int k=1; k<=KMAX_MAX; k++)
      {
	wght_phi_xl->set(m,k, wgts2d->get(m)*phi_xl->get(m,k));
	wght_phi_xr->set(m,k, wgts2d->get(m)*phi_xr->get(m,k));
	wght_phi_yl->set(m,k, wgts2d->get(m)*phi_yl->get(m,k));
	wght_phi_yr->set(m,k, wgts2d->get(m)*phi_yr->get(m,k));
	wght_phi_zl->set(m,k, wgts2d->get(m)*phi_zl->get(m,k));
	wght_phi_zr->set(m,k, wgts2d->get(m)*phi_zr->get(m,k));
      }
  
}
