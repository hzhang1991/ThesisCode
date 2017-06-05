#include "Legendre4d.h"

Legendre4d::Legendre4d():
    faceData(0)
{
    faceData = new FaceData;
    faceData->init();
};

Legendre4d::~Legendre4d()
{
}

Legendre4d& Legendre4d::instance()
{
    static Legendre4d* legendre4d = new Legendre4d;
    return *legendre4d;
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

    const int md2 = mpoints1d/2;
    dTensor1 w1d(mpoints1d),x1d(mpoints1d);

    switch (mpoints1d)
    {
        case 1:
            w1d.set(1, 2.0 );

            x1d.set(1, 0.0 );
            break;

        case 2:
            w1d.set(1,  1.0 );
            w1d.set(2,  1.0 );

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );
            break;

        case 3:
            w1d.set(1,  5.0/9.0 );
            w1d.set(2,  8.0/9.0 );
            w1d.set(3,  5.0/9.0 );

            x1d.set(1,  sq3/sq5 );
            x1d.set(2,  0.0 );
            x1d.set(3, -sq3/sq5 );
            break;

        case 4:
            w1d.set(1, (18.0 - sq3*sq10)/36.0 );
            w1d.set(2, (18.0 + sq3*sq10)/36.0 );
            w1d.set(3, w1d.get(2) );
            w1d.set(4, w1d.get(1) );

            x1d.set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            x1d.set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            x1d.set(3, -x1d.get(2) );
            x1d.set(4, -x1d.get(1) );       
            break;

        case 5:      
            w1d.set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            w1d.set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            w1d.set(3, 128.0/225.0 );
            w1d.set(4, w1d.get(2) );
            w1d.set(5, w1d.get(1) );

            x1d.set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            x1d.set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );
            break;

        case 6:
            w1d.set(1, 0.1713244923791703450402961 );
            w1d.set(2, 0.3607615730481386075698335 );
            w1d.set(3, 0.4679139345726910473898703 );      

            x1d.set(1, 0.9324695142031520278123016 );
            x1d.set(2, 0.6612093864662645136613996 );
            x1d.set(3, 0.2386191860831969086305017 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 8:
            w1d.set(1, 0.1012285362903762591525314 );
            w1d.set(2, 0.2223810344533744705443560 );
            w1d.set(3, 0.3137066458778872873379622 );
            w1d.set(4, 0.3626837833783619829651504 );

            x1d.set(1, 0.9602898564975362316835609 );
            x1d.set(2, 0.7966664774136267395915539 );
            x1d.set(3, 0.5255324099163289858177390 );
            x1d.set(4, 0.1834346424956498049394761 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 10:
            w1d.set(1, 0.0666713443086881375935688 );
            w1d.set(2, 0.1494513491505805931457763 );
            w1d.set(3, 0.2190863625159820439955349 );
            w1d.set(4, 0.2692667193099963550912269 );
            w1d.set(5, 0.2955242247147528701738930 );

            x1d.set(1, 0.9739065285171717200779640 );
            x1d.set(2, 0.8650633666889845107320967 );
            x1d.set(3, 0.6794095682990244062343274 );
            x1d.set(4, 0.4333953941292471907992659 );
            x1d.set(5, 0.1488743389816312108848260 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 12:
            w1d.set(1, 0.0471753363865118271946160 );
            w1d.set(2, 0.1069393259953184309602547 );
            w1d.set(3, 0.1600783285433462263346525 );
            w1d.set(4, 0.2031674267230659217490645 );
            w1d.set(5, 0.2334925365383548087608499 );
            w1d.set(6, 0.2491470458134027850005624 );

            x1d.set(1, 0.9815606342467192506905491 );
            x1d.set(2, 0.9041172563704748566784659 );
            x1d.set(3, 0.7699026741943046870368938 );
            x1d.set(4, 0.5873179542866174472967024 );
            x1d.set(5, 0.3678314989981801937526915 );
            x1d.set(6, 0.1252334085114689154724414 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 14:
            w1d.set(1, 0.0351194603317518630318329 );
            w1d.set(2, 0.0801580871597602098056333 );
            w1d.set(3, 0.1215185706879031846894148 );
            w1d.set(4, 0.1572031671581935345696019 );
            w1d.set(5, 0.1855383974779378137417166 );
            w1d.set(6, 0.2051984637212956039659241 );
            w1d.set(7, 0.2152638534631577901958764 );

            x1d.set(1, 0.9862838086968123388415973 );
            x1d.set(2, 0.9284348836635735173363911 );
            x1d.set(3, 0.8272013150697649931897947 );
            x1d.set(4, 0.6872929048116854701480198 );
            x1d.set(5, 0.5152486363581540919652907 );
            x1d.set(6, 0.3191123689278897604356718 );
            x1d.set(7, 0.1080549487073436620662447 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 16:
            w1d.set(1, 0.0271524594117540948517806 );
            w1d.set(2, 0.0622535239386478928628438 );
            w1d.set(3, 0.0951585116824927848099251 );
            w1d.set(4, 0.1246289712555338720524763 );
            w1d.set(5, 0.1495959888165767320815017 );
            w1d.set(6, 0.1691565193950025381893121 );
            w1d.set(7, 0.1826034150449235888667637 );
            w1d.set(8, 0.1894506104550684962853967 );

            x1d.set(1, 0.9894009349916499325961542 );
            x1d.set(2, 0.9445750230732325760779884 );
            x1d.set(3, 0.8656312023878317438804679 );
            x1d.set(4, 0.7554044083550030338951012 );
            x1d.set(5, 0.6178762444026437484466718 );
            x1d.set(6, 0.4580167776572273863424194 );
            x1d.set(7, 0.2816035507792589132304605 );
            x1d.set(8, 0.0950125098376374401853193 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 18:
            w1d.set(1, 0.0216160135264833103133427 );
            w1d.set(2, 0.0497145488949697964533349 );
            w1d.set(3, 0.0764257302548890565291297 );
            w1d.set(4, 0.1009420441062871655628140 );
            w1d.set(5, 0.1225552067114784601845191 );
            w1d.set(6, 0.1406429146706506512047313 );
            w1d.set(7, 0.1546846751262652449254180 );
            w1d.set(8, 0.1642764837458327229860538 );
            w1d.set(9, 0.1691423829631435918406565 );

            x1d.set(1, 0.9915651684209309467300160 );
            x1d.set(2, 0.9558239495713977551811959 );
            x1d.set(3, 0.8926024664975557392060606 );
            x1d.set(4, 0.8037049589725231156824175 );
            x1d.set(5, 0.6916870430603532078748911 );
            x1d.set(6, 0.5597708310739475346078715 );
            x1d.set(7, 0.4117511614628426460359318 );
            x1d.set(8, 0.2518862256915055095889729 );
            x1d.set(9, 0.0847750130417353012422619 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 20:
            w1d.set(1,  0.0176140071391521183118620 );
            w1d.set(2,  0.0406014298003869413310400 );
            w1d.set(3,  0.0626720483341090635695065 );
            w1d.set(4,  0.0832767415767047487247581 );
            w1d.set(5,  0.1019301198172404350367501 );
            w1d.set(6,  0.1181945319615184173123774 );
            w1d.set(7,  0.1316886384491766268984945 );
            w1d.set(8,  0.1420961093183820513292983 );
            w1d.set(9,  0.1491729864726037467878287 );
            w1d.set(10, 0.1527533871307258506980843 );

            x1d.set(1,  0.9931285991850949247861224 );
            x1d.set(2,  0.9639719272779137912676661 );
            x1d.set(3,  0.9122344282513259058677524 );
            x1d.set(4,  0.8391169718222188233945291 );
            x1d.set(5,  0.7463319064601507926143051 );
            x1d.set(6,  0.6360536807265150254528367 );
            x1d.set(7,  0.5108670019508270980043641 );
            x1d.set(8,  0.3737060887154195606725482 );
            x1d.set(9,  0.2277858511416450780804962 );
            x1d.set(10, 0.0765265211334973337546404 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        default:
            unsupported_value_error(mpoints1d);
    }

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
