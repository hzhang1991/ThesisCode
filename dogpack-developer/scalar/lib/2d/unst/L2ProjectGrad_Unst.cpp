#include "dog_math.h"
#include "dogdefs.h"
#include "mesh.h"
#include "MonomialsToLegendre.h"

// All-purpose routine for computing the L2-projection
// of various functions onto the gradient of the Legendre basis
//     (Unstructured grid version)
//
void L2ProjectGrad_Unst(
    const dTensor2* vel_vec,
    const int istart, 
    const int iend, 
    const int QuadOrder, 
    const int BasisOrder_qin,
    const int BasisOrder_auxin,
    const int BasisOrder_fout,
    const mesh& Mesh, 
    const dTensor3* qin, 
    const dTensor3* auxin, 
    dTensor3* fout, 
    void (*Func)(const dTensor2* vel_vec,
        const dTensor2&,const dTensor2&,
        const dTensor2&,dTensor3&))
{
    // starting and ending indeces
    const int   NumElems = Mesh.get_NumElems();
    assert_ge(istart,1);
    assert_le(iend,NumElems);

    // qin variable
    assert_eq(NumElems,qin->getsize(1));
    const int     meqn = qin->getsize(2);
    const int kmax_qin = qin->getsize(3);
    assert_eq(kmax_qin,(BasisOrder_qin*(BasisOrder_qin+1))/2);

    // auxin variable
    assert_eq(NumElems,auxin->getsize(1));
    const int       maux = auxin->getsize(2);
    const int kmax_auxin = auxin->getsize(3);
    assert_eq(kmax_auxin,(BasisOrder_auxin*(BasisOrder_auxin+1))/2);

    // fout variables
    assert_eq(NumElems,fout->getsize(1));
    const int mcomps_out = fout->getsize(2);
    const int  kmax_fout = fout->getsize(3);
    assert_eq(kmax_fout,(BasisOrder_fout*(BasisOrder_fout+1))/2);

    // number of quadrature points
    assert_ge(QuadOrder,1);
    assert_le(QuadOrder,5);
    int mpoints;
    switch ( QuadOrder )
    {
        case 1:
            mpoints = 0;
            break;

        case 2:
            mpoints = 1;
            break;

        case 3:
            mpoints = 6;
            break;

        case 4:
            mpoints = 7;
            break;

        case 5:	     
            mpoints = 16;
            break;
    }

    // trivial case
    if ( QuadOrder==1 )
    {
        for (int i=istart; i<=iend; i++)
        for (int m=1; m<=mcomps_out; m++) 
        for (int k=1; k<=kmax_fout; k++) 
        {  fout->set(i,m,k, 0.0 );  }
    }
    else
    {
        const int kmax = iMax(iMax(kmax_qin,kmax_auxin),kmax_fout);
        dTensor2    spts(mpoints,2);
        dTensor1    wgts(mpoints);
        dTensor2    xpts(mpoints,2);
        dTensor2   qvals(mpoints,meqn);
        dTensor2 auxvals(mpoints,maux);
        dTensor3   fvals(mpoints,mcomps_out,2);
        dTensor2      mu(mpoints,kmax); // monomial basis (non-orthogonal)
        dTensor2     phi(mpoints,kmax); // Legendre basis (orthogonal)
        dTensor2   mu_xi(mpoints,kmax_fout);   //  xi-derivative of monomial basis (non-orthogonal)
        dTensor2  mu_eta(mpoints,kmax_fout);   // eta-derivative of monomial basis (non-orthogonal)
        dTensor2  phi_xi(mpoints,kmax_fout);   //  xi-derivative of Legendre basis (orthogonal)
        dTensor2 phi_eta(mpoints,kmax_fout);   // eta-derivative of Legendre basis (orthogonal)
        dTensor2   phi_x(mpoints,kmax_fout);   //   x-derivative of Legendre basis (orthogonal)
        dTensor2   phi_y(mpoints,kmax_fout);   //   y-derivative of Legendre basis (orthogonal)

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

            // monomial basis functions at each gaussian quadrature point
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

            // Loop over each quadrature point and construct Legendre polys
            for (int i=1; i<=kmax; i++)
            {
                double tmp = 0.0;
                for (int j=1; j<=i; j++)
                {  tmp = tmp + Mmat[i-1][j-1]*mu.get(m,j);  }

                phi.set(m,i, tmp );
            }	

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

        // -------------------------------------------------------------
        // Loop over every grid cell indexed by user supplied parameters
        // described by istart...iend
        // -------------------------------------------------------------
#pragma omp parallel for
        for (int i=istart; i<=iend; i++)
        {	  
            // Find center of current cell
            const int i1 = Mesh.get_tnode(i,1);
            const int i2 = Mesh.get_tnode(i,2);
            const int i3 = Mesh.get_tnode(i,3);
            const double x1 = Mesh.get_node(i1,1);
            const double y1 = Mesh.get_node(i1,2);
            const double x2 = Mesh.get_node(i2,1);
            const double y2 = Mesh.get_node(i2,2);
            const double x3 = Mesh.get_node(i3,1);
            const double y3 = Mesh.get_node(i3,2);

            const double xc = (x1+x2+x3)/3.0;
            const double yc = (y1+y2+y3)/3.0;

            // Compute q, aux and fvals at each Gaussian Quadrature point
            // for this current cell indexed by (i,j)
            // Save results into dTensor2 qvals, auxvals and fvals.
            for (int m=1; m<=mpoints; m++)
            {
                // convert phi_xi and phi_eta derivatives
                // to phi_x and phi_y derivatives through Jacobian
                for (int k=1; k<=kmax_fout; k++)
                {
                    phi_x.set(m,k, Mesh.get_jmat(i,1,1)*phi_xi.get(m,k)
                            + Mesh.get_jmat(i,1,2)*phi_eta.get(m,k) );
                    phi_y.set(m,k, Mesh.get_jmat(i,2,1)*phi_xi.get(m,k)
                            + Mesh.get_jmat(i,2,2)*phi_eta.get(m,k) );
                }

                // point on the unit triangle
                const double s = spts.get(m,1);
                const double t = spts.get(m,2);

                // point on the physical triangle
                xpts.set(m,1, xc + (x2-x1)*s + (x3-x1)*t );
                xpts.set(m,2, yc + (y2-y1)*s + (y3-y1)*t );

                // Solution values (q) at each grid point
                for (int me=1; me<=meqn; me++)
                {
                    qvals.set(m,me, 0.0 );

                    for (int k=1; k<=kmax_qin; k++)
                    {
                        qvals.set(m,me, qvals.get(m,me) 
                                + phi.get(m,k) * qin->get(i,me,k) );
                    }
                }

                // Auxiliary values (aux) at each grid point
                for (int ma=1; ma<=maux; ma++)
                {
                    auxvals.set(m,ma, 0.0 );

                    for (int k=1; k<=kmax_auxin; k++)
                    {
                        auxvals.set(m,ma, auxvals.get(m,ma) 
                                + phi.get(m,k) * auxin->get(i,ma,k) );
                    }
                } 
            }

            // Call user-supplied function to set fvals
            Func(vel_vec, xpts, qvals, auxvals, fvals);

            // Evaluate integral on current cell (project onto Legendre basis) 
            // using Gaussian Quadrature for the integration
            for (int m1=1; m1<=mcomps_out; m1++)		
            for (int m2=1; m2<=kmax_fout; m2++)
            {
                double tmp = 0.0;
                for (int k=1; k<=mpoints; k++)
                {
                    tmp = tmp + wgts.get(k)*
                        ( fvals.get(k,m1,1)*phi_x.get(k,m2) +
                          fvals.get(k,m1,2)*phi_y.get(k,m2) );
                }
                fout->set(i, m1, m2,  2.0*tmp );
            }

        }
    }
}
