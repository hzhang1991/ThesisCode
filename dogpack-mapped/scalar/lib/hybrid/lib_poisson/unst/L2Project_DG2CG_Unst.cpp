#include "dog_math.h"
#include "dogdefs.h"
#include "mesh.h"
#include "MonomialsToLegendre.h"

// -------------------------------------------------------------
// Routine for computing the L2-projection of a DG function
// onto onto a CG basis
// -------------------------------------------------------------
void L2Project_DG2CG_Unst(const int istart, 
        const int iend, 
        const int QuadOrder,            
        const int cg_order,
        const int dg_comp,
        const mesh& Mesh,
        const dTensor3& fdg,
        dTensor1& fcg)
{
    const int     NumElems = Mesh.get_NumElems();
    const int NumPhysElems = Mesh.get_NumPhysElems();
    const int         kmax = fdg.getsize(3);

    // initialize fcg
    fcg.setall(0.);
//  for (int i=1; i<=fcg.getsize(); i++)
//  {  fcg.set(i, 0.0 );  }

    // number of quadrature points
    assert_ge(QuadOrder,1);
    assert_le(QuadOrder,5);
    int mpoints;
    switch ( QuadOrder )
    {
        case 1:
            mpoints = 1;
            break;

        case 2:
            mpoints = 3;
            break;

        case 3:
            mpoints = 6;
            break;

        case 4:
            mpoints = 12;
            break;

        case 5:         
            mpoints = 16;
            break;
    }
    int iMax(int,int);
    dTensor2      mu(mpoints,kmax); // monomial basis (non-orthogonal)
    dTensor2     phi(mpoints,kmax); // Legendre basis (orthogonal)
    dTensor2    spts(mpoints,2);
    dTensor1    wgts(mpoints);
    dTensor2    xpts(mpoints,2);

    // ---------------------------------
    // Set quadrature weights and points
    // ---------------------------------  
    switch ( QuadOrder )
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

    // Inverse matrix to compute CG basis functions on each element
    const int cg_num_pts = ((cg_order+1)*(cg_order+2))/2;
    dTensor2 MatInv(cg_num_pts,cg_num_pts);
    switch(cg_order)
    {
        case 1:
            MatInv.set(1,1, onethird );
            MatInv.set(1,2, onethird );
            MatInv.set(1,3, onethird );

            MatInv.set(2,1, -1.0 );
            MatInv.set(2,2,  1.0 );
            MatInv.set(2,3,  0.0 );

            MatInv.set(3,1, -1.0 );
            MatInv.set(3,2,  0.0 );
            MatInv.set(3,3,  1.0 );
            break;

        case 2:
            MatInv.set(1,1, -oneninth     );
            MatInv.set(1,2,  4.0*oneninth );
            MatInv.set(1,3, -oneninth     );
            MatInv.set(1,4,  4.0*oneninth );
            MatInv.set(1,5,  4.0*oneninth );
            MatInv.set(1,6, -oneninth     );

            MatInv.set(2,1, -onethird     );
            MatInv.set(2,2,  0.0          );
            MatInv.set(2,3,  onethird     );
            MatInv.set(2,4, -4.0*onethird );
            MatInv.set(2,5,  4.0*onethird );
            MatInv.set(2,6,  0.0          );

            MatInv.set(3,1, -onethird     );
            MatInv.set(3,2, -4.0*onethird );
            MatInv.set(3,3,  0.0          );
            MatInv.set(3,4,  0.0          );
            MatInv.set(3,5,  4.0*onethird );
            MatInv.set(3,6,  onethird     );

            MatInv.set(4,1,  4.0          );
            MatInv.set(4,2, -4.0          );
            MatInv.set(4,3,  0.0          );
            MatInv.set(4,4, -4.0          );
            MatInv.set(4,5,  4.0          );
            MatInv.set(4,6,  0.0          );

            MatInv.set(5,1,  2.0          );
            MatInv.set(5,2, -4.0          );
            MatInv.set(5,3,  2.0          );
            MatInv.set(5,4,  0.0          );
            MatInv.set(5,5,  0.0          );
            MatInv.set(5,6,  0.0          );

            MatInv.set(6,1,  2.0          );
            MatInv.set(6,2,  0.0          );
            MatInv.set(6,3,  0.0          );
            MatInv.set(6,4, -4.0          );
            MatInv.set(6,5,  0.0          );
            MatInv.set(6,6,  2.0          );
            break;

        case 3:

            MatInv.set(1,1,   0.0          );
            MatInv.set(1,2,   0.0          );
            MatInv.set(1,3,   0.0          );
            MatInv.set(1,4,   0.0          );
            MatInv.set(1,5,   0.0          );
            MatInv.set(1,6,   1.0          );
            MatInv.set(1,7,   0.0          );
            MatInv.set(1,8,   0.0          );
            MatInv.set(1,9,   0.0          );
            MatInv.set(1,10,  0.0          );

            MatInv.set(2,1,   0.5          );
            MatInv.set(2,2,  -1.5          );
            MatInv.set(2,3,   1.5          );
            MatInv.set(2,4,  -0.5          );
            MatInv.set(2,5,  -1.5          );
            MatInv.set(2,6,   0.0          );
            MatInv.set(2,7,   1.5          );
            MatInv.set(2,8,   0.0          );
            MatInv.set(2,9,   0.0          );
            MatInv.set(2,10,  0.0          );

            MatInv.set(3,1,   0.5          );
            MatInv.set(3,2,  -1.5          );
            MatInv.set(3,3,   0.0          );
            MatInv.set(3,4,   0.0          );
            MatInv.set(3,5,  -1.5          );
            MatInv.set(3,6,   0.0          );
            MatInv.set(3,7,   0.0          );
            MatInv.set(3,8,   1.5          );
            MatInv.set(3,9,   1.5          );
            MatInv.set(3,10, -0.5          );

            MatInv.set(4,1,   0.0          );
            MatInv.set(4,2,   4.5          );
            MatInv.set(4,3,  -4.5          );
            MatInv.set(4,4,   0.0          );
            MatInv.set(4,5,   4.5          );
            MatInv.set(4,6,  -9.0          );
            MatInv.set(4,7,   4.5          );
            MatInv.set(4,8,  -4.5          );
            MatInv.set(4,9,   4.5          );
            MatInv.set(4,10,  0.0          );

            MatInv.set(5,1,   0.0          );
            MatInv.set(5,2,   0.0          );
            MatInv.set(5,3,   0.0          );
            MatInv.set(5,4,   0.0          );
            MatInv.set(5,5,   4.5          );
            MatInv.set(5,6,  -9.0          );
            MatInv.set(5,7,   4.5          );
            MatInv.set(5,8,   0.0          );
            MatInv.set(5,9,   0.0          );
            MatInv.set(5,10,  0.0          );

            MatInv.set(6,1,   0.0          );
            MatInv.set(6,2,   4.5          );
            MatInv.set(6,3,   0.0          );
            MatInv.set(6,4,   0.0          );
            MatInv.set(6,5,   0.0          );
            MatInv.set(6,6,  -9.0          );
            MatInv.set(6,7,   0.0          );
            MatInv.set(6,8,   0.0          );
            MatInv.set(6,9,   4.5          );
            MatInv.set(6,10,  0.0          );

            MatInv.set(7,1,  -4.5          );
            MatInv.set(7,2,  13.5          );
            MatInv.set(7,3, -13.5          );
            MatInv.set(7,4,   4.5          );
            MatInv.set(7,5,   0.0          );
            MatInv.set(7,6,   0.0          );
            MatInv.set(7,7,   0.0          );
            MatInv.set(7,8,   0.0          );
            MatInv.set(7,9,   0.0          );
            MatInv.set(7,10,  0.0          );

            MatInv.set(8,1, -13.5          );
            MatInv.set(8,2,  27.0          );
            MatInv.set(8,3, -13.5          );
            MatInv.set(8,4,   0.0          );
            MatInv.set(8,5,  13.5          );
            MatInv.set(8,6, -27.0          );
            MatInv.set(8,7,  13.5          );
            MatInv.set(8,8,   0.0          );
            MatInv.set(8,9,   0.0          );
            MatInv.set(8,10,  0.0          );

            MatInv.set(9,1, -13.5          );
            MatInv.set(9,2,  13.5          );
            MatInv.set(9,3,   0.0          );
            MatInv.set(9,4,   0.0          );
            MatInv.set(9,5,  27.0          );
            MatInv.set(9,6, -27.0          );
            MatInv.set(9,7,   0.0          );
            MatInv.set(9,8, -13.5          );
            MatInv.set(9,9,  13.5          );
            MatInv.set(9,10,  0.0          );

            MatInv.set(10,1,  -4.5          );
            MatInv.set(10,2,   0.0          );
            MatInv.set(10,3,   0.0          );
            MatInv.set(10,4,   0.0          );
            MatInv.set(10,5,  13.5          );
            MatInv.set(10,6,   0.0          );
            MatInv.set(10,7,   0.0          );
            MatInv.set(10,8, -13.5          );
            MatInv.set(10,9,   0.0          );
            MatInv.set(10,10,  4.5          );

            break;

        default:
            printf("\n");
            printf(" ERROR in L2Project_DG2CG_Unst.cpp: cg_order value not supported.\n");
            printf("       cg_order = %i\n",cg_order);
            printf("\n");
            exit(1);
            break;      
    }


    // Evaluate CG basis functions on quadrature points
    dTensor2 psi(mpoints,cg_num_pts);
    switch(cg_order)
    {
        case 1:     
            assert_eq(kmax,1);
            assert_eq(cg_num_pts,3);
            assert_eq(QuadOrder,2);
            assert_eq(mpoints,3);
            for (int m=1; m<=(mpoints); m++)
                for (int k=1; k<=(cg_num_pts); k++)
                {
                    psi.set(m,k, MatInv.get(1,k) 
                            + MatInv.get(2,k)*spts.get(m,1)
                            + MatInv.get(3,k)*spts.get(m,2) );
                }
            break;

        case 2:
            assert_eq(kmax,3);
            assert_eq(cg_num_pts,6);
            assert_eq(QuadOrder,3);
            assert_eq(mpoints,6);
            for (int m=1; m<=(mpoints); m++)
            {
                double s1 = spts.get(m,1);
                double s2 = spts.get(m,2);
                for (int k=1; k<=(cg_num_pts); k++)
                {
                    psi.set(m,k, MatInv.get(1,k) 
                            + MatInv.get(2,k)*s1
                            + MatInv.get(3,k)*s2
                            + MatInv.get(4,k)*s1*s2
                            + MatInv.get(5,k)*s1*s1
                            + MatInv.get(6,k)*s2*s2 );
                }
            }
            break;

        case 3:
            assert_eq(kmax,6);
            assert_eq(cg_num_pts,10);
            assert_eq(QuadOrder,4);
            assert_eq(mpoints,12);
            for (int m=1; m<=(mpoints); m++)
            {
                double s1 = spts.get(m,1);
                double s2 = spts.get(m,2);
                for (int k=1; k<=(cg_num_pts); k++)
                {
                    psi.set(m,k, MatInv.get(1,k) 
                            + MatInv.get(2,k) *s1
                            + MatInv.get(3,k) *s2
                            + MatInv.get(4,k) *s1*s2
                            + MatInv.get(5,k) *s1*s1
                            + MatInv.get(6,k) *s2*s2
                            + MatInv.get(7,k) *s1*s1*s1
                            + MatInv.get(8,k) *s1*s1*s2
                            + MatInv.get(9,k) *s1*s2*s2
                            + MatInv.get(10,k)*s2*s2*s2);
                }
            }
            break;

        default:
            printf("\n");
            printf(" ERROR in L2Project_DG2CG_Unst.cpp: cg_order value not supported.\n");
            printf("       cg_order = %i\n",cg_order);
            printf("\n");
            exit(1);
            break;
    }
    // ---------------------------------------------------------

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

    // Loop over each quadrature point and construct Legendre polys
    for (int m=1; m<=mpoints; m++)
        for (int i=1; i<=kmax; i++)
        {
            double tmp = 0.0;
            for (int j=1; j<=i; j++)
            {  tmp = tmp + Mmat[i-1][j-1]*mu.get(m,j);  }

            phi.set(m,i, tmp );      
        }

    // -------------------------------------------------------------
    // Loop over every element indexed by user supplied parameters
    // described by istart...iend
    // -------------------------------------------------------------
    for (int i=istart; i<=iend; i++)
    {      
        // find nodes attached to current element
        iTensor1 tt(cg_num_pts);
        switch(cg_order)
        {
            case 1:
                for (int k=1; k<=3; k++)
                {  tt.set(k, Mesh.get_tnode(i,k) );  }
                break;

            case 2:
                for (int k=1; k<=6; k++)
                {  tt.set(k, Mesh.get_node_subs(i,k) );  }
                break;

            case 3:
                for (int k=1; k<=10; k++)
                {  tt.set(k, Mesh.get_node_subs(i,k) );  }
                break;

            default:
                printf("\n");
                printf(" ERROR in L2Project_DG2CG_Unst.cpp: cg_order value not supported.\n");
                printf("       cg_order = %i\n",cg_order);
                printf("\n");
                exit(1);
                break;
        }

        // ---------------------------------------------------------
        // Evaluate DG function on quadrature points
        // ---------------------------------------------------------
        dTensor1 fdg_vals(mpoints);
        for (int m=1; m<=(mpoints); m++)
        {
            fdg_vals.set(m, 0.0 );

            for (int k=1; k<=kmax; k++)
            {
                fdg_vals.set(m, fdg_vals.get(m) 
                        + phi.get(m,k) * fdg.get(i,dg_comp,k) );
            }
        }
        // ---------------------------------------------------------

        // Evaluate integral on current cell (project onto CG basis) 
        // using Gaussian Quadrature for the integration    
        double Area = Mesh.get_area_prim(i);
        for (int k=1; k<=cg_num_pts; k++)
        {
            double integral = 0.0;

            for (int m=1; m<=mpoints; m++)
            {
                integral = integral + 2.0*Area*wgts.get(m)*psi.get(m,k)*fdg_vals.get(m);
            }
            fcg.set(tt.get(k), fcg.get(tt.get(k)) + integral );
        }
    }
}
