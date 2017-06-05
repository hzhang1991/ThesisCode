#include <cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "assert.h"
#include "QuadratureRules.h"

// Read in q(i,j,:, 1:meqn=1, 1:kmax4d) into qout(:, 1:mpoints_cart, 1:kmax2d )
void ReadSlice( 
    const QuadratureRules& QuadFuncs, 
    const dTensorBC5& qin, const int i, const int j, dTensor3& qout )
{

    const int NumElems     = qout.getsize(1);
    const int mpoints_cart = qout.getsize(2);
    const int kmax2d       = qout.getsize(3);

    const int kmax = qin.getsize(5);

    // 2D weights, points and polynomials:
    const dTensor1& wgts_unst = QuadFuncs.get_wgt_unst();
    const dTensor2& spts_unst = QuadFuncs.get_spts_unst();
    const dTensor2& phi_unst  = QuadFuncs.get_phi_unst();

    const int mpoints_unst = spts_unst.getsize(1);

    // 4D weights, points and polynomials:
    const dTensor2& wgts = QuadFuncs.get_wgt();
    const dTensor3& spts = QuadFuncs.get_spts();
    const dTensor3& phi  = QuadFuncs.get_phi();

    // This loop is assumed to be 'fast', because it is in 2D, and not 4D
    for( int n=1; n <= NumElems; n++ )
    for( int m2=1; m2 <= mpoints_cart; m2++ )
    {
        for( int k2=1; k2 <= kmax2d; k2++ )
        {
            double tmp = 0.;
            for( int m1=1; m1 <= mpoints_unst; m1++ )
            {
                // evaluate q at this point:
                double q = 0.;
                for( int  k=1; k <= kmax; k++ )
                {
                    q += qin.get(i,j,n, 1, k ) * phi.get(m1,m2,k);
                }
                tmp += q * phi_unst.get(m1,k2) * wgts_unst.get(m1);
            }
            qout.set(n, m2, k2, 2.0*tmp );
        }
    }

}

void WriteSlice( 
    const QuadratureRules& QuadFuncs, 
    const dTensor3& qin, const int i, const int j, dTensorBC5& qout )
{

    const int sorder = dogParams.get_space_order();

    const int mx   = dogParamsCart2.get_mx();
    const int my   = dogParamsCart2.get_my();
    const int kmax = qout.getsize(5);

    const int NumElems = qin.getsize(1);

    assert_eq( kmax, QuadFuncs.get_kmax() );

    const int& kmax2d = QuadFuncs.get_kmax2d();

    // 4D weights, points and polynomials:
    const dTensor2& wgts = QuadFuncs.get_wgt();
    const dTensor3& spts = QuadFuncs.get_spts();
    const dTensor3& phi  = QuadFuncs.get_phi();

    const int& mpoints_unst = wgts.getsize(1);   
    const int& mpoints_cart = wgts.getsize(2);   
    assert_eq( mpoints_cart, qin.getsize(2) );

    const dTensor2& phi_unst  = QuadFuncs.get_phi_unst();
    for( int n=1; n <= NumElems; n++ )
    {

        // loop over each polynomial:
        for( int k=1; k <= kmax; k++ )
        {
            double tmp = 0.;
            for( int m1=1; m1 <= mpoints_unst; m1++ )
            {
                for( int m2=1; m2 <= mpoints_cart; m2++ )
                {
                    // Evaluate q at this point:
                    double q = 0.;
                    for( int k2=1; k2 <= kmax2d; k2++ )
                    {
                        q += phi_unst.get(m1,k2) * qin.get(n, m2, k2 );
                    }

                    // Compute phi*q, add this to the running total:
                    tmp += q * phi.get(m1,m2,k)*wgts.get(m1,m2);

                }
            }
            qout.set(i, j, n, 1, k, 0.5*tmp );
        }
    }

}

// Read in q(:,:,n, 1:meqn=1, 1:kmax4d) into qout(:, :, 1:npts_unst, 1:kmax2d )
void ReadSlice( 
    const QuadratureRules& QuadFuncs, 
    const dTensorBC5& qin, const int n, dTensorBC4& qout )
{

    const int mx            = qout.getsize(1);  assert_eq( mx, dogParamsCart2.get_mx() );
    const int my            = qout.getsize(2);  assert_eq( my, dogParamsCart2.get_my() );
    const int mpoints_unst  = qout.getsize(3);
    const int kmax2d        = qout.getsize(4);

    const int kmax = qin.getsize(5);

    // 2D weights, points and polynomials:
    const dTensor1& wgts_cart = QuadFuncs.get_wgt_cart();
    const dTensor2& spts_cart = QuadFuncs.get_spts_cart();
    const dTensor2& phi_cart  = QuadFuncs.get_phi_cart();

    const int mpoints_cart = spts_cart.getsize(1);

    // 4D weights, points and polynomials:
    const dTensor2& wgts = QuadFuncs.get_wgt();
    const dTensor3& spts = QuadFuncs.get_spts();
    const dTensor3& phi  = QuadFuncs.get_phi();

    // This loop is assumed to be 'fast', because it is in 2D, and not 4D
    for( int i=1; i <= mx; i++ )
    for( int j=1; j <= my; j++ )
    for( int m1=1; m1 <= mpoints_unst; m1++ )
    {
        for( int k2=1; k2 <= kmax2d; k2++ )
        {
            double tmp = 0.;
            for( int m2=1; m2 <= mpoints_cart; m2++ )
            {
                // evaluate q at this point:
                double q = 0.;
                for( int  k=1; k <= kmax; k++ )
                {
                    q += qin.get(i,j,n, 1, k ) * phi.get(m1,m2,k);
                }
                tmp += q * phi_cart.get(m2,k2) * wgts_cart.get(m2);
            }
            qout.set(i, j, m1, k2, 0.25*tmp );
        }
    }

}

void WriteSlice( 
    const QuadratureRules& QuadFuncs, 
    const dTensorBC4& qin, const int n, dTensorBC5& qout )
{

    const int mx = qin.getsize(1);  assert_eq( mx, dogParamsCart2.get_mx() );
    const int my = qin.getsize(2);  assert_eq( my, dogParamsCart2.get_my() );
    const int kmax2d = qin.getsize(4);
    const int kmax = qout.getsize(5);

    // 4D weights, points and polynomials:
    const dTensor2& wgts = QuadFuncs.get_wgt();
    const dTensor3& spts = QuadFuncs.get_spts();
    const dTensor3& phi  = QuadFuncs.get_phi();

    const int& mpoints_unst = wgts.getsize(1);   
    const int& mpoints_cart = wgts.getsize(2);   

    const dTensor2& phi_cart  = QuadFuncs.get_phi_cart();

    for( int i=1; i <= mx; i++ )
    for( int j=1; j <= my; j++ )
    {

        // loop over each polynomial:
        for( int k=1; k <= kmax; k++ )
        {
            double tmp = 0.;
            for( int m2=1; m2 <= mpoints_cart; m2++ )
            {
                for( int m1=1; m1 <= mpoints_unst; m1++ )
                {
                    // Evaluate q at this point:
                    double q = 0.;
                    for( int k2=1; k2 <= kmax2d; k2++ )
                    {
                        q += phi_cart.get(m2,k2) * qin.get(i,j, m1, k2 );
                    }

                    // Compute phi*q, add this to the running total:
                    tmp += q * phi.get(m1,m2,k)*wgts.get(m1,m2);
                }
            }
            qout.set(i, j, n, 1, k, 0.5*tmp );
        }
    }
}
