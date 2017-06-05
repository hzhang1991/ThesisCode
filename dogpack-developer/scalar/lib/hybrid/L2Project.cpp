#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "QuadratureRules.h"

#include "assert.h"

void L2Project_Unst(
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
    void (*Func)(
        const dTensor2* vel_vec,
        const dTensor2&,const dTensor2&,
        const dTensor2&,dTensor2&));

void AuxFuncWrapper(
    const dTensor2* vel_vec,
    const dTensor2& xpts,
    const dTensor2& NOT_USED_1,
    const dTensor2& NOT_USED_2,
    dTensor2& auxvals);

void QinitFuncWrapper(
    const dTensor2* vel_vec,
    const dTensor2& xpts,
    const dTensor2& NOT_USED_1,
    const dTensor2& NOT_USED_2,
    dTensor2& qvals);

void WriteSlice( 
    const QuadratureRules& QuadFunc, 
    const dTensor3&   qin, const int i, const int j, dTensorBC5& qout );

// -------------------------------------------------------------
// Routine for computing the L2-projection of an input function
// onto an orthonormal Legendre basis
// -------------------------------------------------------------
void L2Project( const mesh& Mesh, dTensorBC5& q)
{
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    const int space_order       = dogParams.get_space_order();
    const int num_quad_pts_cart = space_order*space_order;

    const int mx  = dogParamsCart2.get_mx();    assert_eq( mx, q.getsize(1) );
    const int my  = dogParamsCart2.get_my();    assert_eq( my, q.getsize(2) );
    const int kmax = q.getsize(5);

    const int NumElems = q.getsize(3);

    // Figure out what all of the quadrature weights and points are:
    QuadratureRules QuadFuncs;
    QuadFuncs.init( dogParams.get_space_order() );

    const dTensor2& phi_unst  = QuadFuncs.get_phi_unst();
    const dTensor2& spts_unst = QuadFuncs.get_spts_unst();

    const dTensor2& phi_cart  = QuadFuncs.get_phi_cart();
    const dTensor2& spts_cart = QuadFuncs.get_spts_cart();

    const int kmax2d = phi_unst.getsize(2);

// TODO - uncommenting this statement produces a bug.  Run
// test_initial_conditions with space_order == 2.  There is a single triangle
// that produces garbage for each moment.
// #pragma omp parallel for
    for( int i=1; i <= mx; i++ )
    {

        dTensor3 qslice(NumElems, num_quad_pts_cart, kmax2d );
        dTensor3    aux(NumElems,                 0, kmax2d );

        dTensor2* vel_vec = new dTensor2(num_quad_pts_cart,2);

        for( int j=1; j <= my; j++ )
        {


            // Evaluate the speeds, vx, and vy throughout the current cell at
            // each required quadrature point:
            const double xc = dogParamsCart2.get_xc(i);
            const double yc = dogParamsCart2.get_yc(j);
            for( int s=1; s <= num_quad_pts_cart; s++ )
            {  
                vel_vec->set(s,1, xc + 0.5*dx*spts_cart.get(s, 1 ) );
                vel_vec->set(s,2, yc + 0.5*dy*spts_cart.get(s, 2 ) );
            }
            L2Project_Unst(vel_vec, 1, NumElems,
                    space_order,space_order,space_order,space_order,		       
                    Mesh,&qslice,&aux,&qslice,&QinitFuncWrapper );  
            WriteSlice( QuadFuncs, qslice, i, j , q );

        }

        delete vel_vec;

    }

}
