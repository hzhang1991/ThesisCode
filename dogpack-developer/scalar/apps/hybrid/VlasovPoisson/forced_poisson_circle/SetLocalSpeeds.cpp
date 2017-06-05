#include "dogdefs.h"
#include "QuadratureRules.h"

// Default routine for setting the advection speeds
//
// This routine evaluates
//     E1(n,1:kmax2d), and E2(n,1:kmax2d), 
//
// at the set of quadrature points, and produces
//
//    speeds( 1:num_quad_pts, 1:2 )
//
//
void SetLocalSpeeds( const QuadratureRules& QuadFuncs, 
    const int n, 
    const dTensor2& E1, const dTensor2& E2, 
    dTensor2& speeds )
{

    const dTensor2& phi_unst = QuadFuncs.get_phi_unst();

    const int num_quad_pts   = phi_unst.getsize(1);
    const int kmax2d         = phi_unst.getsize(2);
    for( int me=1; me <= num_quad_pts; me++ )
    {

        // Evaluate e at this point:
        double e1 = 0.;
        double e2 = 0.;
        for( int k2=1; k2 <= kmax2d; k2++ )
        {
            e1 += phi_unst.get(me,k2) * E1.get(n, k2 );
            e2 += phi_unst.get(me,k2) * E2.get(n, k2 );
        }

        // f_t + v \cdot f_x - E \cdot f_v = 0 - we take the negative E here.
        speeds.set(me,1, -e1 );
        speeds.set(me,2, 0.  );

    }

}
