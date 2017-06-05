#include "dogdefs.h"
#include "QuadratureRules.h"

// No electric field for this problem
//
//    speeds( 1:num_quad_pts, 1:2 )
//
void SetLocalSpeeds( const QuadratureRules& QuadFuncs, 
    const int n, 
    const dTensor2& E1, const dTensor2& E2, 
    dTensor2& speeds )
{

    const int meqn = speeds.getsize(1);
    for( int me=1; me <= meqn; me++ )
    {
        // Turn off the electric field for this problem only
        speeds.set(me,1, 0.);
        speeds.set(me,2, 0.);
    }

}
