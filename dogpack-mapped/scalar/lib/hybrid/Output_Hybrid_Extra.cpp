#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "MonomialsToLegendre.h"
#include "mesh.h"
#include <sstream>
#include <string>
#include "assert.h"

/* 
 * Extra call to output.  This can be swapped out without changing
 * Output_Hybrid.
 * 
 * This is a 'do nothing' file, that needs to be swapped out ot be useful.
 */
void Output_Hybrid_Extra(const mesh& Mesh, const dTensorBC5& q, 
        double t, int nframe, string outputdir)
{ 

    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);
    const int     meqn = q.getsize(4);
    const int     kmax = q.getsize(5);
  
    // TODO - of course these need to be modified!
    const int NumMoments = 1;  

    int kmax2d;
    switch( dogParams.get_space_order() )
    {

        case 3:
            kmax2d = 6;
            break;

        case 2:
            kmax2d = 3;
            break;

        case 1:
            kmax2d = 1;
            break;
        default:
            unsupported_value_error( dogParams.get_space_order() );
    }

}
