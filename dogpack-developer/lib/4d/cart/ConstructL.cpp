//#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
//#include "FaceData.h"
//#include "DogParams.h"
//#include "DogParamsCart4.h"
//#include "L2Project.h"
//#include "Legendre4d.h"
//#include "RiemannSolve.h"

// Right-hand side for hyperbolic PDE in divergence form
//
void ConstructL(
    const dTensorBC6& aux, const dTensorBC6& q,
    dTensorBC6& Lstar, dTensorBC5& smax)
{
    printf("\n");
    printf("   ConstructL has not yet been implemented in 4D \n");
    printf("\n");
    exit(1);
}
