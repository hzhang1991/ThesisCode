#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "RiemannSolve.h"

// to use this RiemannSolver the user must implement the following
// "linking callbacks":
// * FluxFunc
// * SetWaveSpd

void SetBndFlux(const dTensor1& nvec, const dTensor1& xedge,
		    dTensor1& Ql, dTensor1& Qr,
		    const dTensor1& Auxl, const dTensor1& Auxr,
		    int side)
{

}

