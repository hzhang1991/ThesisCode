#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "tensors.h"
#include "dog_math.h"
#include "assert.h"

using namespace std;

// Maximal principle preserving flux limiter.  This routine should be defined
// separately for each set of equations it's intended to solve.
//
// Advantage: provable positivity, provided that Lax-Friedrichs flux is
// positive for the problem of interest.
//
// Disadvantage: high-order accuracy is only preserved numerically.  We
// currently do not have a proof of this for a multi-D problem.
//
// References: http://arxiv.org/abs/1411.0328
//
// This limiter has been implemented for 1D advection, and 1D Euler is a work
// in progress.
//
void ApplyPosMPPLimiter(const double dt, const int method[], const dTensor2& node,
        const dTensorBC1& smax,
        const dTensorBC3& aux, const dTensorBC3& q, 
        dTensorBC2& Fm, dTensorBC2& Fp )
{

    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // Lax-Friedrich's fluxes
    dTensorBC2 FmLF(melems, meqn, mbc );
    dTensorBC2 FpLF(melems, meqn, mbc );

    // Limiting parameter
    dTensorBC2 Theta(melems, meqn, mbc );

    // Range for limiting parameter, Lambda_{+1/2}, Lambda_{-1/2}.
    dTensorBC2 LambdaP(melems, meqn, mbc );
    dTensorBC2 LambdaM(melems, meqn, mbc );

    printf("Warning: $DOGPACK/1d/lib/ApplyPosMPPLimiter.cpp : no limiter defined\n");

}
