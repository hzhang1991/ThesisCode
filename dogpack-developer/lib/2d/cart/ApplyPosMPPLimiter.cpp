#include "stdio.h"
#include <cmath>
#include<iostream>
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "constants.h"
#include "tensors.h"
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
void ApplyPosMPPLimiter( const double dt, const dTensorBC3& smax,
        const dTensorBC3& aux, const dTensorBC3& q, 
        dTensorBC3& fHat, dTensorBC3& gHat )
{
    // Parameters for the current grid
    const int   meqn = dogParams.get_meqn();
    const int   maux = dogParams.get_maux();
    const int     mx = dogParamsCart2.get_mx();
    const int     my = dogParamsCart2.get_my();
    const int   mbc  = q.getmbc();

    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();

    //Nothing for Now
    printf("Warning: $DOGPACK/2d/lib/ApplyPosMPPLimiter.cpp : no limiter defined\n");

}
