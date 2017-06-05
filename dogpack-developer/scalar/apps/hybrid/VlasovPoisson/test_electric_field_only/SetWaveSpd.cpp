#include "dog_math.h"
#include "tensors.h"
#include "assert.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Advection equation
//
void SetWaveSpd(const dTensor2* vel_vec, 
    const dTensor1& nvec, const dTensor1& xedge,
    const dTensor1& Ql, const dTensor1& Qr,
    const dTensor1& Auxl, const dTensor1& Auxr,
    double& s1,double& s2)
{

    s1 = 0.;
    s2 = 0.;

}
