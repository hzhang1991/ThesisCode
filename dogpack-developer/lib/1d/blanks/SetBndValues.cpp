#include "tensors.h"

// *REQUIRED*
//
// Boundary conditions specified the the application.
//
// Examples include: periodic, zeroth-order extrapolation, 
//                   reflective (Euler eqns.), ...
//
// Inputs:
//
// Outputs:
//
//     Modified boundary data for the following arrays:
//
//     aux(1-mbc:mx+mbc, 1:maux, 1:kmax ) - The auxiliary function
//       q(1-mbc:mx+mbc, 1:meqn, 1:kmax ) - The vector of conserved variables
//
// See also: 
void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

}
