#include "dogdefs.h"

// Optional call to function that adds extra piece into
// solution right-hand side (i.e., q_t = Lstar).  This routine is called towards
// the end of ConstructL.  That is, the boundary integral (with Riemann solves),
// and interiour integrals have already been added into Lstar when this function
// gets called.
//
// Input:
// ------
//
//      q( 1:mx, 1:my, 1:meqn, 1:kmax ) - solution that was used to ConstructL
//    aux( 1:mx, 1:my, 1:maux, 1:kmax ) - Auxiliary array
//
// Returns:
// --------
//
//  Lstar( 1:mx, 1:my, 1:meqn, 1:kmax ) - Right hand side of MOL formulation
//
// See also: ConstructL
void LstarExtra(const dTensorBC4* q, const dTensorBC4* aux, dTensorBC4* Lstar)
{

}
