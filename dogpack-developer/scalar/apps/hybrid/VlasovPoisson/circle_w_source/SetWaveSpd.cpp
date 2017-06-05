#include "dog_math.h"
#include "tensors.h"

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

//  double un_left,un_right,un_av;

    // Normal component of velocity
//  un_left  = nvec.get(1)*Auxl.get(1) + nvec.get(2)*Auxl.get(2);
//  un_right = nvec.get(1)*Auxr.get(1) + nvec.get(2)*Auxr.get(2);
//  un_av = 0.5*(un_left + un_right);

    // Minimum speed
//  s1 = Min(un_left, un_av);

    // Maximum speed
//  s2 = Max(un_right, un_av);

    // vel_vec( 1:meqn , 1:2 ) defines the wave speed for each equation in this
    // set.  Here, we need to define the fastest and slowest wave speed for each
    // equation:
    s1 = nvec.get(1)*vel_vec->get(1,1) + nvec.get(2)*vel_vec->get(1,2); 
    s2 = s1;
    for( int n=2; n <= vel_vec->getsize(1); n++ )
    {
        s1 = Min( s1, nvec.get(1)*vel_vec->get(n,1) + nvec.get(2)*vel_vec->get(n,2) );
        s2 = Max( s2, nvec.get(1)*vel_vec->get(n,1) + nvec.get(2)*vel_vec->get(n,2) );
    }

}
