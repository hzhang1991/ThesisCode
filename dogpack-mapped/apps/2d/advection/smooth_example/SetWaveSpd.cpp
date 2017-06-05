#include "dog_math.h"
#include "tensors.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Simple advection equation
//
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
		const dTensor1& Ql, const dTensor1& Qr,
		const dTensor1& Auxl, const dTensor1& Auxr,
		double& s1,double& s2)
{
  double un_left,un_right,un_av;
  
  // Normal component of velocity
  un_left  = nvec.get(1)*Auxl.get(1) + nvec.get(2)*Auxl.get(2);
  un_right = nvec.get(1)*Auxr.get(1) + nvec.get(2)*Auxr.get(2);
  un_av    = 0.5*(un_left + un_right);
  
  // Minimum speed
  s1 = Min(un_left, un_av);
  
  // Maximum speed
  s2 = Max(un_right, un_av);
  
}
