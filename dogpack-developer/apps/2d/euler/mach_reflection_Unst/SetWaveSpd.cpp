#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "EulerParams.h"
#include <iostream>
using namespace std;

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Euler equations for gas dynamics
//
void SetWaveSpd(const dTensor1& nvec, 
		const dTensor1& xedge,
		const dTensor1& Ql, 
		const dTensor1& Qr,
		const dTensor1& Auxl, 
		const dTensor1& Auxr,
		double& s1,double& s2)
{
  // Gas constant
  const double gamma = 1.4;//eulerParams.gamma;
  
  // Left states
  const double rhol    = Ql.get(1);
  const double u1l     = Ql.get(2)/rhol;
  const double u2l     = Ql.get(3)/rhol;
  const double u3l     = Ql.get(4)/rhol;
  const double energyl = Ql.get(5);
  const double pressl  = (gamma-1.0e0)*(energyl-0.5e0*rhol*(u1l*u1l+u2l*u2l+u3l*u3l));
  
  // Right states
  const double rhor    = Qr.get(1);
  const double u1r     = Qr.get(2)/rhor;
  const double u2r     = Qr.get(3)/rhor;
  const double u3r     = Qr.get(4)/rhor;
  const double energyr = Qr.get(5);
  const double pressr  = (gamma-1.0e0)*(energyr-0.5e0*rhor*(u1r*u1r+u2r*u2r+u3r*u3r));
  
  // Average states
  const double rho     = 0.5e0*(rhol+rhor);
  const double u1      = 0.5e0*(u1l+u1r);
  const double u2      = 0.5e0*(u2l+u2r);
  const double u3      = 0.5e0*(u3l+u3r);
  const double press   = 0.5e0*(pressl+pressr);
  
  // Sound speeds
  const double nmag = sqrt(pow(nvec.get(1),2) + pow(nvec.get(2),2));
  const double cl = nmag*sqrt((gamma*pressl/rhol));
  const double cr = nmag*sqrt((gamma*pressr/rhor));
  const double c  = nmag*sqrt((gamma*press/rho));
  // normal velocities
  const double un  = nvec.get(1)*u1  + nvec.get(2)*u2;
  const double unl = nvec.get(1)*u1l + nvec.get(2)*u2l;
  const double unr = nvec.get(1)*u1r + nvec.get(2)*u2r;

  // Minimum speed
  s1 = Min(unl-cl, un-c);
  
  // Maximum speed
  s2 = Max(unr+cr, un+c);
}
