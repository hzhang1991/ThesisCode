#include <cmath>
#include "dogdefs.h"
#include "EulerParams.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(int ixy, 
		     const dTensor1& Aux_ave,
		     const dTensor1& Q_ave, 
		     const dTensor2& Wvals, 
		     dTensor2& Qvals)
{    
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2)+1;
  
  // Direction
  int mu,mv,mw;
  if (ixy==1)
    {  
      mu = 2;
      mv = 3;
      mw = 4;
    }
  else
    {
      mu = 3;
      mv = 2;
      mw = 4;
    }
  
  // Average states
  const double gamma = eulerParams.gamma;
  const double rho    = Q_ave.get(1);
  const double u1     = Q_ave.get(mu)/rho;
  const double u2     = Q_ave.get(mv)/rho;
  const double u3     = Q_ave.get(mw)/rho;
  const double energy = Q_ave.get(5);
  const double umag2  = (u1*u1 + u2*u2 + u3*u3);
  const double press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
  const double c      = sqrt(fabs(gamma*press/rho));
  const double H      = (energy+press)/rho;  
  
  // Project onto right eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Qvals.set(1, k, Wvals.get(1,k) + Wvals.get(mu,k) + Wvals.get(5,k)  );
      
      Qvals.set(mu,k, (u1-c)*Wvals.get(1,k) + u1*Wvals.get(mu,k) 
		+ (u1+c)*Wvals.get(5,k) );
      
      Qvals.set(mv,k, u2*(Wvals.get(1,k) + Wvals.get(mu,k) + Wvals.get(5,k))
		+ Wvals.get(mv,k) );
      
      Qvals.set(mw,k, u3*(Wvals.get(1,k) + Wvals.get(mu,k) + Wvals.get(5,k))
		+ Wvals.get(mw,k) );
      
      Qvals.set(5, k, (H-u1*c)*Wvals.get(1,k) + umag2/2.0*Wvals.get(mu,k)
		+ u2*Wvals.get(mv,k) + u3*Wvals.get(mw,k) 
		+ (H+u1*c)*Wvals.get(5,k) );
    }
}
