#include <cmath>
#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, 
		     const dTensor2& Wvals,
		     dTensor2& Qvals)
{    
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2)+1;
  
  // Average states
  const double gamma  = eulerParams.gamma;
  const double rho    = Q_ave.get(1);
  const double u1     = Q_ave.get(2)/rho;
  const double u2     = Q_ave.get(3)/rho;
  const double u3     = Q_ave.get(4)/rho;
  const double energy = Q_ave.get(5);
  const double umag2  = (u1*u1 + u2*u2 + u3*u3);
  const double press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
  const double c      = sqrt(fabs(gamma*press/rho));
  const double H      = (energy+press)/rho;  
  
  // Project onto right eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Qvals.set(1,k, Wvals.get(1,k) + Wvals.get(2,k) + Wvals.get(5,k)  );
      
      Qvals.set(2,k, (u1-c)*Wvals.get(1,k) + u1*Wvals.get(2,k) 
		+ (u1+c)*Wvals.get(5,k) );
      
      Qvals.set(3,k, u2*(Wvals.get(1,k) + Wvals.get(2,k) + Wvals.get(5,k))
		+ Wvals.get(3,k) );
      
      Qvals.set(4,k, u3*(Wvals.get(1,k) + Wvals.get(2,k) + Wvals.get(5,k))
		+ Wvals.get(4,k) );
      
      Qvals.set(5,k, (H-u1*c)*Wvals.get(1,k) + umag2/2.0*Wvals.get(2,k)
		+ u2*Wvals.get(3,k) + u3*Wvals.get(4,k) 
		+ (H+u1*c)*Wvals.get(5,k) );
    }
}
