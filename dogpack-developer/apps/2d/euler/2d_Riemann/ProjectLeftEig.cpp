#include <cmath>
#include "dogdefs.h"
#include "EulerParams.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(int ixy, 
		    const dTensor1& Aux_ave,
		    const dTensor1& Q_ave, 
		    const dTensor2& Qvals, 
		    dTensor2& Wvals)
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
  double const gamma  = eulerParams.gamma;
  double const rho    = Q_ave.get(1);
  double const u1     = Q_ave.get(mu)/rho;
  double const u2     = Q_ave.get(mv)/rho;
  double const u3     = Q_ave.get(mw)/rho;
  double const energy = Q_ave.get(5);
  double const umag2  = (u1*u1 + u2*u2 + u3*u3);
  double const press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
  double const c      = sqrt(gamma*press/rho);
  double const H      = (energy+press)/rho; 
  
  // Project onto left eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Wvals.set(1, k, ((umag2/2.0-H-u1*c)*Qvals.get(mu,k) 
		       + (umag2/2.0*(c-u1)+H*u1)*Qvals.get(1,k) 
		       + c*(Qvals.get(5,k)-u2*Qvals.get(mv,k)
			    - u3*Qvals.get(mw,k)))/(c*(2.0*H-umag2)) );
      
      Wvals.set(mu,k, 2.0*((H-umag2)*Qvals.get(1,k) + u1*Qvals.get(mu,k) 
			   + u2*Qvals.get(mv,k) + u3*Qvals.get(mw,k) 
			   - Qvals.get(5,k))/(2.0*H-umag2) );
      
      Wvals.set(mv,k, Qvals.get(mv,k)-u2*Qvals.get(1,k) );
      
      Wvals.set(mw,k, Qvals.get(mw,k)-u3*Qvals.get(1,k) );
      
      Wvals.set(5, k, ((H-umag2/2.0-u1*c)*Qvals.get(mu,k) 
		       + (umag2/2.0*(c+u1)-H*u1)*Qvals.get(1,k)
		       + c*(Qvals.get(5,k)-u2*Qvals.get(mv,k) 
			    - u3*Qvals.get(mw,k)))/(c*(2.0*H-umag2)) );
    }
}
