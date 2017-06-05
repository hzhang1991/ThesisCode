#include <cmath>
#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(const dTensor1& Aux_ave, 
        const dTensor1& Q_ave,
        const dTensor2& Qvals,
        dTensor2& Wvals)
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

    // Project onto left eigenvectors
    for (int k=1; k<=(kmax-1); k++)
    {
        Wvals.set(1,k, ((umag2/2.0-H-u1*c)*Qvals.get(2,k) 
                    + (umag2/2.0*(c-u1)+H*u1)*Qvals.get(1,k) 
                    + c*(Qvals.get(5,k)-u2*Qvals.get(3,k)
                        - u3*Qvals.get(4,k)))/(c*(2.0*H-umag2)) );

        Wvals.set(2,k, 2.0*((H-umag2)*Qvals.get(1,k) + u1*Qvals.get(2,k) 
                    + u2*Qvals.get(3,k) + u3*Qvals.get(4,k) 
                    - Qvals.get(5,k))/(2.0*H-umag2) );

        Wvals.set(3,k, Qvals.get(3,k)-u2*Qvals.get(1,k) );

        Wvals.set(4,k, Qvals.get(4,k)-u3*Qvals.get(1,k) );

        Wvals.set(5,k, ((H-umag2/2.0-u1*c)*Qvals.get(2,k) 
                    + (umag2/2.0*(c+u1)-H*u1)*Qvals.get(1,k)
                    + c*(Qvals.get(5,k)-u2*Qvals.get(3,k) 
                        - u3*Qvals.get(4,k)))/(c*(2.0*H-umag2)) );
    }
}
