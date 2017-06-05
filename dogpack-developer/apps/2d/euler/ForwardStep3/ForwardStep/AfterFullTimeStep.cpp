#include "dogdefs.h"
#include "DogSolverCart2.h"

// Function that is called after a full time step
// (i.e., after all stages are complete)
//
// This routine will zero out all of the elements located outside of the
// computational domain.  The next call to SetBndValues will override the bad
// elements.
void AfterFullTimeStep(DogSolverCart2& solver)
{
    const double t = solver.get_state().get_time();
    dTensorBC4& aux = solver.fetch_state().fetch_aux();
    dTensorBC4& q = solver.fetch_state().fetch_q();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    // Compute index where the step is located.  
    // q(istep,:) is inside the domain.
    const int istep = (mx / 5)+1;
    const int jstep = ((my / 5))+1;


    double gamma=1.4;

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;

    double energy = press/(gamma-1.0e0)
       + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);



    for( int i = istep; i <= mx+mbc; i++)
    for( int j = 1-mbc; j <= jstep-1; j++)
    {
       for( int m = 1; m <= meqn; m++ )
       for( int k = 1; k <= kmax; k++ )
       {
          q.set(i,j,m,k, 0.0 );
       }
       q.set(i,j,1,1, rho );
       q.set(i,j,2,1, u1 );
       q.set(i,j,3,1, u2 );
       q.set(i,j,4,1, u3 );
       q.set(i,j,5,1, energy );
    }
}   
