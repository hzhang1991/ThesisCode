#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig( int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& Wvals, dTensor2& Qvals)
{    

    const int meqn = Qvals.getsize(1);
    const int kmax = Qvals.getsize(2)+1;


    // Direction
    int mu,mv;
    if (ixy==1)
    {  
        mu = 2;
        mv = 3;
    }
    else
    {
        mu = 3;
        mv = 2;
    }
   
    double c = Aux_ave.get(1)*Aux_ave.get(2);

    // Project onto right eigenvectors
    for (int k=1; k<=(kmax-1); k++)
    {
        Qvals.set(1, k, ( -c ) * Wvals.get(1,k) + ( c ) * Wvals.get(3,k)  );
        Qvals.set(mu,k,          Wvals.get(1,k) +         Wvals.get(3,k) );
        Qvals.set(mv,k,          Wvals.get(2,k) );
    }
  // Project onto right eigenvectors


/*
 for (int m=1;m<=meqn;m++)
  for (int k=1; k<=kmax; k++)
    {
      Qvals.set(m,k, Wvals.get(m,k) );
    }*/

}
