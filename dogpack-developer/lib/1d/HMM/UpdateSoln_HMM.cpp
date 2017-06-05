#include "../defs.h"

// Update the solution using the constructed Lstar
void UpdateSoln_HMM(double alpha1,double alpha2,double beta,double dt,
		    dTensor2 node,dTensorBC3 qn,dTensorBC3 qstar,
		    dTensorBC3 Lstar,dTensorBC3& qnew)
{
    int j,m,k;
    double tmp;
    int melems = qnew.getsize(1);
    int   meqn = qnew.getsize(2);
    int   kmax = qnew.getsize(3);
    int mbc = qnew.getmbc();
   
    for (k=1; k<=kmax; k++)
    {
        for (m=1; m<=meqn; m++)
        {
            for (j=(2-mbc); j<=(melems+mbc-1); j++)
            {
                tmp = alpha1*qn.get(j,m,k) + alpha2*qstar.get(j,m,k)
                      + beta*dt*Lstar.get(j,m,k);
                
                qnew.set(j,m,k, tmp );
            }
        }
    }

}
