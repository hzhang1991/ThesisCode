#include <cmath>
#include "tensors.h"
#include "stdlib.h"

void EvalRiemmanData(const int& i, const dTensorBC3& q, const dTensorBC3& aux,  
		     const dTensorBC3& LxWF, dTensor1& Ql, dTensor1& Qr, 
		     dTensor1& LxFFluxl, dTensor1& LxFFluxr,
		     dTensor1& Auxl, dTensor1& Auxr)
{

    int k,m;
    double tmp;
    int melems = q.getsize(1);
    int   meqn = q.getsize(2);
    int   kmax = q.getsize(3);
    int   maux = aux.getsize(2);
    int    mbc = q.getmbc();

    for (m=1; m<=meqn; m++)
    {
        // Riemann data - q
        Ql.set(m, 0.0 );
        Qr.set(m, 0.0 );
        LxFFluxl.set(m, 0.0 );
        LxFFluxr.set(m, 0.0 );
        for (k=1; k<=kmax; k++)
        {
            Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                    *q.get(i-1,m,k) );
            Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                    *sqrt(2.0*double(k)-1.0)*q.get(i,m,k) ); 
            LxFFluxl.set(m, LxFFluxl.get(m) + sqrt(2.0*double(k)-1.0)
                    *LxWF.get(i-1,m,k) );
            LxFFluxr.set(m, LxFFluxr.get(m) + pow(-1.0,k+1)
                    *sqrt(2.0*double(k)-1.0)*LxWF.get(i,m,k) ); 
        }

    }

    for(m=1; m<=maux; m++)
    {
        // Riemann data - aux
        Auxl.set(m, 0.0 );
        Auxr.set(m, 0.0 );
        for (k=1; k<=kmax; k++)
        {
            Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                    *aux.get(i-1,m,k) );
            Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                    *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) ); 
        }
    }
}
