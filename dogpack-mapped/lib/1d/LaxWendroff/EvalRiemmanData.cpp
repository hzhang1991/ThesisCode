#include <cmath>
#include "tensors.h"
#include "stdlib.h"
#include "stdio.h"

// Evaluate the Riemman data
//
// This access q(i) q(i-1), aux(i) aux(i-1) and LxWF(i), and LxWF(i-1).
// 
// This routine simply evaluates the basis elements provided by 
//     q(i, 1:meqn, 1:kmax) and q(i-1, 1:meqn, 1:kmax )
// 
// to produce two constant values F_l and F_r to be fed into the Riemann
// solver.
//
// For the Lax-Wendroff scheme, we do NOT want to call the flux function here,
// because we need to use the time derivatives and the LxW flux function in
// place of f(ql) and f(qr).  In this case, we assume that the use has already saved
// \tilde{f} = LxWF, the Lax-Wendroff flux function that contains extra terms
// from time derivatives.
//
void EvalRiemmanData(const int& i, const dTensorBC3& q, const dTensorBC3& aux,  
        const dTensorBC3& LxWF, dTensor1& Ql, dTensor1& Qr, 
        dTensor1& LxFFluxl, dTensor1& LxFFluxr,
        dTensor1& Auxl, dTensor1& Auxr)
{

    // Static variables
    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    for (int m=1; m<=meqn; m++)
    {
        // Riemann data - q
        Ql.set(m, 0.0 );
        Qr.set(m, 0.0 );
        LxFFluxl.set(m, 0.0 );
        LxFFluxr.set(m, 0.0 );
        for (int k=1; k<=kmax; k++)
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

    for(int m=1; m<=maux; m++)
    {
        // Riemann data - aux
        Auxl.set(m, 0.0 );
        Auxr.set(m, 0.0 );
        for (int k=1; k<=kmax; k++)
        {
            Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                    *aux.get(i-1,m,k) );
            Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                    *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) ); 
        }
    }
}
