#include <cmath>
#include "DogParams.h"
#include "tensors.h"
#include "dog_math.h"

#include "assert.h"         // TODO - remove this later!

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
// This function constructs the Lax-Friedrichs flux
void ConstructL_LLF(const int method[],
        const dTensor2& node,
        const dTensorBC1& smax,
        const dTensorBC3& aux,
        const dTensorBC3& q,      // setbndy conditions modifies q
        dTensorBC2& Fm,
        dTensorBC2& Fp)
{

    double RiemannSolve(const dTensor1& xedge,
            const dTensor1& Ql,
            const dTensor1& Qr,
            const dTensor1& Auxl,
            const dTensor1& Auxr,
            dTensor1& Fl,
            dTensor1& Fr,
            void (*FluxFunc)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&),
            void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
                const dTensor1&,double&,double&));
    void L2Project(int,int,int,const dTensor2&,const dTensorBC3&,const dTensorBC3&,dTensorBC3&,
            void (*Func)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&));
    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void SetWaveSpd(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void LstarExtra(const dTensor2&,dTensorBC3&,dTensorBC3&,dTensorBC3&);
    void ArtificialViscosity(const dTensor2&,dTensorBC3&,dTensorBC3&,dTensorBC3&);

    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double xlower = node.get(1,1);
    const double     dx = node.get(2,1) - node.get(1,1);

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    // TODO - do we want to assume these have been set before calling this
    // routine!?
//  SetBndValues(node,aux,q);

    // Loop over interior edges and solve Riemann problems
    //#pragma omp parallel for
    //
    // This sets both smax(i) and smax(i-1), so can't be parallelized!
    //
    for (int i=(2-mbc); i<=(melems+mbc); i++)
    {
        dTensor1   Ql(meqn),   Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {
            // Only consider the constant value for Q
            Ql.set(m, q.get(i-1,m,1) );
            Qr.set(m, q.get(i,  m,1) ); 
        }

        // Riemann data - aux
        for (int m=1; m<=maux; m++)
        {
            Auxl.set(m, 0.0 );
            Auxr.set(m, 0.0 );

            const int k = 1;
            Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                    *aux.get(i-1,m,k) );
            Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                    *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) ); 
        }

        // Solve Riemann problem
        dTensor1 xedge(1);
        xedge.set(1, xlower + (double(i)-1.0)*dx );

        dTensor1 Fl(meqn);
        dTensor1 Fr(meqn);
        double smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
                &FluxFunc,&SetWaveSpd);

        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {
            Fm.set(i,  m, Fr.get(m) );
            Fp.set(i-1,m, Fl.get(m) );
        }
    }
    // ---------------------------------------------------------

}
