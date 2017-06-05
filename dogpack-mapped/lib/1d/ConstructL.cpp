#include <cmath>
#include "DogParams.h"
#include "tensors.h"
#include "dog_math.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
void ConstructL(const int method[],
        const dTensor2& node,
        dTensorBC3& aux,
        dTensorBC3& q,      // setbndy conditions modifies q
        dTensorBC3& Lstar,
        dTensorBC1& smax)
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
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
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

    // Flux values, interior integrals and sourse term, respectively.
    dTensorBC2 Fm(melems,meqn,mbc);
    dTensorBC2 Fp(melems,meqn,mbc);
    dTensorBC3  N(melems,meqn,kmax,mbc);
    dTensorBC3 Psi(melems,meqn,kmax,mbc);

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double xlower = node.get(1,1);
    const double     dx = node.get(2,1) - node.get(1,1);

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    SetBndValues(node,aux,q);

    // Loop over interior edges and solve Riemann problems
    //#pragma omp parallel for
    //
    // This sets both smax(i) and smax(i-1), so can't be parallelized!
    //
    for (int i=(2-mbc); i<=(melems+mbc); i++)
    {
        dTensor1 Ql(meqn);
        dTensor1 Qr(meqn);
        dTensor1 Auxl(meqn);
        dTensor1 Auxr(meqn);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {
            Ql.set(m, 0.0 );
            Qr.set(m, 0.0 );

            for (int k=1; k<=method[1]; k++)
            {
                Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                        *q.get(i-1,m,k) );
                Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*q.get(i,m,k) ); 
            }
        }

        // Riemann data - aux
        for (int m=1; m<=maux; m++)
        {

            Auxl.set(m, 0.0 );
            Auxr.set(m, 0.0 );

            for (int k=1; k<=method[1]; k++)
            {
                Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                        *aux.get(i-1,m,k) );
                Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) ); 
            }
        }

        // Solve Riemann problem
        dTensor1 xedge(1);
        xedge.set(1, xlower + (double(i)-1.0)*dx );

        dTensor1 Fl(meqn);
        dTensor1 Fr(meqn);
        double smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
                &FluxFunc,&SetWaveSpd);

        // This is a problem for the pragma statements! (-DS)
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i)) );

        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {
            Fm.set(i,  m, Fr.get(m) );
            Fp.set(i-1,m, Fl.get(m) );
        }
    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part II: compute intra-element contributions
    // ---------------------------------------------------------
    //
    //   N = int( f(q,x,t) * phi_x, x )/dx
    //
    // Compute ``N'' by projecting flux function onto the 
    // gradient of Legendre polynomials
    if (method[1]>1)
    { 
        L2Project(1,1-mbc,melems+mbc,node,q,aux,N,&FluxFunc);
    }
    else
    {
#pragma omp parallel for
        for (int i=1-mbc; i<=(melems+mbc); i++)
            for (int m=1; m<=meqn; m++)
            {
                N.set(i,m,1, 0.0 );   
            }	
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part III: compute source term
    // --------------------------------------------------------- 
    if ( method[7]>0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        L2Project(0,1-mbc,melems+mbc,node,q,aux,Psi,&SourceTermFunc);
    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part IV: construct Lstar
    // ---------------------------------------------------------
    if (method[7]==0)  // Without Source Term
    {
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
            for (int m=1; m<=meqn; m++)
                for (int k=1; k<=method[1]; k++)
                {
                    double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                        ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx;

                    Lstar.set(i,m,k, tmp );	      
                }
    }
    else  // With Source Term
    {
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
            for (int m=1; m<=meqn; m++)	    
                for (int k=1; k<=method[1]; k++)
                {
                    double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                        ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx
                        + Psi.get(i,m,k);

                    Lstar.set(i,m,k, tmp );
                }
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part V: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra(node,aux,q,Lstar);
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part VI: artificial viscosity limiter
    // ---------------------------------------------------------
    if (dogParams.get_space_order()>1)
    {
        if (dogParams.using_viscosity_limiter())
        {  ArtificialViscosity(node,aux,q,Lstar);  }
    }
    // ---------------------------------------------------------
}
