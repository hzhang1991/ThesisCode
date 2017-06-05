#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "tensors.h"
#include "dog_math.h"
#include "LaxWendroffTD.h"
#include "dogdefs.h"

using namespace std;

// Right-hand side of the Lax-Wendroff discretization:
//
//      ( q(t+dt) - q(t) )/dt = -F_{x}.
//
// This routine constructs the function F.
//
// The function F is created by looking at derivatives of q.  There are two
// inputs here to handle combinations of each input value:
//
//    q^{n+1} = q^n + dt * (alpha1 q1_t + dt * beta1 q1_tt ) + 
//                    dt * (alpha2 q2_t + dt * beta2 q2_tt ),
//
// Which means that
//
//    ( q^{n+1} - q^n ) / dt =   (alpha1 q1_t + dt * beta1 q1_tt ) 
//                             + (alpha2 q2_t + dt * beta2 q2_tt )
//
// Each time derivative leaves a single x derivative that can be pulled out,
// so we write
//
//   q_{t}  = ( -f(q) )_x, 
//   q_{tt} = ( f'(q) * f_x )_x,
//
// And then work with the weak formulation for F:
//
//    F = ( f - dt * f'(q)*f_x )
//
// For the analogous RK time integrator, see ConstructL
//
void LaxWendroffTD(double dt, const int method[], const dTensor2& node,
    double alpha1,  double beta1,
    dTensorBC3& aux1, dTensorBC3& q1,   // set bndy values modifies these
    double alpha2, double beta2,
    dTensorBC3& aux2, dTensorBC3& q2,   // set bndy values modifies these
    dTensorBC3& Lstar,dTensorBC1& smax)
{

    const int melems = q1.getsize(1);
    const int   meqn = q1.getsize(2);
    const int   kmax = q1.getsize(3);
    const int   maux = aux1.getsize(2);
    const int    mbc = q1.getmbc();

    // -----------------
    // Quick error check
    // -----------------
    if (meqn<1 || maux < 0 || kmax < 1 || kmax>4 ) 
    {
        cout << " Error in LaxWendroffTD.cpp ... " << endl;
        cout << "      meqn    = " << meqn << endl;
        cout << "      kmax    = " << kmax << endl;
        cout << "      maux    = " << maux << endl;
        cout << endl;
        exit(1);
    }

    // Storage (for \tilde{f})
    dTensorBC2 Fm   (melems,meqn,mbc);
    dTensorBC2 Fp   (melems,meqn,mbc);
    dTensorBC3 N    (melems,meqn,kmax,mbc);
    dTensorBC3 Psi  (melems,meqn,kmax,mbc);

    // Weighted average value of the two input values.  
    //
    // The LLF Riemman solver needs values of Q to determine flux values,
    // we'll use a weighted average of q1 and q2 to determine those values.
    //
    // Specifically, we'll define:
    //     qavg = alpha1 * q1 + alpha2 * q2
    //
    // We need to assume here that alpha1 + alpha2 = 1.
    //
    // These values are passed into the LLF solver, which needs left and right
    // values for q.
    //
    // Note: for the unique fourth order time integrator, this will always set 
    // qavg = q1.
    dTensorBC3 qavg    (melems,meqn,kmax,mbc);
    dTensorBC3 auxavg  (melems,maux,kmax,mbc);

    // lax-wendroff flux function
    dTensorBC3  LxWF(melems,meqn,kmax,mbc);

    // storage (TODO - isn't this a bit of overkill ... ? )
    dTensorBC3 F1_1(melems,meqn,kmax,mbc);
    dTensorBC3 F1_2(melems,meqn,kmax,mbc);
    dTensorBC3 N1_1(melems,meqn,kmax,mbc);
    dTensorBC3 N1_2(melems,meqn,kmax,mbc);

    dTensorBC3 F2_1(melems,meqn,kmax,mbc);
    dTensorBC3 F2_2(melems,meqn,kmax,mbc);
    dTensorBC3 N2_1(melems,meqn,kmax,mbc);
    dTensorBC3 N2_2(melems,meqn,kmax,mbc);

    // quick error check
    if( method[7] > 0 )
    {
        cout << "    error: have not implemented source term for LxW solver "
            << endl;
        exit(1);
    }

    // Grid spacing
    const double xlower = node.get(1,1);
    const double dx     = node.get(2,1) - node.get(1,1);

    // double check to make sure these were initialized to zero...
    F1_1.setall(0.); F1_2.setall(0.);
    F2_1.setall(0.); F2_2.setall(0.);

    N1_1.setall(0.); N1_2.setall(0.);
    N2_1.setall(0.); N2_2.setall(0.);

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    SetBndValues(node,aux1,q1);
    SetBndValues(node,aux2,q2);

// TODO - for now I'm simply assuming that we're taking a single time step on
// q1 rather than some combination of q1 and q2.  I think this is a reasonable
// assumption to be making, for the time being ... (-DS)
CopyQ( q1,   qavg   );
CopyQ( aux1, auxavg );


    // compute projections onto basis functions, and derivatives of basis
    // functions of f, f'*f_x:
    L2ProjectLxWTD(1-mbc, melems+mbc, node, q1, aux1,  
            F1_1, F1_2, N1_1,  N1_2,
            &FluxFunc, &DFluxFunc );
    L2ProjectLxWTD(1-mbc, melems+mbc, node, q2, aux2,  
            F2_1, F2_2, N2_1,  N2_2,
            &FluxFunc, &DFluxFunc );

    // evaluate coefficients of LaxWendroff Flux Function 
#pragma omp parallel for
    for(int i=(1-mbc); i <= (melems+mbc); i++ )
    for(int k=1; k<= kmax; k++)	    
    for(int m=1; m<=meqn; m++)
    {

        double tmp  = alpha1 * F1_1.get(i,m,k) + ( beta1 * dt ) * F1_2.get(i,m,k);
        tmp        += alpha2 * F2_1.get(i,m,k) + ( beta2 * dt ) * F2_2.get(i,m,k);
        LxWF.set(i,m,k, tmp );
    }

    // Loop over interior edges and solve Riemann problems
#pragma omp parallel for
    for (int i=(2-mbc); i<=(melems+mbc); i++)
    {

        // Local storage
        dTensor1 Ql(meqn),Qr(meqn),Auxl(maux),Auxr(maux);
        dTensor1 LxFFluxl(meqn),  LxFFluxr(meqn);
        dTensor1 Fl(meqn),Fr(meqn);

        // Grid information
        dTensor1 xedge(1);
        xedge.set(1, xlower + (double(i)-1.0)*dx );

        // evaluate the Riemman data
        // EvalRiemmanData accesses q(i) and q(i-1), as well as LxW(i) and LxW(i-1).
        EvalRiemmanData(i, qavg, auxavg, LxWF, Ql, Qr, LxFFluxl, LxFFluxr, Auxl, Auxr); 
        double smax_edge = RiemannSolveLxW(xedge,Ql,Qr,Auxl,Auxr,
                LxFFluxl, LxFFluxr, Fl, Fr, &SetWaveSpd);

        // TODO: this *should* be a problem if the openmp flags are turned
        // on.
        // 
        // For some reason, we don't see this as a problem in the 2D code ...
        // (-DS)
        //
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i  )) );

        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {
            Fm.set(i,   m, Fr.get(m) );
            Fp.set(i-1, m, Fl.get(m) );
        }

    }

    // clean up boundary values
    SetBndValues(node,aux1,q1);
    SetBndValues(node,aux2,q2);
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part II: compute intra-element contributions
    // ---------------------------------------------------------
    //
    //   N = int( F(q,x,t) * phi_x, x )/dx
    //
    // Compute ``N'' by projecting flux function onto the 
    // gradient of Legendre polynomials
    // 
    // In other sections of the code, this is calling the function,
    // L2ProjectGrad.  Here, we already have access to these interior values,
    // so no need to perform that again.
    //
    if (method[1]>1)
    {

#pragma omp parallel for
        for(int i=1-mbc; i<=(melems+mbc); i++)           
        for(int k=1; k<=kmax; k++)
        for(int m=1; m<=meqn; m++)
        {

            double tmp  = alpha1 * N1_1.get(i,m,k) + ( beta1 * dt ) * N1_2.get(i,m,k);
            tmp        += alpha2 * N2_1.get(i,m,k) + ( beta2 * dt ) * N2_2.get(i,m,k);
            N.set( i, m, k, tmp );

        }

    }
    else // case kmax == 1 - there is no integral term in this case
    { N.setall(0.); }

    SetBndValues(node,aux1,q1);   // clean up boundary values
    SetBndValues(node,aux2,q2);   // clean up boundary values
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part III: compute source term
    // --------------------------------------------------------- 
    if ( method[7]>0 )
    // TODO this needs some work - needs to be expanded in taylor series as
    // well....
    {   
        cout << "    Warning: LxW solver needs to expand source term." << endl;
        cout << " method will be first order for source term function" << endl;

        // Set source term on computational grid
        // Set values and apply L2-projection
        L2Project(0,1-mbc,melems+mbc,node,qavg,auxavg,Psi,&SourceTermFunc);
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
    LstarExtra(node,aux1,q1,Lstar);
    // ---------------------------------------------------------

}
