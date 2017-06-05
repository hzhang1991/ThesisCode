#include <cmath>
#include <iostream>
#include "stdio.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "constants.h"
#include "tensors.h"
#include "EulerParams.h"        // for gas constant, gamma

using namespace std;

// TODO - document this file.
//
// MPP limiter for 1D Euler equations.
//
// This is based off a Harten & Zwas flux limiter that considers linear
// combinations of low- and high-order fluxes.
//
// 
void ApplyPosMPPLimiter( 
        const double dt, const int method[], const dTensor2& node,
        const dTensorBC1& smax,
        const dTensorBC3& aux, const dTensorBC3& q, 
        dTensorBC2& Fm, dTensorBC2& Fp )
//  const dTensorBC2& fLF, dTensorBC2& fHat )
{


    // Parameters for the current grid
    const int     mx = q.getsize(1);
    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    const double dx      = dogParamsCart1.get_dx();
    const double eps     = 1.0e-13;
    const double dtx     = dt/dx;

    const double rho_min = eps;
    const double gamma   = eulerParams.gamma;
    const double gm1     = gamma - 1.0;

    // Storage
    double gmin[mx],amin[mx][2],thex[mx+1],ff[2],aa[2];
    double qh[meqn],qtmp[meqn],qlf[meqn],fhat_local[2][meqn],flf_local[2][meqn];

    // local variables
    double rho,u1,u2,u3,energy,plow,phigh,ffsum,ftmp;
    double ac,bc,cc,delta,root1,root2,rate;
    int isgn[2];

    // Lax-Friedrich's fluxes
    dTensorBC2 FmLF(melems, meqn, mbc );
    dTensorBC2 FpLF(melems, meqn, mbc );

    // Limiting parameter
    dTensorBC2 Theta(melems, meqn, mbc );

    // Range for limiting parameter, Lambda_{+1/2}, Lambda_{-1/2}.
    dTensorBC2 LambdaP(melems, meqn, mbc );
    dTensorBC2 LambdaM(melems, meqn, mbc );

    // Construct the flux for a Lax-Friedrichs solver
    void ConstructL_LLF(const int method[],
            const dTensor2& node,
            const dTensorBC1& smax,
            const dTensorBC3& aux,
            const dTensorBC3& q,
            dTensorBC2& Fm,
            dTensorBC2& Fp);
    ConstructL_LLF(method, node, smax, aux, q, FmLF, FpLF);
    
    // limiting on rho
    for (int i=1; i<=mx; i++)
        gmin[i-1] = rho_min-( q.get(i,1,1) + dtx*( FpLF.get(i-1,1) - FpLF.get(i,1)) );

    for (int i=1; i<=mx; i++)
    {
        ff[0] = dtx*(Fm.get(i,1) - FmLF.get(i,1)  );
        ff[1] =-dtx*(Fp.get(i,1) - FpLF.get(i,1));

        for (int k=0; k<=1; k++)
        {
            if (ff[k]<0.0)
                isgn[k] = 1;
            else
                isgn[k] = 0;
        }

        ffsum = isgn[0]*ff[0]+isgn[1]*ff[1];

        for (int k=0; k<=1; k++)
        {
            if (isgn[k]==1)
                amin[i-1][k] = Min(gmin[i-1]/(ffsum-eps),1.0);
            else
                amin[i-1][k] = 1.0;
        }
    }

    // limiting on pressure
    for (int i=1; i<=mx; i++)
    {
        for( int m=1; m <= meqn; m++ )  
            qlf[m-1]=q.get(i,m,1) + dtx*(FpLF.get(i-1,m)-FmLF.get(i,m));

        rho    = qlf[0];
        u1     = qlf[1]/rho;
        u2     = qlf[2]/rho;
        u3     = qlf[3]/rho;
        energy = qlf[4];
        plow   = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

        if (plow < 0.0)
            cout << "Negative solution in Lax-Fridrichs flux" << endl;

        for( int m=1; m <= meqn; m++ )
        {
            fhat_local[0][m-1] = Fm.get(i,m);
            fhat_local[1][m-1] = Fp.get(i,m);
            flf_local[0][m-1]  =  FmLF.get(i,m);
            flf_local[1][m-1]  =  FpLF.get(i,m);
            qtmp[m-1] = q.get(i,m,1);
        }

        plow = min(eps,max(eps,plow));

        double rescale[2][2] = {{1.0,1.0},{1.0,1.0}};

        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        {
            aa[0] = i1*amin[i-1][0];
            aa[1] = i2*amin[i-1][1];

            for( int m=1; m <= meqn; m++ )
            {
                ff[0] = aa[0]*fhat_local[0][m-1]+(1.0-aa[0])*flf_local[0][m-1];
                ff[1] = aa[1]*fhat_local[1][m-1]+(1.0-aa[1])*flf_local[1][m-1];
                qh[m-1] = qtmp[m-1] + dtx*(ff[0]-ff[1]);
            }

            rho    = qh[0];
            u1     = qh[1]/rho;
            u2     = qh[2]/rho;
            u3     = qh[3]/rho;
            energy = qh[4];
            phigh  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

            if (phigh < 0.0)
            {
                ac = (qh[1]-qlf[1])*(qh[5]-qlf[5])-0.5*( pow(qh[2]-qlf[2],2)+pow(qh[3]-qlf[3],2)+pow(qh[4]-qlf[4],2) );
                bc = (qh[5]-qlf[5])*qlf[1]+(qh[1]-qlf[1])*qlf[5]-(qh[2]-qlf[2])*qlf[2]-(qh[3]-qlf[3])*qlf[3]-(qh[4]-qlf[4])*qlf[4]-plow/gm1*(qh[1]-qlf[1]);
                cc = qlf[1]*qlf[5]-0.5*(pow(qlf[2],2)+pow(qlf[3],2)+pow(qlf[4],2)) - plow/gm1*qlf[1];

                if (fabs(ac) >= eps)
                {
                    delta = sqrt(fabs(bc*bc-4.0*ac*cc));
                    root1 = (-bc-delta)/(2.0*ac);
                    root2 = (-bc+delta)/(2.0*ac);
                    if (ac > 0.0) {
                    rate = 1.0;
                    if (root1 >=0.0)
                        rate = Min(root1,rate);
                    if (root2 >=0.0)
                        rate = Min(root2,rate); }
                    else {
                        rate = 0.0;
                        if (root1 <=1.0)
                            rate = Max(rate,root1);
                        if (root2 <=1.0)
                            rate = Max(rate,root2); }
                }
                else 
                {
                    rate = 0.0;
                    if ( fabs(bc)>0.0 )
                    {
                        root1 = -cc/bc; 
                        if ((root1>=0.0) && (root1<=1.0))
                            rate = root1;
                    }
                }
                rescale[i1][i2] = rate;
            }

        }   //end of for i1,i2 loop
        
        double rescale2[2] = {1.0, 1.0};
        for(int i1=0; i1<2; i1++){
            rescale2[0] = Min(rescale2[0],rescale[1][i1]);
            rescale2[1] = Min(rescale2[1],rescale[i1][1]);
        }
            
        for(int i1=0; i1<2; i1++)
            amin[i-1][i1] = rescale2[i1]*amin[i-1][i1];    
    }


    for (int i=1; i<=mx-1; i++)
    {
        thex[i] = Min(amin[i][0],amin[i-1][1]);
    }
    thex[0] = amin[0][0];
    thex[mx]= amin[mx-1][1];

    for (int i=1; i<=mx+1; i++)
        for( int m=1; m <= meqn; m++ )
        {
            ftmp = thex[i-1]*(Fm.get(i,m)-FmLF.get(i,m))+FmLF.get(i,m);
            Fm.set(i,m,ftmp);
            Fp.set(i-1,m,ftmp);
        }


}
