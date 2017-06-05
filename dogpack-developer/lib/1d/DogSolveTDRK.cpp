// Two-Derivative Runge-Kutta Solver.
//
// This solver is a hybrid between the full Lax-Wendroff method and a
// Runge-Kutta solver.
//
// See the fourth order scheme from
// http://www.springerlink.com/content/412506004227h393/
//
// "On Two Derivative Runge-Kutta Methods."
//
// This method is fourth order accurate, and uses two stage values, and two
// derivatives of q.  Hence, any application wishing to use this scheme need
// to implement a Dflux, but not a D2Flux as is needed for the 3rd order
// Lax-Wendroff scheme.

#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "RKinfo.h"

using namespace std;

void DogSolveUser(const dTensor2& node, const dTensor1& prim_vol, dTensorBC3& aux, 
        dTensorBC3& qold, dTensorBC3& qnew, dTensorBC1& smax,
        double tstart, double tend,int nv, const int method[],
        double dtv[], const double cflv[],string outputdir)
{

    // ------------------------------------------------------------
    // Function definitions
    void CopyQ(const dTensorBC3&,dTensorBC3&);
    void ConSoln(const int[],const dTensor2&,const dTensorBC3&,const dTensorBC3&,
            double,string);
    void BeforeStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
    void AfterStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
    void ApplyLimiter(const dTensor2&,const dTensorBC3&,dTensorBC3&,
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,
                const dTensor2&,dTensor2&),
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,
                const dTensor2&,dTensor2&));
    void ProjectRightEig(const dTensor1&,const dTensor1&,
            const dTensor2&,dTensor2&);
    void ProjectLeftEig(const dTensor1&,const dTensor1&,
            const dTensor2&,dTensor2&);
    void UpdateSoln(double,double,double,double,const dTensor2&,const dTensorBC3&,
            const dTensorBC3&, const dTensorBC3&,dTensorBC3&);
    double GetCFL(double,double,const dTensor1&,const int[],
            const dTensorBC3&,const dTensorBC1&);
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
    void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
            dTensorBC3& auxstar, dTensorBC3& aux, 
            dTensorBC3& qold, dTensorBC3& q);
    void LaxWendroffTD(double dt, const int method[], const dTensor2& node,
            double alpha1, double beta1,
            dTensorBC3& aux1, dTensorBC3& q1, 
            double alpha2, double beta2,
            dTensorBC3& aux2, dTensorBC3& q2, 
            dTensorBC3& Lstar,dTensorBC1& smax);
    void EulerStep( double dt, const dTensorBC3& qold, 
        const dTensorBC3& Lstar, dTensorBC3& qnew );
    // ------------------------------------------------------------

    // define local variables
    double t  = tstart;
    double dt = dtv[1];
    const double CFL_max    = cflv[1];
    const double CFL_target = cflv[2];
    double dtmin = dt;
    double dtmax = dt;

    // grid and solution parameters:
    const int melems = qold.getsize(1);
    const int meqn   = qold.getsize(2);
    const int kmax   = qold.getsize(3);
    const int maux   =  aux.getsize(2);
    const int mbc    = qold.getmbc();

    // Intermediate Storage
    dTensorBC3   qstar(melems,meqn,method[1],mbc);
    dTensorBC3 auxstar(melems,maux,method[1],mbc);
    dTensorBC3   Lstar(melems,meqn,method[1],mbc);    

    // Set initialize qstar and auxstar values
    CopyQ(qold,qstar);
    CopyQ(aux,auxstar);

    // ------------------------- //
    // Main time stepping loop   //
    // ------------------------- //

    int n_step = 0;
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveLxWTD.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }    

        // Save old value in the case of a restart
        CopyQ(qnew,qold);

        // The integrator will fail if rejecting first step without this here:
        CopyQ(qnew,qstar);  

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // --------------------------------------------- //
            // ---- Take a Full time step of length dt ----- //
            // --------------------------------------------- //

            BeforeStep(dt,node,aux,qnew);

            // -- First Stage -- //

            // Construct Right hand side for a single time step of length
            // dt/2:
            LaxWendroffTD(0.5*dt, method, node, 
                    1.0e0, 0.5e0, aux, qnew, 
                    0.0, 0.0, auxstar, qstar, 
                    Lstar, smax );

            // qstar = qold + (0.5*dt)*Lstar:
            EulerStep( 0.5*dt, qold, Lstar, qstar );

            if( dogParams.using_moment_limiter() )
            {  ApplyLimiter( node, auxstar, qstar, &ProjectRightEig, &ProjectLeftEig);  }

            AfterStep(dt, node, aux, qstar);

            // -- Second Stage -- //

            LaxWendroffTD(dt, method, node, 
                    1.0e0, (1.0/6.0), aux, qnew, 
                    0.0,   (1.0/3.0), auxstar, qstar, 
                    Lstar, smax );

            // qnew = qold + (dt)*Lstar:
            EulerStep( dt, qold, Lstar, qnew );

            if( dogParams.using_moment_limiter() )
            {  ApplyLimiter( node, aux, qnew, &ProjectRightEig, &ProjectLeftEig );  }

            AfterStep(dt, node, aux, qnew);

            // --------------------------------------------------- //
            // ------- Finished taking a full time step ---------- //
            // --------------------------------------------------- //

            // compute cfl number
            double cfl = GetCFL(dt,dtv[2],prim_vol,method,aux,smax);

            // output time step information
            if (method[4]>0) 
            {
                cout << setprecision(3);
                cout << "DogSolve1D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max) // accept
            { m_accept = 1; }
            else 
            {   
                //reject
                t = told;
                if (method[4]>0)
                {
                    cout<<"DogSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                CopyQ(qold,qnew);
            }

        }

        // compute conservation and print to file
        ConSoln(method,node,aux,qnew,t,outputdir);
    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

}

void EulerStep( 
    double dt, 
    const dTensorBC3& qold, 
    const dTensorBC3& Lstar,
    dTensorBC3& qnew )
{

    // grid and solution parameters:
    const int melems = qold.getsize(1);
    const int meqn   = qold.getsize(2);
    const int kmax   = qold.getsize(3);
    const int mbc    = qold.getmbc();

#pragma omp parallel for
    for( int i=1-mbc; i <= (melems+mbc); i++ )
    for( int me=1; me <= meqn; me++ )
    for( int k=1; k <= kmax; k++ )
    { qnew.set(i, me, k, qold.get(i, me, k) + dt*Lstar.get(i, me, k) ); }

}
