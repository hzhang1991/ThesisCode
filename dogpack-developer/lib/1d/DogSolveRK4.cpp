#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"

// If we want to use DogSolver from the top-level library, this needs to be
// written:
// #include "DogState1d.h"  

using namespace std;

// Classical Runge-Kutta 4 Time stepping
//
void DogSolveUser(
    const dTensor2& node, 
        const dTensor1& prim_vol, 
        dTensorBC3& aux, 
        dTensorBC3& qold,
        dTensorBC3& qnew,
        dTensorBC1& smax,
        double tstart, double tend,int nv, const int method[],
        double dtv[], const double cflv[],string outputdir)
{  
    void DogSolveRK4(
        const dTensor2& node, const dTensor1& prim_vol, 
            dTensorBC3& aux, dTensorBC3& qold, dTensorBC3& qnew,
            dTensorBC1& smax, 
            double tstart, double tend,int nv, const int method[],
            double dtv[], const double cflv[],string outputdir);

    DogSolveRK4( node, prim_vol, aux, qold, qnew, smax, tstart, tend,
        nv, method, dtv, cflv, outputdir );
}

void DogSolveRK4(
    const dTensor2& node, 
        const dTensor1& prim_vol, 
        dTensorBC3& aux, 
        dTensorBC3& qold,
        dTensorBC3& qnew,
        dTensorBC1& smax,
        double tstart, double tend,int nv, const int method[],
        double dtv[], const double cflv[],string outputdir)
{

    // ------------------------------------------------------------
    // Function definitions
    void CopyQ(const dTensorBC3&,dTensorBC3&);
    void ConSoln(const int[],const dTensor2&,const dTensorBC3&,const dTensorBC3&,double,string);
    void BeforeStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
    void AfterStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
    void ApplyLimiter(const dTensor2&,const dTensorBC3&,dTensorBC3&,
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    void ViscousLimiter(const dTensor2& node, const dTensorBC3& aux, 
            const dTensorBC3& qold, const dTensorBC3& q, dTensorBC3& Lstar,
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    void ProjectRightEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);
    void ProjectLeftEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);
    void ConstructL(const int method[], const dTensor2& node, 
                    dTensorBC3& aux, dTensorBC3& q, dTensorBC3& Lstar, dTensorBC1& smax);
    void UpdateSoln(double,double,double,double,const dTensor2&,const dTensorBC3&,
            const dTensorBC3&, const dTensorBC3&,dTensorBC3&);
    double GetCFL(double,double,const dTensor1&,const int[],const dTensorBC3&,const dTensorBC1&);
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
    void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
            dTensorBC3& auxstar, dTensorBC3& aux, 
            dTensorBC3& qold, dTensorBC3& q);
    double EulerStep( double t, double dt, 
        const dTensorBC3& qold, const dTensorBC3& Lstar, dTensorBC3& qnew );
    // ------------------------------------------------------------

    // define local variables
    double t  = tstart;
    double dt = dtv[1];
    const double CFL_max = cflv[1];
    const double CFL_target = cflv[2];
    double cfl = 0.0;
    double dtmin = dt;
    double dtmax = dt;
    const int melems = qold.getsize(1);
    const int meqn   = qold.getsize(2);
    const int maux   = aux.getsize(2);
    const int kmax   = qold.getsize(3);
    const int mbc = qnew.getmbc();

    // intermediate storage:
    dTensorBC3   qstar(melems,meqn,kmax,mbc);
    dTensorBC3 auxstar(melems,maux,kmax,mbc);

    // Runge-Kutta stages:
    dTensorBC3   k1(melems,meqn,kmax,mbc);    
    dTensorBC3   k2(melems,meqn,kmax,mbc);    
    dTensorBC3   k3(melems,meqn,kmax,mbc);    
    dTensorBC3   k4(melems,meqn,kmax,mbc);    

    // Set initialize qstar and auxstar values
    CopyQ(qold,qstar);
    CopyQ(aux,auxstar);

    int n_step = 0;  // counter for number of steps taken
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveUser.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold in case we need to reject a step
        CopyQ(qnew,qold);


        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            const double told = t;
            const double tn   = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            for (int j=1-mbc; j<=(melems+mbc); j++)
            { smax.set(j, 0.0e0 ); }

            // ------------------------------------------------------------- //
            //
            // Single RK4 Time step forced here:
            //
            // ------------------------------------------------------------- //

            // We'll use qstar as q^n (in case we reject steps):
            CopyQ(qold,qstar);

            // Stage 1:
            BeforeStep(0.5*dt,node,aux,qstar);
            ConstructL( method, node, aux, qstar, k1, smax );
            t = EulerStep( tn, 0.5*dt, qstar, k1, qnew );
            if (dogParams.using_moment_limiter())
            {  ApplyLimiter(node,aux,qnew,&ProjectRightEig,&ProjectLeftEig);  }
            AfterStep(0.5*dt,node,aux,qnew);

            // Stage 2:
            BeforeStep(0.5*dt,node,aux,qstar);
            ConstructL( method, node, aux, qnew, k2, smax );
            t = EulerStep( tn, 0.5*dt, qstar, k2, qnew );
            if (dogParams.using_moment_limiter())
            {  ApplyLimiter(node,aux,qnew,&ProjectRightEig,&ProjectLeftEig);  }
            AfterStep(0.5*dt,node,aux,qnew);

            // Stage 3:
            BeforeStep(dt,node,aux,qstar);
            ConstructL( method, node, aux, qnew, k3, smax );
            t = EulerStep( tn, dt, qstar, k3, qnew );
            if (dogParams.using_moment_limiter())
            {  ApplyLimiter(node,aux,qnew,&ProjectRightEig,&ProjectLeftEig);  }
            AfterStep(dt,node,aux,qnew);

            // Stage 4:
            BeforeStep(dt,node,aux,qnew);
            ConstructL( method, node, aux, qnew, k4, smax );

            // Construct Full Update ...
#pragma omp parallel for
            for( int i=1; i <= melems; i++ )
            for( int m=1; m <= meqn; m++ )
            for( int k=1; k <= kmax; k++ )
            {
                double tmp = 0.;
                tmp +=     k1.get(i,m,k);
                tmp += 2.0*k2.get(i,m,k);
                tmp += 2.0*k3.get(i,m,k);
                tmp +=     k4.get(i,m,k);
                qnew.set(i,m,k, qstar.get(i,m,k) + dt*tmp/6.0 );
            }
            if (dogParams.using_moment_limiter())
            {  ApplyLimiter(node,aux,qnew,&ProjectRightEig,&ProjectLeftEig);  }

            // ------------------------------------------------------------- //
            // Finished taking a single RK4 time step
            // ------------------------------------------------------------- //

            // do any extra work      
            AfterFullTimeStep(dt,node,prim_vol,auxstar,aux,qold,qnew);

            // compute cfl number
            cfl = GetCFL(dt,dtv[2],prim_vol,method,aux,smax);

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
            if (cfl<=CFL_max)
                // accept
            { m_accept = 1; }
            else 
                //reject
            {   
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

// Single Eulerian Step
double EulerStep( double t, double dt, 
    const dTensorBC3& qold, const dTensorBC3& Lstar,
    dTensorBC3& qnew )
{

    const int mx      = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells

#pragma omp parallel for
    for( int i=1; i <= mx ; i++ )
    for( int m=1; m <= meqn; m++ )
    for( int k=1; k <= kmax; k++ )
    {
        double tmp = qold.get(i,m,k) + dt*Lstar.get(i,m,k);
        qnew.set(i,m,k, tmp);
    }

    return t+dt;

}
