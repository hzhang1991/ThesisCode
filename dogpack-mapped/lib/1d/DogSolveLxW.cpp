#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "RKinfo.h"

// If we want to use DogSolver from the top-level library, this needs to be
// written:
// #include "DogState1d.h"  

using namespace std;

void DogSolveLxW(const dTensor2& node, const dTensor1& prim_vol, dTensorBC3& aux, 
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
    void ConstructL(const int[],const dTensor2&,const dTensorBC3&,
            const dTensorBC3&,dTensorBC3&,dTensorBC1&);
    void UpdateSoln(double,double,double,double,const dTensor2&,const dTensorBC3&,
            const dTensorBC3&, const dTensorBC3&,dTensorBC3&);
    double GetCFL(double,double,const dTensor1&,const int[],
            const dTensorBC3&,const dTensorBC1&);
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
    void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
            dTensorBC3& auxstar, dTensorBC3& aux, 
            dTensorBC3& qold, dTensorBC3& q);
    void LaxWendroff(const double dt, const int method[], const dTensor2& node,
            dTensorBC3& aux, dTensorBC3& q, 
            dTensorBC3& Lstar,dTensorBC1& smax);
    // ------------------------------------------------------------


    // define local variables
    int j,n_step,m_accept,mtmp;
    double t,dt,CFL_max,CFL_target,dtmin,dtmax;
    double told,cfl;
    n_step = 0;
    t = tstart;
    dt = dtv[1];
    CFL_max = cflv[1];
    CFL_target = cflv[2];
    cfl = 0.0;
    dtmin = dt;
    dtmax = dt;

    const int melems = qold.getsize(1);
    const int meqn   = qold.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc = qnew.getmbc();

    dTensorBC3   qstar(melems,meqn,method[1],mbc);
    dTensorBC3 auxstar(melems,maux,method[1],mbc);
    dTensorBC3   Lstar(melems,meqn,method[1],mbc);    

    // Set initialize qstar and auxstar values
    CopyQ(qold,qstar);
    CopyQ(aux,auxstar);

    while (t<tend)
    {
        // initialize time step
        m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveLxW.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }    

        // copy qnew into qold
        CopyQ(qnew,qold);

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            for (j=1-mbc; j<=(melems+mbc); j++)
            { smax.set(j, 0.0e0 ); }

            // Take a full time step of size dt
            BeforeStep(dt,node,aux,qnew);

            // Construct right hand side function:  
            // (q(t+dt)-q(t)) / dt = Lstar
            LaxWendroff(dt,method,node,aux,qnew,Lstar,smax);

            // qnew = qold + dt*Lstar
            UpdateSoln(1.0e0,0.0e0,1.0e0,dt,node,aux,qnew,Lstar,qnew);
            if (dogParams.using_moment_limiter())
            {  ApplyLimiter( node,aux,qnew,&ProjectRightEig,&ProjectLeftEig);  }

            AfterStep(dt,node,aux,qnew);

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
            if (cfl<=CFL_max)     // accept
            { m_accept = 1; }
            else                  //reject
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
