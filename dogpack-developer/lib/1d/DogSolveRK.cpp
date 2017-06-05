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

void DogSolveRK(const dTensor2& node, const dTensor1& prim_vol, dTensorBC3& aux, 
        dTensorBC3& qold, dTensorBC3& qnew, dTensorBC1& smax,
        double tstart, double tend,int nv, const int method[],
        double dtv[], const double cflv[],string outputdir)
{
    // ------------------------------------------------------------
    // Function definitions
    void CopyQ(const dTensorBC3&,dTensorBC3&);
    void ConSoln(const int[],const dTensor2&,const dTensorBC3&,
            const dTensorBC3&,double,string);
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
    void ConstructL(const int[],const dTensor2&,dTensorBC3&,
            dTensorBC3&,dTensorBC3&,dTensorBC1&);
    void UpdateSoln(double,double,double,double,const dTensor2&,const dTensorBC3&,
            const dTensorBC3&, const dTensorBC3&,dTensorBC3&);
    void UpdateSoln(
            double g1,double g2, double g3, double delta, double beta,double dt,
            const dTensor2& node,const dTensorBC3& aux,
            const dTensorBC3& qold, const dTensorBC3& Lstar,
            dTensorBC3& qstar, dTensorBC3& qnew);
    double GetCFL(double,double,const dTensor1&,const int[],
            const dTensorBC3&,const dTensorBC1&);
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
    void BeforeFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
            dTensorBC3& auxold, dTensorBC3& aux, 
            dTensorBC3& qold, dTensorBC3& q);
    void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
            dTensorBC3& auxstar, dTensorBC3& aux, 
            dTensorBC3& qold, dTensorBC3& q);
    void SetRKinfo(int method2, RKinfo& rk);
    void DeleteRKInfo(RKinfo& rk);
    // ------------------------------------------------------------

    // define local variables
    RKinfo rk;
    SetRKinfo(method[2],rk);
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

    const int mx   = qold.getsize(1);
    const int meqn = qold.getsize(2);
    const int kmax = qold.getsize(3);
    const int mbc  = qnew.getmbc();

    const int maux = aux.getsize(2);

    // Local storage
    dTensorBC3   qstar(mx,meqn,kmax,mbc);
    dTensorBC3 auxstar(mx,maux,kmax,mbc);
    dTensorBC3   Lstar(mx,meqn,kmax,mbc);
    dTensorBC3    Lold(mx,meqn,kmax,mbc);
    dTensorBC3      q1(mx,meqn,kmax,mbc);
    dTensorBC3      q2(mx,meqn,kmax,mbc);

    // Set initialize qstar and auxstar values
    // Note: it is assumed that qnew is the actual state vector.  Qold is used
    // to save this state in case of a rejected time step.
    qstar.copyfrom(qnew);
    auxstar.copyfrom(aux);

    while (t<tend)
    {
        // initialize time step
        m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveRK.cpp: "<< 
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
            for (j=1-mbc; j<=(mx+mbc); j++)
            { smax.set(j, 0.0e0 ); }

            // do any extra work
            BeforeFullTimeStep(dt,node,prim_vol,auxstar,aux,qold,qnew);

            // Take a full time step of size dt
            switch ( abs(method[2]) )
            {
                case 1:  // First order in time

                    // --------------------------------------------------------
                    // Stage #1 (the only one in this case)
                    dogParams.set_time(told);
                    rk.mstage = 1;
                    BeforeStep(dt,node,aux,qnew);
                    ConstructL(method,node,aux,qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,qnew,Lstar,qnew);
                    AfterStep(dt,node,aux,qnew);
                    // --------------------------------------------------------

                    break;

                case 2:  // Second order in time

                    // ---------------------------------------------------------
                    // Stage #1
                    dogParams.set_time(told);
                    rk.mstage = 1;
                    BeforeStep(dt,node,aux,qnew);
                    ConstructL(method,node,aux,qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,qnew,Lstar,qstar);      
                    if (dogParams.using_moment_limiter())
                    {  ApplyLimiter(node,aux,qstar,&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep(dt,node,auxstar,qstar);

                    // ---------------------------------------------------------
                    // Stage #2
                    dogParams.set_time(told+dt);
                    rk.mstage = 2;
                    BeforeStep(dt,node,auxstar,qstar);
                    ConstructL(method,node,auxstar,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,auxstar,qstar,Lstar,qnew);
                    if (dogParams.using_moment_limiter())
                    {  ApplyLimiter(node,auxstar,qnew,&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep(dt,node,aux,qnew); 
                    // ---------------------------------------------------------

                    break;

                case 3:  // Third order in time

                    // ---------------------------------------------------------
                    // Stage #1
                    dogParams.set_time(told);
                    rk.mstage = 1;
                    BeforeStep(dt,node,aux,qnew);    
                    ConstructL(method,node,aux,qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,qnew,Lstar,qstar);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,qstar,&ProjectRightEig,&ProjectLeftEig); }          
                    AfterStep(dt,node,auxstar,qstar);

                    // ---------------------------------------------------------
                    // Stage #2
                    dogParams.set_time(told+dt);
                    rk.mstage = 2;
                    BeforeStep(dt,node,auxstar,qstar);
                    ConstructL(method,node,auxstar,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,auxstar,qnew,Lstar,qstar);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,auxstar,qstar,&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep(dt,node,auxstar,qstar);

                    // ---------------------------------------------------------
                    // Stage #3
                    dogParams.set_time(told+0.5*dt);
                    rk.mstage = 3;
                    BeforeStep(dt,node,auxstar,qstar);
                    ConstructL(method,node,auxstar,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,auxstar,qstar,Lstar,qnew);   
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,auxstar,qnew,&ProjectRightEig,&ProjectLeftEig); }
                    AfterStep(dt,node,aux,qnew);
                    // ---------------------------------------------------------

                    break;

                case 4:  // Fourth order in time (10-stages)

                    // -----------------------------------------------
                    CopyQ(qnew,q1);
                    CopyQ(q1,q2);

                    // Stage: 1,2,3,4, and 5
                    for (int s=1; s<=5; s++)
                    {          
                        dogParams.set_time(told+0.5*onethird*dt*double(s-1));
                        rk.mstage = s;
                        BeforeStep(dt,node,aux,q1);
                        ConstructL(method,node,aux,q1,Lstar,smax);
                        if (s==1)
                        {  CopyQ(Lstar,Lold);  }
                        UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,node,aux,q1,Lstar,q1);
                        if (dogParams.using_moment_limiter())
                        {  ApplyLimiter(node,aux,q1,
                                &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep(dt,node,aux,q1);          
                    }

                    // Temporary storage
                    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {
                                double tmp = (q2.get(i,m,k) + 9.0*q1.get(i,m,k))/25.0;
                                q2.set(i,m,k, tmp );
                                q1.set(i,m,k, 15.0*tmp - 5.0*q1.get(i,m,k) );
                            }

                    // Stage: 6,7,8, and 9
                    for (int s=6; s<=9; s++)
                    {
                        dogParams.set_time(told+0.5*onethird*dt*double(s-4));
                        rk.mstage = s;
                        BeforeStep(dt,node,aux,q1);
                        ConstructL(method,node,aux,q1,Lstar,smax);
                        UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,node,aux,q1,Lstar,q1);
                        if (dogParams.using_moment_limiter())
                        {  ApplyLimiter(node,aux,q1,
                                &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep(dt,node,aux,q1);
                    }

                    // Stage: 10
                    dogParams.set_time(told+dt);
                    rk.mstage = 10;
                    BeforeStep(dt,node,aux,q1);
                    ConstructL(method,node,aux,q1,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,q2,Lstar,q1);
                    if (dogParams.using_moment_limiter())
                    {  ApplyLimiter(node,aux,q1,
                            &ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep(dt,node,aux,q1);

                    CopyQ(q1,qnew);
                    // -----------------------------------------------          
                    break;

                case 5:  // Fifth order in time (8-stages)

                    // -----------------------------------------------
                    CopyQ(qnew,q1);
                    q2.setall(0.);

                    for (int s=1; s<=8; s++)
                    {
                        rk.mstage = s;
                        BeforeStep(dt,node,aux,q1);
                        ConstructL(method,node,aux,q1,Lstar,smax);
                        if (s==1)
                        {  CopyQ(Lstar,Lold);  }

                        UpdateSoln(
                                rk.gamma->get(1,s), 
                                rk.gamma->get(2,s), 
                                rk.gamma->get(3,s), 
                                rk.delta->get(s), rk.beta->get(s),
                                dt, node, aux, qold, Lstar, q1, q2);

                        if (dogParams.using_moment_limiter())
                        {  ApplyLimiter(node,aux,q1,
                                &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep(dt,node,aux,q1);
                    }

                    CopyQ(q1,qnew);
                    // -----------------------------------------------          
                    break;

            }

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
            if (cfl<=CFL_max)       // accept
            { m_accept = 1; }
            else                    //reject
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

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    DeleteRKInfo(rk);

}
