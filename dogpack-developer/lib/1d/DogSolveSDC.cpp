#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "RKinfo.h"
#include "DogState.h"

// If we want to use DogSolver from the top-level library, this needs to be
// written:
// #include "DogState1d.h"  

using namespace std;

void DogSolveSDC(const dTensor2& node, const dTensor1& prim_vol, dTensorBC3& aux, 
        dTensorBC3& qold, dTensorBC3& qnew, dTensorBC1& smax, 
        double tstart, double tend, int nv, const int method[], 
        double dtv[], const double cflv[], string outputdir)
{
    // ------------------------------------------------------------
    // Function definitions
    void CopyQ(const dTensorBC3&,dTensorBC3&);
    void ConSoln(const int[],const dTensor2&,const dTensorBC3&,const dTensorBC3&,double,string);
    void ProjectRightEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);
    void ProjectLeftEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);
    void ConstructL(const int[],const dTensor2&,dTensorBC3&,dTensorBC3&,
            dTensorBC3&,dTensorBC1&);
    void EulerStepSDC(const double& dt, const dTensorBC3& aux, 
            const dTensorBC3& qold, const dTensorBC3& Lrhs, dTensorBC3& qnew);
    void ResInt(double,const dTensorBC3&,const dTensorBC3&,const dTensorBC3&,
            const dTensorBC3&,const dTensorBC3&,const dTensorBC3&,dTensorBC4&);
    double GetCFL(double,double,const dTensor1&,const int[],const dTensorBC3&,const dTensorBC1&);
    void TimeStepSDC(int,double,double,dTensor1&,dTensor1&);
    void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
            dTensorBC3& auxold, dTensorBC3& aux, 
            dTensorBC3& qold, dTensorBC3& q);
    void RelaxLimiter(const dTensor2& node,dTensorBC3& aux,dTensorBC3& q);
    void ApplyLimiter(const dTensor2&,const dTensorBC3&,dTensorBC3&,
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    void BeforeStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
    void AfterStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);

    // This is a bit cleaner than putting all the pre and post steps in here
    void StepSDC(const double& dt, const int method[], const dTensor2& node,
            dTensorBC1& smax, dTensorBC3& Lrhs,
            dTensorBC3& aux, dTensorBC3& qin, 
            dTensorBC3& qnew);
    void StepSDCRK2(const double& dt, const int method[], const dTensor2& node,
            dTensorBC1& smax, dTensorBC3& Lrhs, dTensorBC3& Lstar,
            dTensorBC3& aux, dTensorBC3& qin, dTensorBC3& qstar,
            dTensorBC3& qnew);
    void StepSDCdeltaRK2(const double& dt, const int method[], const dTensor2& node,
            dTensorBC1& smax, dTensorBC3& aux, dTensorBC3& qstar, dTensorBC3& Lstar,
            dTensorBC3& L1, dTensorBC3& L1new, dTensorBC3& L2, dTensorBC3& q1, 
            dTensorBC3& q2, int num, dTensorBC4& IL);

    // ------------------------------------------------------------

    // Define local variables
    int i,j,k,m,N,n_step,m_accept,mtmp;
    double t,dt,CFL_max,CFL_target,dtmin,dtmax;
    double told,cfl,tmp,err1,err2,err3,err4,err5;
    n_step = 0;
    t = tstart;
    dt = dtv[1];
    CFL_max = cflv[1];
    CFL_target = cflv[2];
    cfl = 0.0;
    dtmin = dt;
    dtmax = dt;

    // problem dimensions and parameters
    const int melems = qnew.getsize(1);
    const int meqn   = qnew.getsize(2);
    const int kmax   = qnew.getsize(3);
    const int mbc    = qnew.getmbc();
    const int maux   = aux.getsize(2);    

    // --------------------------------------------------------------
    // Create helper arrays
    dTensor1   dtvec(method[2]-1);
    dTensor1    tvec(method[2]);
    dTensorBC3    q1(melems,meqn,method[1],mbc);
    dTensorBC3    q2(melems,meqn,method[1],mbc);
    dTensorBC3    q3(melems,meqn,method[1],mbc);
    dTensorBC3    q4(melems,meqn,method[1],mbc);
    dTensorBC3    q5(melems,meqn,method[1],mbc);
    dTensorBC3    qs(melems,meqn,method[1],mbc);

    dTensorBC3    L0(melems,meqn,method[1],mbc);
    dTensorBC3    L1(melems,meqn,method[1],mbc);
    dTensorBC3    L2(melems,meqn,method[1],mbc);
    dTensorBC3    L3(melems,meqn,method[1],mbc);
    dTensorBC3    L4(melems,meqn,method[1],mbc);
    dTensorBC3    L5(melems,meqn,method[1],mbc);
    dTensorBC3    Ls(melems,meqn,method[1],mbc);

    dTensorBC3 L1new(melems,meqn,method[1],mbc);
    dTensorBC3 L2new(melems,meqn,method[1],mbc);
    dTensorBC3 L3new(melems,meqn,method[1],mbc);
    dTensorBC3 L4new(melems,meqn,method[1],mbc);

    dTensorBC4 IL(melems,meqn,method[1],method[2]-1,mbc);
    // --------------------------------------------------------------

    // -----------------------
    // Construct initial L
    // -----------------------
    if (method[2]>2)
    { 
        // create current time step vector
        TimeStepSDC(method[2],told,dt,dtvec,tvec);

        // Construct initial right-hand side for SDC
        BeforeStep(dtvec.get(1),node,aux,qnew);
        ConstructL(method,node,aux,qnew,L0,smax);
    }


    // -----------------------
    // Main time-stepping loop
    // -----------------------
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

        // Copy qnew into qold (in order to save data)
        CopyQ(qnew,qold);

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // create current time step vector
            TimeStepSDC(method[2],told,dt,dtvec,tvec);

            // Set initial maximum wave speed to zero
            for (j=1-mbc; j<=(melems+mbc); j++)
            { smax.set(j, 0.0e0 ); }        

            // Choose the SDC order of accuracy
            switch ( method[2] )
            {
                case 2:  // 2nd order in time (this is identical to RK2)

                    // --------------------------------------------------------
                    CopyQ( qnew, q1 );
                    StepSDCRK2(dt, method, node, smax, L0, Ls, aux, q1, qs, qnew);
                    // --------------------------------------------------------

                case 3:  // 3rd order in time

                    // --------------------------------------------------------
                    // Take 1st Euler time step
                    EulerStepSDC(dtvec.get(1),aux,qnew,L0,q1);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q1);  }
                    AfterStep(dtvec.get(1),node,aux,q1);

                    // Take 2nd Euler time step
                    //dogState->set_dt(dtvec.get(2));
                    BeforeStep(dtvec.get(2),node,aux,q1);
                    ConstructL(method,node,aux,q1,L1,smax);
                    EulerStepSDC(dtvec.get(2),aux,q1,L1,q2);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q2);  }
                    AfterStep(dtvec.get(2),node,aux,q2);

                    // Construct new right-hand side
                    BeforeStep(dtvec.get(1),node,aux,q2);
                    ConstructL(method,node,aux,q2,L2,smax);

                    // Iterate to construct 
                    for (N=1; N<=(method[2]-1); N++)
                    {
                        // Integrate residual over time step
                        ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                        // First sub-interval    
                        //dogState->set_dt(dtvec.get(1));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {    
                                    // Get error in first sub-interval
                                    err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                                        + IL.get(i,m,k,1);

                                    // Correct solution
                                    q1.set(i,m,k, q1.get(i,m,k) + err1 );
                                }
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q1);  }
                        AfterStep(dtvec.get(1),node,aux,q1);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(2),node,aux,q1);
                        ConstructL(method,node,aux,q1,L1new,smax);

                        // Second sub-interval            
                        //dogState->set_dt(dtvec.get(2));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {
                                    // Get error in second sub-interval
                                    err2 = dtvec.get(2)*
                                        (L1new.get(i,m,k) - L1.get(i,m,k)) 
                                        - (q2.get(i,m,k)-q1.get(i,m,k)) 
                                        + IL.get(i,m,k,2);

                                    // Correct solution
                                    q2.set(i,m,k, q2.get(i,m,k) + err2 );
                                    L1.set(i,m,k, L1new.get(i,m,k) );
                                }        
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q2);  }
                        AfterStep(dtvec.get(2),node,aux,q2);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(1),node,aux,q2);
                        ConstructL(method,node,aux,q2,L2,smax);
                    }  

                    // Update solution
                    CopyQ(q2,qnew);
                    // --------------------------------------------------------

                    break;

                case 4:  // 4th order in time

                    // RK2 Time Steps on Q
                    StepSDCRK2(dtvec.get(1), method, node, smax, L0, Ls, aux, qnew, qs, q1);
                    StepSDCRK2(dtvec.get(2), method, node, smax, L1, Ls, aux, q1, qs, q2);
                    StepSDCRK2(dtvec.get(3), method, node, smax, L2, Ls, aux, q2, qs, q3);

                    // Construct new right-hand side at final time point
                    BeforeStep(dtvec.get(1),node,aux,q3);
                    ConstructL(method,node,aux,q3,L3,smax);

                    // Add in single correction, with RK2 time steps:
                    // Integrate residual over time step
                    ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                    //                  CopyQ(L0,L5);
                    //                  StepSDCdeltaRK2(dtvec.get(1), method, node, smax, aux,
                    //                      qs, Ls, L0, L5, L1, qnew, q1, 1, IL);

                    StepSDCdeltaRK2(dtvec.get(1), method, node, smax, aux,
                            qs, Ls, L0, L0, L1, qnew, q1, 1, IL);

                    // Second Interval
                    BeforeStep(dt, node, aux, q1);
                    ConstructL(method,node,aux,q1,L1new,smax);
                    StepSDCdeltaRK2(dtvec.get(2), method, node, smax, aux, 
                            qs, Ls, L1, L1new, L2, q1, q2, 2, IL );

                    // Third Interval
                    BeforeStep(dt, node, aux, q2);
                    ConstructL(method,node,aux,q2,L2new,smax);
                    StepSDCdeltaRK2(dtvec.get(3), method, node, smax, aux, 
                            qs, Ls, L2, L2new, L3, q2, q3, 3, IL );

                    // Update solution
                    CopyQ(q3,qnew);

                    break;


                    /****************************************************
                     ** Start of Euler Updates ***************************
                     *****************************************************
                     for(N=1; N<= 3; N++ )
                     {

                    // First sub-interval    
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                    for (m=1; m<=meqn; m++)
                    for (k=1; k<=method[1]; k++)
                    {    
                    // Get error in first sub-interval
                    err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                    + IL.get(i,m,k,1);

                    // Correct solution
                    q1.set(i,m,k, q1.get(i,m,k) + err1 );
                    }
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q1);  }
                    AfterStep(dtvec.get(1),node,aux,q1);

                    // Construct new right-hand side
                    BeforeStep(dtvec.get(2),node,aux,q1);
                    ConstructL(method,node,aux,q1,L1new,smax);

                    // Second sub-interval    
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                    for (m=1; m<=meqn; m++)
                    for (k=1; k<=method[1]; k++)
                    {
                    // Get error in second sub-interval
                    err2 = dtvec.get(2)*
                    (L1new.get(i,m,k) - L1.get(i,m,k)) 
                    - (q2.get(i,m,k)-q1.get(i,m,k)) 
                    + IL.get(i,m,k,2);

                    // Correct solution
                    q2.set(i,m,k, q2.get(i,m,k) + err2 );
                    L1.set(i,m,k, L1new.get(i,m,k) );
                    }        
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q2);  }
                    AfterStep(dtvec.get(2),node,aux,q2);

                    // Construct new right-hand side
                    BeforeStep(dtvec.get(3),node,aux,q2);
                    ConstructL(method,node,aux,q2,L2new,smax);

                    // Third sub-interval    
                    //dogState->set_dt(dtvec.get(3));
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                    for (m=1; m<=meqn; m++)
                    for (k=1; k<=method[1]; k++)
                    {
                    // Get error in third sub-interval
                    err3 = dtvec.get(3)*
                    (L2new.get(i,m,k) - L2.get(i,m,k)) 
                    - (q3.get(i,m,k)-q2.get(i,m,k)) 
                    + IL.get(i,m,k,3);

                    // Correct solution
                    q3.set(i,m,k, q3.get(i,m,k) + err3 );
                    L2.set(i,m,k, L2new.get(i,m,k) );
                    }    
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q3);  }
                    AfterStep(dtvec.get(3),node,aux,q3);

                    // Construct new right-hand side
                    BeforeStep(dtvec.get(4),node,aux,q3);
                    ConstructL(method,node,aux,q3,L3new,smax);

                    // Fourth sub-interval
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                        for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                                // Get error in fourth sub-interval
                                err4 = dtvec.get(4)*
                                    (L3new.get(i,m,k) - L3.get(i,m,k)) 
                                    - (q4.get(i,m,k)-q3.get(i,m,k)) 
                                    + IL.get(i,m,k,4);

                                // Correct solution
                                q4.set(i,m,k, q4.get(i,m,k) + err4 );
                                L3.set(i,m,k, L3new.get(i,m,k) );
                            }        
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q4,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q4);  }
                    AfterStep(dtvec.get(4),node,aux,q4);

                    // Construct new right-hand side
                    BeforeStep(dtvec.get(5),node,aux,q4);
                    ConstructL(method,node,aux,q4,L4new,smax);

                    // Fifth sub-interval
                    //dogState->set_dt(dtvec.get(5));
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                        for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                                // Get error in fourth sub-interval
                                err5 = dtvec.get(5)*
                                    (L4new.get(i,m,k) - L4.get(i,m,k)) 
                                    - (q5.get(i,m,k)-q4.get(i,m,k)) 
                                    + IL.get(i,m,k,5);

                                // Correct solution
                                q5.set(i,m,k, q5.get(i,m,k) + err5 );
                                L4.set(i,m,k, L4new.get(i,m,k) );
                            }        
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q5,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q5);  }
                    AfterStep(dtvec.get(5),node,aux,q5);

                    // Construct new right-hand side
                    BeforeStep(dtvec.get(1),node,aux,q5);
                    ConstructL(method,node,aux,q5,L5,smax);

                    // --------------------------------------------------------
                    // Take 1st Euler time step            
                    //dogState->set_dt(dtvec.get(1));
                    EulerStepSDC(dtvec.get(1),aux,qnew,L0,q1);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q1);  }
                    AfterStep(dtvec.get(1),node,aux,q1);

                    // Take 2nd Euler time step
                    //dogState->set_dt(dtvec.get(2));
                    BeforeStep(dtvec.get(2),node,aux,q1);
                    ConstructL(method,node,aux,q1,L1,smax);
                    EulerStepSDC(dtvec.get(2),aux,q1,L1,q2);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q2);  }
                    AfterStep(dtvec.get(2),node,aux,q2);

                    // Take 3rd Euler time step
                    //dogState->set_dt(dtvec.get(3));
                    BeforeStep(dtvec.get(3),node,aux,q2);
                    ConstructL(method,node,aux,q2,L2,smax);
                    EulerStepSDC(dtvec.get(3),aux,q2,L2,q3);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q3);  }
                    AfterStep(dtvec.get(3),node,aux,q3);

                    // Construct new right-hand side
                    BeforeStep(dtvec.get(1),node,aux,q3);
                    ConstructL(method,node,aux,q3,L3,smax);

                    // Iterate to construct 
                    for (N=1; N<=(method[2]-1); N++)
                    {
                        // Integrate residual over time step
                        ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                        // First sub-interval    
                        //dogState->set_dt(dtvec.get(1));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {    
                                    // Get error in first sub-interval
                                    err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                                        + IL.get(i,m,k,1);

                                    // Correct solution
                                    q1.set(i,m,k, q1.get(i,m,k) + err1 );
                                }
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q1);  }
                        AfterStep(dtvec.get(1),node,aux,q1);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(2),node,aux,q1);
                        ConstructL(method,node,aux,q1,L1new,smax);

                        // Second sub-interval    
                        //dogState->set_dt(dtvec.get(2));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {
                                    // Get error in second sub-interval
                                    err2 = dtvec.get(2)*
                                        (L1new.get(i,m,k) - L1.get(i,m,k)) 
                                        - (q2.get(i,m,k)-q1.get(i,m,k)) 
                                        + IL.get(i,m,k,2);

                                    // Correct solution
                                    q2.set(i,m,k, q2.get(i,m,k) + err2 );
                                    L1.set(i,m,k, L1new.get(i,m,k) );
                                }        
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q2);  }
                        AfterStep(dtvec.get(2),node,aux,q2);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(3),node,aux,q2);
                        ConstructL(method,node,aux,q2,L2new,smax);

                        // Third sub-interval    
                        //dogState->set_dt(dtvec.get(3));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {
                                    // Get error in third sub-interval
                                    err3 = dtvec.get(3)*
                                        (L2new.get(i,m,k) - L2.get(i,m,k)) 
                                        - (q3.get(i,m,k)-q2.get(i,m,k)) 
                                        + IL.get(i,m,k,3);

                                    // Correct solution
                                    q3.set(i,m,k, q3.get(i,m,k) + err3 );
                                    L2.set(i,m,k, L2new.get(i,m,k) );
                                }        
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q3);  }
                        AfterStep(dtvec.get(3),node,aux,q3);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(1),node,aux,q3);
                        ConstructL(method,node,aux,q3,L3,smax);
                    }  

                    // Update solution
                    CopyQ(q3,qnew);
                    // --------------------------------------------------------

                    break;

                    ***********************************
                        ******* End of Euler updates ******
                        **********************************/

                case 5:  // 5th order in time

                        // --------------------------------------------------------
                        // Take 1st Euler time step
                        //dogState->set_dt(dtvec.get(1));
                        EulerStepSDC(dtvec.get(1),aux,qnew,L0,q1);
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q1);  }
                        AfterStep(dtvec.get(1),node,aux,q1);

                        // Take 2nd Euler time step
                        StepSDC(dtvec.get(2), method, node, smax, L1, aux, q1, q2);

                        // Take 3rd Euler time step
                        StepSDC(dtvec.get(3), method, node, smax, L2, aux, q2, q3);

                        // Take 4th Euler time step
                        StepSDC(dtvec.get(4), method, node, smax, L3, aux, q3, q4);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(1),node,aux,q4);
                        ConstructL(method,node,aux,q4,L4,smax);

                        // Iterate to construct 
                        for (N=1; N<=(method[2]-1); N++)
                        {
                            // Integrate residual over time step
                            ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                            // First sub-interval    
                            //dogState->set_dt(dtvec.get(1));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {    
                                        // Get error in first sub-interval
                                        err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                                            + IL.get(i,m,k,1);

                                        // Correct solution
                                        q1.set(i,m,k, q1.get(i,m,k) + err1 );
                                    }
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q1);  }
                            AfterStep(dtvec.get(1),node,aux,q1);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(2),node,aux,q1);
                            ConstructL(method,node,aux,q1,L1new,smax);

                            // Second sub-interval    
                            //dogState->set_dt(dtvec.get(2));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in second sub-interval
                                        err2 = dtvec.get(2)*
                                            (L1new.get(i,m,k) - L1.get(i,m,k)) 
                                            - (q2.get(i,m,k)-q1.get(i,m,k)) 
                                            + IL.get(i,m,k,2);

                                        // Correct solution
                                        q2.set(i,m,k, q2.get(i,m,k) + err2 );
                                        L1.set(i,m,k, L1new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q2);  }
                            AfterStep(dtvec.get(2),node,aux,q2);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(3),node,aux,q2);
                            ConstructL(method,node,aux,q2,L2new,smax);

                            // Third sub-interval    
                            //dogState->set_dt(dtvec.get(3));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in third sub-interval
                                        err3 = dtvec.get(3)*
                                            (L2new.get(i,m,k) - L2.get(i,m,k)) 
                                            - (q3.get(i,m,k)-q2.get(i,m,k)) 
                                            + IL.get(i,m,k,3);

                                        // Correct solution
                                        q3.set(i,m,k, q3.get(i,m,k) + err3 );
                                        L2.set(i,m,k, L2new.get(i,m,k) );
                                    }    
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q3);  }
                            AfterStep(dtvec.get(3),node,aux,q3);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(4),node,aux,q3);
                            ConstructL(method,node,aux,q3,L3new,smax);

                            // Fourth sub-interval
                            //dogState->set_dt(dtvec.get(4));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in fourth sub-interval
                                        err4 = dtvec.get(4)*
                                            (L3new.get(i,m,k) - L3.get(i,m,k)) 
                                            - (q4.get(i,m,k)-q3.get(i,m,k)) 
                                            + IL.get(i,m,k,4);

                                        // Correct solution
                                        q4.set(i,m,k, q4.get(i,m,k) + err4 );
                                        L3.set(i,m,k, L3new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q4,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q4);  }
                            AfterStep(dtvec.get(4),node,aux,q4);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(1),node,aux,q4);
                            ConstructL(method,node,aux,q4,L4,smax);
                        }  

                        // Update solution
                        CopyQ(q4,qnew);
                        // --------------------------------------------------------

                        break;            

                case 6:  // 6th order in time

                        // RK2 Time Steps on Q
                        StepSDCRK2(dtvec.get(1), method, node, smax, L0, Ls, aux, qnew, qs, q1);
                        StepSDCRK2(dtvec.get(2), method, node, smax, L1, Ls, aux, q1, qs, q2);
                        StepSDCRK2(dtvec.get(3), method, node, smax, L2, Ls, aux, q2, qs, q3);
                        StepSDCRK2(dtvec.get(4), method, node, smax, L3, Ls, aux, q3, qs, q4);
                        StepSDCRK2(dtvec.get(5), method, node, smax, L4, Ls, aux, q4, qs, q5);

                        // Construct new right-hand side at final time point
                        BeforeStep(dtvec.get(1),node,aux,q5);
                        ConstructL(method,node,aux,q5,L5,smax);

                        /*
                           StepSDC(dtvec.get(1), method, node, smax, L0, aux, qnew, q1);

                        // Take 2nd Euler time step
                        StepSDC(dtvec.get(2), method, node, smax, L1, aux, q1, q2);

                        // Take 3rd Euler time step
                        StepSDC(dtvec.get(3), method, node, smax, L2, aux, q2, q3);

                        // Take 4th Euler time step
                        StepSDC(dtvec.get(4), method, node, smax, L3, aux, q3, q4);

                        // Take 5th Euler time step
                        StepSDC(dtvec.get(5), method, node, smax, L4, aux, q4, q5);

                        // Construct new right-hand side at final time point
                        BeforeStep(dtvec.get(1),node,aux,q5);
                        ConstructL(method,node,aux,q5,L5,smax);
                         */

                        for (N=1; N<=3; N++)
                        {
                            // Integrate residual over time step
                            ResInt(dt,L0,L1,L2,L3,L4,L5,IL);            

                            StepSDCdeltaRK2(dtvec.get(1), method, node, smax, aux,
                                    qs, Ls, L0, L0, L1, qnew, q1, 1, IL);

                            // Second Interval
                            BeforeStep(dt, node, aux, q1);
                            ConstructL(method,node,aux,q1,L1new,smax);
                            StepSDCdeltaRK2(dtvec.get(2), method, node, smax, aux, 
                                    qs, Ls, L1, L1new, L2, q1, q2, 2, IL );

                            // Third Interval
                            BeforeStep(dt, node, aux, q2);
                            ConstructL(method,node,aux,q2,L2new,smax);
                            StepSDCdeltaRK2(dtvec.get(3), method, node, smax, aux, 
                                    qs, Ls, L2, L2new, L3, q2, q3, 3, IL );

                            // Fourth Interval
                            BeforeStep(dt, node, aux, q3);
                            ConstructL(method,node,aux,q3,L3new,smax);
                            StepSDCdeltaRK2(dtvec.get(4), method, node, smax, aux, 
                                    qs, Ls, L3, L3new, L4, q3, q4, 4, IL );

                            // Final Interval
                            BeforeStep(dt, node, aux, q4);
                            ConstructL(method,node,aux,q4,L4new,smax);
                            StepSDCdeltaRK2(dtvec.get(5), method, node, smax, aux, 
                                    qs, Ls, L4, L4new, L5, q4, q5, 5, IL );

                            // Save the L1s for the next integration
                            CopyQ( L1new, L1 );
                            CopyQ( L2new, L2 );
                            CopyQ( L3new, L3 );
                            CopyQ( L4new, L4 );

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(1),node,aux,q5);
                            ConstructL(method,node,aux,q5,L5,smax);

                            /****************************************************
                             ** Start of Euler Updates ***************************
                             *****************************************************

                            // First sub-interval    
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {    
                            // Get error in first sub-interval
                            err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                            + IL.get(i,m,k,1);

                            // Correct solution
                            q1.set(i,m,k, q1.get(i,m,k) + err1 );
                            }
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q1);  }
                            AfterStep(dtvec.get(1),node,aux,q1);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(2),node,aux,q1);
                            ConstructL(method,node,aux,q1,L1new,smax);

                            // Second sub-interval    
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                            // Get error in second sub-interval
                            err2 = dtvec.get(2)*
                            (L1new.get(i,m,k) - L1.get(i,m,k)) 
                            - (q2.get(i,m,k)-q1.get(i,m,k)) 
                            + IL.get(i,m,k,2);

                            // Correct solution
                            q2.set(i,m,k, q2.get(i,m,k) + err2 );
                            L1.set(i,m,k, L1new.get(i,m,k) );
                            }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q2);  }
                            AfterStep(dtvec.get(2),node,aux,q2);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(3),node,aux,q2);
                            ConstructL(method,node,aux,q2,L2new,smax);

                            // Third sub-interval    
                            //dogState->set_dt(dtvec.get(3));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                            // Get error in third sub-interval
                            err3 = dtvec.get(3)*
                            (L2new.get(i,m,k) - L2.get(i,m,k)) 
                            - (q3.get(i,m,k)-q2.get(i,m,k)) 
                            + IL.get(i,m,k,3);

                            // Correct solution
                            q3.set(i,m,k, q3.get(i,m,k) + err3 );
                            L2.set(i,m,k, L2new.get(i,m,k) );
                            }    
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q3);  }
                            AfterStep(dtvec.get(3),node,aux,q3);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(4),node,aux,q3);
                            ConstructL(method,node,aux,q3,L3new,smax);

                            // Fourth sub-interval
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in fourth sub-interval
                                        err4 = dtvec.get(4)*
                                            (L3new.get(i,m,k) - L3.get(i,m,k)) 
                                            - (q4.get(i,m,k)-q3.get(i,m,k)) 
                                            + IL.get(i,m,k,4);

                                        // Correct solution
                                        q4.set(i,m,k, q4.get(i,m,k) + err4 );
                                        L3.set(i,m,k, L3new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q4,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q4);  }
                            AfterStep(dtvec.get(4),node,aux,q4);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(5),node,aux,q4);
                            ConstructL(method,node,aux,q4,L4new,smax);

                            // Fifth sub-interval
                            //dogState->set_dt(dtvec.get(5));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in fourth sub-interval
                                        err5 = dtvec.get(5)*
                                            (L4new.get(i,m,k) - L4.get(i,m,k)) 
                                            - (q5.get(i,m,k)-q4.get(i,m,k)) 
                                            + IL.get(i,m,k,5);

                                        // Correct solution
                                        q5.set(i,m,k, q5.get(i,m,k) + err5 );
                                        L4.set(i,m,k, L4new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q5,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q5);  }
                            AfterStep(dtvec.get(5),node,aux,q5);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(1),node,aux,q5);
                            ConstructL(method,node,aux,q5,L5,smax);

                            ***********************************
                                ******* End of Euler updates ******
                                **********************************/

                        }

                        // Update solution
                        CopyQ(q5,qnew);
                        // --------------------------------------------------------

                        break;            

            }

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

            // If anything needs to be done after the full time step
            AfterFullTimeStep(dt,node,prim_vol,aux,aux,qnew,qnew);

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
            {  
                m_accept = 1;  

                // Copy old L into new L0
                if (method[2]==3)
                {  CopyQ(L2,L0);  }

                if (method[2]==4)
                {  CopyQ(L3,L0);  }

                if (method[2]==5)
                {  CopyQ(L4,L0);  }        
            }
            else //reject
            {   
                t = told;
                if (method[4]>0)
                {
                    cout<<"DogSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // Copy old data back into qnew
                CopyQ(qold,qnew);
            }

        }

        // compute conservation and print to file
        ConSoln(method,node,aux,qnew,t,outputdir);

    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

}
