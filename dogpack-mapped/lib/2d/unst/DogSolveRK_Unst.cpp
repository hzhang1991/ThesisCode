#include "DogSolveRK_Unst.h"

void DogSolveRK_Unst(const mesh& Mesh, const edge_data_Unst& EdgeData,
        dTensor3& aux, dTensor3& qold, dTensor3& qnew, 
        const double tstart, const double tend, 
        const string outputdir)
{

    const int meqn = qnew.getsize(2);
    const int kmax = qnew.getsize(3);
    const int maux = aux.getsize(2);
    const double* cflv     = dogParams.get_cflv();
    const int     nv       = dogParams.get_nv();

    RKinfo rk;
    void SetRKinfo(int time_order, RKinfo& rk);
    SetRKinfo(dogParams.get_time_order(),rk);

    // define local variables
    int n_step = 0;      // counter for the number of steps
    double t  = tstart;
    double dt = dogStateUnst2.get_initial_dt();
    double CFL_max    = cflv[1];
    double CFL_target = cflv[2];
    double cfl = 0.0;
    double dtmin = dt;
    double dtmax = dt;
    const int NumElems = Mesh.get_NumElems(); // Number of total elements in mesh
    const int NumNodes = Mesh.get_NumNodes(); // Number of nodes in mesh
    const int NumEdges = Mesh.get_NumEdges(); // Number of edges in mesh 

    // Storage
    dTensor3   qstar(NumElems,meqn,kmax);
    dTensor3      q1(NumElems,meqn,kmax);
    dTensor3      q2(NumElems,meqn,kmax);
    dTensor3 auxstar(NumElems,maux,kmax);
    dTensor3   Lstar(NumElems,meqn,kmax);
    dTensor3    Lold(NumElems,meqn,kmax);
    dTensor3  auxold(NumElems,maux,kmax);
    dTensor1    smax(NumEdges);

    // Set initialize qstar and auxstar values
    CopyQ_Unst(qold,qstar);
    CopyQ_Unst(aux,auxstar);

    // Runge-Kutta time stepping
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            eprintf(" Error in DogSolveRK_Unst.cpp: "
                    " Exceeded allowed # of time steps \n"
                    "    n_step = %d\n"
                    "        nv = %d\n\n",
                    n_step,nv);
        }

        // copy qnew into qold
        CopyQ_Unst(qnew,qold);
        CopyQ_Unst(aux,auxold);

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;
            dogStateUnst2.set_time(told);
            dogStateUnst2.set_dt(dt);

            // Set initial maximum wave speed to zero
            for (int i=1; i<=NumEdges; i++)
            {  smax.set(i, 0.0e0 );  }

            // Take a full time step of size dt
            switch ( dogParams.get_time_order() )
            {
                case 1:  // First order in time (1-stage)

                    // -----------------------------------------------
                    // Stage #1 (the only one in this case)  
                    rk.mstage = 1;
                    BeforeStep_Unst(dt,Mesh,aux,qnew);
                    ConstructL_Unst(Mesh,EdgeData,aux,qnew,Lstar,smax);
                    CopyQ_Unst(Lstar,Lold);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qnew,Lstar,qnew);
                    AfterStep_Unst(dt,Mesh,aux,qnew);
                    // -----------------------------------------------
                    break;

                case 2:  // Second order in time (2-stages)

                    // -----------------------------------------------
                    // Stage #1  	        
                    rk.mstage = 1;
                    //dogStateUnst2.set_time(told);
                    dogStateUnst2.set_time(told+dt);
                    BeforeStep_Unst(dt,Mesh,aux,qnew);
                    ConstructL_Unst(Mesh,EdgeData,aux,qnew,Lstar,smax);
                    CopyQ_Unst(Lstar,Lold);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qnew,Lstar,qstar);
                    //if (dogParams.get_use_limiter()==1)
                    //{  ApplyLimiter_Unst(Mesh,aux,qstar,
                    //		&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep_Unst(dt,Mesh,auxstar,qstar);
                    // ------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    //dogStateUnst2.set_time(told+dt);
                    dogStateUnst2.set_time(told+dt);
                    BeforeStep_Unst(dt,Mesh,auxstar,qstar);
                    ConstructL_Unst(Mesh,EdgeData,aux,qstar,Lstar,smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,auxstar,qstar,Lstar,qnew);
                    //if (dogParams.get_use_limiter()==1)
                    //{  ApplyLimiter_Unst(Mesh,auxstar,qnew,
                    //		&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep_Unst(dt,Mesh,aux,qnew);
                    // ------------------------------------------------
                    break;

                case 3:  // Third order in time (3-stages)

                    // ------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    //dogStateUnst2.set_time(told);
                    dogStateUnst2.set_time(told+dt);
                    BeforeStep_Unst(dt,Mesh,aux,qnew);	      
                    ConstructL_Unst(Mesh,EdgeData,aux,qnew,Lstar,smax);
                    CopyQ_Unst(Lstar,Lold);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qnew,Lstar,qstar);
                    //if (dogParams.get_use_limiter()==1)
                    //	{  ApplyLimiter_Unst(Mesh,aux,qstar,
                    //		&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep_Unst(dt,Mesh,aux,qstar);
                    // -------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    //dogStateUnst2.set_time(told+dt);
                    dogStateUnst2.set_time(told+0.5*dt);
                    BeforeStep_Unst(dt,Mesh,aux,qstar);
                    ConstructL_Unst(Mesh,EdgeData,aux,qstar,Lstar,smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qnew,Lstar,qstar);   
                    //if (dogParams.get_use_limiter()==1)
                    //	{  ApplyLimiter_Unst(Mesh,auxstar,qstar,
                    //		&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep_Unst(dt,Mesh,aux,qstar);
                    // --------------------------------------------------
                    // Stage #3
                    rk.mstage = 3;
                    //dogStateUnst2.set_time(told+0.5*dt);
                    dogStateUnst2.set_time(told+dt);
                    BeforeStep_Unst(dt,Mesh,auxstar,qstar);
                    ConstructL_Unst(Mesh,EdgeData,auxstar,qstar,Lstar,smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qstar,Lstar,qnew);   
                    //if (dogParams.get_use_limiter()==1)
                    //	{  ApplyLimiter_Unst(Mesh,auxstar,qnew,
                    //		&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep_Unst(dt,Mesh,aux,qnew);
                    // --------------------------------------------------   

                    break;

                case 4:  // Fourth order in time (10-stages)

                    // -----------------------------------------------
                    CopyQ_Unst(qnew,q1);
                    CopyQ_Unst(q1,q2);

                    // Stage: 1,2,3,4, and 5
                    for (int s=1; s<=5; s++)
                    {
                        rk.mstage = s;
                        BeforeStep_Unst(dt,Mesh,aux,q1);
                        ConstructL_Unst(Mesh,EdgeData,aux,q1,Lstar,smax);
                        if (s==1)
                        {  CopyQ_Unst(Lstar,Lold);  }
                        UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,Mesh,aux,q1,Lstar,q1);
                        //if (dogParams.get_use_limiter()==1)
                        //  {  ApplyLimiter_Unst(Mesh,aux,q1,
                        //		    &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep_Unst(dt,Mesh,aux,q1);
                    }

                    // Temporary storage
                    for (int i=1; i<=NumElems; i++)
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
                        rk.mstage = s;
                        BeforeStep_Unst(dt,Mesh,aux,q1);
                        ConstructL_Unst(Mesh,EdgeData,aux,q1,Lstar,smax);
                        UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,Mesh,aux,q1,Lstar,q1);
                        //if (dogParams.get_use_limiter()==1)
                        // {  ApplyLimiter_Unst(Mesh,aux,q1,
                        //		    &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep_Unst(dt,Mesh,aux,q1);
                    }

                    // Stage: 10
                    rk.mstage = 10;
                    BeforeStep_Unst(dt,Mesh,aux,q1);
                    ConstructL_Unst(Mesh,EdgeData,aux,q1,Lstar,smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,q2,Lstar,q1);
                    //if (dogParams.get_use_limiter()==1)
                    //{  ApplyLimiter_Unst(Mesh,aux,q1,
                    //		&ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep_Unst(dt,Mesh,aux,q1);

                    CopyQ_Unst(q1,qnew);
                    // -----------------------------------------------		  
                    break;

                case 5:  // Fifth order in time (8-stages)

                    // -----------------------------------------------
                    CopyQ_Unst(qnew,q1);
                    q2.setall(0.);

                    for (int s=1; s<=8; s++)
                    {
                        rk.mstage = s;
                        BeforeStep_Unst(dt,Mesh,aux,q1);
                        ConstructL_Unst(Mesh,EdgeData,aux,q1,Lstar,smax);
                        if (s==1)
                        {  CopyQ_Unst(Lstar,Lold);  }

                        rk.gamma->get(1,s); 
                        rk.gamma->get(2,s); 
                        rk.gamma->get(3,s);
                        rk.delta->get(s);
                        rk.beta->get(s);

                        UpdateSoln_Unst(rk.gamma->get(1,s), 
                                rk.gamma->get(2,s), 
                                rk.gamma->get(3,s), 
                                rk.delta->get(s), rk.beta->get(s),
                                dt, Mesh, aux, qold, Lstar, q1, q2);

                        //if (dogParams.using_moment_limiter())
                        //  {  ApplyLimiter_Unst(aux,q1,
                        //			 &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep_Unst(dt,Mesh,aux,q1);
                    }

                    CopyQ_Unst(q1,qnew);
                    // -----------------------------------------------          
                    break;

                default: unsupported_value_error(dogParams.get_time_order());

            }

            // compute cfl number
            cfl = GetCFL_Unst(dt,Mesh,aux,smax);

            // output time step information
            if (dogParams.get_verbosity()>0) 
            {
                printf("DogSolve2D ... Step %5d"
                        "   CFL =%6.3f"
                        "   dt =%11.3e"
                        "   t =%11.3e\n",
                        n_step,cfl,dt,t);
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dogParams.get_max_dt(),dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dogParams.get_max_dt();
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
                // accept
            { 
                m_accept = 1; 
                dogStateUnst2.set_time(t);

                // do any extra work
                AfterFullTimeStep_Unst(dogStateUnst2.get_dt(),Mesh,
                        auxold,qold,Lold,aux,qnew);
            }
            else 
                //reject
            {   
                t = told;
                dogStateUnst2.set_time(told);
                if (dogParams.get_verbosity()>0)
                {
                    printf("DogSolve2D rejecting step..."
                            "CFL number too large\n");
                }

                // copy qold into qnew
                CopyQ_Unst(qold,qnew);
                CopyQ_Unst(auxold,aux);

                // after reject function	      
                AfterReject_Unst(Mesh,dt,aux,qnew);
            }      
        }

        // compute conservation and print to file
        ConSoln_Unst(Mesh,aux,qnew,t,outputdir);

    }

    // set initial time step for next call to DogSolveRK
    dogStateUnst2.set_initial_dt(dt);

    void DeleteRKInfo(RKinfo& rk);
    DeleteRKInfo(rk);

}
