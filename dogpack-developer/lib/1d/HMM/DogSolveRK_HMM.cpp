#include "../defs.h"
#include "dog_math.h"

void DogSolveRK_HMM(double tstart, double tend,int nv,double dtv[], double dtv_sub[],
		    double cflv[], iTensor1 map, iTensor2 leftBC, iTensor2 rightBC, 
		    int method[], dTensor2 node, dTensorBC3& aux, dTensorBC3& q, 
		    dTensorBC1& smax, dTensor1 prim_vol, int method_sub[], 
		    dTensor2 node_sub, dTensorBC3& aux_sub, dTensorBC3& q_sub, 
		    dTensorBC1& smax_sub, dTensor1 prim_vol_sub, string outputdir)
{  
    // ------------------------------------------------------------
    // Function definitions
    void CopyQ(dTensorBC3,dTensorBC3&);
    void ConSoln(dTensor2,int[],dTensorBC3,dTensorBC3,double,string);
    void AdvanceOneStage_RK(double,double,double,double,int[],dTensor2, 
			    dTensorBC3,dTensorBC3,dTensorBC3,dTensorBC3&,
			    dTensorBC3&,dTensorBC1&,iTensor1,iTensor2,
			    iTensor2,int);
    double GetCFL(double,double,dTensor1,int[],dTensorBC3,dTensorBC1);
    void GetHermiteCoeffs(int[],int,int,double,double,dTensorBC3,dTensorBC3,
			  dTensorBC3,dTensorBC3,dTensor2,dTensorBC4&);
    void SetBC_SubGrid(iTensor2,iTensor2,double,double,dTensor2,dTensorBC3,
		       dTensorBC4,dTensorBC3&);
    void SetValues_MainGrid(iTensor1,dTensor2,dTensorBC3,dTensorBC3,
			    dTensor2,dTensorBC3,dTensorBC3&);
    // ------------------------------------------------------------


    // define local variables
    int j,n_step,m_accept,mtmp,m,msubsteps;
    double t,dt,CFL_max,CFL_target,dtmin,dtmax,dt_sub,dt_store;
    double told,cfl,t_sub,cfl_sub,dt_old;
    n_step = 0;
    t = tstart;
    dt     = dtv[1];
    dt_sub = dtv_sub[1];
    CFL_max    = cflv[1];
    CFL_target = cflv[2];
    cfl = 0.0;
    dtmin = dt;
    dtmax = dt;
    int mx   = q.getsize(1);
    int meqn = q.getsize(2);
    int maux = aux.getsize(2);
    int mbc = q.getmbc();
    dTensorBC3    qold(mx,meqn,method[1],mbc);
    dTensorBC3  auxold(mx,maux,method[1],mbc);
    dTensorBC3   qstar(mx,meqn,method[1],mbc);
    dTensorBC3 auxstar(mx,maux,method[1],mbc);
    int mx_sub   = q_sub.getsize(1);
    int meqn_sub = q_sub.getsize(2);
    int maux_sub = aux_sub.getsize(2);
    int mbc_sub  = q_sub.getmbc();
    dTensorBC3    qold_sub(mx_sub,meqn_sub,method_sub[1],mbc_sub);
    dTensorBC3  auxold_sub(mx_sub,maux_sub,method_sub[1],mbc_sub);
    dTensorBC3   qstar_sub(mx_sub,meqn_sub,method_sub[1],mbc_sub);
    dTensorBC3 auxstar_sub(mx_sub,maux_sub,method_sub[1],mbc_sub);
    dTensorBC3   qtmp_sub(mx_sub,meqn_sub,method_sub[1],mbc_sub);
    dTensorBC3 auxtmp_sub(mx_sub,maux_sub,method_sub[1],mbc_sub);
    dTensorBC4 HermiteCoeffs(mx,meqn,method[1],4,mbc);

    // Main loop
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

        // copy q into qold and aux into auxold
        CopyQ(q,qold);
	CopyQ(q,qstar);
	CopyQ(aux,auxold);
	CopyQ(aux,auxstar);

	// -------------------------------
	//     MAIN GRID
	// -------------------------------

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

	    // Set initial maximum wave speed to zero
            for (j=1-mbc; j<=(mx_sub+mbc_sub); j++)
            { smax_sub.set(j, 0.0e0 ); }

            // Take a full time step of size dt
            switch ( abs(method[2]) )
            {
                case 1:  // First order in time

                    // --------------------------------------------------------
                    // MainGrid: Stage #1 (the only one in this case)
		    AdvanceOneStage_RK(1.0e0,0.0e0,1.0e0,dt,method,node,qold,
				       aux,q,aux,q,smax,map,leftBC,rightBC,1);
                    // --------------------------------------------------------    

                    break;

                case 2:  // Second order in time
		    
                    // ---------------------------------------------------------
                    // MainGrid: Stage #1
		    AdvanceOneStage_RK(1.0e0,0.0e0,1.0e0,dt,method,node,qold,
		    		       aux,q,auxstar,qstar,smax,map,leftBC,
		    		       rightBC,1);
		    // --------------------------------------------------------
		    
		    // ---------------------------------------------------------
                    // MainGrid: Stage #2
		    AdvanceOneStage_RK(0.5e0,0.5e0,0.5e0,dt,method,node,qold,
				       auxstar,qstar,aux,q,smax,map,leftBC,
				       rightBC,1);
		    // --------------------------------------------------------

                    break;
                    		    
                case 3:  // Third order in time
                    
                    // ---------------------------------------------------------
                    // MainGrid: Stage #1
  		    AdvanceOneStage_RK(1.0e0,0.0e0,1.0e0,dt,method,node,qold,
				       aux,q,auxstar,qstar,smax,map,leftBC,
				       rightBC,1);
		    // --------------------------------------------------------

		    // ---------------------------------------------------------
                    // MainGrid: Stage #2
		    AdvanceOneStage_RK(0.75e0,0.25e0,0.25e0,dt,method,node,qold,
				       auxstar,qstar,auxstar,qstar,smax,map,leftBC,
				       rightBC,1);
		    // --------------------------------------------------------
		    
		    // ---------------------------------------------------------
                    // MainGrid: Stage #3
		    AdvanceOneStage_RK(1.0e0/3.0e0,2.0e0/3.0e0,2.0e0/3.0e0,
				       dt,method,node,qold,auxstar,qstar,aux,q,
				       smax,map,leftBC,rightBC,1);
		    // --------------------------------------------------------
                    
                    break;
            }

            // compute cfl number -- MainGrid
            cfl = GetCFL(dt,dtv[2],prim_vol,method,aux,smax);	    
          
            // output time step information
            if (method[4]>0) 
            {
                cout << setprecision(3);
                cout << "MainGrid ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
	    dt_old = dt;
            if (cfl>0.0)
            {   	        
                dt = Min(dtv[2],dt*CFL_target/cfl);                
            }
            else
            {
  	        dt = Min(dtv[2],dt_old);
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
                    cout<<"MainGrid rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }
                
                // copy qold into q
                CopyQ(qold,q);
                CopyQ(auxold,aux);
	    }
      
        }

	// -------------------------------------------
	//     Hermite polynomial interpolation in
	//     time -- needed to get BCs for SubGrid
	// -------------------------------------------
	int i,me,k;
	int kmax = method[1];
	int istart,iend;
	if (leftBC.get(mbc_sub,2)==1)
	{
	    istart = leftBC.get(mbc_sub,1);
	    iend   = leftBC.get(1,1);
	    
	    GetHermiteCoeffs(method,istart,iend,dt_old,told,qold,auxold,
			     q,aux,node,HermiteCoeffs);
	}

	if (rightBC.get(1,2)==1)
	{
	    istart = rightBC.get(1,1);
	    iend   = rightBC.get(mbc_sub,1);
	    
	    GetHermiteCoeffs(method,istart,iend,dt_old,told,qold,auxold,
			     q,aux,node,HermiteCoeffs);
	}

	// -------------------------------
	//     SUB GRID
	// -------------------------------

	// choose timestep for SubGrid
	if (dt_sub > dt_old)
	  {  dt_sub = dt_old; }
	msubsteps = floor(dt_old/dt_sub + 0.5);
	dt_sub = dt_old/double(msubsteps);
	m_accept = 0;
	double cfl_tally = 0.0;
	double t_sub_old;

	// keep track of starting values in case we
	// need to reduce the timestep and try again
	CopyQ(q_sub,qtmp_sub);
	CopyQ(aux_sub,auxtmp_sub);

	// loop to compute msubsteps
	while (m_accept==0)
        {
	    m = 0;
	    cfl_tally = 0.0;
	    while(m<msubsteps && cfl_tally<=CFL_max)
	    {
	        // ending time
 	        m = m+1;
	        t_sub = told + double(m)*dt_sub;
		t_sub_old = t_sub - dt_sub;

		// Set boundary conditions on "SubGrid"
		SetBC_SubGrid(leftBC,rightBC,told,t_sub_old,node,aux,
			      HermiteCoeffs,q_sub);		

		// copy q into qold and aux into auxold
	        CopyQ(q_sub,qold_sub);
		CopyQ(aux_sub,auxold_sub);
		CopyQ(q_sub,qstar_sub);
		CopyQ(aux_sub,auxstar_sub);

		// Take a full time step of size dt
		switch ( abs(method_sub[2]) )
		{
                    case 1:  // First order in time
		      
			// --------------------------------------------------------
			// SubGrid: Stage #1 (the only one in this case)
			AdvanceOneStage_RK(1.0e0,0.0e0,1.0e0,dt_sub,method_sub,node_sub,
					   qold_sub,aux_sub,q_sub,aux_sub,q_sub,
					   smax_sub,map,leftBC,rightBC,2);
			// --------------------------------------------------------

			break;

		    case 2:  // Second order in time
			
			// --------------------------------------------------------
			// SubGrid: Stage #1
			AdvanceOneStage_RK(1.0e0,0.0e0,1.0e0,dt_sub,method_sub,node_sub,
					   qold_sub,aux_sub,q_sub,auxstar_sub,
					   qstar_sub,smax_sub,map,leftBC,rightBC,2);
			// --------------------------------------------------------
			
			// --------------------------------------------------------
			// SubGrid: Stage #2
			AdvanceOneStage_RK(0.5e0,0.5e0,0.5e0,dt_sub,method_sub,node_sub,
					   qold_sub,auxstar_sub,qstar_sub,aux_sub,
					   q_sub,smax_sub,map,leftBC,rightBC,2);
			// --------------------------------------------------------

			break;

		    case 3:  // Third order in time
			
			// --------------------------------------------------------
			// SubGrid: Stage #1
			AdvanceOneStage_RK(1.0e0,0.0e0,1.0e0,dt_sub,method_sub,node_sub,
					   qold_sub,aux_sub,q_sub,auxstar_sub,qstar_sub,
					   smax_sub,map,leftBC,rightBC,2);
			// --------------------------------------------------------
			
			// --------------------------------------------------------
			// SubGrid: Stage #2
			AdvanceOneStage_RK(0.75e0,0.25e0,0.25e0,dt_sub,method_sub,node_sub,
					   qold_sub,auxstar_sub,qstar_sub,auxstar_sub,
					   qstar_sub,smax_sub,map,leftBC,rightBC,2);
			// --------------------------------------------------------
			
			// --------------------------------------------------------
			// SubGrid: Stage #3
			AdvanceOneStage_RK(1.0e0/3.0e0,2.0e0/3.0e0,2.0e0/3.0e0,
					   dt_sub,method_sub,node_sub,qold_sub,
					   auxstar_sub,qstar_sub,aux_sub,q_sub,
					   smax_sub,map,leftBC,rightBC,2);
			// --------------------------------------------------------
                    
			break;
		}	    
	    
		// compute cfl number -- SubGrid
		cfl_sub = GetCFL(dt_sub,dtv_sub[2],prim_vol_sub,
				 method_sub,aux_sub,smax_sub);

		cfl_tally = Max(cfl_tally,cfl_sub);
	    
		// output time step information
		if (method_sub[4]>0) 
		{
		    cout << setprecision(3);
		    cout << "    SubGrid ... Step" << setw(5) << m;
		    cout << "   CFL =" << setw(6) << fixed << cfl_sub;
		    cout << "   dt =" << setw(11) << scientific << dt_sub;
		    cout << "   t =" << setw(11) << scientific << t_sub <<endl;
		}

	    }

	    // see whether to accept or reject this step
	    // for the SubGrid we only reject the first substep,
	    // after that we accept all steps until we reach the
	    // same output time as the MainGrid (told+dt)
	    if (cfl_tally<=CFL_max && m==msubsteps)
	    // accept
	    { 
	        m_accept = 1; 

		// choose new time step
                dt_sub = dt_sub*CFL_target/cfl_sub;
	    }
	    else
	    //reject
	    {  	        
	        dt_sub = dt_sub*CFL_max/cfl_tally;	       
		msubsteps = ceil(dt_old/dt_sub);		
		dt_sub = dt_old/double(msubsteps);
		cfl_tally = 0.0;
 
		if (method_sub[4]>0)
		{
		    cout<<"    SubGrid rejecting step...";
		    cout<<"CFL number too large";
		    cout<<endl;
		}
                
		// copy qold into q
		CopyQ(qtmp_sub,q_sub);
		CopyQ(auxtmp_sub,aux_sub);
	    }
	 
	}

	// --------------------------------------------------------
	// Map the values on SubGrid onto MainGrid
	SetValues_MainGrid(map,node_sub,aux_sub,q_sub,node,aux,q);
	// --------------------------------------------------------

        // compute conservation and print to file
        ConSoln(node,method,aux,q,t,outputdir);
    }

    // set initial time step for next call to DogSolve
    dtv[1]     = dt;
    dtv_sub[1] = dt_sub;

}
