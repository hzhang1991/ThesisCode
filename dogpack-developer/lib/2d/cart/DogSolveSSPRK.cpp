#include "dogdefs.h"
#include "dog_math.h"

void DogSolveSSPRK(dTensorBC4& aux, dTensorBC4& qold, dTensorBC4& qnew,
		   dTensorBC3& smax, double tstart, double tend, int nv, const int method[], 
		   double dtv[], const double cflv[], string outputdir)
{
    int i,j,m,k,n_step,m_accept,mtmp,s;
    double t,dt,CFL_max,CFL_target,dtmin,dtmax;
    double told,cfl;
    double tmp;
    int mx   = qnew.getsize(1);
    int my   = qnew.getsize(2);
    int meqn = qnew.getsize(3);
    int kmax = qnew.getsize(4);
    int mbc  = qnew.getmbc();
    int maux = aux.getsize(3);
  
    // ------------------------------------------------------------
    // Function definitions
    void CopyQ(const dTensorBC4&,dTensorBC4&);
    void ConSoln(const int method[], 
		 const dTensorBC4& aux,
		 const dTensorBC4& q, double t, string outputdir);
    void BeforeStep(double,const dTensor3&,dTensorBC4&,dTensorBC4&);
    void UpdateSoln(double alpha1, double alpha2, double beta, double dt,
		    const dTensorBC4& aux, 
		    const dTensorBC4& qstar, const dTensorBC4& Lstar, 
		    dTensorBC4& qnew);
    void  AfterStep(double,const dTensor3&,dTensorBC4&,dTensorBC4&);
    void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
		      void (*ProjectRightEig)(int,const dTensor1&,
					      const dTensor1&,const dTensor2&,
					      dTensor2&),
		      void (*ProjectLeftEig)(int,const dTensor1&,
					     const dTensor1&,const dTensor2&,
					     dTensor2&));
    void ConstructL(const int method[],
		    dTensorBC4& aux, // SetBndValues modifies ghost cells
		    dTensorBC4& q,   // SetBndValues modifies ghost cells
		    dTensorBC4& Lstar, dTensorBC3& smax);
    double GetCFL(const int method[], double dt,
		  const dTensorBC4& aux, const dTensorBC3& smax);
    void ProjectRightEig(int,const dTensor1&,const dTensor1&,
			 const dTensor2&,dTensor2&);
    void ProjectLeftEig(int,const dTensor1&,const dTensor1&,
			const dTensor2&,dTensor2&);
    // ------------------------------------------------------------

    // define local variables
    n_step = 0;
    t = tstart;
    dt = dtv[1];
    CFL_max = cflv[1];
    CFL_target = cflv[2];
    cfl = 0.0;
    dtmin = dt;
    dtmax = dt;
    dTensorBC4   q1(mx,my,meqn,kmax,mbc);
    dTensorBC4   q2(mx,my,meqn,kmax,mbc);
    mtmp = maux;
    dTensorBC4   Lstar(mx,my,meqn,kmax,mbc);  
    dTensorBC4 auxold(mx,my,mtmp,kmax,mbc);
	
    // Main loop
    while (t<tend)
    {
        // initialize time step
        m_accept = 0;      
        n_step = n_step + 1;

	// check if max number of time steps exceeded
	if (n_step>nv)
	{
            eprintf(" Error in DogSolveSSPRK.cpp: "
	        " Exceeded allowed # of time steps \n"
	        "    n_step = %d\n"
	        "        nv = %d\n\n",
                n_step,nv);
	}

        // copy qnew into qold
        CopyQ(qnew,qold);
	CopyQ(aux,auxold);

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {
            // set current time
            told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;
            
            // Set initial maximum wave speed to zero
            for (j=1-mbc; j<=(my+mbc); j++)
	      for (i=1-mbc; i<=(mx+mbc); i++)
                {
		  smax.set(i,j,1, 0.0e0 );
		  smax.set(i,j,2, 0.0e0 );
                }
	    
            // Take a full time step of size dt
            switch ( method[2] )
            {
                case 1:  // First order in time  (1-stage)          		  
                    // -----------------------------------------------
                    BeforeStep(dt,aux,qnew);		    
                    ConstructL(method,aux,qnew,Lstar,smax);
                    UpdateSoln(1.0,0.0,1.0,dt,aux,qnew,Lstar,qnew);
                    AfterStep(dt,aux,qnew);
                    // -----------------------------------------------
                    break;
                    
                case 2:  // Second order in time (3-stages)
                    // ----------------------------------------------- 
                    CopyQ(qnew,q1);

		    // Stage: 1 and 2
		    for (s=1; s<=2; s++)
		      {
			BeforeStep(dt,aux,q1);
			ConstructL(method,aux,q1,Lstar,smax);
			UpdateSoln(1.0,0.0,0.5,dt,aux,q1,Lstar,q1);
			if (method[3]==1)
			  {  ApplyLimiter(aux,q1,
					  &ProjectRightEig,&ProjectLeftEig);  }
			AfterStep(dt,aux,q1);
		      }

		    // Stage: 3
		    BeforeStep(dt,aux,q1);
		    ConstructL(method,aux,q1,Lstar,smax);
		    UpdateSoln(2.0/3.0,1.0/3.0,1.0/3.0,dt,aux,q1,Lstar,qnew);
		    if (method[3]==1)
		      {  ApplyLimiter(aux,qnew,
				      &ProjectRightEig,&ProjectLeftEig);  }
		    AfterStep(dt,aux,qnew);
                    // -----------------------------------------------
                    break;
                    
                case 3:  // Third order in time (9-stages)
		    // -----------------------------------------------
		    CopyQ(qnew,q1);

		    // Stage: 1
		    BeforeStep(dt,aux,q1);
		    ConstructL(method,aux,q1,Lstar,smax);
		    UpdateSoln(1.0,0.0,1.0/6.0,dt,aux,q1,Lstar,q1);
		    if (method[3]==1)
		      {  ApplyLimiter(aux,q1,
				      &ProjectRightEig,&ProjectLeftEig);  }
		    AfterStep(dt,aux,q1);		     

		    CopyQ(q1,q2);

		    // Stage: 2, 3, 4, and 5
		    for (s=2; s<=5; s++)
		      {
			BeforeStep(dt,aux,q1);
			ConstructL(method,aux,q1,Lstar,smax);
			UpdateSoln(1.0,0.0,1.0/6.0,dt,aux,q1,Lstar,q1);
			if (method[3]==1)
			  {  ApplyLimiter(aux,q1,
					  &ProjectRightEig,&ProjectLeftEig);  }
			AfterStep(dt,aux,q1);
		      }

		    // Stage: 6
		    BeforeStep(dt,aux,q1);
		    ConstructL(method,aux,q1,Lstar,smax);
		    UpdateSoln(3.0/5.0,2.0/5.0,1.0/15.0,dt,aux,q2,Lstar,q1);
		    if (method[3]==1)
		      {  ApplyLimiter(aux,q1,
				      &ProjectRightEig,&ProjectLeftEig);  }
		    AfterStep(dt,aux,q1);

		    // Stage: 7,8, and 9
		    for (s=7; s<=9; s++)
		      {
			BeforeStep(dt,aux,q1);
			ConstructL(method,aux,q1,Lstar,smax);
			UpdateSoln(1.0,0.0,1.0/6.0,dt,aux,q1,Lstar,q1);
			if (method[3]==1)
			  {  ApplyLimiter(aux,q1,
					  &ProjectRightEig,&ProjectLeftEig);  }
			AfterStep(dt,aux,q1);
		      }
		    
		    CopyQ(q1,qnew);		  
		    // -----------------------------------------------
                    break;

	        case 4:  // Fourth order in time (10-stages)		  
		    // -----------------------------------------------
  		    CopyQ(qnew,q1);
		    CopyQ(q1,q2);

		    // Stage: 1,2,3,4, and 5
		    for (s=1; s<=5; s++)
		      {
			BeforeStep(dt,aux,q1);
			ConstructL(method,aux,q1,Lstar,smax);
			UpdateSoln(1.0,0.0,1.0/6.0,dt,aux,q1,Lstar,q1);
			if (method[3]==1)
			  {  ApplyLimiter(aux,q1,
					  &ProjectRightEig,&ProjectLeftEig);  }
			AfterStep(dt,aux,q1);
		      }
		    
		    // Temporary storage
		    for (i=(2-mbc); i<=(mx+mbc-1); i++)
		      for (j=(2-mbc); j<=(my+mbc-1); j++)
			for (m=1; m<=meqn; m++)
			  for (k=1; k<=kmax; k++)
			    {
			      tmp = (q2.get(i,j,m,k) + 9.0*q1.get(i,j,m,k))/25.0;
			      q2.set(i,j,m,k, tmp );
			      q1.set(i,j,m,k, 15.0*tmp - 5.0*q1.get(i,j,m,k) );
			    }

		    // Stage: 6,7,8, and 9
		    for (s=6; s<=9; s++)
		      {
			BeforeStep(dt,aux,q1);
			ConstructL(method,aux,q1,Lstar,smax);
			UpdateSoln(1.0,0.0,1.0/6.0,dt,aux,q1,Lstar,q1);
			if (method[3]==1)
			  {  ApplyLimiter(aux,q1,
					  &ProjectRightEig,&ProjectLeftEig);  }
			AfterStep(dt,aux,q1);
		      }
        
		    // Stage: 10		    
		    BeforeStep(dt,aux,q1);
		    ConstructL(method,aux,q1,Lstar,smax);
		    UpdateSoln(1.0,3.0/5.0,1.0/10.0,dt,aux,q2,Lstar,q1);
		    if (method[3]==1)
		      {  ApplyLimiter(aux,q1,
				      &ProjectRightEig,&ProjectLeftEig);  }
		    AfterStep(dt,aux,q1);
		      
		    CopyQ(q1,qnew);
  		    // -----------------------------------------------		  
		    break;
            }

            // compute cfl number
            cfl = GetCFL(method,dt,aux,smax);
	  
            // output time step information
            if (method[4]>0) 
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
                    printf("DogSolve2D rejecting step..."
		           "CFL number too large\n");
                }
                
                // copy qold into qnew
                CopyQ(qold,qnew);
		CopyQ(auxold,aux);
            }      
        }

        // compute conservation and print to file
        ConSoln(method,aux,qnew,t,outputdir);
    }
    
    // set initial time step for next call to DogSolveRK
    dtv[1] = dt;

}
