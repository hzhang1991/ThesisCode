#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include "DogParams.h"
#include "DogStateUnst2.h"

void DogSolveUser_Unst(const mesh& Mesh, const edge_data_Unst& EdgeData,
		       dTensor3& aux, dTensor3& qold, dTensor3& qnew,                       
		       const double tstart, const double tend, 
		       const string outputdir)
{
  int i,n_step,m_accept,mtmp;
  double t,dt,CFL_max,CFL_target,dtmin,dtmax;
  double told,cfl;
  int mx   = qnew.getsize(1);
  int meqn = qnew.getsize(2);
  int kmax = qnew.getsize(3);
  int maux = aux.getsize(2);
  const double* cflv     = dogParams.get_cflv();
  int           nv       = dogParams.get_nv();

  // ------------------------------------------------------------
  // Function definitions
  void CopyQ_Unst(const dTensor3&,dTensor3&);
  void ConSoln_Unst(const mesh& Mesh, 
		    const dTensor3& aux, const dTensor3& q, 
		    double t, string outputdir);
  void UpdateSoln_Unst(const double alpha1, const double alpha2, 
		       const double beta, const double dt, const mesh& Mesh,
		       const dTensor3& aux, const dTensor3& qstar, 
		       const dTensor3& Lstar, dTensor3& qnew);
  void BeforeStep_Unst(const double,const mesh&,dTensor3&,dTensor3&);
  void  AfterStep_Unst(const double,const mesh&,dTensor3&,dTensor3&);
  void AfterFullTimeStep_Unst(const double dt, const mesh& Mesh,
			      const dTensor3& auxold, const dTensor3& qold,
			      const dTensor3& Lold, dTensor3& aux, dTensor3& q);
  void UpdateSoln_Unst(const double alpha1, const double alpha2, 
		       const double beta, const double dt, const mesh& Mesh,
		       dTensor3& aux, const dTensor3& qstar, 
		       const dTensor3& Lstar, dTensor3& qnew);
  void ConstructL_Unst(const mesh& Mesh,
		       const edge_data_Unst& EdgeData,
		       dTensor3& aux, // SetBndValues modifies ghost cells
		       dTensor3& q,   // SetBndValues modifies ghost cells
		       dTensor3& Lstar,
		       dTensor1& smax);
  double GetCFL_Unst(double dt, const mesh& Mesh,
		     const dTensor3& aux, const dTensor1& smax);
  // ------------------------------------------------------------

  // define local variables
  int s,m,k;
  double tmp;
  n_step = 0;
  t = tstart;
  dt = dogStateUnst2.get_initial_dt();
  CFL_max = cflv[1];
  CFL_target = cflv[2];
  cfl = 0.0;
  dtmin = dt;
  dtmax = dt;
  int NumElems = Mesh.get_NumElems(); // Number of total elements in mesh
  int NumNodes = Mesh.get_NumNodes(); // Number of nodes in mesh
  int NumEdges = Mesh.get_NumEdges(); // Number of edges in mesh 
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
      m_accept = 0;      
      n_step = n_step + 1;
      
      // check if max number of time steps exceeded
      if (n_step>nv)
	{
          eprintf(" Error in DogSolveUser_Unst.cpp: "
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
	  told = t;
	  if (told+dt > tend)
            { dt = tend - told; }
	  t = told + dt;
	  dogStateUnst2.set_time(told);
          dogStateUnst2.set_dt(dt);
	  
	  // Set initial maximum wave speed to zero
	  for (i=1; i<=NumEdges; i++)
	    {  smax.set(i, 0.0e0 );  }
	    
	  // ----------------------------------------------------------------
	  //
	  //    THIS IS WHERE THE USER-DEFINED TIME-STEPPING SCHEME
	  //    SHOULD BE ADDED. IN THE DEFAULT FILE: DogSolveUser_Unst.cpp,
	  //    THE PROGRAM WILL NOW RETURN AN ERROR MESSAGE.
	  // 
	  // ----------------------------------------------------------------
          eprintf("\n"
	          " No user-defined time-stepping scheme has been defined yet. \n"
	          " Copy $DOGPACK/lib/2d/cart/DogSolveUser.cpp into the current \n"
	          " directory and modify as needed.\n\n");
	  // ----------------------------------------------------------------
	  
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
            }      
        }

      // compute conservation and print to file
      ConSoln_Unst(Mesh,aux,qnew,t,outputdir);
    }
    
  // set initial time step for next call to DogSolveRK
  dogStateUnst2.set_initial_dt(dt);
}
