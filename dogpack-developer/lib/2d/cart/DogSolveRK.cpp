/*
 * This method is implemented in lib/DogSolver
 *
 */

//  #include "dogdefs.h"
//  #include "dog_math.h"
//  #include "edge_data.h"
//  #include "RKinfo.h"
//  #include "DogParams.h"
//  #include "DogStateCart2.h"
//  #include "DogSolverCart2.h"
//  #include <string>

//  // This method could be implemented in
//  // the DogSolver class, since the ODE solver we choose
//  // generally does not involve anything specific to the grid.
//  //
//  void DogSolverCart2::DogSolveRK(double tstart, double tend)
//  {

//      assert_eq(tstart,get_state().get_time());
//      const int nv = dogParams.get_nv();
//      const double* cflv = dogParams.get_cflv();
//      dTensorBC4& qnew = fetch_state().fetch_q();
//      dTensorBC4& aux  = fetch_state().fetch_aux();
//      dTensorBC3& smax = fetch_smax();

//      // Define local variables
//      int n_step = 0;
//      double t = tstart;
//      assert_eq(tstart,get_state().get_time());
//      double dt = get_dt();
//      const double CFL_max    = cflv[1];
//      const double CFL_target = cflv[2];

//      while (t<tend)
//      {
//          // initialize time step
//          int m_accept = 0;      
//          n_step = n_step + 1;

//          // check if max number of time steps exceeded
//          if (n_step>nv)
//          {
//              eprintf(" Error in DogSolveRK.cpp: "
//                      " Exceeded allowed # of time steps \n"
//                      "    n_step = %d\n"
//                      "        nv = %d\n\n",
//                      n_step,nv);
//          }

//          // Copy qnew into qold (in order to save data)
//          saveState();

//          // keep trying until we get a dt that does not violate CFL condition
//          while (m_accept==0)
//          {
//              // set current time
//              double told = t;
//              if (told+dt > tend)
//              { dt = tend - told; }

//              // what the time will be at the conclusion of this step
//              t = told + dt;

//              assert_eq(get_state().get_time(), told);
//              set_dt(dt);

//              // Set initial maximum wave speed to zero
//              smax.setall(0.);

//              // Take a single time step (TODO: why is memory allocated here and
//              // not above?)
//              advanceTimeStepRK(dt);

//              // compute cfl number
//              double cfl = GetCFL(dt);

//              // output time step information
//              if (dogParams.get_verbosity()>0) 
//              {
//                  printf("DogSolve2D ... Step %5d"
//                          "   CFL =%6.3f"
//                          "   dt =%11.3e"
//                          "   t =%11.3e\n",
//                          n_step,cfl,dt,t);
//              }

//              // choose new time step
//              if (cfl>0.0)
//              {   
//                  dt = Min(dogParams.get_max_dt(),dt*CFL_target/cfl);
//              }
//              else
//              {
//                  dt = dogParams.get_max_dt();
//                  printf("setting dt equal to max dt = %f\n", dt );
//              }

//              // see whether to accept or reject this step
//              if (cfl<=CFL_max)
//              // accept
//              { 
//                  m_accept = 1; 
//                  fetch_state().set_time(t);

//                  void AfterFullTimeStep(DogSolverCart2& solver);
//                  AfterFullTimeStep(*this);
//              }
//              else 
//              //reject
//              {   
//                  if (dogParams.get_verbosity()>0)
//                  {
//                      printf("DogSolve2D rejecting step..."
//                              "CFL number too large\n");
//                  }

//                  // revert to old state
//                  //
//                  revertToSavedState();
//                  t = get_state().get_time();
//                  assert_eq(t, told);

//              }      
//          }

//          // compute conservation and print to file
//          // The arguments here should be changed to *this
//          // (received as const)
//          //
//          //       -- Why is it better to pass *this than aux, q and time?  how do
//          //       you access the time for printing conserved quantities if you
//          //       don't explicitely pass it in here?
//          //                                              -DS
//          ConSoln(aux,qnew,t);
//      }

//      // set initial time step for next call to DogSolve
//      set_dt(dt);

//  }

//  // Take a full time step of size dt
//  //
//  // This function needs to be changed to maintain
//  // dogState->get_time() at the current time
//  // at each time stage (in order to agree with the
//  // rest of the state data in q and aux).
//  // I think that this could be done as part of
//  // implementing an AdvanceTimeStage() method.
//  // This will be important for time-dependent
//  // source terms as needed for manufactured
//  // solutions. -eaj
//  //
//  int DogSolverCart2::advanceTimeStepRK(double dt)
//  {
//      void CopyQ(const dTensorBC4& qin,dTensorBC4& qout);

//      // access variable state data
//      dTensorBC4& qnew = fetch_state().fetch_q();
//      dTensorBC4& aux  = fetch_state().fetch_aux();
//      dTensorBC3& smax = fetch_smax();

//      const int mx   = qnew.getsize(1);
//      const int my   = qnew.getsize(2);
//      const int meqn = qnew.getsize(3);
//      const int kmax = qnew.getsize(4);
//      const int mbc  = qnew.getmbc();
//      const int maux = aux.getsize(3);

//      RKinfo rk;
//      void SetRKinfo(int time_order, RKinfo& rk);
//      SetRKinfo(dogParams.get_time_order(),rk);
//      dTensorBC4 Lold(mx,my,meqn,kmax,mbc);

//      // should use gprof to see if memory allocation here is expensive;
//      // if so, put it back in DogSolveRK
//      //
//      // value of inf or nan would cause problems in qstar so we initialize to 0
//      //
//      dTensorBC4   qstar(mx,my,meqn,kmax,mbc); qstar.setall(0.);

//  // -------------------------------------------------------------------------- //
//      // !!! BUG: auxstar needs to be set at some point during the simulation!!!
//      // Where's the correct spot to do it? (-DS)
//      //
//      //dTensorBC4 auxstar(mx,my,maux,kmax,mbc); auxstar.setall(0.);
//      dTensorBC4 auxstar(mx,my,maux,kmax,mbc); auxstar.copyfrom( aux );
//  // -------------------------------------------------------------------------- //

//      dTensorBC4      q1(mx,my,meqn,kmax,mbc);
//      dTensorBC4      q2(mx,my,meqn,kmax,mbc);
//      dTensorBC4   Lstar(mx,my,meqn,kmax,mbc);

//      void UpdateSoln(double alpha1, double alpha2, double beta, double dt,
//              const dTensorBC4& aux, const dTensorBC4& qstar,
//              const dTensorBC4& Lstar, dTensorBC4& qnew);
//      void UpdateSoln( double g1,double g2, double g3, double delta, 
//              double beta,double dt, const dTensorBC4& aux,
//              const dTensorBC4& qold, const dTensorBC4& Lstar,
//              dTensorBC4& q1, dTensorBC4& q2);

//      // I think we should create an AdvanceTimeStage() method. -eaj
//      switch ( dogParams.get_time_order() )
//      {
//          case 1:  // First order in time (1-stage)

//              // -----------------------------------------------
//              // Stage #1 (the only one in this case)  
//              rk.mstage = 1;
//              BeforeStep(dt,aux,qnew,*this);
//              ::ConstructL(aux,qnew,Lstar,smax);
//              CopyQ(Lstar,Lold);
//              UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                      rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qnew);
//              AfterStep(dt,aux,qnew,*this);            
//              // -----------------------------------------------
//              break;

//          case 2:  // Second order in time (2-stages)

//              // -----------------------------------------------
//              // Stage #1  
//              rk.mstage = 1;
//              BeforeStep(dt,aux,qnew,*this);
//              ::ConstructL(aux,qnew,Lstar,smax);
//              CopyQ(Lstar,Lold);
//              UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                      rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qstar);
//              if(dogParams.using_moment_limiter())
//              {  ApplyLimiter(aux,qstar,&ProjectRightEig,&ProjectLeftEig);  }
//              AfterStep(dt,auxstar,qstar,*this);

//              // ------------------------------------------------
//              // Stage #2
//              rk.mstage = 2;
//  // -------------------------------------------------------------------------- //
//  // BUG: AUXSTAR IS NEVER SAVED BEFORE THIS CALL!
//  // -------------------------------------------------------------------------- //
//              BeforeStep(dt,auxstar,qstar,*this);
//              ::ConstructL(aux,qstar,Lstar,smax);
//              UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                      rk.beta->get(rk.mstage),dt,auxstar,qstar,Lstar,qnew);
//              if(dogParams.using_moment_limiter())
//              {  ApplyLimiter(auxstar,qnew,&ProjectRightEig,&ProjectLeftEig);  }
//              AfterStep(dt,aux,qnew,*this);
//              // ------------------------------------------------
//              break;

//          case 3:  // Third order in time (3-stages)

//              // ------------------------------------------------
//              // Stage #1
//              rk.mstage = 1;
//              BeforeStep(dt,aux,qnew,*this);
//              ::ConstructL(aux,qnew,Lstar,smax);
//              CopyQ(Lstar,Lold);
//              UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                      rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qstar);
//              if(dogParams.using_moment_limiter())
//              {  ApplyLimiter(aux,qstar,&ProjectRightEig,&ProjectLeftEig);  }
//              AfterStep(dt,aux,qstar,*this);
//              // -------------------------------------------------
//              // Stage #2
//              rk.mstage = 2;
//              BeforeStep(dt,aux,qstar,*this);
//              ::ConstructL(aux,qstar,Lstar,smax);
//              UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                      rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qstar);   
//              if(dogParams.using_moment_limiter())
//              {  ApplyLimiter(auxstar,qstar,&ProjectRightEig,&ProjectLeftEig);  }
//              AfterStep(dt,aux,qstar,*this);
//              // --------------------------------------------------
//              // Stage #3
//              rk.mstage = 3;
//              BeforeStep(dt,auxstar,qstar,*this);
//              ::ConstructL(auxstar,qstar,Lstar,smax);
//              UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                      rk.beta->get(rk.mstage),dt,aux,qstar,Lstar,qnew);   
//              if(dogParams.using_moment_limiter())
//              {  ApplyLimiter(auxstar,qnew,&ProjectRightEig,&ProjectLeftEig);  }
//              AfterStep(dt,aux,qnew,*this);
//              // --------------------------------------------------   
//              break;

//          case 4:  // Fourth order in time (10-stages)

//              // -----------------------------------------------
//              CopyQ(qnew,q1);
//              CopyQ(q1,q2);

//              // Stage: 1,2,3,4, and 5
//              for (int s=1; s<=5; s++)
//              {
//                  rk.mstage = s;
//                  BeforeStep(dt,aux,q1,*this);
//                  ::ConstructL(aux,q1,Lstar,smax);
//                  if (s==1)
//                  {  CopyQ(Lstar,Lold);  }
//                  UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                          rk.beta->get(rk.mstage),dt,aux,q1,Lstar,q1);
//                  if(dogParams.using_moment_limiter())
//                  {  ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig);  }
//                  AfterStep(dt,aux,q1,*this);
//              }

//              // Temporary storage
//  #pragma omp parallel for
//              for (int i=(2-mbc); i<=(mx+mbc-1); i++)
//                  for (int j=(2-mbc); j<=(my+mbc-1); j++)
//                      for (int m=1; m<=meqn; m++)
//                          for (int k=1; k<=kmax; k++)
//                          {
//                              const double tmp = (q2.get(i,j,m,k) + 9.0*q1.get(i,j,m,k))/25.0;
//                              q2.set(i,j,m,k, tmp );
//                              q1.set(i,j,m,k, 15.0*tmp - 5.0*q1.get(i,j,m,k) );
//                          }

//              // Stage: 6,7,8, and 9
//              for (int s=6; s<=9; s++)
//              {      
//                  rk.mstage = s;

//                  BeforeStep(dt,aux,q1,*this);
//                  ::ConstructL(aux,q1,Lstar,smax);
//                  UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                          rk.beta->get(rk.mstage),dt,aux,q1,Lstar,q1);
//                  if(dogParams.using_moment_limiter())
//                  {  ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig);  }
//                  AfterStep(dt,aux,q1,*this);
//              }

//              // Stage: 10
//              rk.mstage = 10;

//              BeforeStep(dt,aux,q1,*this);
//              ::ConstructL(aux,q1,Lstar,smax);
//              UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
//                      rk.beta->get(rk.mstage),dt,aux,q2,Lstar,q1);
//              if(dogParams.using_moment_limiter())
//              {  ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig);  }
//              AfterStep(dt,aux,q1,*this);

//              CopyQ(q1,qnew);
//              // -----------------------------------------------          
//              break;

//          case 5:  // Fifth order in time (8-stages)

//              // -----------------------------------------------
//              CopyQ(qnew,q1);
//              q2.setall(0.);

//              for (int s=1; s<=8; s++)
//              {
//                  rk.mstage = s;
//                  BeforeStep(dt,aux,q1,*this);
//                  ::ConstructL(aux,q1,Lstar,smax);
//                  if (s==1)
//                  {  CopyQ(Lstar,Lold);  }

//                  rk.gamma->get(1,s); 
//                  rk.gamma->get(2,s); 
//                  rk.gamma->get(3,s);
//                  rk.delta->get(s);
//                  rk.beta->get(s);

//                  UpdateSoln(rk.gamma->get(1,s), 
//                          rk.gamma->get(2,s), 
//                          rk.gamma->get(3,s), 
//                          rk.delta->get(s), rk.beta->get(s),
//                          //dt, aux, qold, Lstar, q1, q2);
//                          dt, aux, qnew, Lstar, q1, q2);

//                  if (dogParams.using_moment_limiter())
//                  {  ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig);  }
//                  AfterStep(dt,aux,q1,*this);
//              }

//              CopyQ(q1,qnew);
//              // -----------------------------------------------          
//              break;

//          default: unsupported_value_error(dogParams.get_time_order());
//      }

//      // free space allocated for RKInfo :
//      void DeleteRKInfo(RKinfo& rk);    
//      DeleteRKInfo( rk );

//      return 0;

//  }

//
