//#undef CHECK_BOUNDS // this file is okay, so omit tensor bounds check
#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "FaceData.h"
#include "DogParams.h"
#include "DogParamsCart3.h"
#include "L2Project.h"
#include "Legendre3d.h"
#include "RiemannSolve.h"
#include "float.h" // for debugging

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t = RHS = -[ f(q,x,y,z,t)_x + g(q,x,y,z,t)_y + h(q,x,y,z,t)_z ] + Psi(q,x,y,z,t)
//
void ConstructL(const dTensorBC5& aux, 
		const dTensorBC5& q,
		dTensorBC5& Lstar, 
		dTensorBC4& smax)
{
  const FaceData& faceData = Legendre3d::get_faceData();
  const int space_order = dogParams.get_space_order();
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mz   = q.getsize(3);
  const int meqn = q.getsize(4);
  const int kmax = q.getsize(5);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(4);
  const int mpoints_face = space_order*space_order;
  
  dTensorBC5 Fm(mx,my,mz,meqn,mpoints_face,mbc,3);
  dTensorBC5 Fp(mx,my,mz,meqn,mpoints_face,mbc,3);
  dTensorBC5 Gm(mx,my,mz,meqn,mpoints_face,mbc,3);
  dTensorBC5 Gp(mx,my,mz,meqn,mpoints_face,mbc,3);
  dTensorBC5 Hm(mx,my,mz,meqn,mpoints_face,mbc,3);
  dTensorBC5 Hp(mx,my,mz,meqn,mpoints_face,mbc,3);
  
  // If you need access elsewhere to one of these declarations
  // move it into DogSolverCart3.h rather than copying the
  // declaration. (By using a single declaration we can append
  // arguments with default values without modifying existing
  // calls.) -eaj
  //
  void FluxFunc(const dTensor2& xpts,
		const dTensor2& Q,
		const dTensor2& Aux,
		dTensor3& flux);
  void LstarExtra(const dTensorBC5*,
		  const dTensorBC5*,
		  dTensorBC5*);
  
  // Grid information
  const double xlower = dogParamsCart3.get_xlow();
  const double ylower = dogParamsCart3.get_ylow();
  const double zlower = dogParamsCart3.get_zlow();
  const double dx = dogParamsCart3.get_dx();
  const double dy = dogParamsCart3.get_dy();
  const double dz = dogParamsCart3.get_dz();

  // ---------------------------------------------------------
  // Part I: compute source term
  // --------------------------------------------------------- 
  if ( dogParams.get_source_term()>0 )
    {        
      // Set source term on computational grid
      // Set values and apply L2-projection
      void SourceTermFunc(const dTensor2& xpts, 
			  const dTensor2& qvals,
			  const dTensor2& auxvals, 
			  dTensor2& source);
      L2Project(1-mbc,mx+mbc,
		1-mbc,my+mbc,
		1-mbc,mz+mbc,
                space_order,
		space_order,
		space_order,
		space_order,
                &q,
		&aux,
		&Lstar,
		&SourceTermFunc);
    }
  else
    {
      Lstar.setall(0.);
    }

  if ( !dogParams.get_flux_term() )
    {  return;  }
  // ---------------------------------------------------------


  // ---------------------------------------------------------
  // Part II: compute inter-element interaction fluxes
  // ---------------------------------------------------------
  
  // 1-direction: loop over interior faces and solve Riemann problems
  dTensor1 nvec(3);
  nvec.set(1, 1.0e0 );
  nvec.set(2, 0.0e0 );
  nvec.set(3, 0.0e0 );

#pragma omp parallel for
  for (int j=(2-mbc); j<=(my+mbc-1); j++)
    for (int k=(2-mbc); k<=(mz+mbc-1); k++)	  
      {
	dTensor1 Ql(meqn),Qr(meqn);
	dTensor1 Auxl(maux),Auxr(maux);
	dTensor1 Fl(meqn),Fr(meqn);
	dTensor1 xface(3);
	RiemannSolver riemannSolver(meqn,maux);

	for (int i=(2-mbc); i<=(mx+mbc); i++)
	  for (int ell=1; ell<=mpoints_face; ell++)
	    {	      
	      // Riemann data - q	      	      
	      for (int m=1; m<=meqn; m++)
		{
		  Ql.set(m, 0.0 );
		  Qr.set(m, 0.0 );
		  
		  for (int km=1; km<=kmax; km++)
		    {
		      Ql.fetch(m) += faceData.phi_xl->get(ell,km)*q.get(i-1,j,k,m,km);
		      Qr.fetch(m) += faceData.phi_xr->get(ell,km)*q.get(i,j,k,m,km);
		    }
		}
	      
	      // Riemann data - aux
	      for (int m=1; m<=maux; m++)
		{
		  Auxl.set(m, 0.0 );
		  Auxr.set(m, 0.0 );
		  
		  for (int km=1; km<=kmax; km++)
		    {
		      Auxl.fetch(m) += faceData.phi_xl->get(ell,km)*aux.get(i-1,j,k,m,km);
		      Auxr.fetch(m) += faceData.phi_xr->get(ell,km)*aux.get(i,j,k,m,km);
		    }
		}

	      // Solve Riemann problem
	      xface.set(1, xlower + (double(i)-1.0)*dx );
	      xface.set(2, ylower + (double(j)-0.5)*dy );
	      xface.set(3, zlower + (double(k)-0.5)*dz );
	      
	      const double smax_face
		= riemannSolver.solve(nvec,xface,Ql,Qr,Auxl,Auxr,Fl,Fr);
	      
	      smax.set(i-1, j, k, 1, Max(dy*dz*smax_face, smax.get(i-1, j, k, 1) ) );
	      smax.set(i,   j, k, 1, Max(dy*dz*smax_face, smax.get(i,   j, k, 1) ) );
	    
	      // Construct fluxes
	      for (int m=1; m<=meqn; m++)
		{
		  Fm.set(i  ,j,k,m,ell, Fr.get(m) );
		  Fp.set(i-1,j,k,m,ell, Fl.get(m) );
		}
	    }
      }

  // 2-direction: loop over interior faces and solve Riemann problems
  nvec.set(1, 0.0e0 );
  nvec.set(2, 1.0e0 );
  nvec.set(3, 0.0e0 );
  
#pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    for (int k=(2-mbc); k<=(mz+mbc-1); k++)	  
      {
	dTensor1 Ql(meqn),Qr(meqn);
	dTensor1 Auxl(maux),Auxr(maux);
	dTensor1 Fl(meqn),Fr(meqn);
	dTensor1 xface(3);
	RiemannSolver riemannSolver(meqn,maux);

	for (int j=(2-mbc); j<=(my+mbc); j++)
	  for (int ell=1; ell<=mpoints_face; ell++)
	    {
	      // Riemann data - q	      	      
	      for (int m=1; m<=meqn; m++)
		{
		  Ql.set(m, 0.0 );
		  Qr.set(m, 0.0 );
		  
		  for (int km=1; km<=kmax; km++)
		    {
		      Ql.fetch(m) += faceData.phi_yl->get(ell,km)*q.get(i,j-1,k,m,km);
		      Qr.fetch(m) += faceData.phi_yr->get(ell,km)*q.get(i,j,k,m,km);
		    }
		}
	      
	      // Riemann data - aux
	      for (int m=1; m<=maux; m++)
		{
		  Auxl.set(m, 0.0 );
		  Auxr.set(m, 0.0 );
		  
		  for (int km=1; km<=kmax; km++)
		    {
		      Auxl.fetch(m) += faceData.phi_yl->get(ell,km)*aux.get(i,j-1,k,m,km);
		      Auxr.fetch(m) += faceData.phi_yr->get(ell,km)*aux.get(i,j,k,m,km);
		    }
		}

	      // Solve Riemann problem
	      xface.set(1, xlower + (double(i)-0.5)*dx );
	      xface.set(2, ylower + (double(j)-1.0)*dy );
	      xface.set(3, zlower + (double(k)-0.5)*dz );
	      
	      const double smax_face
		= riemannSolver.solve(nvec,xface,Ql,Qr,Auxl,Auxr,Fl,Fr);
	      
	      smax.set(i, j-1, k, 2, Max(dx*dz*smax_face, smax.get(i, j-1, k, 1) ) );
	      smax.set(i,   j, k, 2, Max(dx*dz*smax_face, smax.get(i,   j, k, 1) ) );
	    
	      // Construct fluxes
	      for (int m=1; m<=meqn; m++)
		{
		  Gm.set(i,j  ,k,m,ell, Fr.get(m) );
		  Gp.set(i,j-1,k,m,ell, Fl.get(m) );
		}
	    }
      }

  // 3-direction: loop over interior faces and solve Riemann problems
  nvec.set(1, 0.0e0 );
  nvec.set(2, 0.0e0 );
  nvec.set(3, 1.0e0 );
  
#pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc); i++)
    for (int j=(2-mbc); j<=(my+mbc); j++)	  
      {
	dTensor1 Ql(meqn),Qr(meqn);
	dTensor1 Auxl(maux),Auxr(maux);
	dTensor1 Fl(meqn),Fr(meqn);
	dTensor1 xface(3);
	RiemannSolver riemannSolver(meqn,maux);

	for (int k=(2-mbc); k<=(mz+mbc-1); k++)
	  for (int ell=1; ell<=mpoints_face; ell++)
	    {
	      // Riemann data - q	      	      
	      for (int m=1; m<=meqn; m++)
		{
		  Ql.set(m, 0.0 );
		  Qr.set(m, 0.0 );
		  
		  for (int km=1; km<=kmax; km++)
		    {
		      Ql.fetch(m) += faceData.phi_zl->get(ell,km)*q.get(i,j,k-1,m,km);
		      Qr.fetch(m) += faceData.phi_zr->get(ell,km)*q.get(i,j,k,m,km);
		    }
		}
	      
	      // Riemann data - aux
	      for (int m=1; m<=maux; m++)
		{
		  Auxl.set(m, 0.0 );
		  Auxr.set(m, 0.0 );
		  
		  for (int km=1; km<=kmax; km++)
		    {
		      Auxl.fetch(m) += faceData.phi_zl->get(ell,km)*aux.get(i,j,k-1,m,km);
		      Auxr.fetch(m) += faceData.phi_zr->get(ell,km)*aux.get(i,j,k,m,km);
		    }
		}

	      // Solve Riemann problem
	      xface.set(1, xlower + (double(i)-0.5)*dx );
	      xface.set(2, ylower + (double(j)-0.5)*dy );
	      xface.set(3, zlower + (double(k)-1.0)*dz );
	      
	      const double smax_face
		= riemannSolver.solve(nvec,xface,Ql,Qr,Auxl,Auxr,Fl,Fr);
	      
	      smax.set(i, j, k-1, 3, Max(dx*dy*smax_face, smax.get(i, j, k-1, 1) ) );
	      smax.set(i, j, k  , 3, Max(dx*dy*smax_face, smax.get(i, j, k,   1) ) );
	    
	      // Construct fluxes
	      for (int m=1; m<=meqn; m++)
		{
		  Hm.set(i,j,k  ,m,ell, Fr.get(m) );
		  Hp.set(i,j,k-1,m,ell, Fl.get(m) );
		}
	    }
      }

    // Compute ``flux differences'' dF, dG, dH
    const double quarter_dx_inv = 0.25/dx;
    const double quarter_dy_inv = 0.25/dy;
    const double quarter_dz_inv = 0.25/dz;
    const int mlength = Lstar.getsize(4);

#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)    
      for (int j=(2-mbc); j<=(my+mbc-1); j++)
	for (int k=(2-mbc); k<=(mz+mbc-1); k++)
	  for (int m=1; m<=mlength; m++)
	    for (int km=1; km<=kmax; km++)
	      {
		// 1-direction: dF
		double F1 = 0.0;
		double F2 = 0.0;
		for (int ell=1; ell<=mpoints_face; ell++)
		  {
		    F1 = F1 + faceData.wght_phi_xr->get(ell,km)*Fm.get(i,j,k,m,ell);
		    F2 = F2 + faceData.wght_phi_xl->get(ell,km)*Fp.get(i,j,k,m,ell);
		  }

		// 2-direction: dG
		double G1 = 0.0;
		double G2 = 0.0;
		for (int ell=1; ell<=mpoints_face; ell++)
		  {
		    G1 = G1 + faceData.wght_phi_yr->get(ell,km)*Gm.get(i,j,k,m,ell);
		    G2 = G2 + faceData.wght_phi_yl->get(ell,km)*Gp.get(i,j,k,m,ell);
		  }

		// 3-direction: dH
		double H1 = 0.0;
		double H2 = 0.0;
		for (int ell=1; ell<=mpoints_face; ell++)
		  {
		    H1 = H1 + faceData.wght_phi_zr->get(ell,km)*Hm.get(i,j,k,m,ell);
		    H2 = H2 + faceData.wght_phi_zl->get(ell,km)*Hp.get(i,j,k,m,ell);
		  }

		// Update
		Lstar.fetch(i,j,k,m,km) -= ( quarter_dx_inv*(F2-F1) +
					     quarter_dy_inv*(G2-G1) +
					     quarter_dz_inv*(H2-H1) );
	      }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part III: compute intra-element contributions
    // ---------------------------------------------------------
    // No need to call this if first-order in space    
    if(dogParams.get_space_order()>1)
      {
        L2ProjectGradAdd(1-mbc,mx+mbc,
			 1-mbc,my+mbc,
			 1-mbc,mz+mbc,
			 space_order,
			 space_order,
			 space_order,
			 space_order,
			 &q,
			 &aux,
			 &Lstar,
			 &FluxFunc);
      }
    // ---------------------------------------------------------  


    // ---------------------------------------------------------
    // Part IV: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra(&q,&aux,&Lstar);
    // ---------------------------------------------------------
}
