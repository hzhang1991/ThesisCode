//#undef CHECK_BOUNDS // this file is okay, so omit tensor bounds check
#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "L2Project.h"
#include "RiemannSolve.h"
#include "Legendre2d.h"
#include "float.h" // for debugging
#include <iostream>
using namespace std;
// ------------------------------------------------------------------------- //
// LLF fluxes for hyperbolic PDE in divergence form
//
//       q_t = RHS =  -[ f(q,x,y,t)_x + g(q,x,y,t)_y ] + Psi(q,x,y,t)
//
// It is expected that the user sets the correct boundary conditions before
// calling this routine.
//
// Input:
// ------
//
//     aux(1:mx, 1:my, 1:maux, 1:kmax ) - auxiliary function
//       q(1:mx, 1:my, 1:meqn, 1:kmax ) - vector of conserved variables
//
// Returns:
// --------
//
//   Lstar(1:mx, 1:my, 1:meqn, 1:kmax ) - right hand side of MOL formulation
//    smax(1:mx, 1:my, 1:2 )            - maximum speed observed 
//
// ------------------------------------------------------------------------- //
void ConstructL_LLF(
        dTensorBC3& smax,
        dTensorBC3& aux,  
		dTensorBC3& q,
		dTensorBC3& Fm,
        dTensorBC3& Fp,
		dTensorBC3& Gm,
        dTensorBC3& Gp)
{

  const edge_data& edgeData = Legendre2d::get_edgeData();
  const int space_order = 1;// We only need a single flux value//dogParams.get_space_order();
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);

  // If you need access elsewhere to one of these declarations
  // move it into DogSolverCart2.h rather than copying the
  // declaration. (By using a single declaration we can append
  // arguments with default values without modifying existing
  // calls.) -eaj
  //
  void FluxFunc(const dTensor2& xpts, 
		const dTensor2& Q, const dTensor2& Aux, dTensor3& flux);
  void LstarExtra(const dTensorBC4*, const dTensorBC4*, dTensorBC4*);
  void ArtificialViscosity(const dTensorBC4* aux, const dTensorBC4* q, 
			   dTensorBC4* Lstar);

  // Grid information
  const double xlower = dogParamsCart2.get_xlow();
  const double ylower = dogParamsCart2.get_ylow();
  const double dx     = dogParamsCart2.get_dx();
  const double dy     = dogParamsCart2.get_dy();


    void SetBndValues(dTensorBC3& q, dTensorBC3& aux);
    void SetBndValuesX(dTensorBC3& q, dTensorBC3& aux);
    void SetBndValuesY(dTensorBC3& q, dTensorBC3& aux);
    SetBndValues(q, aux);

  // ---------------------------------------------------------
  // Part II: compute inter-element interaction fluxes
  // ---------------------------------------------------------

  // 1-direction: loop over interior edges and solve Riemann problems
  dTensor1 nvec(2);
  nvec.set(1, 1.0e0 );
  nvec.set(2, 0.0e0 );
      RiemannSolver riemannSolver(meqn,maux);
int k;

SetBndValuesX(q, aux);
#pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc); i++)
    {
      dTensor1 Ql(meqn),Qr(meqn);
      dTensor1 Auxl(maux),Auxr(maux);
      dTensor1 Fl(meqn),Fr(meqn);
      dTensor1 xedge(2);


      for (int j=(2-mbc); j<=(my+mbc-1); j++)
        {
	  // "ell" indexes Riemann point along the edge
	  //for (int ell=1; ell<=space_order; ell++)
            {
	      // Riemann data - q
	      for (int m=1; m<=meqn; m++)
                {
		  Ql.set(m, 0.0 );
		  Qr.set(m, 0.0 );

		  k=1;
                    {
		      Ql.fetch(m) += q.get(i-1,j,m);
		      Qr.fetch(m) += q.get(i,j,m);
                    }
                }

	      // Riemann data - aux
	      for (int m=1; m<=maux; m++)
                {
		  Auxl.set(m, 0.0 );
		  Auxr.set(m, 0.0 );

		  k=1;
                    {
		      Auxl.fetch(m) += aux.get(i-1,j,m);
		      Auxr.fetch(m) += aux.get(i,j,m);
                    }
                }

	      // Solve Riemann problem
	      xedge.set(1, xlower + (double(i)-1.0)*dx );
	      xedge.set(2, ylower + (double(j)-0.5)*dy );

	      const double smax_edge
		= riemannSolver.solve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);

	      // TODO: See comments in all ConstructL methods.  This *should* be a 
	      // bug when OpenMP loops are turned on. (-DS)
	      //
	      // If we reverse the role of i and j in the above for loop,
	      // then it will be OK, but the tensors don't run as fast when
	      // that's the case.
	      //
	      //smax.set(i-1, j, 1, Max(dy*smax_edge, smax.get(i-1, j, 1) ) );
	      //smax.set(i,   j, 1, Max(dy*smax_edge, smax.get(i,   j, 1) ) );

	      // Construct fluxes
	      for (int m=1; m<=meqn; m++)
                {
		  Fm.set(i  ,j,m, Fr.get(m) );
		  Fp.set(i-1,j,m, Fl.get(m) );
          if (i==1 && j==20 && isnan(Fm.get(i,j,  m)) ){cout<<i<<" "<<j<<" "<<m<<" "<<Fm.get(i,j,1)<<" "<<Ql.get(m)<<" "<<Qr.get(m)<<" "<<i<<" "<<j<<" Bad pressure21 "<<endl; exit(1);}
                }
            }
        }
    }


  // ---------------------------------------------------------------------- //
  // 2-direction: loop over interior edges and solve Riemann problems
  // ---------------------------------------------------------------------- //
  nvec.set(1, 0.0e0 );
  nvec.set(2, 1.0e0 );
SetBndValuesY(q, aux);
#pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    {

      dTensor1 Ql(meqn),Qr(meqn);
      dTensor1 Auxl(maux),Auxr(maux);
      dTensor1 xedge(2);
      dTensor1 Fl(meqn),Fr(meqn);

      for (int j=(2-mbc); j<=(my+mbc); j++)
	  {
	    // Riemann data - q
	    for (int m=1; m<=meqn; m++)
	      {
		Ql.set(m, 0.0 );
		Qr.set(m, 0.0 );

		k=1;
		  {
		    Ql.fetch(m) += q.get(i,j-1,m);
		    Qr.fetch(m) += q.get(i,j,m);
		  }
	      }

	    // Riemann data - aux
	    for (int m=1; m<=maux; m++)
	      {
		Auxl.set(m, 0.0 );
		Auxr.set(m, 0.0 );

		k=1;
		  {
		    Auxl.fetch(m) += aux.get(i,j-1,m);
		    Auxr.fetch(m) += aux.get(i,j,m);
		  }
	      }

	    // Solve Riemann problem
	    xedge.set(1, xlower + (double(i)-0.5)*dx );
	    xedge.set(2, ylower + (double(j)-1.0)*dy );

	    // TODO: why is this different than what was used earlier in
	    // the code? (-DS)
	    const double smax_edge = riemannSolver.solve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);

	    // TODO: this looks like a potential problem for the upper
	    // pragma statement.  See above comment. (-DS)
	    //
	    //smax.set(i, j-1, 2, Max(dx*smax_edge,smax.get(i, j-1,2) ) );
	    //smax.set(i, j,   2, Max(dx*smax_edge,smax.get(i, j,  2) ) );

	    // Construct fluxes
	    for (int m=1; m<=meqn; m++)
	      {
		Gm.set(i,j,  m, Fr.get(m) );
		Gp.set(i,j-1,m, Fl.get(m) );
        if (i==4 && j==20 && isnan(Gm.get(i,j,  m)) ){cout<<Fm.get(i,j,1)<<" "<<Gm.get(i,j,1)<<" "<<Gm.get(i,j,1)<<" "<<i<<" "<<j<<" Bad pressure "<<endl; exit(1);}
	      }
	  }
    }
    SetBndValues(q, aux);
  // ---------------------------------------------------------------------- //
  // Modify boundary fluxes if necessary
  // ---------------------------------------------------------------------- //
  /*void SetBndFluxes(const dTensorBC4& q, 
		    const dTensorBC4& aux,
		    dTensorBC4& Fm,
		    dTensorBC4& Fp,
		    dTensorBC4& Gm,
		    dTensorBC4& Gp);
  SetBndFluxes(q,aux,Fm,Fp,Gm,Gp);*/

}
