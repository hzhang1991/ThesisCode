#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "tensors.h"
#include "dog_math.h"
#include "mesh.h"
using namespace std;

// Right-hand side of the Lax-Wendroff discretization:
//
//      ( q(t+dt) - q(t) )/dt = F
//
// where the the flux is given by
//
//      F = F1 + dt/2 * F2 + dt^2/6 * F3 + dt^3/24 * F4 + dt^4/120 * F5
//
//        F1 = -u   * q_{x}
//        F2 =  u^2 * q_{xx}
//        F3 = -u^3 * q_{xxx}
//        F4 =  u^4 * q_{xxxx} 
//        F5 = -u^5 * q_{xxxxx}
//
void LaxWendroff(const double dt, const mesh& Mesh,
		 dTensorBC3& aux, dTensorBC3& q, 
		 dTensorBC3& Lstar, dTensorBC1& smax)
{
  /*
  void EvalRiemmanData(const int& i, const dTensorBC3& q, const dTensorBC3& aux,  
		       const dTensorBC3& LxWF, dTensor1& Ql, dTensor1& Qr, 
		       dTensor1& LxFFluxl, dTensor1& LxFFluxr,
		       dTensor1& Auxl, dTensor1& Auxr);
  int i,k,m;
  double tmp,z2,z3,xc,xl,dx,xlower;
  double smax_edge = 0.0e0;
  double nu = 0.0;
  double s1,s2;
  int melems = q.getsize(1);
  int   meqn = q.getsize(2);
  int   kmax = q.getsize(3);
  int   maux = aux.getsize(2);
  int    mbc = q.getmbc();
  dTensor1 xedge(1);
  dTensorBC1 nu_vec(melems,mbc);
  dTensorBC3 Fm(melems,meqn,3,mbc);
  dTensorBC3 Fp(melems,meqn,3,mbc);
  dTensorBC3 N(melems,meqn,kmax,mbc);
  dTensorBC3 Psi(melems,meqn,kmax,mbc);
  dTensor1 Ql(meqn),Qr(meqn),Auxl(maux),Auxr(maux);
  dTensor1 LxFFluxl1(meqn), LxFFluxr1(meqn);
  dTensor1 LxFFluxl2(meqn), LxFFluxr2(meqn);
  dTensor1 LxFFluxl3(meqn), LxFFluxr3(meqn);
  dTensor1 LxFFluxl(meqn),  LxFFluxr(meqn);
  dTensor1 Fl(meqn),Fr(meqn);
  dTensor1 Fl1(meqn),Fr1(meqn);
  dTensor1 Fl2(meqn),Fr2(meqn);
  dTensor1 Fl3(meqn),Fr3(meqn);
 
  // lax-wendroff flux function
  dTensorBC3  LxWF(melems,meqn,kmax,mbc);
  
  // these are necessary for stability corrections.  they are the
  // lax-wendroff flux function broken up into dt^0, dt^1, and dt^2
  // components.
  dTensorBC3 LxWF1(melems,meqn,kmax,mbc);
  dTensorBC3 LxWF2(melems,meqn,kmax,mbc);
  dTensorBC3 LxWF3(melems,meqn,kmax,mbc);
  dTensorBC3 IntF1(melems,meqn,kmax,mbc);
  dTensorBC3 IntF2(melems,meqn,kmax,mbc);
  dTensorBC3 IntF3(melems,meqn,kmax,mbc);
  
  
  double RiemannSolveLxW(const dTensor1 xedge,
			 const dTensor1& Ql, const dTensor1& Qr, const dTensor1& Auxl,
			 const dTensor1& Auxr, const dTensor1& ffl, const dTensor1& ffr,
			 dTensor1& Fl, dTensor1& Fr,
			 void (*SetWaveSpd)(dTensor1,dTensor1,dTensor1,dTensor1,
					    dTensor1,double&,double&));
  double RiemannSolveLxWupwind(const dTensor1 xedge,
			       const dTensor1& Ql, const dTensor1& Qr, const dTensor1& Auxl,
			       const dTensor1& Auxr, const dTensor1& ffl, const dTensor1& ffr,
			       dTensor1& Fl, dTensor1& Fr,
			       void (*SetWaveSpd)(dTensor1,dTensor1,dTensor1,dTensor1,
						  dTensor1,double&,double&));
  void SetBndValues(dTensor2,dTensorBC3&,dTensorBC3&);
  void L2Project(int,int,int,dTensor2,dTensorBC3,dTensorBC3,dTensorBC3&,
		 void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&));
  void L2ProjectLxW(const int stab_cor, const int istart, const int iend, 
		    const dTensor2 node, const dTensorBC3 qin, const dTensorBC3 auxin,  
		    dTensorBC3& F1, dTensorBC3& F2, dTensorBC3& F3,
		    dTensorBC3& LxWFlux1, dTensorBC3& LxWFlux2, dTensorBC3& LxWFlux3,
		    void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&),
		    void (*DFunc)(dTensor1, dTensor2, dTensor2, dTensor3&),
		    void (*D2Func)(dTensor1, dTensor2, dTensor2, dTensor4&));
  void FluxFunc(dTensor1,dTensor2,dTensor2,dTensor2&);
  void DFluxFunc (dTensor1,dTensor2,dTensor2,dTensor3&);
  void D2FluxFunc(dTensor1,dTensor2,dTensor2,dTensor4&);
  void SetWaveSpd(dTensor1,dTensor1,dTensor1,dTensor1,
		  dTensor1,double&,double&);
  void SourceTermFunc(dTensor1,dTensor2,dTensor2,dTensor2&);
  void LstarExtra(const dTensor2&,dTensorBC3&,dTensorBC3&,dTensorBC3&);
  
  // quick error check
  if( method[7] > 0 )
    {
      cout << "    error: have not implemented source term for LxW solver "
	   << endl;
      exit(1);
    }
  
  // Grid spacing
  xlower = node.get(1,1);
  dx = node.get(2,1) - node.get(1,1);

  // Number of ghost cells
  mbc = q.getmbc();
  
  // double check to make sure these were initialized to zero...
  for(k=1; k<= kmax; k++)
    for(i=1-mbc; i<=(melems+mbc); i++)      
      for(m=1; m<=meqn; m++)
	{
	  LxWF1.set(i,m,k,0.0e0);
	  LxWF2.set(i,m,k,0.0e0);
	  LxWF3.set(i,m,k,0.0e0);
	  IntF1.set(i,m,k,0.0e0);
	  IntF2.set(i,m,k,0.0e0);
	  IntF3.set(i,m,k,0.0e0);
	}

  // Boundary conditions
  SetBndValues(node,aux,q);

  // compute necessary coefficients
  //L2Project(1,1-mbc,melems+mbc,node,q,aux,N,&FluxFunc);
  L2ProjectLxW(stab_cor, 1-mbc, melems+mbc, node, q, aux,  
	       IntF1, IntF2, IntF3,  LxWF1,  LxWF2,  LxWF3,
	       &FluxFunc, &DFluxFunc, &D2FluxFunc);
  
  // ---------------------------------------------------------
  // Part I: compute inter-element interaction fluxes
  // ---------------------------------------------------------
  
  // Boundary conditions
  SetBndValues(node,aux,q);
  
  // evaluate the Riemman data for the very first and very last cell
  if( stab_cor == 0 )
    {
      nu_vec.set(1-mbc, 0.0);
    } 
  else
    {
      EvalRiemmanData(2-mbc, q, aux, LxWF, Ql, Qr, LxFFluxl, LxFFluxr, Auxl, Auxr); 
      SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1, s2);
      nu = fabs(s1) * dt / dx;  // can get away with just using s1 when meqn == 1
      nu_vec.set(1-mbc, nu);
      // cout << "   nu = " << nu << "   for i = " << 1-mbc << endl;
      //        EvalRiemmanData(melems+mbc, q, aux, LxWF, Ql, Qr, LxFFluxl, LxFFluxr, Auxl, Auxr); 
      //        SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1, s2);
      //        nu = fabs(s1) * dt / dx;  // can get away with just using s1 when meqn == 1
      //        nu_vec.set(1+mbc, nu);      
    }
  // Loop over interior edges and solve Riemann problems
  for (i=(2-mbc); i<=(melems+mbc); i++)
    {
      xedge.set(1, xlower + (double(i)-1.0)*dx );
      if( meqn == 1 && stab_cor != 0)  // stab corrections only work when meqn == 0
        {	  
	  // evaluate the Riemann Data for edge x_{i-1/2} and
	  // solve the riemann problem.
	  switch(kmax)  // use cascading switch statement to reduce number of riemann solves
            {
	    case 3:
	      EvalRiemmanData(i, q, aux, LxWF3, Ql, Qr, LxFFluxl3, LxFFluxr3, Auxl, Auxr);
	      smax_edge = RiemannSolveLxWupwind(xedge,Ql,Qr,Auxl,Auxr,
						LxFFluxl3, LxFFluxr3, Fl3, Fr3, &SetWaveSpd);
	    case 2:
	      EvalRiemmanData(i, q, aux, LxWF2, Ql, Qr, LxFFluxl2, LxFFluxr2, Auxl, Auxr);
	      tmp = RiemannSolveLxWupwind(xedge,Ql,Qr,Auxl,Auxr,
					  LxFFluxl2, LxFFluxr2, Fl2, Fr2, &SetWaveSpd);
	      smax_edge = Max( smax_edge, tmp );
	    case 1:
	      EvalRiemmanData(i, q, aux, LxWF1, Ql, Qr, LxFFluxl1, LxFFluxr1, Auxl, Auxr);
	      tmp = RiemannSolveLxWupwind(xedge,Ql,Qr,Auxl,Auxr,
					  LxFFluxl1, LxFFluxr1, Fl1, Fr1, &SetWaveSpd);
	      smax_edge = Max( smax_edge, tmp );
            }
        } 
      else // meqn > 1 or no stability corrections
        {

	  // evaluate coefficients of LaxWendroff Flux Function 
	  for(k=1; k<= kmax; k++)	    
	    for(m=1; m<=meqn; m++)
	      {
		tmp = LxWF1.get(i,m,k)
		  + LxWF2.get(i,m,k) * dt / 2.0
		  + LxWF3.get(i,m,k) * pow(dt,2) / 6.0;
		LxWF.set(i,m,k, tmp );
	      }

	  // evaluate the Riemman data
	  EvalRiemmanData(i, q, aux, LxWF, Ql, Qr, LxFFluxl, LxFFluxr, Auxl, Auxr); 
	  smax_edge = RiemannSolveLxW(xedge,Ql,Qr,Auxl,Auxr,
				      LxFFluxl, LxFFluxr, Fl, Fr, &SetWaveSpd);
	  //smax_edge = RiemannSolveLxWupwind(xedge,Ql,Qr,Auxl,Auxr,
	  //        LxFFluxl, LxFFluxr, Fl, Fr, &SetWaveSpd);
        }

      smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
      smax.set(i,   Max(smax_edge,smax.get(i)) );
      
      // stability correction term.  (this is used for integrals too, so 
      // need to save as a vector)
      // Calculate minimum and maximum HLLE speeds
      SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1, s2);
      nu = fabs(s1) * dt / dx;  // can get away with just using one value when meqn == 1

      // term for integration corrections
      nu_vec.set(i, 0.5 * (nu + nu_vec.get(i-1)) );
      //        nu_vec.set(i, nu );
      
      // Construct fluxes
      if( meqn == 1 && stab_cor != 0)
        {
	  for (m=1; m<=meqn; m++)
            {
	      
	      switch( kmax )
                {
		case 3:
		  
		  // k == 1 polynomial
		  tmp = Fr1.get(m) 
		    + Fr2.get(m) * dt / 2.0 
		    + Fr3.get(m) *  pow(dt,2) / 6.0;
		  Fm.set(i, m, 1, tmp );
		  
		  // k == 2 polynomial
		  if( s1 > 0 )  // positive velocity
		    {
		      tmp = (1.0e0 - nu)         * Fl1.get(m) 
                                + (1.0e0 - 2.0/3.0*nu) * Fl2.get(m) * dt / 2.0
			+ (1.0e0 -0.5*nu)      * Fl3.get(m) * pow(dt,2) / 6.0;
		    }
		  else // negative velocity
		    {
		      tmp = (1.0e0 - nu)             * Fr1.get(m) 
			+ (1.0e0 - 2.0/3.0*nu)     * Fr2.get(m) * dt / 2.0
			+ (1.0e0 -0.5*nu - 1.0/nu) * Fr3.get(m) * pow(dt,2) / 6.0;
		    }
		  Fm.set(i, m, 2, tmp );
		  
		  // k == 2 polynomial
		  if( s1 <=  0 )  // negative velocity
		    {
		      tmp = (1.0e0 - nu)         * Fr1.get(m) 
			+ (1.0e0 - 2.0/3.0*nu) * Fr2.get(m) * dt / 2.0
			+ (1.0e0 -0.5*nu)      * Fr3.get(m) * pow(dt,2) / 6.0;
		    }
		  else // positive velocity
		    {
		      tmp = (1.0e0 - nu)             * Fl1.get(m) 
			+ (1.0e0 - 2.0/3.0*nu)     * Fl2.get(m) * dt / 2.0
			+ (1.0e0 -0.5*nu - 1.0/nu) * Fl3.get(m) * pow(dt,2) / 6.0;
		    }
		  Fp.set(i-1, m, 2, tmp );
		  
		  // k == 3 polynomial
		  tmp = (1.0e0 - 3.0*nu + 2.0*pow(nu,2) ) * Fr1.get(m) 
		    + (1.0e0 - 2.0*nu + pow(nu,2))      * Fr2.get(m) * dt / 2.0
		    + (1.0e0 - 3.0/2.0*nu+3.0/5.0*pow(nu,2)) 
		    * Fr3.get(m) * pow(dt,2) / 6.0;
		  Fm.set(i, m, 3, tmp );
		  
		  // k == 3 polynomial
		  tmp = (1.0e0 - 3.0*nu + 2.0*pow(nu,2) ) * Fl1.get(m) 
		    + (1.0e0 - 2.0*nu + pow(nu,2))      * Fl2.get(m) * dt / 2.0
		    + (1.0e0 - 3.0/2.0*nu+3.0/5.0*pow(nu,2)) 
		    * Fl3.get(m) * pow(dt,2) / 6.0;
		  Fp.set(i-1, m, 3, tmp );
		  
		  // k == 1 polynomial
		  tmp = Fl1.get(m) 
		    + Fl2.get(m) * dt / 2.0
		    + Fl3.get(m) * pow(dt,2) / 6.0;
		  Fp.set(i-1, m, 1, tmp );

		  break;
		  
		case 2:
		  // k == 1 polynomial
		  tmp = Fr1.get(m) 
		    + Fr2.get(m) * dt / 2.0;
		  Fm.set(i, m, 1, tmp );
		  
		  // k == 2 polynomial
		  tmp = (1.0e0 - nu)       * Fr1.get(m) 
		    + (1.0e0-2.0/3.0*nu) * Fr2.get(m) * dt / 2.0;
		  Fm.set(i, m, 2, tmp );
		  
		  // k == 1 polynomial
		  tmp = Fl1.get(m) 
		    + Fl2.get(m) * dt / 2.0;
		  Fp.set(i-1, m, 1, tmp );
		  
		  // k == 2 polynomial
		  tmp = (1.0e0 - nu)       * Fl1.get(m) 
		    + (1.0e0-2.0/3.0*nu) * Fl2.get(m) * dt / 2.0;
		  Fp.set(i-1, m, 2, tmp );
		  
		  break;

		case 1:
		  
		  // k == 1 polynomial minus
		  tmp = Fr1.get(m) ;
		  Fm.set(i, m, 1, tmp );
		  
		  // k == 1 polynomial plus
		  tmp = Fl1.get(m); 
		  Fp.set(i-1, m, 1, tmp );

                }

            } // end of looping over each eqn
        } 
      else // meqn > 1 ( have not tried to implement corrections! )
        {
	  for (m=1; m<=meqn; m++)
	    for (k=1; k<=kmax; k++)
	      {
		Fm.set(i,  m, k, Fr.get(m) );
		Fp.set(i-1,m, k, Fl.get(m) );
	      }
	} // end of looping over each meqn
        
    }
  SetBndValues(node,aux,q);  // clean up boundary values
  // ---------------------------------------------------------
  
  // ---------------------------------------------------------
  // Part II: compute intra-element contributions
  // ---------------------------------------------------------
  //
  //   N = int( F(q,x,t) * phi_x, x )/dx
  //
  // Compute ``N'' by projecting flux function onto the 
  // gradient of Legendre polynomials
  if (method[1]>1)
    {
      if( meqn == 1 && stab_cor != 0)  
        { // stability corrections turned on
	  for (i=1-mbc; i<=(melems+mbc); i++)
            {
	      
	      for (m=1; m<=meqn; m++)
                {	  
		  nu = nu_vec.get(i);
		  switch(kmax)
                    {
		    case 3:
		      
		      // k == 1 polynomial
		      //                           tmp = IntF1.get(i,m,1)
		      //                                + IntF2.get(i,m,1) * dt / 2.0 
		      //                                + IntF3.get(i,m,1) *  pow(dt,2) / 6.0;
		      //                            N.set(i, m, 1, tmp );
		      N.set(i, m, 1, 0.0 );
		      
		      // k == 2 polynomial
		      tmp = (1.0e0 - nu)                * IntF1.get(i,m,2) 
			+ (1.0e0 - 2.0/3.0*nu)        * IntF2.get(i,m,2) * dt / 2.0
			+ (1.0e0 - 0.5*nu - 1.0 / nu) * IntF3.get(i,m,2) * pow(dt,2) / 6.0;
		      N.set(i, m, 2, tmp );
		      
		      // k == 3 polynomial
		      tmp = (1.0e0 - 3.0*nu + 2.0*pow(nu,2))
			* IntF1.get(i,m,3) 
			+ (1.0e0 - 2.0*nu + pow(nu,2) )
			* IntF2.get(i,m,3) * dt / 2.0
			+ (1.0e0 - 3.0/2.0 * nu + 3.0/5.0 * pow(nu,2) )
			* IntF3.get(i,m,3) * pow(dt,2) / 6.0;
		      N.set(i, m, 3, tmp );
		      break;
		      
		    case 2:
		      
		      // k == 1 polynomial
		      //                            tmp = IntF1.get(i,m,1)
		      //                                + IntF2.get(i,m,1) * dt / 2.0;
		      N.set(i, m, 1, 0.0e0 );
		      
		      // k == 2 polynomial
		      tmp = (1.0e0 - nu)        * IntF1.get(i,m,2) 
			+ (1.0e0-2.0/3.0*nu)  * IntF2.get(i,m,2) * dt / 2.0;
		      N.set(i, m, 2, tmp );
		      
		      break;

		    case 1:

		      // k == 1 polynomial
		      tmp = IntF1.get(i,m,1);
		      N.set(i, m, 1, 0.0e0 );
		      break;
                    }

                } // end of looping over each eqn
            }
        }
      else // meqn > 1 or stability corrections turned off
        {
	  for (i=1-mbc; i<=(melems+mbc); i++)           
	    for(k=1; k<=kmax; k++)
	      for (m=1; m<=meqn; m++)
		{
		  tmp = IntF1.get(i,m,k) 
		    + IntF2.get(i,m,k) * dt/2.0
		    + IntF3.get(i,m,k) * pow(dt,2)/6.0;
		  N.set(i,m,k, tmp );
		}
        }    
    }
  else // case kmax == 1 - there is no integral term in this case
    {
      for (i=1-mbc; i<=(melems+mbc); i++)        
	for (m=1; m<=meqn; m++)
	  {
	    N.set(i,m,1, 0.0 );   
	  }    
    }
  SetBndValues(node,aux,q);   // clean up boundary values
  // ---------------------------------------------------------
  
    
  // ---------------------------------------------------------
  // Part III: compute source term
  // --------------------------------------------------------- 
  if ( method[7]>0 )
    // TODO this needs some work - needs to be expanded in taylor series as
    // well....
    {   
      cout << "    Warning: LxW solver needs to expand source term." << endl;
      cout << " method will be first order for source term function" << endl;
      
      // Set source term on computational grid
      // Set values and apply L2-projection
      L2Project(0,1-mbc,melems+mbc,node,q,aux,Psi,&SourceTermFunc);
    }
  // ---------------------------------------------------------


  // ---------------------------------------------------------
  // Part IV: construct Lstar
  // ---------------------------------------------------------
  if (method[7]==0)  // Without Source Term
    { 
      for (i=(2-mbc); i<=(melems+mbc-1); i++)
        for (m=1; m<=meqn; m++)
	  for (k=1; k<=method[1]; k++)
	    {
	      tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                ( Fp.get(i,m,k) + pow(-1.0,k)*Fm.get(i,m,k) )/dx;
	      
	      Lstar.set(i,m,k, tmp );
	    }
    }
  else  // With Source Term
    {
      for (i=(2-mbc); i<=(melems+mbc-1); i++)
        for (m=1; m<=meqn; m++)
	  for (k=1; k<=method[1]; k++)
	    {
	      tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                ( Fp.get(i,m,k) + pow(-1.0,k)*Fm.get(i,m,k) )/dx
                + Psi.get(i,m,k);
	      
	      Lstar.set(i,m,k, tmp );
	    }
    }
  // ---------------------------------------------------------

  // ---------------------------------------------------------
  // Part V: add extra contributions to Lstar
  // ---------------------------------------------------------
  // Call LstarExtra
  LstarExtra(node,aux,q,Lstar);
  // ---------------------------------------------------------
  */
}
