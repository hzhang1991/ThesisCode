#include "../defs.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
// This special routine is only used to set up values for 
// Hermite interpolation in time for SubGrid boundary conditions
//
void ConstructL_HMM_short(int istart, int iend, int method[], dTensor2 node,
			  dTensorBC3 aux, dTensorBC3 q, dTensorBC3& Lstar,
			  void (*FluxFunc)(dTensor1,dTensor2,dTensor2,dTensor2&),
			  void (*SetWaveSpd)(dTensor1,dTensor2,dTensor2,dTensor2,
					     dTensor2,double&,double&),
			  void (*SourceTermFunc)(dTensor1,dTensor2,dTensor2,dTensor2&))
{
    int i,k,m;
    double tmp,z2,z3,xc,xl,dx,xlower;
    double smax_edge = 0.0e0;
    int melems = q.getsize(1);
    int   meqn = q.getsize(2);
    int   kmax = q.getsize(3);
    int   maux = aux.getsize(2);
    int    mbc = q.getmbc();
    dTensor1 xedge(1);
    dTensorBC2 Fm(melems,meqn,mbc);
    dTensorBC2 Fp(melems,meqn,mbc);
    dTensorBC3 N(melems,meqn,kmax,mbc);
    dTensorBC3 Psi(melems,meqn,kmax,mbc);
    dTensor2 Ql(1,meqn),Qr(1,meqn),Auxl(1,maux),Auxr(1,maux);
    dTensor2 Fl(1,meqn),Fr(1,meqn);
    double Max(double,double);
    double RiemannSolve(dTensor1,dTensor2,dTensor2,dTensor2,dTensor2,
			dTensor2&,dTensor2&,
			void (*Func1)(dTensor1,dTensor2,dTensor2,dTensor2&),
			void (*Func2)(dTensor1,dTensor2,dTensor2,dTensor2,
				      dTensor2,double&,double&));
    void SetBndValues(dTensor2,dTensorBC3&,dTensorBC3&);
    void L2Project(int,int,int,dTensor2,dTensorBC3,dTensorBC3,dTensorBC3&,
		   void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&));

    // Grid spacing
    xlower = node.get(1,1);
    dx = node.get(2,1) - node.get(1,1);

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    SetBndValues(node,aux,q);
        
    // Loop over interior edges and solve Riemann problems
    for (i=istart; i<=(iend+1); i++)
    {
        // Riemann data - q
        for (m=1; m<=meqn; m++)
        {
	    Ql.set(1,m, 0.0 );
	    Qr.set(1,m, 0.0 );
	    
	    for (k=1; k<=method[1]; k++)
	    {
	        Ql.set(1,m, Ql.get(1,m) + sqrt(2.0*double(k)-1.0)
		       *q.get(i-1,m,k) );
		Qr.set(1,m, Qr.get(1,m) + pow(-1.0,k+1)
		       *sqrt(2.0*double(k)-1.0)*q.get(i,m,k) ); 
	    }
        }
        
        // Riemann data - aux
        for (m=1; m<=maux; m++)
        {
	    Auxl.set(1,m, 0.0 );
	    Auxr.set(1,m, 0.0 );

	    for (k=1; k<=method[1]; k++)
	    {
	      Auxl.set(1,m, Auxl.get(1,m) + sqrt(2.0*double(k)-1.0)
		       *aux.get(i-1,m,k) );
	      Auxr.set(1,m, Auxr.get(1,m) + pow(-1.0,k+1)
		       *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) ); 
	    }
        }
            
        // Solve Riemann problem
	xedge.set(1, xlower + (double(i)-1.0)*dx );
	smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
				 FluxFunc,SetWaveSpd);
        
        // Construct fluxes
        for (m=1; m<=meqn; m++)
        {
	    Fm.set(i,  m, Fr.get(1,m) );
            Fp.set(i-1,m, Fl.get(1,m) );
        }
    }
    // ---------------------------------------------------------
    
    
    // ---------------------------------------------------------
    // Part II: compute intra-element contributions
    // ---------------------------------------------------------
    //
    //   N = int( f(q,x,t) * phi_x, x )/dx
    //
    // Compute ``N'' by projecting flux function onto the 
    // gradient of Legendre polynomials
    if (method[1]>1)
      {  L2Project(1,istart,iend,node,q,aux,N,FluxFunc);  }
    else
    {
	for (i=istart; i<=iend; i++)
	{   
	    for (m=1; m<=meqn; m++)
	    {
		N.set(i,m,1, 0.0 );   
	    }
	}
    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part III: compute source term
    // --------------------------------------------------------- 
    if ( method[7]>0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        L2Project(0,istart,iend,node,q,aux,Psi,SourceTermFunc);
    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part IV: construct Lstar
    // ---------------------------------------------------------
    if (method[7]==0)  // Without Source Term
    { 
	for (i=istart; i<=iend; i++)
	{
	    for (m=1; m<=meqn; m++)
	    {
		for (k=1; k<=method[1]; k++)
		{
		    tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
		      ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx;
		    
		    Lstar.set(i,m,k, tmp );
		}
	    }
	}
    }
    else  // With Source Term
    {
	for (i=istart; i<=iend; i++)
	{
	    for (m=1; m<=meqn; m++)
	    {
		for (k=1; k<=method[1]; k++)
		{
		    tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
		      ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx
		      + Psi.get(i,m,k);
		    
		    Lstar.set(i,m,k, tmp );
		}
	    }
	}
    }
    // ---------------------------------------------------------

}
