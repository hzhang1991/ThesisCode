#include "dogdefs.h"
#include "stdlib.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "dog_math.h"

// All-purpose routine for computing the L2-projection
// of various functions onto:
//     mopt==0:   the Legendre basis
//     mopt==1:   the derivatives of Legendre basis
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensor2    node(mnodes,1)
//           dTensorBC3 auxin(1-mbc:mnodes+mbc,maux,mpoints)
//           dTensorBC3   qin(1-mbc:mnodes+mbc,meqn,mpoints)
//           dTensorBC3  Fout(1-mbc:mnodes+mbc,mlength,mpoints)
// ---------------------------------------------------------------------

void L2Project(int mopt, int istart, int iend,
	       const dTensor2& node,
	       const dTensorBC3& qin, 
	       const dTensorBC3& auxin,  
	       dTensorBC3& Fout,
	       void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&))
{    
  const int kmax = dogParams.get_space_order();
  const int meqn = qin.getsize(2);
  const int maux = auxin.getsize(2);
  const int mlength = Fout.getsize(2);
  const int mpoints = Fout.getsize(3);  

  int mtmp = iMax(1,mpoints-mopt);
  dTensor1 wgt(mtmp), spts(mtmp);
  dTensor2 phi(mtmp,kmax), phi_x(mtmp,kmax);

  // -----------------
  // Quick error check
  // -----------------
  if (meqn<1 || maux <0 || mpoints<1 || mpoints>6 || mlength<1 
      || mopt < 0 || mopt > 1)
    {
      printf(" Error in L2project.cpp ... \n");
      printf("         meqn = %i\n",meqn);
      printf("         maux = %i\n",maux);
      printf("      mpoints = %i\n",mpoints);
      printf("      mlength = %i\n",mlength);
      printf("       istart = %i\n",istart);
      printf("         iend = %i\n",iend);
      printf("        mopts = %i\n",mopt);
      printf("\n");
      exit(1);
    }

  // ---------------------------------------------
  // Check for trivial case in the case of mopt==1
  // ---------------------------------------------
  if ( mpoints == mopt )
    { Fout.setall(0.); }
  else
    {
      // Set quadrature weights and points
      void SetQuadPts(dTensor1&,dTensor1&);
      SetQuadPts(wgt,spts);

      // Sample basis function at quadrature points
      void SampleBasis(const dTensor1&,dTensor2&);
      SampleBasis(spts,phi);

      // Sample gradient of basis function at quadrature points
      void SampleBasisGrad(const dTensor1&,dTensor2&);
      SampleBasisGrad(spts,phi_x);

      // ----------------------------------
      // Loop over all elements of interest
      // ----------------------------------    
      const double xlow = dogParamsCart1.get_xlow();
      const double dx = dogParamsCart1.get_dx();

#pragma omp parallel for
      for (int i=istart; i<=iend; i++)
        {
	  double xc = xlow + (double(i)-0.5)*dx;

	  // Each of these three items needs to be private to each thread ..
	  dTensor1 xpts(mtmp);
	  dTensor2 qvals(mtmp,meqn);
	  dTensor2 auxvals(mtmp,maux);
	  dTensor2 fvals(mtmp,mlength);

	  // Loop over each quadrature point
	  for (int m=1; m<=mtmp; m++)
            {
	      // grid point x
	      xpts.set( m, xc + 0.5*dx*spts.get(m) );

	      // Solution values (q) at each grid point
	      for (int me=1; me<=meqn; me++)
                {
		  qvals.set(m,me, 0.0 );

		  for (int k=1; k<=mpoints; k++)
                    {
		      qvals.set(m,me, qvals.get(m,me) 
                                + phi.get(m,k) * qin.get(i,me,k) );
                    }
                }

	      // Auxiliary values (aux) at each grid point
	      for (int ma=1; ma<=maux; ma++)
                {
		  auxvals.set(m,ma, 0.0 );

		  for (int k=1; k<=mpoints; k++)
                    {
		      auxvals.set(m,ma, auxvals.get(m,ma) 
				  + phi.get(m,k) * auxin.get(i,ma,k) );
                    }
                }
            }

	  // Call user-supplied function to set fvals
	  Func(xpts,qvals,auxvals,fvals);

	  // Evaluate integrals
	  if (mopt==0) // project onto Legendre basis
            {
	      for (int m1=1; m1<=mlength; m1++)
		for (int m2=1; m2<=mpoints; m2++)
		  {
		    double tmp = 0.0;
		    for (int k=1; k<=mpoints; k++)
		      {
			tmp += wgt.get(k)*fvals.get(k,m1)*phi.get(k,m2);
		      }
		    Fout.set(i,m1,m2, 0.5*tmp );
		  }

            }
	  else // project onto derivatives of Legendre basis
            {
	      for (int m1=1; m1<=mlength; m1++)             
		for (int m2=1; m2<=mpoints; m2++)
		  {
		    double tmp = 0.0;
		    for (int k=1; k<=(mpoints-mopt); k++)
		      {
			tmp += wgt.get(k)*fvals.get(k,m1)*phi_x.get(k,m2);
		      }
		    Fout.set(i,m1,m2, 0.5*tmp );
		  }
            }
        }
    }

}

void SetQuadPts(dTensor1& wgt,
		dTensor1& spts)
{
  const int mpoints = wgt.getsize();

  switch ( mpoints )
    {
    case 1:
      wgt.set(1,  2.0e0 );

      spts.set(1, 0.0e0 );

      break;

    case 2:
      wgt.set(1,   1.0 );
      wgt.set(2,   1.0 );

      spts.set(1, -1.0/sq3 );
      spts.set(2,  1.0/sq3 );

      break;

    case 3:
      wgt.set(1, 5.0e0/9.0e0 );
      wgt.set(2, 8.0e0/9.0e0 );
      wgt.set(3, 5.0e0/9.0e0 );

      spts.set(1, -sq3/sq5 );
      spts.set(2,  0.0e0 );
      spts.set(3,  sq3/sq5 );

      break;

    case 4:

      wgt.set(1, (18.0 - sq30)/36.0 );
      wgt.set(2, (18.0 + sq30)/36.0 );
      wgt.set(3, (18.0 + sq30)/36.0 );
      wgt.set(4, (18.0 - sq30)/36.0 );

      spts.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
      spts.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
      spts.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
      spts.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );

      break;

    case 5:
      wgt.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
      wgt.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
      wgt.set(3,  128.0/225.0 );
      wgt.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
      wgt.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );

      spts.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
      spts.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
      spts.set(3,  0.0 );
      spts.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
      spts.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );

      break;

    case 6:

      wgt.set(1,  0.171324492379170);
      wgt.set(2,  0.360761573048139);
      wgt.set(3,  0.467913934572691);
      wgt.set(4,  0.467913934572691);
      wgt.set(5,  0.360761573048139);
      wgt.set(6,  0.171324492379170);


      spts.set(1, -0.932469514203152);
      spts.set(2, -0.661209386466265);
      spts.set(3, -0.238619186083197);
      spts.set(4, -spts.get(3) );
      spts.set(5, -spts.get(2) );
      spts.set(6, -spts.get(1) );

      break;

    }
}

void SampleBasis(const dTensor1& spts,
		 dTensor2& phi)
{
  const int mpoints = phi.getsize(1);
  const int    kmax = phi.getsize(2);

  switch (kmax)
    {
    case 1:      
      for (int m=1; m<=mpoints; m++)
	{
	  phi.set( m,1, 1.0 );
	}
      break;

    case 2:      
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);

	  phi.set( m,1, 1.0 );
	  phi.set( m,2, sq3*xi );
	}
      break;

    case 3:      
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;

	  phi.set( m,1, 1.0 );
	  phi.set( m,2, sq3*xi );
	  phi.set( m,3, 0.5*sq5*( 3.0*xi2 - 1.0 ) );
	}
      break;

    case 4:      
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;

	  // Legendre basis functions at each grid point
	  phi.set( m,1, 1.0 );
	  phi.set( m,2, sq3*xi );
	  phi.set( m,3, 0.5*sq5*( 3.0*xi2 - 1.0 ) );
	  phi.set( m,4, 0.5*sq7*xi*(5.0*xi2 - 3.0) );
	}
      break;

    case 5:      
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;
	  const double xi4 = xi3*xi;

	  // Legendre basis functions at each grid point
	  phi.set( m,1, 1.0 );
	  phi.set( m,2, sq3*xi );
	  phi.set( m,3, 0.5*sq5*( 3.0*xi2 - 1.0 ) );
	  phi.set( m,4, 0.5*sq7*xi*(5.0*xi2 - 3.0) );
	  phi.set( m,5, 13.125*xi4- 11.25*xi2 + 1.125 );
	}
      break;

    case 6:      
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;
	  const double xi4 = xi3*xi;
	  const double xi5 = xi4*xi;

	  // Legendre basis functions at each grid point
	  phi.set( m,1, 1.0 );
	  phi.set( m,2, sq3*xi );
	  phi.set( m,3, 0.5*sq5*( 3.0*xi2 - 1.0 ) );
	  phi.set( m,4, 0.5*sq7*xi*(5.0*xi2 - 3.0) );
	  phi.set( m,5, 13.125*xi4- 11.25*xi2 + 1.125 );
	  phi.set( m, 6, 7.875*sq11 * (xi5 - 10.0*oneninth*xi3+ 5.0*onethird*oneseventh*xi) );
	}
      break;

    }
}

void SampleBasisGrad(const dTensor1& spts,
		     dTensor2& phi_x)
{
  const double odx = 1.0/dogParamsCart1.get_dx();
  const int mpoints = phi_x.getsize(1);
  const int    kmax = phi_x.getsize(2);

  switch( kmax)
    {
    case 1:
      for (int m=1; m<=mpoints; m++)
	{
	  phi_x.set( m, 1, 0.0 );
	}
      break;

    case 2:
      for (int m=1; m<=mpoints; m++)
	{
	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	}
      break;

    case 3:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	}
      break;

    case 4:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	  phi_x.set( m, 4, 3.0*sq7*(5.0*xi2-1.0)*odx );
	}
      break;

    case 5:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	  phi_x.set( m, 4, 3.0*sq7*(5.0*xi2-1.0)*odx );
	  phi_x.set( m, 5, 15.0*xi* (7.0*xi2-3.0)*odx );
	}
      break;

    case 6:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;
	  const double xi4 = xi3*xi;

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	  phi_x.set( m, 4, 3.0*sq7*(5.0*xi2-1.0)*odx );
	  phi_x.set( m, 5, 15.0*xi* (7.0*xi2-3.0)*odx );
	  phi_x.set( m, 6, (2.0*odx)*7.875*sq11*(5.0*xi4-10.0*onethird*xi2+5.0*onethird*oneseventh) );
	}
      break;

    }
}
