#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"
#include "DogParamsCart1.h"
using namespace std;

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
  const int KMAX = 6;
  const int meqn = qin.getsize(2);
  const int maux = auxin.getsize(2);
  const int mlength = Fout.getsize(2);
  const int mpoints = Fout.getsize(3);  
  dTensor1 wgt(mpoints), spts(mpoints);
  dTensor2 phi(mpoints,KMAX), phi_x(mpoints,KMAX);
  
  // -----------------
  // Quick error check
  // -----------------
  if (meqn<1 || maux <0 || mpoints<1 || mpoints>KMAX || mlength<1 
      || mopt < 0 || mopt > 1)
    {
      cout << " Error in L2project.cpp ... " << endl;
      cout << "         meqn = " << meqn << endl;
      cout << "         maux = " << maux << endl;
      cout << "      mpoints = " << mpoints << endl;
      cout << "      mlength = " << mlength << endl;
      cout << "       istart = " << istart << endl;
      cout << "         iend = " << iend << endl;
      cout << "        mopts = " << mopt << endl;
      cout << endl;
      exit(1);
    }
  
  // ---------------------------------------------
  // Check for trivial case in the case of mopt==1
  // ---------------------------------------------
  if ( mpoints == mopt )
    { Fout.setall(0.); }
  else
    {
      // ---------------------------------
      // Set quadrature weights and points
      // ---------------------------------
      switch ( (mpoints-mopt) )
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
      
      // Loop over each quadrature point to construct Legendre polys
      const double dx = dogParamsCart1.get_dx();//node.get(2,1) - node.get(1,1);
      for (int m=1; m<=(mpoints-mopt); m++)
        {
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;
	  const double xi4 = xi3*xi;
	  const double xi5 = xi4*xi;
	  
	  // Legendre basis functions at each grid point
	  phi.set( m,1, 1.0 );
	  phi.set( m,2, sq3*xi );
	  phi.set( m,3, 0.5*sq5*( 3.0*(xi2) - 1.0 ) );
	  phi.set( m,4, 0.5*sq7*xi
		   *(5.0*(xi2) - 3.0) );
	  phi.set( m,5, (105.0/8.0)*(xi4) 
		   - (45.0/4.0)*(xi2) + (9.0/8.0) );
	  phi.set( m, 6, (63.0/8.0)*sq11 * ( (xi5) - (10.0/9.0)*(xi3)+ (5.0/21.0)*xi ) );
	  
	  // 1st derivative of Legendre basis functions at each grid point
	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3/dx );
	  phi_x.set( m, 3, 6.0*sq5*xi/dx );
	  phi_x.set( m, 4, 3.0*sq7*(5.0*(xi2)-1.0)/dx );
	  phi_x.set( m, 5, 15.0*xi* (7.0*(xi2)-3.0)/dx );
	  phi_x.set( m, 6, (2.0/dx)*(63.0/8.0)*sq11*(5.0*(xi4)-(10.0/3.0)*(xi2)+5.0/21.0) );
	  
        }
      
      // ----------------------------------
      // Loop over all elements of interest
      // ----------------------------------    
      const double s_area = 2.0;
      const double xlow = dogParamsCart1.get_xlow();//node.get(1,1);
      
#pragma omp parallel for
      for (int i=istart; i<=iend; i++)
        {
	  double xc = xlow + (double(i)-0.5)*dx;
	  
	  // each of these three items needs to be private to each thread ..
	  dTensor1 xpts(mpoints-mopt);
	  dTensor2 qvals(mpoints-mopt,meqn), auxvals(mpoints-mopt,maux);
	  dTensor2 fvals(mpoints-mopt,mlength);
	  
	  // Loop over each quadrature point
	  for (int m=1; m<=(mpoints-mopt); m++)
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
			tmp = tmp + wgt.get(k)*fvals.get(k,m1)
			  *phi.get(k,m2);
		      }
		    Fout.set(i,m1,m2, tmp/s_area );         
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
			tmp = tmp + wgt.get(k)*fvals.get(k,m1)
			  *phi_x.get(k,m2);
		      }
		    Fout.set(i,m1,m2, tmp/s_area );
		  }
            }
        }
    }
}
