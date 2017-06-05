#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"
#include "DogAsynch1d.h"
using namespace std;

// All-purpose routine for computing the L2-projection
// of various functions onto:
//     mopt==0:   the Legendre basis
//     mopt==1:   the derivatives of Legendre basis
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensor2    dogAsynch1d.node(mx+1,1)
//           dTensorBC3  dogAsynch1d.aux(1-mbc:mx+mbc,maux,mpoints)
//           dTensorBC3  dogAsynch1d.q(1-mbc:mx+mbc,meqn,mpoints)
//           dTensorBC3  Fout(1-mbc:mx+mbc,mlength,mpoints)
// ---------------------------------------------------------------------

void L2Project(int mopt, int istart, int iend, const DogAsynch1d& dogAsynch1d, dTensorBC3* Fout,
	       void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&))
{    
  int i,k,m,me,ma,m1,m2;
  int meqn = dogAsynch1d.get_meqn();
  int maux = dogAsynch1d.get_maux();
  int mx = Fout->getsize(1);  assert(mx==dogAsynch1d.get_mx());
  int mlength = Fout->getsize(2);
  int mpoints = Fout->getsize(3);  
  double xc,s_area,tmp;
  dTensor1 wgt(mpoints), spts(mpoints), xpts(mpoints);
  dTensor2 phi(mpoints,5), phi_x(mpoints,5);
  dTensor2 qvals(mpoints,meqn), auxvals(mpoints,maux);
  dTensor2 fvals(mpoints,mlength);
  
  // -----------------
  // Quick error check
  // -----------------
  if (meqn<1 || maux <1 || mpoints<1 || mpoints>5 || mlength<1 
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
    {
      for (i=istart; i<=iend; i++)
        {
	  for (m=1; m<=mlength; m++)
            {
	      Fout->set(i,m,1, 0.0 );
            }
        }
    }
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
	    wgt.set(1, (18.0 - sqrt(30.0))/36.0 );
	    wgt.set(2, (18.0 + sqrt(30.0))/36.0 );
	    wgt.set(3, (18.0 + sqrt(30.0))/36.0 );
	    wgt.set(4, (18.0 - sqrt(30.0))/36.0 );
	    
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
      }
      
      // Loop over each quadrature point to construct Legendre polys
      for (m=1; m<=(mpoints-mopt); m++)
	{
	  // Legendre basis functions at each grid point
	  phi.set( m,1, 1.0 );
	  phi.set( m,2, sq3*spts.get(m) );
	  phi.set( m,3, 0.5*sq5*( 3.0*pow(spts.get(m),2) - 1.0 ) );
	  phi.set( m,4, 0.5*sq7*spts.get(m)
		   *(5.0*pow(spts.get(m),2) - 3.0) );
	  phi.set( m,5, (105.0/8.0)*pow(spts.get(m),4) 
		   - (45.0/4.0)*pow(spts.get(m),2) + (9.0/8.0) );
	  
	  // 1st derivative of Legendre basis functions at each grid point
	  // (These are derivatives with respect to xi, the variable that
	  //  runs from -1 to 1. In order to properly compute an x-derivative,
	  //  will need to multiply by 2/dx(i) for each element i.)
	  phi_x.set( m,1, 0.0 );
	  phi_x.set( m,2, sq3 );
	  phi_x.set( m,3, 3.0*sq5*spts.get(m) );
	  phi_x.set( m,4, 1.5*sq7*(5.0*pow(spts.get(m),2)-1.0) );
	  phi_x.set( m,5, 7.5*spts.get(m)*(7.0*pow(spts.get(m),2)-3.0) );
	}
      
      // ----------------------------------
      // Loop over all elements of interest
      // ----------------------------------	
      s_area = 2.0;
      
      double xleft = dogAsynch1d.node->get(1);
      double xright = xleft;
      for (i=istart; i<=iend; i++)
	{
	  xleft  = xright;
	  xright = dogAsynch1d.node->get(i+1);
	  xc = 0.5*(xleft+xright);
	  
	  // Loop over each quadrature point
	  for (m=1; m<=(mpoints-mopt); m++)
	    {
	      // grid point x
	      xpts.set( m, xc + 0.5*dogAsynch1d.dx->get(i)*spts.get(m) );
	      
	      // Solution values (q) at each grid point
	      for (me=1; me<=meqn; me++)
		{
		  qvals.set(m,me, 0.0 );
		  
		  for (k=1; k<=mpoints; k++)
		    {
		      qvals.set(m,me, qvals.get(m,me) 
				+ phi.get(m,k) * dogAsynch1d.q->get(i,me,k) );
		    }
		}
	      
	      // Auxiliary values (aux) at each grid point
	      for (ma=1; ma<=maux; ma++)
		{
		  auxvals.set(m,ma, 0.0 );
		  
		  for (k=1; k<=mpoints; k++)
		    {
		      auxvals.set(m,ma, auxvals.get(m,ma) 
				  + phi.get(m,k) * dogAsynch1d.aux->get(i,ma,k) );
		    }
		}
	      
	      // Call user-supplied function to set fvals
	      Func(xpts,qvals,auxvals,fvals);
	    }

	  // Evaluate integrals
	  if (mopt==0) // project onto Legendre basis
	    {
	      for (m1=1; m1<=mlength; m1++)
		for (m2=1; m2<=mpoints; m2++)
		  {		
		    tmp = 0.0;
		    for (k=1; k<=mpoints; k++)
		      {
			tmp = tmp + wgt.get(k)*fvals.get(k,m1)
			  *phi.get(k,m2);
		      }
		    Fout->set(i,m1,m2, tmp/s_area );
		  }
	    }
	  else // project onto derivatives of Legendre basis
	    {
	      for (m1=1; m1<=mlength; m1++)		
		for (m2=1; m2<=mpoints; m2++)
		  {
		    tmp = 0.0;
		    for (k=1; k<=(mpoints-mopt); k++)
		      {
			tmp = tmp + wgt.get(k)*fvals.get(k,m1)
			  *phi_x.get(k,m2)*(2.0/dogAsynch1d.dx->get(i));
		      }
		    Fout->set(i,m1,m2, tmp/s_area );
		  }		
	    }
	}
    }

}
