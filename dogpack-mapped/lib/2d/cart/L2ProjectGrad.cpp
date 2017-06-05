#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Legendre2d.h"
#include <iostream>
using namespace std;
// -------------------------------------------------------------------------- //
// All-purpose routine for computing the L2-projection
// of various functions onto Gradient (e.g divergence) of the Legendre basis.
//
// The routine L2ProjectGradAdd assumes that there is information stored in 
// fout that the user doesn't want overwritten.  In this case, start by 
// projecting Func onto the basis functions, and then "add" the moments onto
// fout.
//
// This is the case in "Part III: compute intra-element contributions" of
// ConstructL, where we already have inter-element constributions saved.
//
// The form of Func is inconsistent with what is currently in L2Project.  The
// reason being that this routine was only ever used to project FluxFunc onto
// the basis functions.  The idea is the following:
//
// A divergence of phi produces two components: div( phi ) = < phi_x, phi_y >.
//
// Therefore, to "project" this onto the basis functions, we need to have
// access to something with two components: f = < f1, f2 >.
//
// The fourth argument to Func sets fvals( 1:numpts, 1:MEQN, 1:NDIMS ), which
// supplies us with the necessary second argument for f.
//
// At the end of the day, we end up computing the following quantity:
//
//    fout( istart:iend, jstart:jend, 1:length(f), k ) = 
//
//              1 / dA \iint \div( phi^{(k)} ) \cdot f\, dx dy
//
// -------------------------------------------------------------------------- //

void mapc2p(double& xc,double& yc);
    double mapx(double xi,double eta,double xp1,double xp2,double xp3,double xp4,double yp1,double yp2,double yp3,double yp4);
    double mapy(double xi,double eta,double xp1,double xp2,double xp3,double xp4,double yp1,double yp2,double yp3,double yp4);
double dxidx(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4);
double dxidy(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4);
double detadx(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4);
double detady(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4);
double jacobian(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4);

void L2ProjectGradAdd(const int istart, 
		      const int iend, 
		      const int jstart, 
		      const int jend,
		      const int QuadOrder, 
		      const int BasisOrder_qin,
		      const int BasisOrder_auxin,
		      const int BasisOrder_fout,
		      const dTensorBC4* qin,
		      const dTensorBC4* auxin, 
		      dTensorBC4* fout,
		      void Func( const dTensor2&, const dTensor2&,
				 const dTensor2&, dTensor3&) )
{
  // dx and dy
  const double dx   = dogParamsCart2.get_dx();
  const double dy   = dogParamsCart2.get_dy();
  const double xlow = dogParamsCart2.get_xlow();
  const double ylow = dogParamsCart2.get_ylow();
  
  // mbc
  const int mbc = qin->getmbc();
  assert_eq(mbc, auxin->getmbc());
  assert_eq(mbc, fout->getmbc());
  
  // qin variable
  const int       mx = qin->getsize(1);
  const int       my = qin->getsize(2);
  const int     meqn = qin->getsize(3);
  const int kmax_qin = qin->getsize(4);
  assert_eq(kmax_qin, (BasisOrder_qin*BasisOrder_qin));
  
  // auxin variable
  assert_eq(mx, auxin->getsize(1));
  assert_eq(my, auxin->getsize(2));
  const int       maux = auxin->getsize(3);
  const int kmax_auxin = auxin->getsize(4);
  assert_eq(kmax_auxin, (BasisOrder_auxin*BasisOrder_auxin));
  
  // fout variables
  assert_eq(mx, fout->getsize(1));
  assert_eq(my, fout->getsize(2));
  const int mcomps_out = fout->getsize(3);
  const int  kmax_fout = fout->getsize(4);
  assert_eq(kmax_fout, (BasisOrder_fout*BasisOrder_fout));
  
  // starting and ending indeces
  assert_ge(istart, 1-mbc);
  assert_le(iend, mx+mbc);
  assert_ge(jstart, 1-mbc);
  assert_le(jend, my+mbc);
  
  // number of quadrature points
  assert_ge(QuadOrder, 1);
  assert_le(QuadOrder, 6);
  const int mpoints = iMax((QuadOrder+1)*(QuadOrder+1), 1);
  
  // trivial case
  if ( BasisOrder_fout==1 )
    {
      fout->setall(0.);
      return;
    }
  
  // set quadrature point and weight information
  void SetQuadWgtsPts(const int, dTensor1&, dTensor2&);
  dTensor1 wgt(mpoints);
  dTensor2 spts(mpoints, 2);
  SetQuadWgtsPts(QuadOrder+1, wgt, spts);
  
  // Loop over each quadrature point to construct Legendre polys
  const int kmax = iMax(iMax(kmax_qin, kmax_auxin), kmax_fout);
  dTensor2 phi(mpoints, kmax);
  dTensor2 phi_x(mpoints, kmax_fout);
  dTensor2 phi_y(mpoints, kmax_fout);
  
  void SetLegendrePolys(const int, const int, const dTensor2&, dTensor2&);
  SetLegendrePolys(mpoints, kmax, spts, phi);
  
  void SetLegendrePolysGrad(const double, 
			    const double, 
			    const int,
			    const int,
			    const dTensor2&,
			    dTensor2&,
			    dTensor2&);
  SetLegendrePolysGrad(dx, dy, mpoints, kmax_fout, spts, phi_x, phi_y);
  
  dTensor2 wgt_phi_x_tr(kmax_fout, mpoints);
  dTensor2 wgt_phi_y_tr(kmax_fout, mpoints);
  iTensorBase nonzero_flags(kmax+1);
  nonzero_flags.setall(0);
  for(int k=1; k<=kmax_fout; k++)
    for(int mp=1; mp<=mpoints; mp++)
      {
	double xval = wgt.get(mp)*phi_x.get(mp, k);
	double yval = wgt.get(mp)*phi_y.get(mp, k);
	wgt_phi_x_tr.set(k, mp, xval );
	wgt_phi_y_tr.set(k, mp, yval );
	if(xval) nonzero_flags.vfetch(k) |= 1;
	if(yval) nonzero_flags.vfetch(k) |= 2;
      }
  
  // ----------------------------------
  // Loop over all elements of interest
  // ----------------------------------

#pragma omp parallel
{
  const double s_area = 1.0;
  const double s_area_inv = 1.0;//1./s_area;

        dTensor2 xpts(mpoints, 2);
        dTensor1 Jacob(mpoints);
        dTensor2 ddx(mpoints, 2);
        dTensor2 ddy(mpoints, 2);
        dTensor2 qvals(mpoints, meqn);
        dTensor2 auxvals(mpoints, maux);
        dTensor3 fvalso(mpoints, mcomps_out, 2);
        dTensor3 fvals(mpoints, mcomps_out, 2); 

 
#pragma omp for
  for (int i=istart; i<=iend; i++)
    for (int j=jstart; j<=jend; j++)
      {
	// find center of current cell
	const double xc = xlow + (double(i)-0.5)*dx;
	const double yc = ylow + (double(j)-0.5)*dy;
	
        double xp1 = xlow + (double(i)-1.0)*dx;
        double yp1 = ylow + (double(j)-1.0)*dy;
        double xp2 = xlow + (double(i))*dx;
        double yp2 = ylow + (double(j-1))*dy; 
        double xp3 = xlow + (double(i))*dx;
        double yp3 = ylow + (double(j))*dy;
        double xp4 = xlow + (double(i)-1)*dx;
        double yp4 = ylow + (double(j))*dy;
        mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);

	
	// Loop over each quadrature point
	for (int m=1; m<=(mpoints); m++)
	  {
	    //assert_eq(phi.get(m, 1), 1.);
	    
	    // grid point (x, y)
	    const double xi  = spts.get(m, 1);
	    const double eta = spts.get(m, 2);          
            double x1=mapx(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4);
            double y1=mapy(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4);
            xpts.set( m, 1, x1  );
            xpts.set( m, 2, y1 );
	    ddx.set( m, 1, dxidx(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4));
	    ddx.set( m, 2, detadx(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4)); 
	    ddy.set( m, 1, dxidy(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4));
	    ddy.set( m, 2, detady(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4)); 
          Jacob.set(m, jacobian(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4));    

	    // Solution values (q) at each grid point
	    for (int me=1; me<=meqn; me++)
	      {
		double val = qin->get(i, j, me, 1);
		
		for (int k=2; k<=kmax_qin; k++)
		  {
		    val += phi.get(m, k) * qin->get(i, j, me, k);
		  }
		qvals.set(m, me, val );
	      }
	    
	    // Auxiliary values (aux) at each grid point
	    for (int ma=1; ma<=maux; ma++)
	      {
		double val = auxin->get(i, j, ma, 1);
		
		for (int k=2; k<=kmax_auxin; k++)
		  {
		    val += phi.get(m, k) * auxin->get(i, j, ma, k);
		  }
		auxvals.set(m, ma, val );
	      }
	  }
	// Call user-supplied function to set fvals
	Func(xpts, qvals, auxvals, fvalso);

	// Loop over each quadrature point
	for (int m=1; m<=(mpoints); m++)
	  {
	    //assert_eq(phi.get(m, 1), 1.);
	    
	    // grid point (x, y)
	    const double xi  = spts.get(m, 1);
	    const double eta = spts.get(m, 2);                  

	    // Reorder Flux at each point
	    for (int me=1; me<=mcomps_out; me++)
	      {
        //      fvals.set(m,me,1,fvalso.get(m,me,1));
        //     fvals.set(m,me,2,fvalso.get(m,me,2));
           
            fvals.set(m,me,1,ddx.get(m,1)*fvalso.get(m,me,1)+ddy.get(m,1)*fvalso.get(m,me,2));
            fvals.set(m,me,2,ddx.get(m,2)*fvalso.get(m,me,1)+ddy.get(m,2)*fvalso.get(m,me,2));
	      }
	    
	  }
	

	
	for (int mi=1; mi<=mcomps_out; mi++)
	  {
	    // This is the second-most executed block in the code after
	    // the loop in L2Project.
	    // can start with 2 since phi_*(:, 1)=0
	    for (int k2=2; k2<=kmax_fout; k2++)
	      {
		double tmp = 0.;
		// switch statement to avoid unnecessary computations
		switch( 3 )//nonzero_flags.vget(k2) )
		  {
		  case 0:
		    eprintf("this case should only arise when k2=%d equals 1", k2);
		    
		  case 1:
		    for (int mp=1; mp<=mpoints; mp++)
		      {
			tmp += fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp)*Jacob.get(mp);
		      }
		    break;
		    
		  case 2:
		    for (int mp=1; mp<=mpoints; mp++)
		      {
			tmp += fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp)*Jacob.get(mp);
		      }
		    break;
		    
		  case 3:
		  default:
		    for (int mp=1; mp<=mpoints; mp++)
		      {
			tmp += ( fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp) +
				 fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp) )*Jacob.get(mp);
		      }
		    break;
		  }
		fout->fetch(i, j, mi, k2) += tmp*s_area_inv;
	      }
	  }
      }
 }
  
}// end of function L2ProjectGrad

void L2ProjectGradAdd(int istart, int iend, int jstart, int jend,
		      const dTensorBC4& qin,
		      const dTensorBC4& auxin, dTensorBC4& Fout,
		      void Func(const dTensor2&, const dTensor2&,
				const dTensor2&, dTensor3&))
{
  const int space_order = dogParams.get_space_order();
  L2ProjectGradAdd(istart, iend, jstart, jend,
		   space_order, space_order, space_order, space_order,
		   &qin, &auxin, &Fout, Func);
}

// Loop over each quadrature point to construct Legendre polys
void SetLegendrePolysGrad(const double dx,
			  const double dy,
			  const int mpoints, 
			  const int kmax, 
			  const dTensor2& spts, 
			  dTensor2& phi_x,
			  dTensor2& phi_y)
{
  
  {
    assert_eq(mpoints, spts.getsize(1));
    assert_eq(mpoints, phi_x.getsize(1));
    assert_eq(kmax, phi_x.getsize(2));
    assert_eq(mpoints, phi_y.getsize(1));
    assert_eq(kmax, phi_y.getsize(2));
  }
  
  const double tmpx = 1.0;//2.0/dx;
  const double tmpy = 1.0;//2.0/dy;

  for (int m=1; m<=(mpoints); m++)
    {
      // grid point (x, y)
      const double xi  = spts.get(m, 1);
      const double eta = spts.get(m, 2);
      const double xi2 = xi*xi;
      const double xi3 = xi*xi2;
      const double xi4 = xi*xi3;
      const double eta2 = eta*eta;
      const double eta3 = eta*eta2;
      const double eta4 = eta*eta3;     
      
      // Gradient of Legendre basis functions at each gaussian quadrature point
      switch( kmax )
        {
	case 25:  // fifth order
          phi_x.set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi) );
          phi_x.set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*tmpx*sq7*(7.5*xi2 - 1.5) );
          phi_x.set( m,23, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi)*sq7*(2.5*eta3 - 1.5*eta) );
          phi_x.set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*tmpx*sq5*3.0*xi );
          phi_x.set( m,21, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi)*sq5*(1.5*eta2 - 0.5) );
          phi_x.set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*tmpx );
          phi_x.set( m,19, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi)*sq3*eta );
          phi_x.set( m,18, 0.0 );
          phi_x.set( m,17, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi) );


          phi_y.set( m,25, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
          phi_y.set( m,24, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*sq7*(2.5*xi3 - 1.5*xi) );
          phi_y.set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*tmpy*sq7*(7.5*eta2 - 1.5) );
          phi_y.set( m,22, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*sq5*(1.5*xi2 - 0.5) );
          phi_y.set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*tmpy*sq5*3.0*eta );
          phi_y.set( m,20, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta)*sq3*xi );
          phi_y.set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*tmpy );
          phi_y.set( m,18, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta) );
          phi_y.set( m,17, 0.0 );

	case 16:  // fourth order
          phi_x.set( m,16, sq7*(2.5*eta3 - 1.5*eta)*tmpx*sq7*(7.5*xi2 - 1.5)  );
          phi_x.set( m,15, sq7*(2.5*eta3 - 1.5*eta)*tmpx*sq5*3.0*xi );
          phi_x.set( m,14, tmpx*sq7*(7.5*xi2 - 1.5) *sq5*(1.5*eta2 - 0.5) );
          phi_x.set( m,13, sq7*(2.5*eta3 - 1.5*eta)*sq3*tmpx );
          phi_x.set( m,12, tmpx*sq7*(7.5*xi2 - 1.5) *sq3*eta );
          phi_x.set( m,11, 0.0 );
          phi_x.set( m,10,  tmpx*sq7*(7.5*xi2 - 1.5)  );


          phi_y.set( m,16, tmpy*sq7*(7.5*eta2 - 1.5)*sq7*(2.5*xi3 - 1.5*xi) );
          phi_y.set( m,15, tmpy*sq7*(7.5*eta2 - 1.5)*sq5*(1.5*xi2 - 0.5) );
          phi_y.set( m,14, sq7*(2.5*xi3 - 1.5*xi)*tmpy*sq5*3.0*eta );
          phi_y.set( m,13, tmpy*sq7*(7.5*eta2 - 1.5)*sq3*xi );
          phi_y.set( m,12, sq7*(2.5*xi3 - 1.5*xi)*sq3*tmpy );
          phi_y.set( m,11, tmpy*sq7*(7.5*eta2 - 1.5) );
          phi_y.set( m,10,  0.0 );

	case 9:  // third order
          phi_x.set( m,9,  sq5*tmpx*3.0*xi*sq5*(1.5*eta2 - 0.5) );
          phi_x.set( m,8,  sq3*sq5*tmpx*(1.5*eta2 - 0.5) );
          phi_x.set( m,7,  sq3*eta*tmpx*sq5*3.0*xi );
          phi_x.set( m,6,  0.0 );
          phi_x.set( m,5,  tmpx*sq5*3.0*xi );


          phi_y.set( m,9,  sq5*(1.5*xi2 - 0.5)*tmpy*sq5*3.0*eta );
          phi_y.set( m,8,  sq3*xi*tmpy*sq5*3.0*eta );
          phi_y.set( m,7,  sq3*sq5*tmpy*(1.5*xi2 - 0.5) );
          phi_y.set( m,6,  tmpy*sq5*3.0*eta );
          phi_y.set( m,5,  0.0 );

	case 4:  // second order
          phi_x.set( m,4,  3.0*tmpx*eta );
          phi_x.set( m,3, 0.0 );
          phi_x.set( m,2, sq3*tmpx  );

          phi_y.set( m,4,  3.0*xi*tmpy );
          phi_y.set( m,3, sq3*tmpy );
          phi_y.set( m,2, 0.0  );

	case 1:  // first order
	  phi_x.set( m, 1,  0.0 );
	  
	  phi_y.set( m, 1,  0.0 );
	  
	  break;
	  
	default:
	  unsupported_value_error(kmax);
        }
    }
}

// Similar function to L2ProjectGradAdd, but this time in place of using a user
// supplied function, we pass in Legendre weights for a function we already
// have the Legendre expansion of.
//
// Parameters:
// ----------
//
// The function we're projecting onto the basis functions is <F,G>:
//
//    F : F( 1-mbc, mx+mbc, 1-mbc, my+mbc, 1:meqn, 1:kmax )
//    G : G( 1-mbc, mx+mbc, 1-mbc, my+mbc, 1:meqn, 1:kmax )
//
// Returns:
// -------
//
// fout : the projected function of 
//
//        1/dA \iint \div( phi^{(k)} ) \cdot < F, G >.
//
//            
void L2ProjectGradAddLegendre(const int istart, 
			      const int iend, 
			      const int jstart, 
			      const int jend,
			      const int QuadOrder, 
			      const dTensorBC4* F, 
			      const dTensorBC4* G, 
			      dTensorBC4* fout )
{ cout<<"BLAH!"<<endl;
  
  // dx and dy
  const double dx   = dogParamsCart2.get_dx();
  const double dy   = dogParamsCart2.get_dy();
  const double xlow = dogParamsCart2.get_xlow();
  const double ylow = dogParamsCart2.get_ylow();
  
  // mbc
  const int mbc = F->getmbc();
  assert_eq(mbc,  G->getmbc());
  assert_eq(mbc, fout->getmbc());
  
  // fin variable
  const int       mx = F->getsize(1);
  const int       my = F->getsize(2);
  const int     meqn = F->getsize(3);
  const int     kmax = F->getsize(4);
  
  // fout variables
  assert_eq(mx, fout->getsize(1));
  assert_eq(my, fout->getsize(2));
  const int mcomps_out = fout->getsize(3);   assert_eq( meqn, mcomps_out );
  const int  kmax_fout = fout->getsize(4);   assert_eq( kmax, kmax_fout  );
  
  // starting and ending indeces
  assert_ge(istart, 1-mbc);
  assert_le(iend, mx+mbc);
  assert_ge(jstart, 1-mbc);
  assert_le(jend, my+mbc);
  
  // number of quadrature points
  assert_ge(QuadOrder, 1);
  assert_le(QuadOrder, 6);
  const int mpoints = iMax((QuadOrder)*(QuadOrder), 1);
    
  // trivial case
  if ( QuadOrder ==1 )
    {
      fout->setall(0.);
      return;
    }

  // set quadrature point and weight information
  void SetQuadWgtsPts(const int, dTensor1&, dTensor2&);
  dTensor1 wgt(mpoints);
  dTensor2 spts(mpoints, 2);
  SetQuadWgtsPts(QuadOrder-1, wgt, spts);
  
  // Loop over each quadrature point to construct Legendre polys
  dTensor2 phi(mpoints,   kmax     );
  dTensor2 phi_x(mpoints, kmax_fout);
  dTensor2 phi_y(mpoints, kmax_fout);
  
  void SetLegendrePolys(const int, const int, const dTensor2&, dTensor2&);
  SetLegendrePolys(mpoints, kmax, spts, phi);
  
  void SetLegendrePolysGrad(const double, 
			    const double, 
			    const int,
			    const int,
			    const dTensor2&,
			    dTensor2&,
			    dTensor2&);
  SetLegendrePolysGrad(dx, dy, mpoints, kmax_fout, spts, phi_x, phi_y);
  
  dTensor2 wgt_phi_x_tr(kmax_fout, mpoints);
  dTensor2 wgt_phi_y_tr(kmax_fout, mpoints);
  iTensorBase nonzero_flags(kmax+1);
  nonzero_flags.setall(0);
  for(int k=1; k<=kmax_fout; k++)
    for(int mp=1; mp<=mpoints; mp++)
      {
	double xval = wgt.get(mp)*phi_x.get(mp, k);
	double yval = wgt.get(mp)*phi_y.get(mp, k);
	wgt_phi_x_tr.set(k, mp, xval );
	wgt_phi_y_tr.set(k, mp, yval );
	if(xval) nonzero_flags.vfetch(k) |= 1;
	if(yval) nonzero_flags.vfetch(k) |= 2;
      }
  
  // ----------------------------------
  // Loop over all elements of interest
  // ----------------------------------
  const double s_area = 1.0;
  const double s_area_inv = 1./s_area;
  
#pragma omp parallel for
  for (int i=istart; i<=iend; i++)
    for (int j=jstart; j<=jend; j++)
      {
	// find center of current cell
	const double xc = xlow + (double(i)-0.5)*dx;
	const double yc = ylow + (double(j)-0.5)*dy;

	dTensor2  xpts(mpoints, 2);
	dTensor3 fvals(mpoints, mcomps_out, 2);



            double xp1 = xc-0.5*dx;
            double yp1= yc-0.5*dy;
            double xp2 = xc+0.5*dx;
            double yp2= yc-0.5*dy;
            double xp3 = xc+0.5*dx;
            double yp3= yc+0.5*dy;
            double xp4 = xc-0.5*dx;
            double yp4= yc+0.5*dy;


            mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);




	
	// Loop over each quadrature point
	for (int m=1; m<=(mpoints); m++)
	  {

	    // grid point (x, y)
	    const double xi  = spts.get(m, 1);
	    const double eta = spts.get(m, 2);          
            double x1=mapx(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4);
            double y1=mapy(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4);
	    xpts.set( m, 1, x1  );
	    xpts.set( m, 2, y1 );         

	    // Solution values (q) at each grid point
	    for (int me=1; me<=meqn; me++)
	      {
		double val1 = F->get(i, j, me, 1);
		double val2 = G->get(i, j, me, 1);
		for (int k=2; k<=kmax; k++)
		  {
		    val1 += phi.get(m, k) * F->get(i, j, me, k);
		    val2 += phi.get(m, k) * G->get(i, j, me, k);
		  }
		fvals.set(m, me, 1, val1 );
		fvals.set(m, me, 2, val2 );
	      }
	    
	  }
	
	// Call user-supplied function to set fvals
	// 
	// No need to call this, because we use the basis functions to
	// evaluate F and G
	//
	//Func(xpts, qvals, auxvals, fvals);
	
	for (int mi=1; mi<=mcomps_out; mi++)
	  {
	    
	    // can start with 2 since phi_*(:, 1)=0
	    for (int k2=2; k2<=kmax_fout; k2++)
	      {
		
		double tmp = 0.;
		// switch statement to avoid unnecessary computations
		switch( nonzero_flags.vget(k2) )
		  {
		  case 0:
		    eprintf("this case should only arise when k2=%d equals 1", k2);

		  case 1:
		    for (int mp=1; mp<=mpoints; mp++)
		      {
			tmp += fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp);
		      }
		    break;

		  case 2:
		    for (int mp=1; mp<=mpoints; mp++)
		      {
			tmp += fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp);
		      }
		    break;
		    
		  case 3:
		  default:
		    for (int mp=1; mp<=mpoints; mp++)
		      {
			tmp += ( fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp) +
				 fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp) );
		      }
		    break;
		  }
		fout->fetch(i, j, mi, k2) += tmp*s_area_inv;
	      }
	  }
      }
  
}
