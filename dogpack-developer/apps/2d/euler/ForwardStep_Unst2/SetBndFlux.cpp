#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "RiemannSolve.h"
#include "stdio.h"
#include <iostream>
using namespace std;

// to use this RiemannSolver the user must implement the following
// "linking callbacks":
// * FluxFunc
// * SetWaveSpd

void qinitfuncs(double x,double y,double &rho,double &press,double &u1,double &u2,double &u3,double &energy);

void SetBndFlux(const dTensor1& nvec, const dTensor1& xedge,
		    dTensor1& Ql, dTensor1& Qr,
		    const dTensor1& Auxl, const dTensor1& Auxr,
		    int side)
{
    int m,mcase;
    double smax_edge = 0.0e0;
    int meqn = Ql.getsize();
    int maux = Auxl.getsize();
    dTensor1 ffl(meqn), ffr(meqn);
    dTensor3 ftmp(1,meqn,2);
    double s1,s2,Qstar;
    dTensor2 xedge_tmp(1,2);
    dTensor2 Auxl_tmp(1,maux),Auxr_tmp(1,maux);

    double x1=xedge.get(1);
    double y1=xedge.get(2);

    double gamma=1.4;

    double rho   =  1.0;
    double u1    =  3.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0/1.4;
    double energy = press/(gamma-1.0e0)
                + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);


  
    //qinitfuncs(0.0, 7.0, rho, press, u1, u2, u3, energy);

    if(abs(x1-0.0)<1.0e-12)
    {
           if(side==1)
           {
             Ql.set(1, rho);
             Ql.set(2, rho*u1 );
             Ql.set(3, rho*u2 );
             Ql.set(4, rho*u3 );
             Ql.set(5, energy );
           }
           else
           {
             Qr.set(1, rho );
             Qr.set(2, rho*u1 );
             Qr.set(3, rho*u2 );
             Qr.set(4, rho*u3 );
             Qr.set(5, energy );
           }
    }
    else if(x1>3.0-1.0e-12)
    {//Zeroth order extrapolation on right top and bottom edges
           if(side==1)
           {
             Ql.set(1, Qr.get(1));
             Ql.set(2, Qr.get(2) );
             Ql.set(3, Qr.get(3) );
             Ql.set(4, Qr.get(4) );
             Ql.set(5, Qr.get(5) );
           }
           else if (side==2)
           {
             Qr.set(1, Ql.get(1) );
             Qr.set(2, Ql.get(2) );
             Qr.set(3, Ql.get(3) );
             Qr.set(4, Ql.get(4) );
             Qr.set(5, Ql.get(5) );
           }
     }
     else{
           double nx=nvec.get(1);
           double ny=nvec.get(2);
           double nmag=sqrt(nx*nx+ny*ny);
           dTensor1 nhat(2);
           nhat.set(1, nx/nmag );
           nhat.set(2, ny/nmag );
           dTensor1 that(2);
           that.set(1, ny/nmag );
           that.set(2, -nx/nmag );
           double ux,uy;
           if(side==1)
           {
              ux=Qr.get(2)/Qr.get(1);
              uy=Qr.get(3)/Qr.get(1);
           }
           if(side==2)
           {
              ux=Ql.get(2)/Ql.get(1);
              uy=Ql.get(3)/Ql.get(1);          
           }
           double un=ux*nhat.get(1)+uy*nhat.get(2);
           double ut=ux*that.get(1)+uy*that.get(2);
           double uxnew=-un*nhat.get(1)+ut*that.get(1);
           double uynew=-un*nhat.get(2)+ut*that.get(2);
           if(side==1)
           {
             //printf(" here %f %f \n",nhat.get(1),that.get(1));
             double rho=Qr.get(1);
             double energy=Qr.get(5);
             double press1=0.4*(energy-0.5*rho*(ux*ux+uy*uy));
             energy=press1 / (0.4) + 0.5*rho*(uxnew*uxnew + uynew*uynew);
             Ql.set(1, rho);  
             Ql.set(2, rho*uxnew );  
             Ql.set(3, rho*uynew );    
             Ql.set(4, Qr.get(4) );
             Ql.set(5, energy );  
           }
           else
           {
             double rho=Ql.get(1);
             double energy=Ql.get(5);
             double press1=0.4*(energy-0.5*rho*(ux*ux+uy*uy));
             energy=press1/(0.4)+ 0.5e0*rho*(uxnew*uxnew + uynew*uynew);
             Qr.set(1, rho );  
             Qr.set(2, rho*uxnew );  
             Qr.set(3, rho*uynew );    
             Qr.set(4, Ql.get(4) );
             Qr.set(5, energy );  
           }
      }
}

