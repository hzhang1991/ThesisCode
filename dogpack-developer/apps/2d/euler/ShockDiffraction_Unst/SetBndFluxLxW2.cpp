#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "RiemannSolve.h"
#include <iostream>
#include "stdio.h"
using namespace std;

// to use this RiemannSolver the user must implement the following
// "linking callbacks":
// * FluxFunc
// * SetWaveSpd

void qinitfuncs1(double x,double y,double &rho,double &press,double &u1,double &u2,double &u3,double &energy);

void SetBndFluxLxW(const dTensor1& nvec, const dTensor1& xedge,
		    dTensor1& Ql, dTensor1& Qr,
		    dTensor1& Qlt, dTensor1& Qrt,
		    dTensor1& Qltt, dTensor1& Qrtt,
		    dTensor1& Qlttt, dTensor1& Qrttt,
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


  

    if(abs(x1-0.0)<1.0e-14 && y1>6.0)
    {

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;
    double energy = press/(gamma-1.0e0)
                + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

           qinitfuncs1(0.0, 7.0, rho, press, u1, u2, u3, energy);
           if(side==1)
           {
             Ql.set(1, rho);
             Ql.set(2, rho*u1 );
             Ql.set(3, rho*u2 );
             Ql.set(4, rho*u3 );
             Ql.set(5, energy );

             Qlt.set(1, 0.0);
             Qlt.set(2, 0.0 );
             Qlt.set(3, 0.0 );
             Qlt.set(4, 0.0 );
             Qlt.set(5, 0.0 );

             Qltt.set(1, 0.0);
             Qltt.set(2, 0.0 );
             Qltt.set(3, 0.0 );
             Qltt.set(4, 0.0 );
             Qltt.set(5, 0.0 );

             Qlttt.set(1, 0.0);
             Qlttt.set(2, 0.0 );
             Qlttt.set(3, 0.0 );
             Qlttt.set(4, 0.0 );
             Qlttt.set(5, 0.0 );


           }
           else
           {
             Qr.set(1, rho );
             Qr.set(2, rho*u1 );
             Qr.set(3, rho*u2 );
             Qr.set(4, rho*u3 );
             Qr.set(5, energy );

             Qrt.set(1, 0.0);
             Qrt.set(2, 0.0 );
             Qrt.set(3, 0.0 );
             Qrt.set(4, 0.0 );
             Qrt.set(5, 0.0 );

             Qrtt.set(1, 0.0);
             Qrtt.set(2, 0.0 );
             Qrtt.set(3, 0.0 );
             Qrtt.set(4, 0.0 );
             Qrtt.set(5, 0.0 );

             Qrttt.set(1, 0.0);
             Qrttt.set(2, 0.0 );
             Qrttt.set(3, 0.0 );
             Qrttt.set(4, 0.0 );
             Qrttt.set(5, 0.0 );
           }


    }
    else if( x1<0.0 || y1>6.0 || x1>3.4 || abs(y1)<1.0e-14 )
    {//Zeroth order extrapolation on right top and bottom edges
           if(side==1)
           {
             Ql.set(1, Qr.get(1));
             Ql.set(2, Qr.get(2) );
             Ql.set(3, Qr.get(3) );
             Ql.set(4, Qr.get(4) );
             Ql.set(5, Qr.get(5) );

             Qlt.set(1, Qrt.get(1));
             Qlt.set(2, Qrt.get(2) );
             Qlt.set(3, Qrt.get(3) );
             Qlt.set(4, Qrt.get(4) );
             Qlt.set(5, Qrt.get(5) );


             Qltt.set(1, Qrtt.get(1));
             Qltt.set(2, Qrtt.get(2) );
             Qltt.set(3, Qrtt.get(3) );
             Qltt.set(4, Qrtt.get(4) );
             Qltt.set(5, Qrtt.get(5) );

             Qlttt.set(1, Qrttt.get(1));
             Qlttt.set(2, Qrttt.get(2) );
             Qlttt.set(3, Qrttt.get(3) );
             Qlttt.set(4, Qrttt.get(4) );
             Qlttt.set(5, Qrttt.get(5) );
           }
           else if (side==2)
           {
             Qr.set(1, Ql.get(1) );
             Qr.set(2, Ql.get(2) );
             Qr.set(3, Ql.get(3) );
             Qr.set(4, Ql.get(4) );
             Qr.set(5, Ql.get(5) );

             Qrt.set(1, Qlt.get(1) );
             Qrt.set(2, Qlt.get(2) );
             Qrt.set(3, Qlt.get(3) );
             Qrt.set(4, Qlt.get(4) );
             Qrt.set(5, Qlt.get(5) );

             Qrtt.set(1, Qltt.get(1) );
             Qrtt.set(2, Qltt.get(2) );
             Qrtt.set(3, Qltt.get(3) );
             Qrtt.set(4, Qltt.get(4) );
             Qrtt.set(5, Qltt.get(5) );


             Qrttt.set(1, Qlttt.get(1) );
             Qrttt.set(2, Qlttt.get(2) );
             Qrttt.set(3, Qlttt.get(3) );
             Qrttt.set(4, Qlttt.get(4) );
             Qrttt.set(5, Qlttt.get(5) );
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
           double uxt,uyt;
           double uxtt,uytt;
           double uxttt,uyttt;

           if(side==1)
           {
              ux=Qr.get(2);
              uy=Qr.get(3);

              uxt=Qrt.get(2);
              uyt=Qrt.get(3);

              uxtt=Qrtt.get(2);
              uytt=Qrtt.get(3);

              uxttt=Qrttt.get(2);
              uyttt=Qrttt.get(3);
           }
           if(side==2)
           {
              ux=Ql.get(2);
              uy=Ql.get(3);          

              uxt=Qlt.get(2);
              uyt=Qlt.get(3);          

              uxtt=Qltt.get(2);
              uytt=Qltt.get(3);          

              uxttt=Qlttt.get(2);
              uyttt=Qlttt.get(3);          
           }
           double un=ux*nhat.get(1)+uy*nhat.get(2);
           double ut=ux*that.get(1)+uy*that.get(2);

           double un_t=uxt*nhat.get(1)+uyt*nhat.get(2);
           double ut_t=uxt*that.get(1)+uyt*that.get(2);

           double un_tt=uxtt*nhat.get(1)+uytt*nhat.get(2);
           double ut_tt=uxtt*that.get(1)+uytt*that.get(2);

           double un_ttt=uxttt*nhat.get(1)+uyttt*nhat.get(2);
           double ut_ttt=uxttt*that.get(1)+uyttt*that.get(2);

           double uxnew=-un*nhat.get(1)+ut*that.get(1);
           double uynew=-un*nhat.get(2)+ut*that.get(2);

           double uxnewt=-un_t*nhat.get(1)+ut_t*that.get(1);
           double uynewt=-un_t*nhat.get(2)+ut_t*that.get(2);

           double uxnewtt=-un_tt*nhat.get(1)+ut_tt*that.get(1);
           double uynewtt=-un_tt*nhat.get(2)+ut_tt*that.get(2);

           double uxnewttt=-un_ttt*nhat.get(1)+ut_ttt*that.get(1);
           double uynewttt=-un_ttt*nhat.get(2)+ut_ttt*that.get(2);

           //double uxnew=ut*that.get(1);
           //double uynew=ut*that.get(2);

           if(side==1)
           {
             double rho=Qr.get(1);
             double energy=Qr.get(5);

             /*
             if(abs(Qr.get(2)-uxnew)>1.0e-12 || abs(Qr.get(3)-uynew)>1.0e-12)
             {printf("HERE %e %e %e %e \n",Qr.get(2),uxnew,Qr.get(3),uynew);
             printf("HERE %e %e %e %e \n",Qr.get(2),uxnew,Qr.get(3),uynew);
             printf("HERE %e %e %e %e \n",Qr.get(2),uxnew,Qr.get(3),uynew);
             printf("HERE %e %e %e %e \n",Qr.get(2),uxnew,Qr.get(3),uynew);

             printf("HEREt %e %e %e %e \n",Qrt.get(2),uxnewt,Qrt.get(3),uynewt);
             printf("HEREt %e %e %e %e \n",Qrt.get(2),uxnewt,Qrt.get(3),uynewt);
             printf("HEREt %e %e %e %e \n",Qrt.get(2),uxnewt,Qrt.get(3),uynewt);
             printf("HEREt %e %e %e %e \n",Qrt.get(2),uxnewt,Qrt.get(3),uynewt);
             }*/
             /*
             Ql.set(1, Qr.get(1));
             Ql.set(2, Qr.get(2) );
             Ql.set(3, Qr.get(3) );
             Ql.set(4, Qr.get(4) );
             Ql.set(5, Qr.get(5) );

             Qlt.set(1, Qrt.get(1));
             Qlt.set(2, Qrt.get(2) );
             Qlt.set(3, Qrt.get(3) );
             Qlt.set(4, Qrt.get(4) );
             Qlt.set(5, Qrt.get(5) );


             Qltt.set(1, Qrtt.get(1));
             Qltt.set(2, Qrtt.get(2) );
             Qltt.set(3, Qrtt.get(3) );
             Qltt.set(4, Qrtt.get(4) );
             Qltt.set(5, Qrtt.get(5) );

             Qlttt.set(1, Qrttt.get(1));
             Qlttt.set(2, Qrttt.get(2) );
             Qlttt.set(3, Qrttt.get(3) );
             Qlttt.set(4, Qrttt.get(4) );
             Qlttt.set(5, Qrttt.get(5) );
             */
             
             Ql.set(1, rho);  
             Ql.set(2, uxnew );  
             Ql.set(3, uynew );    
             Ql.set(4, Qr.get(4));  
             Ql.set(5, energy );  

             rho=Qrt.get(1);
             energy=Qrt.get(5);
             Qlt.set(1, rho);
             Qlt.set(2, uxnewt );    
             Qlt.set(3, uynewt );    
             Qlt.set(4, Qrt.get(4));  
             Qlt.set(5, energy );

             rho=Qrtt.get(1);
             energy=Qrtt.get(5);
             Qltt.set(1, rho);
             Qltt.set(2, uxnewtt );    
             Qltt.set(3, uynewtt );    
             Qltt.set(4, Qrtt.get(4));  
             Qltt.set(5, energy );

             rho=Qrttt.get(1);
             energy=Qrttt.get(5);
             Qlttt.set(1, rho);
             Qlttt.set(2, uxnewttt );    
             Qlttt.set(3, uynewttt );    
             Qlttt.set(4, Qrttt.get(4));  
             Qlttt.set(5, energy );
             

           }
           else
           {
             /*
             Qr.set(1, Ql.get(1) );
             Qr.set(2, Ql.get(2) );
             Qr.set(3, Ql.get(3) );
             Qr.set(4, Ql.get(4) );
             Qr.set(5, Ql.get(5) );

             Qrt.set(1, Qlt.get(1) );
             Qrt.set(2, Qlt.get(2) );
             Qrt.set(3, Qlt.get(3) );
             Qrt.set(4, Qlt.get(4) );
             Qrt.set(5, Qlt.get(5) );

             Qrtt.set(1, Qltt.get(1) );
             Qrtt.set(2, Qltt.get(2) );
             Qrtt.set(3, Qltt.get(3) );
             Qrtt.set(4, Qltt.get(4) );
             Qrtt.set(5, Qltt.get(5) );


             Qrttt.set(1, Qlttt.get(1) );
             Qrttt.set(2, Qlttt.get(2) );
             Qrttt.set(3, Qlttt.get(3) );
             Qrttt.set(4, Qlttt.get(4) );
             Qrttt.set(5, Qlttt.get(5) );
             */
             
             double rho=Ql.get(1);
             double energy=Ql.get(5);
             Qr.set(1, rho );  
             Qr.set(2, uxnew );  
             Qr.set(3, uynew );    
             Qr.set(4, Ql.get(4) );    
             Qr.set(5, energy );  

             rho=Qlt.get(1);
             energy=Qlt.get(5);
             Qrt.set(1, rho);
             Qrt.set(2, uxnewt );
             Qrt.set(3, uynewt );
             Qrt.set(4, Qlt.get(4) );    
             Qrt.set(5, energy );

             rho=Qltt.get(1);
             energy=Qltt.get(5);
             Qrtt.set(1, rho);
             Qrtt.set(2, uxnewtt );
             Qrtt.set(3, uynewtt );
             Qrtt.set(4, Qltt.get(4) );    
             Qrtt.set(5, energy );

             rho=Qlttt.get(1);
             energy=Qlttt.get(5);
             Qrttt.set(1, rho);
             Qrttt.set(2, uxnewttt );
             Qrttt.set(3, uynewttt );
             Qrttt.set(4, Qlttt.get(4) );    
             Qrttt.set(5, energy );
            
           }


      }


//printf("HEREl %e %e %e %e %e \n",Ql.get(1),Ql.get(2),Ql.get(3),Ql.get(4),Ql.get(5));
//printf("HEREr %e %e %e %e %e \n",Qr.get(1),Qr.get(2),Qr.get(3),Qr.get(4),Qr.get(5));

//printf("HERElt %e %e %e %e %e \n",Qlt.get(1),Qlt.get(2),Qlt.get(3),Qlt.get(4),Qlt.get(5));
//printf("HERErt %e %e %e %e %e \n",Qrt.get(1),Qrt.get(2),Qrt.get(3),Qrt.get(4),Qrt.get(5));

//printf("HEREltt %e %e %e %e %e \n",Qltt.get(1),Qltt.get(2),Qltt.get(3),Qltt.get(4),Qltt.get(5));
//printf("HERErtt %e %e %e %e %e \n",Qrtt.get(1),Qrtt.get(2),Qrtt.get(3),Qrtt.get(4),Qrtt.get(5));

//printf("HERElttt %e %e %e %e %e \n",Qlttt.get(1),Qlttt.get(2),Qlttt.get(3),Qlttt.get(4),Qlttt.get(5));
//printf("HERErttt %e %e %e %e %e \n",Qrttt.get(1),Qrttt.get(2),Qrttt.get(3),Qrttt.get(4),Qrttt.get(5));

}

