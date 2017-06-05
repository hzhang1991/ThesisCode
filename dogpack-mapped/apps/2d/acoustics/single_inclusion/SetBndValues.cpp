#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"
#include <iostream>
using namespace std;    

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void mapc2p(double& xc,double& yc);
double averages(double px1,double py1,double px2,double py2,double px3,double py3,double px4,double py4,double dx1,double dy1,double a[6],int degree);

void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{
    int i,j,ell,m;
    double tmp;
    int mx   = q.getsize(1);
    int my   = q.getsize(2);
    int meqn = q.getsize(3);
    int kmax = q.getsize(4);
    int mbc  = q.getmbc();
    int maux = aux.getsize(3);
    const double t = dogParams.get_time();
    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();
    const double xlow = dogParamsCart2.get_xlow();
    const double ylow = dogParamsCart2.get_ylow();
    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------
 
    double s=0.0;
    if(t<0.25){s=1.0;}
      
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (i=0; i>=(1-mbc); i--)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
               // for (m=1; m<=meqn; m++)
               // {
                //    tmp = q.get(1,j,m,ell);
                //    q.set(i,j,m,ell, tmp );
               // }

            int i1=1-i;
            int j1=j;
      //      cout<<i1<<" "<<i<<endl;
            double xp1 = xlow+(i1-1)*dx;
            double yp1= ylow+(j1-1)*dy;
            double xp2 = xlow+(i1)*dx;
            double yp2= ylow+(j1-1)*dy;
            double xp3 = xlow+(i1)*dx;
            double yp3= ylow+(j1)*dy;
            double xp4 = xlow+(i1-1)*dx;
            double yp4= ylow+(j1)*dy;

            mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);
            double xmax=max(max(xp1,xp2),max(xp3,xp4));double xmin=min(min(xp1,xp2),min(xp3,xp4));
            double ymax=max(max(yp1,yp2),max(yp3,yp4));double ymin=min(min(yp1,yp2),min(yp3,yp4));
            double dx1=xmax-xmin;double dy1=ymax-ymin;
            double xpts[4];double ypts[4];
            xpts[0]=xp1;ypts[0]=yp1;
            xpts[1]=xp2;ypts[1]=yp2;
            xpts[2]=xp3;ypts[2]=yp3;
            xpts[3]=xp4;ypts[3]=yp4;
            double intermediate=0.0;
            double px1,py1,px2,py2,px3,py3,px4,py4=0.0;

             double a[6];double C1_now;

            px1=(xpts[0]-xmin)/dx1*2.0-1.0;py1=(ypts[0]-ymin)/dy1*2.0-1.0;
            px2=(xpts[1]-xmin)/dx1*2.0-1.0;py2=(ypts[1]-ymin)/dy1*2.0-1.0;
            px3=(xpts[2]-xmin)/dx1*2.0-1.0;py3=(ypts[2]-ymin)/dy1*2.0-1.0;
            px4=(xpts[3]-xmin)/dx1*2.0-1.0;py4=(ypts[3]-ymin)/dy1*2.0-1.0;

            a[0]=q.get(1-i,j,1,1);
            a[1]=q.get(1-i,j,1,2);
            a[2]=q.get(1-i,j,1,3);
            a[3]=q.get(1-i,j,1,4);
            a[4]=q.get(1-i,j,1,5);
            a[5]=q.get(1-i,j,1,6);
            C1_now=a[0];//averages(px1,py1,px2,py2,px3,py3,px4,py4,dx1,dy1,a,1);
            q.set(i,j,1,1,C1_now);
                for (ell=2; ell<=kmax; ell++)
    {
            q.set(i,j,1,ell, 0.0 );}
            a[0]=q.get(1-i,j,2,1);
            a[1]=q.get(1-i,j,2,2);
            a[2]=q.get(1-i,j,2,3);
            a[3]=q.get(1-i,j,2,4);
            a[4]=q.get(1-i,j,2,5);
            a[5]=q.get(1-i,j,2,6);
            C1_now=a[0];//averages(px1,py1,px2,py2,px3,py3,px4,py4,dx1,dy1,a,1);
            q.set(i,j,2,1,2.0*s-C1_now);       
            for (ell=2; ell<=kmax; ell++)
            {
               q.set(i,j,2,ell, 0.0 );}
            a[0]=q.get(1-i,j,3,1);
            a[1]=q.get(1-i,j,3,2);
            a[2]=q.get(1-i,j,3,3);
            a[3]=q.get(1-i,j,3,4);
            a[4]=q.get(1-i,j,3,5);
            a[5]=q.get(1-i,j,3,6);
            C1_now=a[0];//averages(px1,py1,px2,py2,px3,py3,px4,py4,dx1,dy1,a,1);
                  
                   q.set(i,j,3,1,C1_now);
                   for (ell=2; ell<=kmax; ell++)
                   { q.set(i,j,3,ell, 0.0 );}
            }
        } 
    for (ell=1; ell<=kmax; ell++)
    {

        // ***********************************************
        
        
        
        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (i=(mx+1); i<=(mx+mbc); i++)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
  		    tmp = q.get(mx,j,m,ell);                    
                    q.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // BOTTOM BOUNDARY
        // ***********************************************
        for (j=0; j>=(1-mbc); j--)
        {
	    for (i=(2-mbc); i<=(mx+mbc-1); i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,1,m,ell);                    
                    q.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // TOP BOUNDARY
        // ***********************************************
        for (j=(my+1); j<=(my+mbc); j++)
        {
	    for (i=(2-mbc); i<=(mx+mbc-1); i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,my,m,ell);                    
                    q.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
    }
    
}
