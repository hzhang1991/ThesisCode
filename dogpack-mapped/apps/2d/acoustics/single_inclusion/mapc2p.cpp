#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector> 
#include "dogdefs.h"
#include "dog_math.h"
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::abs;




void mapc2p(double& xc,double& yc)
{
   if(true)
   {// cout<<"AHOY"<<endl;
    double x0circ=0.6;
    double y0circ=0.5;
    double r1circ=0.1;
    double r2circ=0.15;
    double zin=12.0;
    double cin=0.3; 

    double xp = xc;
    double yp = yc;


      double   x0 = x0circ;
      double   y0 = y0circ;
      double   r1 = r1circ;
      double   r2 = r2circ;
      double   xc0 = fabs(xc-x0);
      double   yc0 = fabs(yc-y0);
      //cout<<xp<<" "<<xc-x0<<" "<<yp<<" "<<yc0<<endl;
         if ((max(xc0,yc0)-r2)<=1.0e-15) 
            {
   //         cout<<"first "<<xp<<" "<<yp<<endl;
            double xc1 = xc0/r2;
            double yc1 = yc0/r2;
            double d = max(max(xc1,yc1), 1.0e-14);
            d = min(d, 0.9999);
            double d1 = d*r2/sqrt(2.0);
            double d2;
            double rad = r1;             // for constant curvature grid lines
            if (d > r1/r2) 
              { d1 = r1/sqrt(2.0) + (d - r1/r2) * (r2 - r1/sqrt(2.0))/ (1.0 - r1/r2);
               d2 = max(1.0 - d, 1.0e-8);
               rad = r1 * pow(((1.0 - r1/r2) / d2) , (r2/r1 + 0.5));
               }
            double xp2 = d1/d * xc1;
            double yp2 = d1/d * yc1;
            double center = d1 - sqrt(rad*rad - d1*d1);
            if (abs(xc1-d)<1.0e-15) {xp2 = center + sqrt(rad*rad - yp2*yp2);}
            if (abs(yc1-d)<1.0e-15) {yp2 = center + sqrt(rad*rad - xp2*xp2);}
            if ((xc-x0)<1.0e-15){xp2=-fabs(xp2);}
            if ((yc-y0)<1.0e-15){yp2=-fabs(yp2);}
            xp = x0 + xp2;
            yp = y0 + yp2;
       //     cout<<"second "<<xc<<" "<<yc<<" "<<xp<<" "<<yp<<endl;
            xc=xp;
            yc=yp;
            }

            }

}


double phin(double xin,double yin,double x1,double y1,double dx1,double dy1,int a)
{double c0,c1,c2,c3,c4,c5;

double x=xin;
double y=yin;
//double x=2.0*(xin-x1)/dx1-1.0;
//double y=2.0*(yin-y1)/dy1-1.0;
double o6=1.0/6.0;double o3=1.0/3.0;
switch (a)
{
  case 1:
    c0=1.0;
    c1=0.0;
    c2=0.0;
    c3=0.0;
    c4=0.0;
    c5=0.0;
    break;
  case 2:
    c0=0.0;
    c1=sqrt(3.0);
    c2=0.0;
    c3=0.0;
    c4=0.0;
    c5=0.0;
    break;
  case 3:
    c0=0.0;
    c1=0.0;
    c2=sqrt(3.0);
    c3=0.0;
    c4=0.0;
    c5=0.0;
    break;
  case 4:
    c0=0.0;
    c1=0.0;
    c2=0.0;
    c3=3.0;
    c4=0.0;
    c5=0.0;
    break;
  case 5:
    c0=-sqrt(5.0)*0.5;
    c1=0.0;
    c2=0.0;
    c3=0.0;
    c4=1.5*sqrt(5.0);
    c5=0.0;
    break;
  case 6:
    c0=-sqrt(5.0)*0.5;
    c1=0.0;
    c2=0.0;
    c3=0.0;
    c5=1.5*sqrt(5.0);
    c4=0.0;
    break;
  }

double   phi1=c0+c1*x+c2*y+c3*x*y+c4*x*x+c5*y*y;
//cout<<"here= "<<(x-x1)/dx1<<" "<<(y-y1)/dx1<<" "<<c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<y1<<" "<<dy1<<endl;
//cout<<"phi1= "<<phi1<<" "<<15.0+440.0*y+3200.0*y*y<<endl;
   return phi1;}

double dphinx(double xin,double yin,double x1,double y1,double dx1,double dy1,int a)
{double c0,c1,c2,c3,c4,c5;
double o6=1.0/6.0;double o3=1.0/3.0;
double x=2.0*(xin-x1)/dx1-1.0;
double y=2.0*(yin-y1)/dy1-1.0;

switch (a)
{
  case 1:
    c0=-o6;
    c1=-o6;
    c2=-o6;
    c3=o3;
    c4=o6;
    c5=o3;
    break;
  case 2:
    c0=1.0;
    c1=0.0;
    c2=-o3;
    c3=0.0;
    c4=-2.0*o3;
    c5=-2.0*o3;
    break;
  case 3:
    c0=-o6;
    c1=o6;
    c2=-o6;
    c3=-o3;
    c4=o6;
    c5=o3;
    break;
  case 4:
    c0=o6;
    c1=o3;
    c2=o6;
    c3=o3;
    c4=o3;
    c5=-o3;
    break;
  case 5:
    c0=0.0;
    c1=0.0;
    c2=o3;
    c3=0.0;
    c4=-o3;
    c5=2.0*o3;
    break;
  case 6:
    c0=o6;
    c1=-o3;
    c2=o6;
    c3=-o3;
    c4=o3;
    c5=-o3;
    break;
  }


double   phi1=c1+c3*y+2.0*c4*x;
//cout<<"here= "<<(x-x1)/dx1<<" "<<(y-y1)/dx1<<" "<<c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<y1<<" "<<dy1<<endl;
//cout<<"phi1= "<<phi1<<" "<<15.0+440.0*y+3200.0*y*y<<endl;
   return 2.0*phi1/dx1;}
double qin(double xin,double yin,double x1,double y1,double dx1,double dy1,double a[6])
{
double qval=0.0;
for(int i1=1;i1<=6;i1++)
{
qval=qval+phin(xin,yin,x1,y1,dx1,dy1,i1)*a[i1-1];
}
return qval;
}

double dphinxx(double xin,double yin,double x1,double y1,double dx1,double dy1,int a)
{double c0,c1,c2,c3,c4,c5;
double o6=1.0/6.0;double o3=1.0/3.0;
double x=2.0*(xin-x1)/dx1-1.0;
double y=2.0*(yin-y1)/dy1-1.0;

switch (a)
{
  case 1:
    c0=-o6;
    c1=-o6;
    c2=-o6;
    c3=o3;
    c4=o6;
    c5=o3;
    break;
  case 2:
    c0=1.0;
    c1=0.0;
    c2=-o3;
    c3=0.0;
    c4=-2.0*o3;
    c5=-2.0*o3;
    break;
  case 3:
    c0=-o6;
    c1=o6;
    c2=-o6;
    c3=-o3;
    c4=o6;
    c5=o3;
    break;
  case 4:
    c0=o6;
    c1=o3;
    c2=o6;
    c3=o3;
    c4=o3;
    c5=-o3;
    break;
  case 5:
    c0=0.0;
    c1=0.0;
    c2=o3;
    c3=0.0;
    c4=-o3;
    c5=2.0*o3;
    break;
  case 6:
    c0=o6;
    c1=-o3;
    c2=o6;
    c3=-o3;
    c4=o3;
    c5=-o3;
    break;
  }
double   phi1=2.0*c4;
//cout<<"here= "<<(x-x1)/dx1<<" "<<(y-y1)/dx1<<" "<<c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<y1<<" "<<dy1<<endl;
//cout<<"phi1= "<<phi1<<" "<<15.0+440.0*y+3200.0*y*y<<endl;
   return 4.0*phi1/dx1/dx1;}


double dphiny(double xin,double yin,double x1,double y1,double dx1,double dy1,int a)
{double c0,c1,c2,c3,c4,c5;
double o6=1.0/6.0;double o3=1.0/3.0;
double x=2.0*(xin-x1)/dx1-1.0;
double y=2.0*(yin-y1)/dy1-1.0;

switch (a)
{
  case 1:
    c0=-o6;
    c1=-o6;
    c2=-o6;
    c3=o3;
    c4=o6;
    c5=o3;
    break;
  case 2:
    c0=1.0;
    c1=0.0;
    c2=-o3;
    c3=0.0;
    c4=-2.0*o3;
    c5=-2.0*o3;
    break;
  case 3:
    c0=-o6;
    c1=o6;
    c2=-o6;
    c3=-o3;
    c4=o6;
    c5=o3;
    break;
  case 4:
    c0=o6;
    c1=o3;
    c2=o6;
    c3=o3;
    c4=o3;
    c5=-o3;
    break;
  case 5:
    c0=0.0;
    c1=0.0;
    c2=o3;
    c3=0.0;
    c4=-o3;
    c5=2.0*o3;
    break;
  case 6:
    c0=o6;
    c1=-o3;
    c2=o6;
    c3=-o3;
    c4=o3;
    c5=-o3;
    break;
  }
double   phi1=c2+c3*x+2.0*c5*y;
//cout<<"here= "<<(x-x1)/dx1<<" "<<(y-y1)/dx1<<" "<<c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<y1<<" "<<dy1<<endl;
//cout<<"phi1= "<<phi1<<" "<<15.0+440.0*y+3200.0*y*y<<endl;
   return 2.0*phi1/dy1;}
double dphinyy(double xin,double yin,double x1,double y1,double dx1,double dy1,int a)
{double c0,c1,c2,c3,c4,c5;
double o6=1.0/6.0;double o3=1.0/3.0;
double x=2.0*(xin-x1)/dx1-1.0;
double y=2.0*(yin-y1)/dy1-1.0;

switch (a)
{
  case 1:
    c0=-o6;
    c1=-o6;
    c2=-o6;
    c3=o3;
    c4=o6;
    c5=o3;
    break;
  case 2:
    c0=1.0;
    c1=0.0;
    c2=-o3;
    c3=0.0;
    c4=-2.0*o3;
    c5=-2.0*o3;
    break;
  case 3:
    c0=-o6;
    c1=o6;
    c2=-o6;
    c3=-o3;
    c4=o6;
    c5=o3;
    break;
  case 4:
    c0=o6;
    c1=o3;
    c2=o6;
    c3=o3;
    c4=o3;
    c5=-o3;
    break;
  case 5:
    c0=0.0;
    c1=0.0;
    c2=o3;
    c3=0.0;
    c4=-o3;
    c5=2.0*o3;
    break;
  case 6:
    c0=o6;
    c1=-o3;
    c2=o6;
    c3=-o3;
    c4=o3;
    c5=-o3;
    break;
  }
double   phi1=2.0*c5;
//cout<<"here= "<<(x-x1)/dx1<<" "<<(y-y1)/dx1<<" "<<c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<y1<<" "<<dy1<<endl;
//cout<<"phi1= "<<phi1<<" "<<15.0+440.0*y+3200.0*y*y<<endl;
   return 4.0*phi1/dy1/dy1;}

double dphinxy(double xin,double yin,double x1,double y1,double dx1,double dy1,int a)
{double c0,c1,c2,c3,c4,c5;
double o6=1.0/6.0;double o3=1.0/3.0;
double x=2.0*(xin-x1)/dx1-1.0;
double y=2.0*(yin-y1)/dy1-1.0;

switch (a)
{
  case 1:
    c0=-o6;
    c1=-o6;
    c2=-o6;
    c3=o3;
    c4=o6;
    c5=o3;
    break;
  case 2:
    c0=1.0;
    c1=0.0;
    c2=-o3;
    c3=0.0;
    c4=-2.0*o3;
    c5=-2.0*o3;
    break;
  case 3:
    c0=-o6;
    c1=o6;
    c2=-o6;
    c3=-o3;
    c4=o6;
    c5=o3;
    break;
  case 4:
    c0=o6;
    c1=o3;
    c2=o6;
    c3=o3;
    c4=o3;
    c5=-o3;
    break;
  case 5:
    c0=0.0;
    c1=0.0;
    c2=o3;
    c3=0.0;
    c4=-o3;
    c5=2.0*o3;
    break;
  case 6:
    c0=o6;
    c1=-o3;
    c2=o6;
    c3=-o3;
    c4=o3;
    c5=-o3;
    break;
  }
double   phi1=c3;
//cout<<"here= "<<(x-x1)/dx1<<" "<<(y-y1)/dx1<<" "<<c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<y1<<" "<<dy1<<endl;
//cout<<"phi1= "<<phi1<<" "<<15.0+440.0*y+3200.0*y*y<<endl;
   return 4.0*phi1/dy1/dx1;}



double jacobian(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4)
{double xi=xiin;//(xiin+1.0)/2.0;
double eta=etain;//(etain+1.0)/2.0;
double Mx=(x1-x2+x3-x4);
double My=(y1-y2+y3-y4);
double area=(0.1250000000*y4*x3-0.1250000000*x4*y3-0.1250000000*y4*x1+0.1250000000*x4*y1+0.1250000000*y3*x2-0.1250000000*y1*x2-0.1250000000*x3*y2+0.1250000000*x1*y2)*4.0;
double J=abs((0.5*(x2-x1)+1.0/4.0*Mx*(eta+1.))*(0.5*(y4-y1)+1.0/4.0*My*(xi+1.0))-(0.5*(x4-x1)+1.0/4.0*Mx*(xi+1.0))*(0.5*(y2-y1)+1.0/4.0*My*(eta+1.0)));
J=J/area;
//cout<<"J= "<<J<<endl;
return J;
}


double dxidx(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4)
{double xi=xiin;//(xiin+1.0)/2.0;
double eta=etain;//(etain+1.0)/2.0;
double Mx=(x1-x2+x3-x4);
double My=(y1-y2+y3-y4);

double J=(x2-x1+Mx*eta)*(y4-y1+My*xi)-(x4-x1+Mx*xi)*(y2-y1+My*eta);
double out=-(2.0*(y4-y1+y1*xi-y2*xi-y2+y3*xi+y3-y4*xi))/(-y4*x3+x4*y3+y4*x1-x4*y1-y3*x2+y1*x2+x3*y2-x1*y2+x1*y3*xi+x3*eta*y1+x2*y4*xi-x2*y3*xi+x4*y1*xi+x4*eta*y3+x1*eta*y2+x3*y2*xi-x3*eta*y4-x3*y1*xi-x4*y2*xi-x1*eta*y3-x4*eta*y2+x2*eta*y4-x2*eta*y1-x1*y4*xi);//1.0/J*(y4-y1+My*xi);
return out;
}


double dxidy(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4)
{double xi=xiin;//(xiin+1.0)/2.0;
double eta=etain;//(etain+1.0)/2.0;
double Mx=(x1-x2+x3-x4);
double My=(y1-y2+y3-y4);

double J=(x2-x1+Mx*eta)*(y4-y1+My*xi)-(x4-x1+Mx*xi)*(y2-y1+My*eta);
double out=(2.0*(x4-x1+x1*xi-x2*xi-x2+x3*xi+x3-x4*xi))/(-y4*x3+x4*y3+y4*x1-x4*y1-y3*x2+y1*x2+x3*y2-x1*y2+x1*y3*xi+x3*eta*y1+x2*y4*xi-x2*y3*xi+x4*y1*xi+x4*eta*y3+x1*eta*y2+x3*y2*xi-x3*eta*y4-x3*y1*xi-x4*y2*xi-x1*eta*y3-x4*eta*y2+x2*eta*y4-x2*eta*y1-x1*y4*xi);//-1.0/J*(x4-x1+Mx*xi);
return out;
}


double detadx(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4)
{double xi=xiin;//(xiin+1.0)/2.0;
double eta=etain;//(etain+1.0)/2.0;
double Mx=(x1-x2+x3-x4);
double My=(y1-y2+y3-y4);

double J=(x2-x1+Mx*eta)*(y4-y1+My*xi)-(x4-x1+Mx*xi)*(y2-y1+My*eta);
double out=(2.0*(y2-y1+y1*eta-y2*eta+y3*eta+y3-y4*eta-y4))/(-y4*x3+x4*y3+y4*x1-x4*y1-y3*x2+y1*x2+x3*y2-x1*y2+x1*y3*xi+x3*eta*y1+x2*y4*xi-x2*y3*xi+x4*y1*xi+x4*eta*y3+x1*eta*y2+x3*y2*xi-x3*eta*y4-x3*y1*xi-x4*y2*xi-x1*eta*y3-x4*eta*y2+x2*eta*y4-x2*eta*y1-x1*y4*xi);//-1.0/J*(y2-y1+My*eta);
return out;
}


double detady(double xiin,double etain,double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4)
{double xi=xiin;//(xiin+1.0)/2.0;
double eta=etain;//(etain+1.0)/2.0;
double Mx=(x1-x2+x3-x4);
double My=(y1-y2+y3-y4);

double J=(x2-x1+Mx*eta)*(y4-y1+My*xi)-(x4-x1+Mx*xi)*(y2-y1+My*eta);
double out=-(2.0*(x2-x1+x1*eta-x2*eta+x3*eta+x3-x4*eta-x4))/(-y4*x3+x4*y3+y4*x1-x4*y1-y3*x2+y1*x2+x3*y2-x1*y2+x1*y3*xi+x3*eta*y1+x2*y4*xi-x2*y3*xi+x4*y1*xi+x4*eta*y3+x1*eta*y2+x3*y2*xi-x3*eta*y4-x3*y1*xi-x4*y2*xi-x1*eta*y3-x4*eta*y2+x2*eta*y4-x2*eta*y1-x1*y4*xi);//1.0/J*(x2-x1+Mx*eta);
return out;
}



std::vector<double> returnright(double x1l,double y1l,double x2l,double y2l,double x3l,double y3l,double x4l,double y4l,int l1)
{

  mapc2p(x1l,y1l);
  mapc2p(x2l,y2l);
  mapc2p(x3l,y3l);
  mapc2p(x4l,y4l);

  double areal1,areal2,arear1,arear2;
  double x1,y1,x2,y2,x3,y3;
  
     
  x1 =x1l;
  y1 =y1l;
  x2 =x2l;
  y2 =y2l;
  x3 =x4l;
  y3 =y4l; 
  areal1=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

  x1 =x3l;
  y1 =y3l;
  x2 =x2l;
  y2 =y2l;
  x3 =x4l;
  y3 =y4l; 
  arear1=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

  x1 =x3l;
  y1 =y3l;
  x2 =x1l;
  y2 =y1l;
  x3 =x4l;
  y3 =y4l; 
  areal2=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

  x1 =x3l;
  y1 =y3l;
  x2 =x2l;
  y2 =y2l;
  x3 =x1l;
  y3 =y1l; 
  arear2=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);
  std::vector<double> x;x.reserve(10);
  if(min(abs(areal1),abs(arear1))>min(abs(areal2),abs(arear2))-1.0e-15)
 {
  x.push_back(areal1);x.push_back(arear1);
  x.push_back(x3l);x.push_back(y3l);
  x.push_back(x2l);x.push_back(y2l);
  x.push_back(x4l);x.push_back(y4l);
  x.push_back(areal2);x.push_back(arear2);
 }
 else
  {
   if(l1==2)
   {
    x.push_back(areal2);x.push_back(arear2);
    x.push_back(x3l);x.push_back(y3l);
    x.push_back(x2l);x.push_back(y2l);
    x.push_back(x1l);x.push_back(y1l);
    x.push_back(areal1);x.push_back(arear1);
   }
   if(l1==1)
    {
    x.push_back(arear2);x.push_back(areal2);
    x.push_back(x3l);x.push_back(y3l);
    x.push_back(x1l);x.push_back(y1l);
    x.push_back(x4l);x.push_back(y4l);
    x.push_back(areal1);x.push_back(arear1);
    }

   }

  //if(min(areal1,arear1)<min(areal2,arear2))
  //{x[
  return x;  
}

std::vector<double> returnleft(double x1l,double y1l,double x2l,double y2l,double x3l,double y3l,double x4l,double y4l,int l1)
{

  mapc2p(x1l,y1l);
  mapc2p(x2l,y2l);
  mapc2p(x3l,y3l);
  mapc2p(x4l,y4l);

  double areal1,areal2,arear1,arear2;
  double x1,y1,x2,y2,x3,y3;
     
  x1 =x1l;
  y1 =y1l;
  x2 =x2l;
  y2 =y2l;
  x3 =x4l;
  y3 =y4l; 
  areal1=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

  x1 =x3l;
  y1 =y3l;
  x2 =x2l;
  y2 =y2l;
  x3 =x4l;
  y3 =y4l; 
  arear1=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

  x1 =x3l;
  y1 =y3l;
  x2 =x1l;
  y2 =y1l;
  x3 =x4l;
  y3 =y4l; 
  areal2=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

  x1 =x3l;
  y1 =y3l;
  x2 =x2l;
  y2 =y2l;
  x3 =x1l;
  y3 =y1l; 
  arear2=1.0/2.0*(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);
  std::vector<double> x;x.reserve(10);
 if(min(abs(areal1),abs(arear1))>=min(abs(areal2),abs(arear2))-1.0e-15)
  {x.push_back(areal1);x.push_back(arear1);
  x.push_back(x1l);x.push_back(y1l);
  x.push_back(x2l);x.push_back(y2l);
  x.push_back(x4l);x.push_back(y4l);
  x.push_back(areal2);x.push_back(arear2);}
 else
  {
   if(l1==2)
   {
    x.push_back(areal2);x.push_back(arear2);
    x.push_back(x3l);x.push_back(y3l);
    x.push_back(x1l);x.push_back(y1l);
    x.push_back(x4l);x.push_back(y4l);
    x.push_back(areal1);x.push_back(arear1);
   }
   if(l1==1)
    {
    x.push_back(arear2);x.push_back(areal2);
    x.push_back(x3l);x.push_back(y3l);
    x.push_back(x2l);x.push_back(y2l);
    x.push_back(x1l);x.push_back(y1l);
    x.push_back(areal1);x.push_back(arear1);
    }

   }
 
  //if(min(areal1,arear1)<min(areal2,arear2))
  //{x[
  return x;  
}

double centermass(double px1,double py1,double px2,double py2,double px3,double py3,double px4,double py4,double dx1,double dy1,double x0,double y0,int degree)
{
double coeff;
 double almostarea= abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) / 0.2e1 +  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) / 0.2e1;
double area=almostarea;
 switch (degree)
 {
 case 1:
  coeff=0.8333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  px4 + 0.8333333337e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  px1 + 0.8333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  px2 + 0.2500000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 + 0.5000000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  x0 + 0.8333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  px4 + 0.8333333337e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  px3 + 0.8333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  px2 + 0.2500000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 + 0.5000000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  x0;

 break;
case 2:
coeff=0.8333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  py4 + 0.8333333337e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  py1 + 0.8333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  py2 + 0.2500000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 + 0.5000000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  y0 + 0.8333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  py4 + 0.8333333337e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  py3 + 0.8333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  py2 + 0.2500000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 + 0.5000000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  y0;

 break;

}

coeff=coeff/area;
return coeff;

}

/*double coefficient(double px1,double py1,double px2,double py2,double px3,double py3,double px4,double py4,double dx1,double dy1,double Ax,double Ay,int degree)
{

 double coeff;
 double almostarea= abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) / 0.2e1 +  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) / 0.2e1;
 double area=almostarea;///4.0*dx1*dy1;
  
 switch (degree)
 {
   case 2:
       coeff =  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1) / 0.6e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px2) / 0.6e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px4) / 0.6e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2) / 0.6e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3) / 0.6e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px4) / 0.6e1;

    break;

   case 3:
       coeff =  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py1) / 0.6e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py2) / 0.6e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py4) / 0.6e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py2) / 0.6e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py3) / 0.6e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py4) / 0.6e1;
       break;

   case 4:
       coeff =   (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ax * Ay) / 0.2e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2 * Ay) / 0.6e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3 * Ay) / 0.6e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px4 * Ay) / 0.6e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ax * Ay) / 0.2e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ax * py1) / 0.6e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ax * py2) / 0.6e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ax * py4) / 0.6e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ax * py2) / 0.6e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ax * py3) / 0.6e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ax * py4) / 0.6e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * Ay) / 0.6e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px2 * Ay) / 0.6e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px4 * Ay) / 0.6e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * py1) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * py2) / 0.24e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * py4) / 0.24e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px2 * py1) / 0.24e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px2 * py2) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px2 * py4) / 0.24e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px4 * py1) / 0.24e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px4 * py2) / 0.24e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px4 * py4) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2 * py2) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2 * py3) / 0.24e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2 * py4) / 0.24e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3 * py2) / 0.24e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3 * py3) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3 * py4) / 0.24e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px4 * py2) / 0.24e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px4 * py3) / 0.24e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px4 * py4) / 0.12e2;
       break;


   case 5:
     coeff =  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ax * Ax) / 0.2e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ax * Ax) / 0.2e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * Ax) / 0.3e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ax * px2) / 0.3e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px4 * Ax) / 0.3e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ax * px2) / 0.3e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3 * Ax) / 0.3e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px4 * Ax) / 0.3e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * px1) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * px2) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px1 * px4) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px2 * px2) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px2 * px4) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * px4 * px4) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2 * px2) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2 * px3) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px2 * px4) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3 * px3) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px3 * px4) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * px4 * px4) / 0.12e2;
     break;



  case 6:
     coeff =  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ay * Ay) / 0.2e1 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ay * Ay) / 0.2e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py1 * Ay) / 0.3e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * Ay * py2) / 0.3e1 -  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py4 * Ay) / 0.3e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * Ay * py2) / 0.3e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py3 * Ay) / 0.3e1 -  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py4 * Ay) / 0.3e1 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py1 * py1) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py1 * py2) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py1 * py4) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py2 * py2) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py2 * py4) / 0.12e2 +  (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * py4 * py4) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py2 * py2) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py2 * py3) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py2 * py4) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py3 * py3) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py3 * py4) / 0.12e2 +  (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * py4 * py4) / 0.12e2;
     break;

 }

coeff=coeff/area;
return coeff;
}*/

double coefficient(double px1,double py1,double px2,double py2,double px3,double py3,double px4,double py4,double dx1,double dy1,double Ax,double Ay,double x0,double y0,int degree)
{

 double coeff;
 double almostarea= abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) / 0.2e1 +  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) / 0.2e1;
 double area=almostarea;///4.0*dx1*dy1;
  
 switch (degree)
 {
   case 2:
       coeff =  0.833333333333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  px4 + 0.833333333333333337e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  px1 + 0.833333333333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 *  px2 + 0.250000000000000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dx1 + 0.500000000000000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  x0 + 0.833333333333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  px4 + 0.833333333333333337e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  px3 + 0.833333333333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 *  px2 + 0.250000000000000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dx1 + 0.500000000000000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  x0;

    break;

   case 3:
       coeff =  0.833333333333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dy1 *  py4 + 0.833333333333333337e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dy1 *  py1 + 0.833333333333333333e-1 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dy1 *  py2 + 0.250000000000000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  dy1 + 0.500000000000000000e0 *  abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) *  y0 + 0.833333333333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dy1 *  py4 + 0.833333333333333337e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dy1 *  py3 + 0.833333333333333333e-1 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dy1 *  py2 + 0.250000000000000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  dy1 + 0.500000000000000000e0 *  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) *  y0;
       break;

   case 4:
       coeff =  0.104166666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px1 * (double) py4 + 0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px1 * (double) py1 + 0.104166666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px1 * (double) py2 + 0.416666666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) dx1 * (double) py2 - 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) Ax * (double) py2 + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px2 * (double) y0 + 0.416666666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px2 * (double) dy1 - 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px2 * (double) Ay + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) x0 * (double) py4 + 0.416666666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) dx1 * (double) py4 - 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) Ax * (double) py4 + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px4 * (double) y0 + 0.416666666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px4 * (double) dy1 - 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px4 * (double) Ay + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px3 * (double) y0 + 0.416666666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px3 * (double) dy1 - 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) px3 * (double) Ay + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) x0 * (double) py2 + 0.416666666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) dx1 * (double) py2 - 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) Ax * (double) py2 + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px2 * (double) y0 + 0.416666666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px2 * (double) dy1 - 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px2 * (double) Ay + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) x0 * (double) py4 + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) x0 * (double) py1 + 0.416666666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) dx1 * (double) py4 + 0.416666666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) dx1 * (double) py1 - 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) Ax * (double) py4 - 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) Ax * (double) py1 + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px4 * (double) y0 + 0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px4 * (double) py4 + 0.104166666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px4 * (double) py1 + 0.104166666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px2 * (double) py4 + 0.104166666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px2 * (double) py1 + 0.104166666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px4 * (double) py2 + 0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 * (double) px2 * (double) py2 + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px4 * (double) py4 + 0.104166666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px3 * (double) py4 + 0.104166666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px2 * (double) py4 + 0.104166666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px4 * (double) py2 + 0.104166666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px3 * (double) py2 + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px2 * (double) py2 + 0.104166666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px4 * (double) py3 + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px3 * (double) py3 + 0.104166666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 * (double) px2 * (double) py3 + 0.416666666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px4 * (double) dy1 - 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px4 * (double) Ay + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) x0 * (double) py2 + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) x0 * (double) py3 + 0.416666666666666667e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) dx1 * (double) py3 - 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) Ax * (double) py3 + 0.125000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) dy1 + 0.125000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) dy1 + 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) x0 * (double) y0 + 0.250000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) x0 * (double) dy1 - 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) x0 * (double) Ay + 0.250000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) y0 - 0.250000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) Ay - 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) Ax * (double) y0 - 0.250000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) Ax * (double) dy1 + 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) x0 * (double) y0 + 0.250000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) x0 * (double) dy1 - 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) x0 * (double) Ay + 0.250000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) y0 - 0.250000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dx1 * (double) Ay - 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) Ax * (double) y0 - 0.250000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) Ax * (double) dy1 + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px1 * (double) y0 + 0.416666666666666667e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px1 * (double) dy1 - 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dx1 * (double) px1 * (double) Ay + 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) Ax * (double) Ay + 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) Ax * (double) Ay;
       break;


   case 5:
     coeff = (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (0.208333333333333336e-1 * (double) (dx1 * dx1) * (double) (px1 * px1) + 0.208333333333333330e-1 * (double) (dx1 * dx1) * (double) px1 * (double) px2 + 0.208333333333333336e-1 * (double) (dx1 * dx1) * (double) px1 * (double) px4 + 0.208333333333333336e-1 * (double) (dx1 * dx1) * (double) (px2 * px2) + 0.208333333333333330e-1 * (double) (dx1 * dx1) * (double) px4 * (double) px2 + 0.208333333333333336e-1 * (double) (dx1 * dx1) * (double) (px4 * px4) - 0.166666666666666667e0 * (double) dx1 * (double) px1 * (double) Ax - 0.166666666666666666e0 * (double) dx1 * (double) Ax * (double) px2 - 0.166666666666666667e0 * (double) dx1 * (double) px4 * (double) Ax + 0.833333333333333335e-1 * (double) (dx1 * dx1) * (double) px1 + 0.833333333333333330e-1 * (double) (dx1 * dx1) * (double) px2 + 0.833333333333333335e-1 * (double) (dx1 * dx1) * (double) px4 + 0.166666666666666667e0 * (double) x0 * (double) dx1 * (double) px1 + 0.166666666666666666e0 * (double) x0 * (double) dx1 * (double) px2 + 0.166666666666666667e0 * (double) x0 * (double) dx1 * (double) px4 + 0.500000000000000000e0 * (double) Ax * (double) Ax - 0.500000000000000000e0 * (double) dx1 * (double) Ax - (double) (x0 * Ax) + 0.125000000000000000e0 * (double) dx1 * (double) dx1 + 0.500000000000000000e0 * (double) x0 * (double) dx1 + 0.500000000000000000e0 * (double) x0 * (double) x0) + (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (0.208333333333333333e-1 * (double) (dx1 * dx1) * (double) (px2 * px2) + 0.125000000000000000e0 * (double) dx1 * (double) dx1 + 0.500000000000000000e0 * (double) x0 * (double) x0 + 0.208333333333333333e-1 * (double) (dx1 * dx1) * (double) (px3 * px3) + 0.833333333333333333e-1 * (double) (dx1 * dx1) * (double) px3 + 0.833333333333333333e-1 * (double) (dx1 * dx1) * (double) px2 + 0.500000000000000000e0 * (double) x0 * (double) dx1 - (double) (x0 * Ax) + 0.208333333333333333e-1 * (double) (dx1 * dx1) * (double) (px4 * px4) + 0.833333333333333333e-1 * (double) (dx1 * dx1) * (double) px4 - 0.500000000000000000e0 * (double) dx1 * (double) Ax + 0.500000000000000000e0 * (double) Ax * (double) Ax - 0.166666666666666667e0 * (double) dx1 * (double) px4 * (double) Ax + 0.166666666666666667e0 * (double) x0 * (double) dx1 * (double) px2 + 0.166666666666666667e0 * (double) x0 * (double) dx1 * (double) px4 + 0.208333333333333333e-1 * (double) (dx1 * dx1) * (double) px4 * (double) px2 - 0.166666666666666667e0 * (double) dx1 * (double) Ax * (double) px2 + 0.166666666666666667e0 * (double) x0 * (double) dx1 * (double) px3 - 0.166666666666666667e0 * (double) dx1 * (double) px3 * (double) Ax + 0.208333333333333333e-1 * (double) (dx1 * dx1) * (double) px3 * (double) px2 + 0.208333333333333334e-1 * (double) (dx1 * dx1) * (double) px4 * (double) px3);
     break;



  case 6:
     coeff =   0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) py1 * (double) py2 + 0.208333333333333334e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) py1 * (double) py4 + 0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) py4 * (double) py2 - 0.166666666666666667e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) py1 * (double) Ay - 0.166666666666666667e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) Ay * (double) py2 - 0.166666666666666667e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) py4 * (double) Ay + 0.166666666666666667e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) y0 * (double) dy1 * (double) py1 + 0.166666666666666667e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) y0 * (double) dy1 * (double) py2 + 0.166666666666666667e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) y0 * (double) dy1 * (double) py4 + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) py3 * (double) py2 + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) py4 * (double) py2 + 0.208333333333333334e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) py3 * (double) py4 - 0.166666666666666667e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) Ay * (double) py2 - 0.166666666666666667e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) py3 * (double) Ay - 0.166666666666666667e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) py4 * (double) Ay + 0.166666666666666667e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) y0 * (double) dy1 * (double) py2 + 0.166666666666666667e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) y0 * (double) dy1 * (double) py3 + 0.166666666666666667e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) y0 * (double) dy1 * (double) py4 + 0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) (py1 * py1) + 0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) (py2 * py2) + 0.208333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) (py4 * py4) + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) py1 + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) py2 + 0.833333333333333333e-1 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) * (double) py4 - 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) dy1 * (double) Ay - (double) (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * y0 * Ay) + 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) y0 * (double) dy1 + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) (py2 * py2) + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) (py3 * py3) + 0.208333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) (py4 * py4) + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) py2 + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) py3 + 0.833333333333333333e-1 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) * (double) py4 - 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) dy1 * (double) Ay - (double) (abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * y0 * Ay) + 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) y0 * (double) dy1 + 0.125000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (dy1 * dy1) + 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (y0 * y0) + 0.125000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (dy1 * dy1) + 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (y0 * y0) + 0.500000000000000000e0 * (double) abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * (double) (Ay * Ay) + 0.500000000000000000e0 * (double) abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * (double) (Ay * Ay);
     break;

 }

coeff=coeff/area;
return coeff;
}

double averages(double px1,double py1,double px2,double py2,double px3,double py3,double px4,double py4,double dx1,double dy1,double a[6],int degree)
{
///calculate the averaged derivative of a given degree...

   double ave=0.0;
    double almostarea= abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) / 0.2e1 +  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) / 0.2e1;
   double area=almostarea;///4.0*dx1*dy1;
 //  cout<<"AREA! "<<area<<endl;
   switch (degree)
   {

   case 1: //average
    
ave =  ((-2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py3 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py4 * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py4 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py4 * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * px1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px1 * px1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * px1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px1 * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * px1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * px2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px2 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px2 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px4 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px4 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 * py1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 * py1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py2 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py4 * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py4 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py4 * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py4 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * px2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px2 * px2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * px2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px2 * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * px2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * px3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px3 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px3 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px4 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px4 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 * py2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 * py2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 * py2 + 36 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] - 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] + 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] + 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] - 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] + 36 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] - 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] + 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] + 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] - 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * px2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px1 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px1 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px1 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px1 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px2 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px2 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py2 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py2 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * px3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px2 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px2 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px2 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px2 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px3 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px3 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py3 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py3 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 * py4) / (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3))) / 0.36e2;
 
break;

   case 2:  ///x deriv
    

ave =  ((-8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px2 - 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px3 - 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px4 - 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px1 - 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px2 - 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py4 + 3 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] + 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] - 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] - 3 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] + 3 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] + 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] - 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] - 3 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0]) / (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3))) / 0.18e2;
ave=ave*2.0/dx1;
break;

  case 3: //y deriv

ave =  ((2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 - 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 + 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 - 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 + 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py4 - 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py4 + 8 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 - 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 + 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 - 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 + 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py4 - 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py4 + 8 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py4 - 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] - 3 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] + 3 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] + 3 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] - 3 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] - 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] - 3 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] + 3 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] + 3 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] - 3 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] + 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] + 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4]) / (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3))) / 0.18e2;
ave=ave*2.0/dy1;
break;

   case 4: //xy deriv

   ave =  a[0] / 0.3e1 - a[2] / 0.3e1 + a[3] / 0.3e1 - a[5] / 0.3e1;
   ave = ave*4.0/dx1/dy1;
  break;

   case 5: //xx deriv

   ave =  a[0] / 0.3e1 - 0.4e1 / 0.3e1 * a[1] + a[2] / 0.3e1 + 0.2e1 / 0.3e1 * a[3] - 0.2e1 / 0.3e1 * a[4] + 0.2e1 / 0.3e1 * a[5];
   ave= ave*4.0/dx1/dx1;
break;

  case 6: //yy deriv

  ave =  0.2e1 / 0.3e1 * a[0] - 0.4e1 / 0.3e1 * a[1] + 0.2e1 / 0.3e1 * a[2] - 0.2e1 / 0.3e1 * a[3] + 0.4e1 / 0.3e1 * a[4] - 0.2e1 / 0.3e1 * a[5];
  ave=ave*4.0/dy1/dy1;

break;
   case 7:
   ave=area/4.0*dx1*dy1;
   }
   //if(a[0]>0.0 && degree==1)
   //{cout<<"as= "<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<" "<<a[4]<<" "<<a[5]<<endl;
   //cout<<"area= "<<area<<" "<<ave<<" "<<ave/area<<endl;}

   //ave=ave/area;

return ave;
}


double integrate(double oxp1,double oyp1,double oxp2,double oyp2,double oxp3,double oyp3,double oxp4,double oyp4,double q[6])
{
int mpoints = 16;int meqn=1;int maux=1;
int numcoeffs=6;int kmax_fout=6;int kmax_qin=6;int kmax=6;
dTensor2    spts(mpoints,2);
dTensor2    spts1(2*mpoints,2);
dTensor1    wgts(mpoints);
dTensor1    wgts1(2*mpoints);int mcomps_out=1;
        dTensor2     phi(2*mpoints,kmax); // Legendre basis (orthogonal)
            double xp1 = oxp1;
            double yp1= oyp1;
            double xp2 = oxp2;
            double yp2= oyp2;
            double xp3 = oxp3;
            double yp3= oyp3;
            double xp4 = oxp4;
            double yp4= oyp4;

            double xr1 = xp1;
            double yr1= yp1;
            double xr2 = xp2;
            double yr2= yp2;
            double xr3 = xp3;
            double yr3= yp3;
            double xr4 = xp4;
            double yr4= yp4;
            double areal,arear;
            mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);

                spts.set(1,1,   0.000000000000000 );
                spts.set(1,2,   0.000000000000000 );

                spts.set(2,1,   0.125959254959390 );
                spts.set(2,2,   0.125959254959390 );

                spts.set(3,1,  -0.251918509918779 );
                spts.set(3,2,   0.125959254959390 );

                spts.set(4,1,   0.125959254959390 );
                spts.set(4,2,  -0.251918509918779 );

                spts.set(5,1,  -0.162764025581573 );
                spts.set(5,2,  -0.162764025581573 );

                spts.set(6,1,   0.325528051163147 );
                spts.set(6,2,  -0.162764025581573 );

                spts.set(7,1,  -0.162764025581573 );
                spts.set(7,2,   0.325528051163147 );

                spts.set(8,1,  -0.282786105016302 );
                spts.set(8,2,  -0.282786105016302 );

                spts.set(9,1,   0.565572210032605 );
                spts.set(9,2,  -0.282786105016302 );

                spts.set(10,1, -0.282786105016302 );
                spts.set(10,2,  0.565572210032605 );

                spts.set(11,1, -0.324938555923375 );
                spts.set(11,2, -0.070220503698695 );

                spts.set(12,1, -0.324938555923375 );
                spts.set(12,2,  0.395159059622071 );

                spts.set(13,1, -0.070220503698695 );
                spts.set(13,2, -0.324938555923375 );

                spts.set(14,1, -0.070220503698695 );
                spts.set(14,2,  0.395159059622071 );

                spts.set(15,1,  0.395159059622071 );
                spts.set(15,2, -0.324938555923375 );

                spts.set(16,1,  0.395159059622071 );
                spts.set(16,2, -0.070220503698695 );

                wgts.set(1,  0.0721578038388935 );
                wgts.set(2,  0.0475458171336425 );
                wgts.set(3,  0.0475458171336425 );
                wgts.set(4,  0.0475458171336425 );
                wgts.set(5,  0.0516086852673590 );
                wgts.set(6,  0.0516086852673590 );
                wgts.set(7,  0.0516086852673590 );
                wgts.set(8,  0.0162292488115990 );
                wgts.set(9,  0.0162292488115990 );
                wgts.set(10, 0.0162292488115990 );
                wgts.set(11, 0.0136151570872175 );
                wgts.set(12, 0.0136151570872175 );
                wgts.set(13, 0.0136151570872175 );
                wgts.set(14, 0.0136151570872175 );
                wgts.set(15, 0.0136151570872175 );
                wgts.set(16, 0.0136151570872175 );

    for(int l1=1;l1<=2;l1++)
    for (int m1=1; m1<=mpoints; m1++)
    {   int j1=0;
        int m;
        double xi,eta;
        if(l1==1){j1=0;m=m1+j1;
        double s1=spts.get(m1,1);
        double t1=spts.get(m1,2);
        xi = 2.0*s1-1.0/3.0;
        eta = 2.0*t1-1.0/3.0;
       }
       if(l1==2){j1=mpoints;m=m1+j1;
        double s1=spts.get(m1,1);
        double t1=spts.get(m1,2);
        xi = -2.0*t1+1.0/3.0;
        eta = -2.0*s1+1.0/3.0;
        }
        // coordinates
        spts1.set(m,1,xi);
        spts1.set(m,2,eta);
        wgts1.set(m,4.0*wgts.get(m1));    




        }
            double xmax=max(max(xp1,xp2),max(xp3,xp4));double xmin=min(min(xp1,xp2),min(xp3,xp4));
            double ymax=max(max(yp1,yp2),max(yp3,yp4));double ymin=min(min(yp1,yp2),min(yp3,yp4));
            double dx1=xmax-xmin;double dy1=ymax-ymin;
            dTensor2    xpts(2*mpoints,2);
            dTensor2   qvals(2*mpoints,meqn);
            vector<double> jacobian;
            double area=0.0;  
         for (int ini1=1;ini1<=2;ini1++)
         {	  



            // These need to be defined locally
            double x1,y1,x2,y2,x3,y3;

           
            vector<double> x;
            if(ini1==1){
    
            x=returnleft(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,1);
            jacobian.push_back(x.at(0));
            areal=x.at(0);
            arear=x.at(1);
            area=area+abs(areal);
            
        //    cout<<"diff1 "<<xi1-x1<<" "<<eta1-y1<<" "<<xi2-x2<<" "<<eta2-y2<<" "<<xi4-x3<<" "<<eta4-y3<<endl;
            }
            if(ini1==2){
 
            x=returnright(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,1);

            jacobian.push_back(x.at(1));
            areal=x.at(0);
            arear=x.at(1);
            area=area+abs(arear);

    }

        
            x1=x.at(2);y1=x.at(3);
            x2=x.at(4);y2=x.at(5);
            x3=x.at(6);y3=x.at(7);

//cout<<i<<" "<<jm<<" kmax_fout "<<kmax_fout<<" "<<mpoints<<endl;
            double xc = (x1+x2+x3)/3.0;
            double yc = (y1+y2+y3)/3.0;
            double jmat[2][2];
            if(ini1==1){
            jmat[0][0]=copysign(1.0,jacobian[ini1-1])*(y3-y1);jmat[0][1]=copysign(1.0,jacobian[ini1-1])*(y1-y2);
            jmat[1][0]=copysign(1.0,jacobian[ini1-1])*(x1-x3);jmat[1][1]=copysign(1.0,jacobian[ini1-1])*(x2-x1);}
            if(ini1==2){
            jmat[0][0]=copysign(1.0,jacobian[ini1-1])*(y2-y1);jmat[0][1]=copysign(1.0,jacobian[ini1-1])*(y1-y3);
            jmat[1][0]=copysign(1.0,jacobian[ini1-1])*(x1-x2);jmat[1][1]=copysign(1.0,jacobian[ini1-1])*(x3-x1);}
           // cout<<i<<" "<<jm<<" "<<ini1<<" jmat "<<jmat[0][0]<<" "<<jmat[0][1]<<" "<<jmat[1][0]<<" "<<jmat[1][1]<<endl;
            jacobian[ini1-1]=abs(jacobian[ini1-1]);
            // Compute q, aux and fvals at each Gaussian Quadrature point
            // for this current cell indexed by (i,j)
            // Save results into dTensor2 qvals, auxvals and fvals.
            for (int m1=1; m1<=mpoints; m1++)
            {    int j1=0;if(ini1==2){j1=mpoints;}
                 int m=m1+j1;

                // point on the unit triangle
                const double s = spts1.get(m,1);
                const double t = spts1.get(m,2);

            if(ini1==1){
            // point on the physical triangle
            xpts.set(m,1, 0.5*(x2+x3) + 0.5*(x2-x1)*s + 0.5*(x3-x1)*t );
            xpts.set(m,2, 0.5*(y2+y3) + 0.5*(y2-y1)*s + 0.5*(y3-y1)*t );}
            else{
            xpts.set(m,1, 0.5*(x2+x3) + 0.5*(x1-x3)*s + 0.5*(x1-x2)*t );
            xpts.set(m,2, 0.5*(y2+y3) + 0.5*(y1-y3)*s + 0.5*(y1-y2)*t );}
            double xo=xpts.get(m,1);double yo=xpts.get(m,2);
      for (int lm=1;lm<=6;lm++){   
        phi.set(m,lm,dphinxx(xpts.get(m,1),xpts.get(m,2),xmin,ymin,dx1,dy1,lm));
       // cout<<"lm= "<<lm<<" "<<dphiny(xpts.get(m,1),xpts.get(m,2),xmin,ymin,dx1,dy1,lm)<<endl;
      }
        

                // Solution values (q) at each grid point

                    qvals.set(m,1, 0.0 );

                    for (int k=1; k<=kmax_qin; k++)
                    {
                        qvals.set(m,1, qvals.get(m,1) 
                                + phi.get(m,k) * q[k-1] );
                    }
//cout<<setprecision(15)<<"pts "<<2.0*(xpts.get(m,1)-xmin)/dx1-1.0<<" "<<2.0*(xpts.get(m,2)-ymin)/dy1-1.0<<" "<<qvals.get(m,1)<<endl;
//cout<<setprecision(15)<<" "<<(xp1-xmin)/dx1*2.0-1.0<<" "<<(yp1-ymin)/dy1*2.0-1.0<<" "<<(xp2-xmin)/dx1*2.0-1.0<<" "<<(yp2-ymin)/dy1*2.0-1.0<<" "<<(xp3-xmin)/dx1*2.0-1.0<<" "<<(yp3-ymin)/dy1*2.0-1.0<<" "<<(xp4-xmin)/dx1*2.0-1.0<<" "<<(yp4-ymin)/dy1*2.0-1.0<<endl;
            }
          }

            // Call user-supplied function to set fvals
          //  cout<<"stuff= "<<fvals.get(1,1,1)<<" "<<fvals.get(1,1,2)<<" "<<auxvals.get(1,1)<<" "<<qvals.get(1,1)<<endl;
            // Evaluate integral on current cell (project onto Legendre basis) 
            // using Gaussian Quadrature for the integration
double tmp=0.0;
               for (int ini1=1;ini1<=2;ini1++)
               {	  
                    
                    //cout<<"tmp= "<<tmp<<endl;
                    for (int k1=1; k1<=mpoints; k1++)
                    { int k=k1;if(ini1==2){k=k1+mpoints;}
                           // if(i==2 && jm==2){cout<<"f= "<<fvals.get(k,1,1)<<" "<<fvals.get(k,1,2)<<endl;}
                        tmp = tmp + wgts1.get(k)*abs(jacobian[ini1-1])/2.0*qvals.get(k,1);
                     //   cout<<"update "<<qvals.get(k,1)<<" "<<tmp<<endl;
                   // cout<<k<<" "<<wgts1.get(k)<<" "<<phi_x.get(k,m2)<<" "<<fvals.get(k,m1,1)<<endl;
                    }
                    //cout<<i<<" "<<jm<<" l1= "<<ini1<<" m1= "<<m1<<" m2= "<<m2<<" "<<fout2.get(i,jm,m1,m2)<<endl;
                }
 //cout<<"areas "<<jacobian[0]<<" "<<jacobian[1]<<" "<<tmp<<endl;




return tmp/area;
}

double integratexc(double oxp1,double oyp1,double oxp2,double oyp2,double oxp3,double oyp3,double oxp4,double oyp4,double q[6])
{
int mpoints = 16;int meqn=1;int maux=1;
int numcoeffs=6;int kmax_fout=6;int kmax_qin=6;int kmax=6;
dTensor2    spts(mpoints,2);
dTensor2    spts1(2*mpoints,2);
dTensor1    wgts(mpoints);
dTensor1    wgts1(2*mpoints);int mcomps_out=1;
        dTensor2     phi(2*mpoints,kmax); // Legendre basis (orthogonal)
            double xp1 = oxp1;
            double yp1= oyp1;
            double xp2 = oxp2;
            double yp2= oyp2;
            double xp3 = oxp3;
            double yp3= oyp3;
            double xp4 = oxp4;
            double yp4= oyp4;

            double xr1 = xp1;
            double yr1= yp1;
            double xr2 = xp2;
            double yr2= yp2;
            double xr3 = xp3;
            double yr3= yp3;
            double xr4 = xp4;
            double yr4= yp4;
            double areal,arear;
            mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);

                spts.set(1,1,   0.000000000000000 );
                spts.set(1,2,   0.000000000000000 );

                spts.set(2,1,   0.125959254959390 );
                spts.set(2,2,   0.125959254959390 );

                spts.set(3,1,  -0.251918509918779 );
                spts.set(3,2,   0.125959254959390 );

                spts.set(4,1,   0.125959254959390 );
                spts.set(4,2,  -0.251918509918779 );

                spts.set(5,1,  -0.162764025581573 );
                spts.set(5,2,  -0.162764025581573 );

                spts.set(6,1,   0.325528051163147 );
                spts.set(6,2,  -0.162764025581573 );

                spts.set(7,1,  -0.162764025581573 );
                spts.set(7,2,   0.325528051163147 );

                spts.set(8,1,  -0.282786105016302 );
                spts.set(8,2,  -0.282786105016302 );

                spts.set(9,1,   0.565572210032605 );
                spts.set(9,2,  -0.282786105016302 );

                spts.set(10,1, -0.282786105016302 );
                spts.set(10,2,  0.565572210032605 );

                spts.set(11,1, -0.324938555923375 );
                spts.set(11,2, -0.070220503698695 );

                spts.set(12,1, -0.324938555923375 );
                spts.set(12,2,  0.395159059622071 );

                spts.set(13,1, -0.070220503698695 );
                spts.set(13,2, -0.324938555923375 );

                spts.set(14,1, -0.070220503698695 );
                spts.set(14,2,  0.395159059622071 );

                spts.set(15,1,  0.395159059622071 );
                spts.set(15,2, -0.324938555923375 );

                spts.set(16,1,  0.395159059622071 );
                spts.set(16,2, -0.070220503698695 );

                wgts.set(1,  0.0721578038388935 );
                wgts.set(2,  0.0475458171336425 );
                wgts.set(3,  0.0475458171336425 );
                wgts.set(4,  0.0475458171336425 );
                wgts.set(5,  0.0516086852673590 );
                wgts.set(6,  0.0516086852673590 );
                wgts.set(7,  0.0516086852673590 );
                wgts.set(8,  0.0162292488115990 );
                wgts.set(9,  0.0162292488115990 );
                wgts.set(10, 0.0162292488115990 );
                wgts.set(11, 0.0136151570872175 );
                wgts.set(12, 0.0136151570872175 );
                wgts.set(13, 0.0136151570872175 );
                wgts.set(14, 0.0136151570872175 );
                wgts.set(15, 0.0136151570872175 );
                wgts.set(16, 0.0136151570872175 );

    for(int l1=1;l1<=2;l1++)
    for (int m1=1; m1<=mpoints; m1++)
    {   int j1=0;
        int m;
        double xi,eta;
        if(l1==1){j1=0;m=m1+j1;
        double s1=spts.get(m1,1);
        double t1=spts.get(m1,2);
        xi = 2.0*s1-1.0/3.0;
        eta = 2.0*t1-1.0/3.0;
       }
       if(l1==2){j1=mpoints;m=m1+j1;
        double s1=spts.get(m1,1);
        double t1=spts.get(m1,2);
        xi = -2.0*t1+1.0/3.0;
        eta = -2.0*s1+1.0/3.0;
        }
        // coordinates
        spts1.set(m,1,xi);
        spts1.set(m,2,eta);
        wgts1.set(m,4.0*wgts.get(m1));    




        }
            double xmax=max(max(xp1,xp2),max(xp3,xp4));double xmin=min(min(xp1,xp2),min(xp3,xp4));
            double ymax=max(max(yp1,yp2),max(yp3,yp4));double ymin=min(min(yp1,yp2),min(yp3,yp4));
            double dx1=xmax-xmin;double dy1=ymax-ymin;
            dTensor2    xpts(2*mpoints,2);
            dTensor2   qvals(2*mpoints,meqn);
            vector<double> jacobian;
            double area=0.0;  
         for (int ini1=1;ini1<=2;ini1++)
         {	  



            // These need to be defined locally
            double x1,y1,x2,y2,x3,y3;

           
            vector<double> x;
            if(ini1==1){
    
            x=returnleft(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,1);
            jacobian.push_back(x.at(0));
            areal=x.at(0);
            arear=x.at(1);
            area=area+abs(areal);
            
        //    cout<<"diff1 "<<xi1-x1<<" "<<eta1-y1<<" "<<xi2-x2<<" "<<eta2-y2<<" "<<xi4-x3<<" "<<eta4-y3<<endl;
            }
            if(ini1==2){
 
            x=returnright(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,1);

            jacobian.push_back(x.at(1));
            areal=x.at(0);
            arear=x.at(1);
            area=area+abs(arear);

    }

        
            x1=x.at(2);y1=x.at(3);
            x2=x.at(4);y2=x.at(5);
            x3=x.at(6);y3=x.at(7);

//cout<<i<<" "<<jm<<" kmax_fout "<<kmax_fout<<" "<<mpoints<<endl;
            double xc = (x1+x2+x3)/3.0;
            double yc = (y1+y2+y3)/3.0;
            double jmat[2][2];
            if(ini1==1){
            jmat[0][0]=copysign(1.0,jacobian[ini1-1])*(y3-y1);jmat[0][1]=copysign(1.0,jacobian[ini1-1])*(y1-y2);
            jmat[1][0]=copysign(1.0,jacobian[ini1-1])*(x1-x3);jmat[1][1]=copysign(1.0,jacobian[ini1-1])*(x2-x1);}
            if(ini1==2){
            jmat[0][0]=copysign(1.0,jacobian[ini1-1])*(y2-y1);jmat[0][1]=copysign(1.0,jacobian[ini1-1])*(y1-y3);
            jmat[1][0]=copysign(1.0,jacobian[ini1-1])*(x1-x2);jmat[1][1]=copysign(1.0,jacobian[ini1-1])*(x3-x1);}
           // cout<<i<<" "<<jm<<" "<<ini1<<" jmat "<<jmat[0][0]<<" "<<jmat[0][1]<<" "<<jmat[1][0]<<" "<<jmat[1][1]<<endl;
            jacobian[ini1-1]=abs(jacobian[ini1-1]);
            // Compute q, aux and fvals at each Gaussian Quadrature point
            // for this current cell indexed by (i,j)
            // Save results into dTensor2 qvals, auxvals and fvals.
            for (int m1=1; m1<=mpoints; m1++)
            {    int j1=0;if(ini1==2){j1=mpoints;}
                 int m=m1+j1;

                // point on the unit triangle
                const double s = spts1.get(m,1);
                const double t = spts1.get(m,2);

            if(ini1==1){
            // point on the physical triangle
            xpts.set(m,1, 0.5*(x2+x3) + 0.5*(x2-x1)*s + 0.5*(x3-x1)*t );
            xpts.set(m,2, 0.5*(y2+y3) + 0.5*(y2-y1)*s + 0.5*(y3-y1)*t );}
            else{
            xpts.set(m,1, 0.5*(x2+x3) + 0.5*(x1-x3)*s + 0.5*(x1-x2)*t );
            xpts.set(m,2, 0.5*(y2+y3) + 0.5*(y1-y3)*s + 0.5*(y1-y2)*t );}
            double xo=xpts.get(m,1);double yo=xpts.get(m,2);

        

                // Solution values (q) at each grid point

                    qvals.set(m,1, xpts.get(m,1) );

//cout<<setprecision(15)<<"pts "<<2.0*(xpts.get(m,1)-xmin)/dx1-1.0<<" "<<2.0*(xpts.get(m,2)-ymin)/dy1-1.0<<" "<<qvals.get(m,1)<<endl;
//cout<<setprecision(15)<<" "<<(xp1-xmin)/dx1*2.0-1.0<<" "<<(yp1-ymin)/dy1*2.0-1.0<<" "<<(xp2-xmin)/dx1*2.0-1.0<<" "<<(yp2-ymin)/dy1*2.0-1.0<<" "<<(xp3-xmin)/dx1*2.0-1.0<<" "<<(yp3-ymin)/dy1*2.0-1.0<<" "<<(xp4-xmin)/dx1*2.0-1.0<<" "<<(yp4-ymin)/dy1*2.0-1.0<<endl;
            }
          }

            // Call user-supplied function to set fvals
          //  cout<<"stuff= "<<fvals.get(1,1,1)<<" "<<fvals.get(1,1,2)<<" "<<auxvals.get(1,1)<<" "<<qvals.get(1,1)<<endl;
            // Evaluate integral on current cell (project onto Legendre basis) 
            // using Gaussian Quadrature for the integration
double tmp=0.0;
               for (int ini1=1;ini1<=2;ini1++)
               {	  
                    
                    //cout<<"tmp= "<<tmp<<endl;
                    for (int k1=1; k1<=mpoints; k1++)
                    { int k=k1;if(ini1==2){k=k1+mpoints;}
                           // if(i==2 && jm==2){cout<<"f= "<<fvals.get(k,1,1)<<" "<<fvals.get(k,1,2)<<endl;}
                        tmp = tmp + wgts1.get(k)*abs(jacobian[ini1-1])/2.0*qvals.get(k,1);
                     //   cout<<"update "<<qvals.get(k,1)<<" "<<tmp<<endl;
                   // cout<<k<<" "<<wgts1.get(k)<<" "<<phi_x.get(k,m2)<<" "<<fvals.get(k,m1,1)<<endl;
                    }
                    //cout<<i<<" "<<jm<<" l1= "<<ini1<<" m1= "<<m1<<" m2= "<<m2<<" "<<fout2.get(i,jm,m1,m2)<<endl;
                }
 //cout<<"areas "<<jacobian[0]<<" "<<jacobian[1]<<" "<<tmp<<endl;




return tmp/area;
}

double integrateyc(double oxp1,double oyp1,double oxp2,double oyp2,double oxp3,double oyp3,double oxp4,double oyp4,double q[6])
{
int mpoints = 16;int meqn=1;int maux=1;
int numcoeffs=6;int kmax_fout=6;int kmax_qin=6;int kmax=6;
dTensor2    spts(mpoints,2);
dTensor2    spts1(2*mpoints,2);
dTensor1    wgts(mpoints);
dTensor1    wgts1(2*mpoints);int mcomps_out=1;
        dTensor2     phi(2*mpoints,kmax); // Legendre basis (orthogonal)
            double xp1 = oxp1;
            double yp1= oyp1;
            double xp2 = oxp2;
            double yp2= oyp2;
            double xp3 = oxp3;
            double yp3= oyp3;
            double xp4 = oxp4;
            double yp4= oyp4;

            double xr1 = xp1;
            double yr1= yp1;
            double xr2 = xp2;
            double yr2= yp2;
            double xr3 = xp3;
            double yr3= yp3;
            double xr4 = xp4;
            double yr4= yp4;
            double areal,arear;
            mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);

                spts.set(1,1,   0.000000000000000 );
                spts.set(1,2,   0.000000000000000 );

                spts.set(2,1,   0.125959254959390 );
                spts.set(2,2,   0.125959254959390 );

                spts.set(3,1,  -0.251918509918779 );
                spts.set(3,2,   0.125959254959390 );

                spts.set(4,1,   0.125959254959390 );
                spts.set(4,2,  -0.251918509918779 );

                spts.set(5,1,  -0.162764025581573 );
                spts.set(5,2,  -0.162764025581573 );

                spts.set(6,1,   0.325528051163147 );
                spts.set(6,2,  -0.162764025581573 );

                spts.set(7,1,  -0.162764025581573 );
                spts.set(7,2,   0.325528051163147 );

                spts.set(8,1,  -0.282786105016302 );
                spts.set(8,2,  -0.282786105016302 );

                spts.set(9,1,   0.565572210032605 );
                spts.set(9,2,  -0.282786105016302 );

                spts.set(10,1, -0.282786105016302 );
                spts.set(10,2,  0.565572210032605 );

                spts.set(11,1, -0.324938555923375 );
                spts.set(11,2, -0.070220503698695 );

                spts.set(12,1, -0.324938555923375 );
                spts.set(12,2,  0.395159059622071 );

                spts.set(13,1, -0.070220503698695 );
                spts.set(13,2, -0.324938555923375 );

                spts.set(14,1, -0.070220503698695 );
                spts.set(14,2,  0.395159059622071 );

                spts.set(15,1,  0.395159059622071 );
                spts.set(15,2, -0.324938555923375 );

                spts.set(16,1,  0.395159059622071 );
                spts.set(16,2, -0.070220503698695 );

                wgts.set(1,  0.0721578038388935 );
                wgts.set(2,  0.0475458171336425 );
                wgts.set(3,  0.0475458171336425 );
                wgts.set(4,  0.0475458171336425 );
                wgts.set(5,  0.0516086852673590 );
                wgts.set(6,  0.0516086852673590 );
                wgts.set(7,  0.0516086852673590 );
                wgts.set(8,  0.0162292488115990 );
                wgts.set(9,  0.0162292488115990 );
                wgts.set(10, 0.0162292488115990 );
                wgts.set(11, 0.0136151570872175 );
                wgts.set(12, 0.0136151570872175 );
                wgts.set(13, 0.0136151570872175 );
                wgts.set(14, 0.0136151570872175 );
                wgts.set(15, 0.0136151570872175 );
                wgts.set(16, 0.0136151570872175 );

    for(int l1=1;l1<=2;l1++)
    for (int m1=1; m1<=mpoints; m1++)
    {   int j1=0;
        int m;
        double xi,eta;
        if(l1==1){j1=0;m=m1+j1;
        double s1=spts.get(m1,1);
        double t1=spts.get(m1,2);
        xi = 2.0*s1-1.0/3.0;
        eta = 2.0*t1-1.0/3.0;
       }
       if(l1==2){j1=mpoints;m=m1+j1;
        double s1=spts.get(m1,1);
        double t1=spts.get(m1,2);
        xi = -2.0*t1+1.0/3.0;
        eta = -2.0*s1+1.0/3.0;
        }
        // coordinates
        spts1.set(m,1,xi);
        spts1.set(m,2,eta);
        wgts1.set(m,4.0*wgts.get(m1));    




        }
            double xmax=max(max(xp1,xp2),max(xp3,xp4));double xmin=min(min(xp1,xp2),min(xp3,xp4));
            double ymax=max(max(yp1,yp2),max(yp3,yp4));double ymin=min(min(yp1,yp2),min(yp3,yp4));
            double dx1=xmax-xmin;double dy1=ymax-ymin;
            dTensor2    xpts(2*mpoints,2);
            dTensor2   qvals(2*mpoints,meqn);
            vector<double> jacobian;
            double area=0.0;  
         for (int ini1=1;ini1<=2;ini1++)
         {	  



            // These need to be defined locally
            double x1,y1,x2,y2,x3,y3;

           
            vector<double> x;
            if(ini1==1){
    
            x=returnleft(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,1);
            jacobian.push_back(x.at(0));
            areal=x.at(0);
            arear=x.at(1);
            area=area+abs(areal);
            
        //    cout<<"diff1 "<<xi1-x1<<" "<<eta1-y1<<" "<<xi2-x2<<" "<<eta2-y2<<" "<<xi4-x3<<" "<<eta4-y3<<endl;
            }
            if(ini1==2){
 
            x=returnright(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,1);

            jacobian.push_back(x.at(1));
            areal=x.at(0);
            arear=x.at(1);
            area=area+abs(arear);

    }

        
            x1=x.at(2);y1=x.at(3);
            x2=x.at(4);y2=x.at(5);
            x3=x.at(6);y3=x.at(7);

//cout<<i<<" "<<jm<<" kmax_fout "<<kmax_fout<<" "<<mpoints<<endl;
            double xc = (x1+x2+x3)/3.0;
            double yc = (y1+y2+y3)/3.0;
            double jmat[2][2];
            if(ini1==1){
            jmat[0][0]=copysign(1.0,jacobian[ini1-1])*(y3-y1);jmat[0][1]=copysign(1.0,jacobian[ini1-1])*(y1-y2);
            jmat[1][0]=copysign(1.0,jacobian[ini1-1])*(x1-x3);jmat[1][1]=copysign(1.0,jacobian[ini1-1])*(x2-x1);}
            if(ini1==2){
            jmat[0][0]=copysign(1.0,jacobian[ini1-1])*(y2-y1);jmat[0][1]=copysign(1.0,jacobian[ini1-1])*(y1-y3);
            jmat[1][0]=copysign(1.0,jacobian[ini1-1])*(x1-x2);jmat[1][1]=copysign(1.0,jacobian[ini1-1])*(x3-x1);}
           // cout<<i<<" "<<jm<<" "<<ini1<<" jmat "<<jmat[0][0]<<" "<<jmat[0][1]<<" "<<jmat[1][0]<<" "<<jmat[1][1]<<endl;
            jacobian[ini1-1]=abs(jacobian[ini1-1]);
            // Compute q, aux and fvals at each Gaussian Quadrature point
            // for this current cell indexed by (i,j)
            // Save results into dTensor2 qvals, auxvals and fvals.
            for (int m1=1; m1<=mpoints; m1++)
            {    int j1=0;if(ini1==2){j1=mpoints;}
                 int m=m1+j1;

                // point on the unit triangle
                const double s = spts1.get(m,1);
                const double t = spts1.get(m,2);

            if(ini1==1){
            // point on the physical triangle
            xpts.set(m,1, 0.5*(x2+x3) + 0.5*(x2-x1)*s + 0.5*(x3-x1)*t );
            xpts.set(m,2, 0.5*(y2+y3) + 0.5*(y2-y1)*s + 0.5*(y3-y1)*t );}
            else{
            xpts.set(m,1, 0.5*(x2+x3) + 0.5*(x1-x3)*s + 0.5*(x1-x2)*t );
            xpts.set(m,2, 0.5*(y2+y3) + 0.5*(y1-y3)*s + 0.5*(y1-y2)*t );}
            double xo=xpts.get(m,1);double yo=xpts.get(m,2);

        

                // Solution values (q) at each grid point

                    qvals.set(m,1, xpts.get(m,2) );

//cout<<setprecision(15)<<"pts "<<2.0*(xpts.get(m,1)-xmin)/dx1-1.0<<" "<<2.0*(xpts.get(m,2)-ymin)/dy1-1.0<<" "<<qvals.get(m,1)<<endl;
//cout<<setprecision(15)<<" "<<(xp1-xmin)/dx1*2.0-1.0<<" "<<(yp1-ymin)/dy1*2.0-1.0<<" "<<(xp2-xmin)/dx1*2.0-1.0<<" "<<(yp2-ymin)/dy1*2.0-1.0<<" "<<(xp3-xmin)/dx1*2.0-1.0<<" "<<(yp3-ymin)/dy1*2.0-1.0<<" "<<(xp4-xmin)/dx1*2.0-1.0<<" "<<(yp4-ymin)/dy1*2.0-1.0<<endl;
            }
          }

            // Call user-supplied function to set fvals
          //  cout<<"stuff= "<<fvals.get(1,1,1)<<" "<<fvals.get(1,1,2)<<" "<<auxvals.get(1,1)<<" "<<qvals.get(1,1)<<endl;
            // Evaluate integral on current cell (project onto Legendre basis) 
            // using Gaussian Quadrature for the integration
double tmp=0.0;
               for (int ini1=1;ini1<=2;ini1++)
               {	  
                    
                    //cout<<"tmp= "<<tmp<<endl;
                    for (int k1=1; k1<=mpoints; k1++)
                    { int k=k1;if(ini1==2){k=k1+mpoints;}
                           // if(i==2 && jm==2){cout<<"f= "<<fvals.get(k,1,1)<<" "<<fvals.get(k,1,2)<<endl;}
                        tmp = tmp + wgts1.get(k)*abs(jacobian[ini1-1])/2.0*qvals.get(k,1);
                     //   cout<<"update "<<qvals.get(k,1)<<" "<<tmp<<endl;
                   // cout<<k<<" "<<wgts1.get(k)<<" "<<phi_x.get(k,m2)<<" "<<fvals.get(k,m1,1)<<endl;
                    }
                    //cout<<i<<" "<<jm<<" l1= "<<ini1<<" m1= "<<m1<<" m2= "<<m2<<" "<<fout2.get(i,jm,m1,m2)<<endl;
                }
 //cout<<"areas "<<jacobian[0]<<" "<<jacobian[1]<<" "<<tmp<<endl;




return tmp/area;
}
