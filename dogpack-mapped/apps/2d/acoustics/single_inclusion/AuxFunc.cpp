#include "dogdefs.h"
#include <iostream>
#include <math.h>
using namespace std;

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
//  --------------------------------------------------------
//  NOTE: Two input values "NOT_USED_1" and "NOT_USED_2"
//        are never used in this function. The reason 
//        these variables are included in the input list
//        is to make "AuxFunc.cpp" conform with the format
//        required by "L2Project.cpp".
//  --------------------------------------------------------
//


double fdisc(double x,double y)
{
double fdisc;

if(x>=0.0)
{fdisc=1.0;}
else{fdisc=-1.0;}
return fdisc;
}



void AuxFunc(const dTensor2& xpts, dTensor2& auxvals)
{
  const int numpts=xpts.getsize(1);
  dTensor1 w1d(pow(numpts,0.5));

  double cl   = 1.0;
  double rhol = 1.0;
  double cr   = 1.5;
  double rhor = 2.0;
  double bulkl = cl*cl * rhol;
  double bulkr = cr*cr * rhor;
  const int md2 = int(pow(numpts,0.5));

    switch (md2)
    {
        case 1:
            w1d.set(1, 2.0 );

            break;

        case 2:
            w1d.set(1,  1.0 );
            w1d.set(2,  1.0 );


            break;

        case 3:
            w1d.set(1,  5.0/9.0 );
            w1d.set(2,  8.0/9.0 );
            w1d.set(3,  5.0/9.0 );


            break;

        case 4:
            w1d.set(1, (18.0 - sq3*sq10)/36.0 );
            w1d.set(2, (18.0 + sq3*sq10)/36.0 );
            w1d.set(3, w1d.get(2) );
            w1d.set(4, w1d.get(1) );

   
            break;

        case 5:      
            w1d.set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            w1d.set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            w1d.set(3, 128.0/225.0 );
            w1d.set(4, w1d.get(2) );
            w1d.set(5, w1d.get(1) );

            break;

        case 6:
            w1d.set(1, 0.1713244923791703450402961 );
            w1d.set(2, 0.3607615730481386075698335 );
            w1d.set(3, 0.4679139345726910473898703 );      


            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      

            break;
        }
  double cav=0.0;
  double rhoav=0.0;
  int indic=0;

  /*for (int i=1;i<=md2;i++)
     {
      for(int j=1;j<=md2;j++)
      {int j1=i+(md2-1)*j;
       if (xpts.get(j1,1)<0.0 && indic==0)
       {indic=1;}
       if (xpts.get(j1,1)<0.0 && indic==-1)
       {indic=1;
         cout<<"BIG PROBLEM!"<<endl;}
       if (xpts.get(j1,1)>=0.0 && indic==0)
       {indic=-1;}
       if (xpts.get(j1,1)>=0.0 && indic==1)
       {indic=-1;
         cout<<"BIG PROBLEM!"<<endl;}
      }
     }


  for (int i=1;i<=md2;i++)
     {
      for(int j=1;j<=md2;j++)
      {int j1=i+(md2-1)*j;
       double f=fdisc(xpts.get(j1,1),xpts.get(j1,2));
       cav=cav+w1d.get(i)*w1d.get(j)/4.0*((f+1.0)/2.0*cr+(1.0-f)/2.0*cl);
       rhoav=rhoav+w1d.get(i)*w1d.get(j)/4.0*((f+1.0)/2.0*rhor+(1.0-f)/2.0*rhol);
      }
     }*/
      double x0circ=1.0;
    double y0circ=0.5;
    double r1circ=0.4;
     double x=xpts.get(1,1);
     double y=xpts.get(1,2);
 
  for (int i=1; i<=numpts; i++)
    {  if((x-x0circ)*(x-x0circ)+(y-y0circ)*(y-y0circ)<r1circ*r1circ)
       {auxvals.set(i,1, cr );
       auxvals.set(i,2, rhor);}
       else
       {auxvals.set(i,1, cl );
       auxvals.set(i,2, rhol);}
    }
     
}
