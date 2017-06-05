#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "RiemannSolve.h"
#include <iostream>
using namespace std;


//

double RiemannSolver::solve(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();
    dTensor3& ftmp  = fetch_flux2 ();
    int mu;int nu;
    //if(nvec.get(1)==1){mu=1;nu=2;}
    //else{mu=2;nu=1;}
     double delta[2];
     double Qvel[2];
     double norm[2];double l=nvec.get(1)*nvec.get(1)+nvec.get(2)*nvec.get(2);
     norm[0]=nvec.get(1)/sqrt(l);
     norm[1]=nvec.get(2)/sqrt(l);
     Qvel[0]=Ql.get(2)*norm[0]+Ql.get(3)*norm[1];
     Qvel[1]=Qr.get(2)*norm[0]+Qr.get(3)*norm[1];
        delta[0] = Qr.get(1) - Ql.get(1);
        delta[1] = Qvel[1]-Qvel[0];
   
   // cout<<"Qr= "<<Qr.get(1)<<" "<<Qr.get(mu+1)<<" "<<Qr.get(nu+1)<<endl;
   // cout<<"Ql= "<<Ql.get(1)<<" "<<Ql.get(mu+1)<<" "<<Ql.get(nu+1)<<endl;
    //        # impedances:
    double   zr = Auxr.get(1)*Auxr.get(2);
    double   zl = Auxl.get(1)*Auxl.get(2);
  //  cout<<" AUX! "<<zr<<" "<<Auxr.get(1)<<" "<<Auxr.get(1)<<endl;
    double    a1 = (-delta[0] + zr*delta[1]) / (zl + zr);
    double    a2 =  (delta[0] + zl*delta[1]) / (zl + zr);
    
    //        # Compute the waves.

   // cout<<a1<<" "<<zl<<endl;
    double Qm[3];
    Qm[0]=Ql.get(1)-a1*zl;
    Qm[1]=Ql.get(2)+a1*norm[0];
    Qm[2]=Ql.get(3)+a1*norm[1];

    ftmp.set(1,1,1, Auxl.get(1)*Auxl.get(1)*Auxl.get(2)*Qm[1]);
    ftmp.set(1,2,1, Qm[0]/Auxl.get(2));
    ftmp.set(1,3,1, 0.0);

    ftmp.set(1,1,2, Auxl.get(1)*Auxl.get(1)*Auxl.get(2)*Qm[2]);
    ftmp.set(1,3,2, Qm[0]/Auxl.get(2));
    ftmp.set(1,2,2, 0.0);

    for (int m=1; m<=meqn; m++)
    {  Fl.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) ); 
       //cout<<"Fl= "<<Fl.get(m)<<endl;
    }
   // std::cout<<mu<<" "<<nu<<" "<<Fl.get(2)<<" "<<Fl.get(3)<<std::endl;
    ftmp.set(1,1,1, Auxr.get(1)*Auxr.get(1)*Auxr.get(2)*Qm[1]);
    ftmp.set(1,2,1, Qm[0]/Auxr.get(2));
    ftmp.set(1,3,1, 0.0);

    ftmp.set(1,1,2, Auxr.get(1)*Auxr.get(1)*Auxr.get(2)*Qm[2]);
    ftmp.set(1,3,2, Qm[0]/Auxr.get(2));
    ftmp.set(1,2,2, 0.0);

    for (int m=1; m<=meqn; m++)
    {  Fr.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );
       //cout<<"Fr= "<<Fr.get(m)<<endl;  
    }

    //for(int m=1;m<=meqn;m++){std::cout<<"Fr "<<m<<" "<<Fr.get(m)-Fl.get(m)<<std::endl;}

   // Fl.set(1, Auxl.get(1)*Auxl.get(1)*Auxl.get(2)*Qm[1]);
   // Fl.set(2, Qm[0]/Auxl.get(2));
   // Fl.set(3, 0.0);


   // Fr.set(1, Auxr.get(1)*Auxr.get(1)*Auxr.get(2)*Qm[1]);
   // Fr.set(2, 0.0);
   // Fr.set(3, Qm[0]/Auxr.get(2));

    // Calculate minimum and maximum speeds
    double s1,s2;
    void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
            const dTensor1& Ql, const dTensor1& Qr,
            const dTensor1& Auxl, const dTensor1& Auxr,
            double& s1,double& s2);
   // SetWaveSpd(nvec,xedge,Ql,Qr,Auxl,Auxr,s1,s2);
   // SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1,s2);
    s1=-Auxl.get(1)*sqrt(l);
    s2=Auxr.get(1)*sqrt(l);
    
    double smax_edge = Max(fabs(s1),fabs(s2));
    return smax_edge;
}


// support the old RiemannSolve interface for backwards compatibility
//
double RiemannSolve(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();
    RiemannSolver riemannSolver(meqn,maux);
    return riemannSolver.solve(nvec, xedge, Ql, Qr, Auxl, Auxr, Fl, Fr);
}
