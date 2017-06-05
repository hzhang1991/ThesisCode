#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "Legendre2d.h"
#include "ApplyPosLimiter.h"
#include <iostream>
using namespace std;

// Routine that applies a modified Barth-Jespersen limiter to higher order moments of
// the conserved variables.  This limiter is locally applied, and conserves
// total mass by not adjusting cell averages.
//
// See: K. Michalak and C. Ollivier-Gooch, "Limiters for Unstructured Higher-Order
// Accurate Solutions of the Euler Equations"
inline double phi_func(double x)
{
   return Min(1.0,x/1.1);//(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}
void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
    void (*ProjectRightEig)(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&),
    void (*ProjectLeftEig)(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&))
{


    double gamma=1.4;
    double Min(double, double);
    const int   mx = q.getsize(1);
    const int   my = q.getsize(2);
    const int meqn = q.getsize(3);
    //const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(2);
    const int space_order = dogParams.get_space_order();

    // ------------------------------------------------ //
    // Number of points where we want to check solution
    //
    // Change this quantity as needed to check more points
    //
    // ------------------------------------------------ //
    const int mpoints   =  4+4*space_order+space_order*space_order;

    double eps=1.0e-10;

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //

    // Quadrature points
    dTensor2 spts(mpoints,2);
    dTensor1 w1d(space_order);
    void SetPositivePointsGauss(const int& space_order, dTensor2& spts,dTensor1& w1d);
    SetPositivePointsGauss(space_order,spts,w1d);

    // Basis functions evaluated at each point
    dTensor2 phi(mpoints,kmax);
    SetLegendreAtPoints(space_order,spts,phi);


    //A matrix to store max and min values on each cell. 
    //We need this to see if a cell's point values fall
    //within the range of the values of its neighbors. 
    //This will determine bounds on each cell.

    dTensorBC3 MaxVal(mx,my,meqn,mbc); MaxVal.setall(-1.0e6);
    dTensorBC3 MinVal(mx,my,meqn,mbc); MinVal.setall(1.0e6);
    dTensorBC3 MaxValC(mx,my,meqn,mbc); MaxVal.setall(-1.0e6);
    dTensorBC3 MinValC(mx,my,meqn,mbc); MinVal.setall(1.0e6);
    dTensorBC3 avecell (mx,my, meqn,mbc);

    dTensorBC2 MTheta(mx,my,mbc); MTheta.setall(1.0);

    double dx=1.0/mx;

    const int istep = (mx / 5)+1;
    const int jstep = (my / 5)+1;


    // -------------------------------------------------------------- //
    //                                                                //
    // Loop over each cell and every quadrature point on each cell    //
    // I do not see a clever way of avoiding this right now 
    // unless we want to compute approximate max and min...           //
    // -------------------------------------------------------------- //
#pragma omp parallel for
    for(int i=1; i <= mx; i++)
    for(int j=1; j <= my; j++)
    //for(int i=1-mbc; i <= mx+mbc; i++)
    //for(int j=1-mbc; j <= my+mbc; j++)
    {

         double minP     = MinVal.get(i,j,5);
         double maxP     = MaxVal.get(i,j,5);
         double minrho   = MinVal.get(i,j,1);
         double maxrho   = MaxVal.get(i,j,1);
         double minu   = MinVal.get(i,j,2);
         double minv   = MinVal.get(i,j,3);
         double minw   = MinVal.get(i,j,4);
         double maxu   = MaxVal.get(i,j,2);
         double maxv   = MaxVal.get(i,j,3);
         double maxw   = MaxVal.get(i,j,4);

         double maxurho=-1.0e6;
         double minurho=1.0e6;
         double maxvrho=-1.0e6;
         double minvrho=1.0e6;
         double maxwrho=-1.0e6;
         double minwrho=1.0e6;
         double maxe=-1.0e6;
         double mine=1.0e6;

         double thetae   = 1.0;
         double thetarho = 1.0;

         double m = eps;
         double ave=q.get(i,j,1,1)*q.get(i,j,5,1)-0.5*(q.get(i,j,2,1)*q.get(i,j,2,1)+q.get(i,j,3,1)*q.get(i,j,3,1)+q.get(i,j,4,1)*q.get(i,j,4,1));
         double thetam = 1.0;


         double mxa  = q.get(i,j,2,1);
         double mya  = q.get(i,j,3,1);
         double mza  = q.get(i,j,4,1);
         double ea   = q.get(i,j,5,1);
         double rhoa = q.get(i,j,1,1);

         int s=1;

         for(int mp=1; mp <= mpoints; mp++)
         {
            // evaluate e and rho at spts(mp)
            double mxnow  = 0.0;
            double mynow  = 0.0;
            double mznow  = 0.0;
            double enow   = 0.0;
            double rhonow = 0.0;

            for( int k=1; k <= kmax; k++ )
            {
               mxnow  += q.get(i,j,2,k) * phi.get(mp,k);
               mynow  += q.get(i,j,3,k) * phi.get(mp,k);
               mznow  += q.get(i,j,4,k) * phi.get(mp,k);
               enow   += q.get(i,j,5,k) * phi.get(mp,k);
               rhonow += q.get(i,j,1,k) * phi.get(mp,k);


            }

            double unow=mxnow/rhonow;
            double vnow=mynow/rhonow;
            double wnow=mznow/rhonow;
            //Pnow=(gamma-1.0e0)*(enowrhonow-0.5*(unow*unow));
            //double Pnow=(gamma-1.0e0)*(enow-0.5/rhonow*(mxnow*mxnow+mynow*mynow));
            double Pnow=(gamma-1.0e0)*(enow*rhonow-0.5*(mxnow*mxnow+mynow*mynow));
            //double Pnow=(gamma-1.0e0)/rhonow*(enow*rhonow-0.5*(mxnow*mxnow+mynow*mynow));

            minrho = Min(minrho,rhonow);
            minurho= Min(minurho,mxnow);
            minvrho= Min(minvrho,mynow);
            minwrho= Min(minwrho,mznow);
            mine   = Min(mine,enow);
            maxrho = Max(maxrho,rhonow);
            maxurho= Max(maxurho,mxnow);
            maxvrho= Max(maxvrho,mynow);
            maxwrho= Max(maxwrho,mznow);
            maxe   = Max(maxe,enow);

            maxu=Max(maxu,unow);
            minu=Min(minu,unow);
            maxv=Max(maxv,vnow);
            minv=Min(minv,vnow);
            maxw=Max(maxw,wnow);
            minw=Min(minw,wnow);
            maxP=Max(maxP,Pnow);
            minP=Min(minP,Pnow);



            double mxc  = mxnow-mxa;
            double myc  = mynow-mya;
            double mzc  = mznow-mza;
            double ec   = enow-ea;
            double rhoc = rhonow-rhoa;
            double slope = ec*rhoa+rhoc*ea-(mxc*mxa+myc*mya+mzc*mza);
            double press = rhonow*enow-0.5*(mxnow*mxnow+mynow*mynow+mznow*mznow);
            //double press = enow-0.5/rhonow*(mxnow*mxnow+mynow*mynow+mznow*mznow);
            if( fabs( slope ) < eps ){ thetam = thetam; }
            else{ thetam = Min( thetam, fabs( (eps-ave) / slope ) ); }

            m = Min(m, press);




         }


    bool nstep=(i<istep+10 || j>jstep+10 );
    bool corner=0;//(i>=istep-10 && i<istep+10 && j<=jstep+10); 
    bool step=~(nstep || corner);

    //nstep=0;step=1;
        /*
        MaxVal.set(i,j,1,maxrho);
        MaxVal.set(i,j,2,nstep*maxu+step*maxurho);
        MaxVal.set(i,j,3,nstep*maxv+step*maxvrho);
        MaxVal.set(i,j,4,nstep*maxw+step*maxwrho);
        MaxVal.set(i,j,5,nstep*maxP+step*maxe);
        MinVal.set(i,j,1,minrho);
        MinVal.set(i,j,2,nstep*minu+step*minurho);
        MinVal.set(i,j,3,nstep*minv+step*minvrho);
        MinVal.set(i,j,4,nstep*minw+step*minwrho);
        MinVal.set(i,j,5,nstep*minP+step*mine);
        avecell.set(i,j,1,q.get(i,j,1,1));
        avecell.set(i,j,2,nstep*q.get(i,j,2,1)/q.get(i,j,1,1)+step*q.get(i,j,2,1));
        avecell.set(i,j,3,nstep*q.get(i,j,3,1)/q.get(i,j,1,1)+step*q.get(i,j,3,1));
        avecell.set(i,j,4,nstep*q.get(i,j,4,1)/q.get(i,j,1,1)+step*q.get(i,j,4,1));
        avecell.set(i,j,5,nstep*(gamma-1.0e0)*(q.get(i,j,5,1)*q.get(i,j,1,1)-0.5*(q.get(i,j,2,1)*q.get(i,j,2,1)+q.get(i,j,3,1)*q.get(i,j,3,1)))+step*q.get(i,j,5,1));
        */

        
        MaxValC.set(i,j,1,maxrho);
        MaxValC.set(i,j,2,maxurho);
        MaxValC.set(i,j,3,maxvrho);
        MaxValC.set(i,j,4,maxwrho);
        MaxValC.set(i,j,5,maxe);
        MinValC.set(i,j,1,minrho);
        MinValC.set(i,j,2,minurho);
        MinValC.set(i,j,3,minvrho);
        MinValC.set(i,j,4,minwrho);
        MinValC.set(i,j,5,mine);
       /* avecell.set(i,j,1,q.get(i,j,1,1));
        avecell.set(i,j,2,q.get(i,j,2,1));
        avecell.set(i,j,3,q.get(i,j,3,1));
        avecell.set(i,j,4,q.get(i,j,4,1));
        avecell.set(i,j,5,q.get(i,j,5,1));
        */
        
        /*
        MaxVal.set(i,j,1,maxrho);
        MaxVal.set(i,j,2,maxurho);
        MaxVal.set(i,j,3,maxvrho);
        MaxVal.set(i,j,4,maxwrho);
        MaxVal.set(i,j,5,maxe);
        MinVal.set(i,j,1,minrho);
        MinVal.set(i,j,2,minurho);
        MinVal.set(i,j,3,minvrho);
        MinVal.set(i,j,4,minwrho);
        MinVal.set(i,j,5,mine);
        avecell.set(i,j,1,q.get(i,j,1,1));
        avecell.set(i,j,2,q.get(i,j,2,1));
        avecell.set(i,j,3,q.get(i,j,3,1));
        avecell.set(i,j,4,q.get(i,j,4,1));
        avecell.set(i,j,5,q.get(i,j,5,1));
        */


        
        MaxVal.set(i,j,1,maxrho);
        MaxVal.set(i,j,2,maxu);
        MaxVal.set(i,j,3,maxv);
        MaxVal.set(i,j,4,maxw);
        MaxVal.set(i,j,5,maxP);
        MinVal.set(i,j,1,minrho);
        MinVal.set(i,j,2,minu);
        MinVal.set(i,j,3,minv);
        MinVal.set(i,j,4,minw);
        MinVal.set(i,j,5,minP);
        avecell.set(i,j,1,q.get(i,j,1,1));
        avecell.set(i,j,2,q.get(i,j,2,1)/q.get(i,j,1,1));
        avecell.set(i,j,3,q.get(i,j,3,1)/q.get(i,j,1,1));
        avecell.set(i,j,4,q.get(i,j,4,1)/q.get(i,j,1,1));
        avecell.set(i,j,5,(gamma-1.0e0)*(q.get(i,j,5,1)*q.get(i,j,1,1)-0.5*(q.get(i,j,2,1)*q.get(i,j,2,1)+q.get(i,j,3,1)*q.get(i,j,3,1))));
       // avecell.set(i,j,5,(gamma-1.0e0)*(q.get(i,j,5,1)-0.5/q.get(i,j,1,1)*(q.get(i,j,2,1)*q.get(i,j,2,1)+q.get(i,j,3,1)*q.get(i,j,3,1))));
               
        double theta = 1.0;
        double Q1 = ave;
        if( fabs( Q1 - m ) < eps ){ theta = 0.0; }
        else if( m < eps){ theta = Min( thetam, fabs( (eps-Q1) / (m-Q1) ) ); }

        Q1 = q.get(i,j,1,1);
        if (Q1 < 0.0 || isnan(Q1))
        {    cout << i<<" "<<j<<" Negative dens "<<Q1 << endl;}
         if( fabs( Q1 - minrho ) < eps || Q1 < eps ){ thetarho = 0.0; }
         else if(  minrho <eps) { thetarho = Min( thetarho, fabs( (eps-Q1) / (minrho-Q1) ) ); }

         theta=Min(Min(theta,thetam),Min(theta,thetarho));        
         MTheta.set(i,j,theta);        
    }

    //double alpha=1000.0*pow(dx,1.0);
    //double alpha=100000.0;
    //double alpha=1.5*pow(dx,1.5);
    double alpha=500.0*pow(dx,1.5);
    double alphat=500.0*pow(dx,1.5);
    double alphac=500.0*pow(dx,1.5);
    double alphas=20.0*pow(dx,1.5);
    double cutoff=10.0;//pow(dx,0.5);
    //double cutoff=0.1;


#pragma omp parallel for
    for(int i=1; i <= mx; i++)
    for(int j=1; j <= my; j++)
    { double thetae=MTheta.get(i,j);
    //bool nstep=(i<istep+10 || j>jstep+10 );
    //bool nstep=(j>my/4 );
    bool nstep=1;
    //bool nstepleft=(j<=12*my/40 && i<=3*mx/20);
    //bool corner=(i>3*mx/20 && i<=2*mx/8 && j<=12*my/40); 
    //bool step=0;//1-(nstep || nstepleft || corner);//~(nstep || corner);
    //bool step=1-(nstep || nstepleft || corner);//~(nstep || corner);
    bool step=1-(nstep );
    double alpha1=nstep*alpha+step*alphas;
    nstep=1-step;
    //double alpha1=nstep*alpha+step*alphas;



    int detector=1;
    //if(detector==1)
    //for(int me=1; me <= meqn; me++)
    for(int me=1; me <= meqn; me++)
    {
        

       //Find the deviation from the max and min values
       //on neighbouring cells from the average value
       //on our current cell.


        //if(detector==1){alpha1=1.5*pow(dx,1.5);}

        //double Q1=q.get(i,j,me,1);
        //double Q1=avecell.get(i,j,me);
        double Q1=nstep*avecell.get(i,j,me)+step*q.get(i,j,me,1);
        /*
        double difflM=MaxVal.get(i-1,j,me)-Q1;
        double difflm=MinVal.get(i-1,j,me)-Q1;

        double diffrM=MaxVal.get(i+1,j,me)-Q1;
        double diffrm=MinVal.get(i+1,j,me)-Q1;

        double diffuM=MaxVal.get(i,j+1,me)-Q1;
        double diffum=MinVal.get(i,j+1,me)-Q1;

        double diffdM=MaxVal.get(i,j-1,me)-Q1;
        double diffdm=MinVal.get(i,j-1,me)-Q1; 
        */ 
        //bool step=~(i<istep+10 || j>jstep+10 );
        //bool nstep=(i<istep+10 || j>jstep+10 );
 
        
        double difflM=nstep*MaxVal.get(i-1,j,me)+step*q.get(i-1,j,me,1)-Q1;
        double difflm=nstep*MinVal.get(i-1,j,me)+step*q.get(i-1,j,me,1)-Q1;

        double diffrM=nstep*MaxVal.get(i+1,j,me)+step*q.get(i+1,j,me,1)-Q1;
        double diffrm=nstep*MinVal.get(i+1,j,me)+step*q.get(i+1,j,me,1)-Q1;

        double diffuM=nstep*MaxVal.get(i,j+1,me)+step*q.get(i,j+1,me,1)-Q1;
        double diffum=nstep*MinVal.get(i,j+1,me)+step*q.get(i,j+1,me,1)-Q1;

        double diffdM=nstep*MaxVal.get(i,j-1,me)+step*q.get(i,j-1,me,1)-Q1;
        double diffdm=nstep*MinVal.get(i,j-1,me)+step*q.get(i,j-1,me,1)-Q1;
        

        /*
        double difflM=nstep*MaxVal.get(i-1,j,me)+step*MaxValC.get(i-1,j,me)-Q1;
        double difflm=nstep*MinVal.get(i-1,j,me)+step*MinValC.get(i-1,j,me)-Q1;

        double diffrM=nstep*MaxVal.get(i+1,j,me)+step*MaxValC.get(i+1,j,me)-Q1;
        double diffrm=nstep*MinVal.get(i+1,j,me)+step*MinValC.get(i+1,j,me)-Q1;

        double diffuM=nstep*MaxVal.get(i,j+1,me)+step*MaxValC.get(i,j+1,me)-Q1;
        double diffum=nstep*MinVal.get(i,j+1,me)+step*MinValC.get(i,j+1,me)-Q1;

        double diffdM=nstep*MaxVal.get(i,j-1,me)+step*MaxValC.get(i,j-1,me)-Q1;
        double diffdm=nstep*MinVal.get(i,j-1,me)+step*MinValC.get(i,j-1,me)-Q1;
        */

         
        //printf("HERE %d %d %d \n",nstep,step,step+nstep);

        /*
        double difflM=MaxVal.get(i-1,j,me)-Q1;
        double difflm=MinVal.get(i-1,j,me)-Q1;

        double diffrM=MaxVal.get(i+1,j,me)-Q1;
        double diffrm=MinVal.get(i+1,j,me)-Q1;

        double diffuM=MaxVal.get(i,j+1,me)-Q1;
        double diffum=MinVal.get(i,j+1,me)-Q1;

        double diffdM=MaxVal.get(i,j-1,me)-Q1;
        double diffdm=MinVal.get(i,j-1,me)-Q1;
        */
    
        double diffM=Max(alpha1,Max(Max(difflM,diffrM),Max(diffuM,diffdM)));
        double diffm=Min(-alpha1,Min(Min(difflm,diffrm),Min(diffum,diffdm)));
         
        //cout<<"here "<<diffM<<endl;
        /*
        double diffdlM=MaxVal.get(i-1,j-1,me)-Q1;
        double diffdlm=MinVal.get(i-1,j-1,me)-Q1;

        double diffdrM=MaxVal.get(i+1,j-1,me)-Q1;
        double diffdrm=MinVal.get(i+1,j-1,me)-Q1;

        double diffurM=MaxVal.get(i+1,j+1,me)-Q1;
        double diffurm=MinVal.get(i+1,j+1,me)-Q1;

        double diffulM=MaxVal.get(i-1,j+1,me)-Q1;
        double diffulm=MinVal.get(i-1,j+1,me)-Q1;       
        */

        double diffdlM=nstep*MaxVal.get(i-1,j-1,me)+step*q.get(i-1,j-1,me,1)-Q1;
        double diffdlm=nstep*MinVal.get(i-1,j-1,me)+step*q.get(i-1,j-1,me,1)-Q1;

        double diffdrM=nstep*MaxVal.get(i+1,j-1,me)+step*q.get(i+1,j-1,me,1)-Q1;
        double diffdrm=nstep*MinVal.get(i+1,j-1,me)+step*q.get(i+1,j-1,me,1)-Q1;

        double diffurM=nstep*MaxVal.get(i+1,j+1,me)+step*q.get(i+1,j+1,me,1)-Q1;
        double diffurm=nstep*MinVal.get(i+1,j+1,me)+step*q.get(i+1,j+1,me,1)-Q1;

        double diffulM=nstep*MaxVal.get(i-1,j+1,me)+step*q.get(i-1,j+1,me,1)-Q1;
        double diffulm=nstep*MinVal.get(i-1,j+1,me)+step*q.get(i-1,j+1,me,1)-Q1;

        diffM=Max(diffM,Max(Max(diffdlM,diffdrM),Max(diffurM,diffulM)));
        diffm=Min(diffm,Min(Min(diffdlm,diffdrm),Min(diffurm,diffulm)));
        //if(me==1){diffm=Max(diffm,-Q1);}
        //cout<<"here2 "<<diffM<<" "<<diffdlM<<" "<<diffdrM<<" "<<diffurM<<" "<<diffulM<<endl;


        double diffcM=nstep*MaxVal.get(i,j,me)+step*MaxValC.get(i,j,me)-Q1;
        double diffcm=nstep*MinVal.get(i,j,me)+step*MinValC.get(i,j,me)-Q1;

        //if(Max(fabs(diffcM),fabs(diffcm))>alpha2){detector=1;}

        //Compute the minimum theta value we need to bound
        //diffcM and diffcm between all of the neighbouring cell
        //differencesl
        double thetam1,thetaM1;
        thetam1=1.0;
        thetaM1=1.0;  
        if (fabs(diffcm)<1.0e-15)
        {
            thetam1=thetam1;
        }
        else
        {
            thetam1=Min(thetam1,phi_func(diffm/diffcm));
        }
 

        if (fabs(diffcM)<1.0e-15)
        {
            thetaM1=thetaM1;
        }
        else
        {
            thetaM1=Min(thetaM1,phi_func(diffM/diffcM));
        }

         

        double theta=Min(thetam1,thetaM1);
        theta=Min(1.0,theta);
        theta=Max(0.0,theta);
        thetae=Min(theta,thetae);
        //if(1)//theta <1.0)
        //{cout<<"theta= "<<theta<<" "<<diffM<<" "<<diffm<<" "<<diffcM<<" "<<diffcm<<" "<<phi_func(diffm/diffcm)<<endl;}
     }
     {  for(int k=2;k<=kmax;k++)
       {
          for(int me=1;me<=meqn;me++)
          {q.set(i,j,me,k,q.get(i,j,me,k)*thetae);}
       }
     }
    }


}


// Gaussian quadrature points.  However, an application will likely want to add in the edge
// points that are used for Riemann solves.
void SetPositivePointsGauss(const int& space_order, dTensor2& spts, dTensor1& w1d)
{

    dTensor1 s1d(space_order);

    // 1D Gaussian quadrature points
    switch ( space_order )
    {
        case 1:
            w1d.set(1, 2.0 );

            s1d.set(1, 0.0e0 );
            break;

        case 2:

            w1d.set(1,  1.0 );
            w1d.set(2,  1.0 );

            s1d.set(1, -1.0/sq3 );
            s1d.set(2,  1.0/sq3 );      
            break;

        case 3:

            w1d.set(1,  5.0/9.0 );
            w1d.set(2,  8.0/9.0 );
            w1d.set(3,  5.0/9.0 );

            s1d.set(1, -sq3/sq5 );
            s1d.set(2,  0.0e0 );
            s1d.set(3,  sq3/sq5 );
            break;

        case 4:

            w1d.set(1, (18.0 - sq3*sq10)/36.0 );
            w1d.set(2, (18.0 + sq3*sq10)/36.0 );
            w1d.set(3, w1d.get(2) );
            w1d.set(4, w1d.get(1) );


            s1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            s1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );      
            break;

        case 5:

            w1d.set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            w1d.set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            w1d.set(3, 128.0/225.0 );
            w1d.set(4, w1d.get(2) );
            w1d.set(5, w1d.get(1) );

            s1d.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            s1d.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(3,  0.0 );
            s1d.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );      
            break;
    }
    // This region has been commented out because we are no longer applying
    // the limiter at the corner points

    // 2D points -- corners
  spts.set(1,1, -1.0e0 );
  spts.set(1,2, -1.0e0 );

  spts.set(2,1,  1.0e0 );
  spts.set(2,2, -1.0e0 );

  spts.set(3,1, -1.0e0 );
  spts.set(3,2,  1.0e0 );

  spts.set(4,1,  1.0e0 );
  spts.set(4,2,  1.0e0 );

    // 2D points -- left, right, bottom and top edges
  for (int m=1; m<=space_order; m++)
  {
      double s = s1d.get(m);

      // left edge
      spts.set(4+m,1, -1.0e0 );
      spts.set(4+m,2,  s     );

      // right edge
      spts.set(4+space_order+m,1,  1.0e0 );
      spts.set(4+space_order+m,2,  s     );

      // bottom edge
      spts.set(4+2*space_order+m,1,  s     );
      spts.set(4+2*space_order+m,2, -1.0e0 );

      // top edge
      spts.set(4+3*space_order+m,1,  s     );
      spts.set(4+3*space_order+m,2,  1.0e0 );
  }

    // 2D points -- all interior points
    int z = 4+4*space_order;
    //z = 0;
    for (int m=1; m<=space_order; m++)
        for (int k=1; k<=space_order; k++)
        {
            double s1 = s1d.get(m);
            double s2 = s1d.get(k);

            z = z+1;

            spts.set(z,1, s1d.get(m) );
            spts.set(z,2, s1d.get(k) );
        }

}
