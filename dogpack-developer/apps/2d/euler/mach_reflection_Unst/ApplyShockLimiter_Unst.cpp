#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"                   // Required for all unstructured stuff
#include "Quadrature.h"             // 1D Quadrature rules
#include "MonomialsToLegendre.h"    // For evaluating basis functions
#include "edge_data_Unst.h"
#include <iostream>
using namespace std;

// Routine that applies a modified Barth-Jespersen limiter to higher order moments of
// the conserved variables.  This limiter is locally applied, and conserves
// total mass by not adjusting cell averages.
//
// See: K. Michalak and C. Ollivier-Gooch, "Limiters for Unstructured Higher-Order
// Accurate Solutions of the Euler Equations"

/*
An example of one of the phi functions we have experimented with. This
phi function penalizes non-monotonic distributions. It is used in place
of what would typically be max(1,x). When we have shocks in several 
neighboring cells we need to use a phi function such that phi(1)<1 and
phi(x)=1 for x>alpha for 2 > alpha >1. The reason for this is
that if neighboring cells have shocks the bounds we would normally
use will be polluted and it is possible that x=1 even while the current
cell is experiencing nonphysical oscillations. However if we use one of these
phi functions over many steps we will penalize cells such that x=1. Eventually
this will eliminate the oscillations because a smooth 
monotonic function should actually satisfy x \ge 2.
*/


inline double phi_func(double x)
{
    //return Min(1.0,x);//(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
   return (x<1.1)*Min(1.0/1.1*x,1.0)+(x>=1.1)*1.0;//(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}


void ApplyShockLimiter_Unst(const mesh& Mesh,const edge_data_Unst& EdgeData, const dTensor3& aux, dTensor3& q)
{

    const int NumElems      = Mesh.get_NumElems();
    const int NumPhysElems  = Mesh.get_NumPhysElems();
    const int NumEdges      = Mesh.get_NumEdges();
    const int NumNodes      = Mesh.get_NumNodes();
    const int meqn          = q.getsize(2);
    const int kmax          = q.getsize(3);
    const int maux          = aux.getsize(2);
    const int space_order   = dogParams.get_space_order();

    // Do nothing in the case of piecewise constants
    if( space_order == 1 )
    { return; }

    // ------------------------------------------------ //
    // number of points where we want to check solution //
    // ------------------------------------------------ //
    const int space_order_sq = space_order*space_order;
    const int mpts_vec[] = {0, 3*space_order_sq, 18, 3*space_order_sq, 3*space_order_sq };  // TODO - FILL IN 2ND-ORDER CASE
    const int mpoints    = mpts_vec[space_order-1];

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //
    dTensor2 spts(mpoints, 2);
    void SetPositivePointsBarJesp_Unst(const int& space_order, dTensor2& spts);
    SetPositivePointsBarJesp_Unst(space_order, spts);

    void SamplePhiAtPositivePointsBarJesp_Unst(const int& space_order, 
            const dTensor2& spts, dTensor2& phi);
    dTensor2 phi(mpoints, kmax);
    SamplePhiAtPositivePointsBarJesp_Unst(space_order, spts, phi);

    dTensor2 edgeval(NumEdges, meqn);edgeval.setall(0.0);
    double eps=1.0e-10;

    dTensor1 MTheta(NumElems); MTheta.setall(1.0);
#pragma omp parallel for
    // Loop over all interior edges
    for (int i=1; i<=NumEdges; i++)
    {
        // Edge coordinates
        double x1 = Mesh.get_edge(i,1);
        double y1 = Mesh.get_edge(i,2);
        double x2 = Mesh.get_edge(i,3);
        double y2 = Mesh.get_edge(i,4);

        // Elements on either side of edge
        int ileft  = Mesh.get_eelem(i,1);
        int iright = Mesh.get_eelem(i,2);  
        double Areal = Mesh.get_area_prim(ileft);
        double Arear = Mesh.get_area_prim(iright);

        // Scaled normal to edge
        dTensor1 nhat(2);      
        nhat.set(1, (y2-y1) );
        nhat.set(2, (x1-x2) );
        double length=sqrt((x1-x2)*(x1-x2)+(y2-y1)*(y2-y1));

        // Variables to store flux integrals along edge
        dTensor2 Fr_tmp(meqn,dogParams.get_space_order());
        dTensor2 Fl_tmp(meqn,dogParams.get_space_order());

        // Loop over number of quadrature points along each edge
        for (int ell=1; ell<=dogParams.get_space_order(); ell++)
        {
            dTensor1 Ql(meqn),Qr(meqn);
            dTensor1 Auxl(maux),Auxr(maux);	  

            // Riemann data - q
            for (int m=1; m<=meqn; m++)
            {
                Ql.set(m, 0.0 );
                Qr.set(m, 0.0 );

                for (int k=1; k<=kmax; k++)
                {
                    Ql.set(m, Ql.get(m) + EdgeData.phi_left->get(i,ell,k) 
                            *q.get(ileft, m,k) );
                    Qr.set(m, Qr.get(m) + EdgeData.phi_right->get(i,ell,k)
                            *q.get(iright,m,k) );
                }
               edgeval.set(i,m,edgeval.get(i,m)+(Qr.get(m)-Ql.get(m))*EdgeData.wgts1d->get(ell));
               //EdgeData.wgts1d->get(ell)
            }
        }
     }

    //--------------------------------------------------------------//
    // We will store the max and min values on each cell....//
    //--------------------------------------------------------------//
    dTensor2 MaxVal(NumElems,meqn);MaxVal.setall(-1.0e6);
    dTensor2 MinVal(NumElems,meqn);MinVal.setall(1.0e6);

    dTensor2 MaxNode(NumNodes,meqn);MaxNode.setall(-1.0e6);
    dTensor2 MinNode(NumNodes,meqn);MinNode.setall(1.0e6);

//#pragma omp parallel for
    for(int  i=1;  i <= NumElems; i++)
    {

         double mine     = MinVal.get(i,5);
         double maxe     = MaxVal.get(i,5);
         double minrho   = MinVal.get(i,1);
         double maxrho   = MaxVal.get(i,1);
         double minurho   = MinVal.get(i,2);
         double minvrho   = MinVal.get(i,3);
         double minwrho   = MinVal.get(i,4);
         double maxurho   = MaxVal.get(i,2);
         double maxvrho   = MaxVal.get(i,3);
         double maxwrho   = MaxVal.get(i,4);
         double thetae   = 1.0;
         double thetarho = 1.0;

         double m = eps;
         double ave=q.get(i,1,1)*q.get(i,5,1)-0.5*(q.get(i,2,1)*q.get(i,2,1)+q.get(i,3,1)*q.get(i,3,1)+q.get(i,4,1)*q.get(i,4,1));
         double thetam = 1.0;


         double mxa  = q.get(i,2,1);
         double mya  = q.get(i,3,1);
         double mza  = q.get(i,4,1);
         double ea   = q.get(i,5,1);
         double rhoa = q.get(i,1,1);

         int s=1;


        double minval = 1.0e16;
        double maxval = -1.0e16;
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double mxnow  = 0.0;
            double mynow  = 0.0;
            double mznow  = 0.0;
            double enow   = 0.0;
            double rhonow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
               mxnow  += q.get(i,2,k) * phi.get(mp,k);
               mynow  += q.get(i,3,k) * phi.get(mp,k);
               mznow  += q.get(i,4,k) * phi.get(mp,k);
               enow   += q.get(i,5,k) * phi.get(mp,k);
               rhonow += q.get(i,1,k) * phi.get(mp,k);
            }
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

            double mxc  = mxnow-mxa;
            double myc  = mynow-mya;
            double mzc  = mznow-mza;
            double ec   = enow-ea;
            double rhoc = rhonow-rhoa;
            double slope = ec*rhoa+rhoc*ea-(mxc*mxa+myc*mya+mzc*mza);
            double press = rhonow*enow-0.5*(mxnow*mxnow+mynow*mynow+mznow*mznow);
            if( fabs( slope ) < eps ){ thetam = thetam; }
            else{ thetam = Min( thetam, fabs( (eps-ave) / slope ) ); }

           
        }

        MaxNode.set(Mesh.get_tnode(i,1),1,Max(maxrho,MaxNode.get(Mesh.get_tnode(i,1),1)));
        MaxNode.set(Mesh.get_tnode(i,2),1,Max(maxrho,MaxNode.get(Mesh.get_tnode(i,2),1)));
        MaxNode.set(Mesh.get_tnode(i,3),1,Max(maxrho,MaxNode.get(Mesh.get_tnode(i,3),1)));
        MinNode.set(Mesh.get_tnode(i,1),1,Min(minrho,MinNode.get(Mesh.get_tnode(i,1),1)));
        MinNode.set(Mesh.get_tnode(i,2),1,Min(minrho,MinNode.get(Mesh.get_tnode(i,2),1)));
        MinNode.set(Mesh.get_tnode(i,3),1,Min(minrho,MinNode.get(Mesh.get_tnode(i,3),1)));

        MaxNode.set(Mesh.get_tnode(i,1),2,Max(maxurho,MaxNode.get(Mesh.get_tnode(i,1),2)));
        MaxNode.set(Mesh.get_tnode(i,2),2,Max(maxurho,MaxNode.get(Mesh.get_tnode(i,2),2)));
        MaxNode.set(Mesh.get_tnode(i,3),2,Max(maxurho,MaxNode.get(Mesh.get_tnode(i,3),2)));
        MinNode.set(Mesh.get_tnode(i,1),2,Min(minurho,MinNode.get(Mesh.get_tnode(i,1),2)));
        MinNode.set(Mesh.get_tnode(i,2),2,Min(minurho,MinNode.get(Mesh.get_tnode(i,2),2)));
        MinNode.set(Mesh.get_tnode(i,3),2,Min(minurho,MinNode.get(Mesh.get_tnode(i,3),2)));

        MaxNode.set(Mesh.get_tnode(i,1),3,Max(maxvrho,MaxNode.get(Mesh.get_tnode(i,1),3)));
        MaxNode.set(Mesh.get_tnode(i,2),3,Max(maxvrho,MaxNode.get(Mesh.get_tnode(i,2),3)));
        MaxNode.set(Mesh.get_tnode(i,3),3,Max(maxvrho,MaxNode.get(Mesh.get_tnode(i,3),3)));
        MinNode.set(Mesh.get_tnode(i,1),3,Min(minvrho,MinNode.get(Mesh.get_tnode(i,1),3)));
        MinNode.set(Mesh.get_tnode(i,2),3,Min(minvrho,MinNode.get(Mesh.get_tnode(i,2),3)));
        MinNode.set(Mesh.get_tnode(i,3),3,Min(minvrho,MinNode.get(Mesh.get_tnode(i,3),3)));

        MaxNode.set(Mesh.get_tnode(i,1),4,Max(maxwrho,MaxNode.get(Mesh.get_tnode(i,1),4)));
        MaxNode.set(Mesh.get_tnode(i,2),4,Max(maxwrho,MaxNode.get(Mesh.get_tnode(i,2),4)));
        MaxNode.set(Mesh.get_tnode(i,3),4,Max(maxwrho,MaxNode.get(Mesh.get_tnode(i,3),4)));
        MinNode.set(Mesh.get_tnode(i,1),4,Min(minwrho,MinNode.get(Mesh.get_tnode(i,1),4)));
        MinNode.set(Mesh.get_tnode(i,2),4,Min(minwrho,MinNode.get(Mesh.get_tnode(i,2),4)));
        MinNode.set(Mesh.get_tnode(i,3),4,Min(minwrho,MinNode.get(Mesh.get_tnode(i,3),4)));

        MaxNode.set(Mesh.get_tnode(i,1),5,Max(maxe,MaxNode.get(Mesh.get_tnode(i,1),5)));
        MaxNode.set(Mesh.get_tnode(i,2),5,Max(maxe,MaxNode.get(Mesh.get_tnode(i,2),5)));
        MaxNode.set(Mesh.get_tnode(i,3),5,Max(maxe,MaxNode.get(Mesh.get_tnode(i,3),5)));
        MinNode.set(Mesh.get_tnode(i,1),5,Min(mine,MinNode.get(Mesh.get_tnode(i,1),5)));
        MinNode.set(Mesh.get_tnode(i,2),5,Min(mine,MinNode.get(Mesh.get_tnode(i,2),5)));
        MinNode.set(Mesh.get_tnode(i,3),5,Min(mine,MinNode.get(Mesh.get_tnode(i,3),5)));


       /* MaxNode.set(Mesh.get_tnode(i,1),me,Max(q.get(i,me,1),MaxNode.get(Mesh.get_tnode(i,1),me)));
        MaxNode.set(Mesh.get_tnode(i,2),me,Max(q.get(i,me,1),MaxNode.get(Mesh.get_tnode(i,2),me)));
        MaxNode.set(Mesh.get_tnode(i,3),me,Max(q.get(i,me,1),MaxNode.get(Mesh.get_tnode(i,3),me)));

        MinNode.set(Mesh.get_tnode(i,1),me,Min(q.get(i,me,1),MinNode.get(Mesh.get_tnode(i,1),me)));
        MinNode.set(Mesh.get_tnode(i,2),me,Min(q.get(i,me,1),MinNode.get(Mesh.get_tnode(i,2),me)));
        MinNode.set(Mesh.get_tnode(i,3),me,Min(q.get(i,me,1),MinNode.get(Mesh.get_tnode(i,3),me)));*/


        MaxVal.set(i,1,maxrho);
        MaxVal.set(i,2,maxurho);
        MaxVal.set(i,3,maxvrho);
        MaxVal.set(i,4,maxwrho);
        MaxVal.set(i,5,maxe);
        MinVal.set(i,1,minrho);
        MinVal.set(i,2,minurho);
        MinVal.set(i,3,minvrho);
        MinVal.set(i,4,minwrho);
        MinVal.set(i,5,mine);

        double theta = 1.0;
        double Q1 = ave;
        if( fabs( Q1 - m ) < eps ){ theta = 0.0; }
        else if( m < eps){ theta = Min( thetam, fabs( (eps-Q1) / (m-Q1) ) ); }

        Q1 = q.get(i,1,1);
        if (Q1 < 0.0 || isnan(Q1))
        {    cout << i<<" Negative dens "<<Q1 << endl;}
         if( fabs( Q1 - minrho ) < eps || Q1 < eps ){ thetarho = 0.0; }
         else if(  minrho <eps) { thetarho = Min( thetarho, fabs( (eps-Q1) / (minrho-Q1) ) ); }

         theta=Min(Min(theta,thetam),Min(theta,thetarho));        
         MTheta.set(i,theta);   
    }

#pragma omp parallel for
    for(int  i=1;  i <= NumPhysElems; i++)
    {
     double thetae=MTheta.get(i);
     int detector=1;
     double Area = Mesh.get_area_prim(i);
     double dx=sqrt(2.0*Area);
     double alpha=1.0*pow(0.1*dx,1.5);
   
     double cutoff=pow(dx,2.5);
    {
    for(int me=1; me <= meqn; me++)
    {
       double Q1=q.get(i,me,1);
       double diffM=-1.0e6;
       double diffm=1.0e6;

       //Find the deviation from the max and min values
       //on neighbouring cells from the average value
       //on our current cell.
       /*for(int e=1;e<=3;e++)
       {
          int edge=Mesh.get_tedge(i,e);
          int index;
          int ileft  = Mesh.get_eelem(edge,1);
          int iright = Mesh.get_eelem(edge,2);
          if(ileft==i){index=iright;}
          else{index=ileft;}
          double diffM=Max(diffM,MaxVal.get(index,me)-Q1);
          double diffm=Min(diffm,MinVal.get(index,me)-Q1);
          if(isnan(diffM1) || isnan(diffm1)){cout<<NumPhysElems<<" problem here "<<index<<endl;exit(1);}
          diffM=diffM1;
       }*/

        for(int e=1;e<=3;e++)
       {
          int nodei=Mesh.get_tnode(i,e);
          diffM=Max(diffM,MaxNode.get(nodei,me)-Q1);
          diffm=Min(diffm,MinNode.get(nodei,me)-Q1);
          if(isnan(diffM) || isnan(diffm)){cout<<NumPhysElems<<" problem here1 "<<endl;}
       }
       diffM=Max(diffM,alpha);
       diffm=Min(diffm,-alpha); 
       double diffcM=MaxVal.get(i,me)-Q1;
       double diffcm=MinVal.get(i,me)-Q1;
        //Compute the minimum theta value we need to bound
        //diffcM and diffcm between all of the neighbouring cell
        //differencesl

       double thetam1,thetaM1;
       if (fabs(diffcm)<1.0e-14)
       {
           thetam1=1.0;
       }
       else
       {
           thetam1=phi_func(diffm/diffcm);
       }


       if (fabs(diffcM)<1.0e-14)
       {
           thetaM1=1.0;
       }
       else
       {
           thetaM1=phi_func(diffM/diffcM);
       }

       double theta=Min(thetam1,thetaM1);
       //cout<<" HERE "<<diffM<<" "<<diffm<<" "<<diffcM<<" "<<diffcm<<" "<<i<<endl;
       //assert_ge(theta,0.0);
       //assert_le(theta,1.0);
       //compute the single theta value to use for all entries in a system
       thetae=Max(0.0,Min(theta,thetae));
       /*for( int k=2;k<=kmax;k++)
       {
            q.set(i,me,k,q.get(i,me,k)*theta);
       }*/
      }
      if(thetae<1.0)//phi_func(1.0))
      {//Limit all of the physical variables with thetae
        for(int me=1; me <= meqn; me++)
        { 
            for( int k=2; k <= kmax; k++ )
            {
                 q.set(i,me,k, q.get(i,me,k) * thetae );
            }  
        }
     }
    }
    }

}



