    #include "stdio.h"
#include <cmath>
#include<iostream>
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "constants.h"
#include "tensors.h"
#include "mesh.h"
#include "edge_data_Unst.h"

using namespace std;

//---------------------------------------
// This code creates a flux that will guarentee that cell average
// values are positive at the next time step.
// This is done by comparing these cell average values to LLF fluxes
// which are known to preserve positivity. This limiter
// is based on the limiter ApplyPosMPPLimiter.cpp which was
// taken from finess. For us these fluxes are integrated high-order
// fluxes along edges.


void ApplyPosMPPLimiter_Unst( const double dt, const mesh& Mesh,const edge_data_Unst& EdgeData, 
        const dTensor1& smax,
        const dTensor3& aux, const dTensor3& q, 
        dTensor3& FeI)
{
    const int NumElems      = Mesh.get_NumElems();
    const int NumPhysElems  = Mesh.get_NumPhysElems();
    const int NumEdges      = Mesh.get_NumEdges();
    const int meqn          = q.getsize(2);
    const int kmax          = q.getsize(3);
    const int maux          = aux.getsize(2);
    const int space_order   = dogParams.get_space_order();
       
    //define some constants (eps will be treated as zero)
    const double gm1 = 1.4;
    const double eps = 1.0e-10;
    const double rho_min = eps;

    //Get the LLF fluxes using cell averages (a finite volume flux)
    dTensor3 FeLLF(NumElems, meqn, 3 );  FeLLF.setall(0.);
    void ConstructL_LLF_Unst(const mesh& Mesh,
        const edge_data_Unst& EdgeData,
        const dTensor3& aux,     
        const dTensor3& q,       
        const dTensor1& smax, dTensor3& FeI);
    ConstructL_LLF_Unst(Mesh,EdgeData,aux,q,smax,FeLLF);
    
    double gmin[NumElems],amin[NumElems][3],thex[NumEdges],ff[3],aa[3];
    double qh[meqn],qtmp[meqn],qtmp2[meqn],qtmp3[meqn],qlf[meqn],fhat_local[3][meqn],flf_local[3][meqn];
    double rho,u1,u2,u3,energy,plow,phigh,ffsum,ftmp,phigh2;
    double ac,bc,cc,delta,root1,root2,rate;
    int isgn[3];
    //Indicator will check if any cell will have negative pressure or density if we do not
    //modify the high order fluxes
    int indicator=0;
    for (int i=1;i<=NumPhysElems;i++)
    {
       double Area = Mesh.get_area_prim(i);
       double rhoa=Area*q.get(i,1,1)-dt*(FeI.get(i,1,1)+FeI.get(i,1,2)+FeI.get(i,1,3));
       double ea=Area*q.get(i,5,1)-dt*(FeI.get(i,5,1)+FeI.get(i,5,2)+FeI.get(i,5,3));
       double momxa=Area*q.get(i,2,1)-dt*(FeI.get(i,2,1)+FeI.get(i,2,2)+FeI.get(i,2,3));
       double momya=Area*q.get(i,3,1)-dt*(FeI.get(i,3,1)+FeI.get(i,3,2)+FeI.get(i,3,3));
       double momza=Area*q.get(i,4,1)-dt*(FeI.get(i,4,1)+FeI.get(i,4,2)+FeI.get(i,4,3));
       double pressa=ea*rhoa-0.5*(momxa*momxa+momya*momya+momza*momza);
       if(rhoa<eps || pressa<eps){indicator=1;}

       gmin[i-1]=-rho_min+Area*q.get(i,1,1);
       double ffsum=0.0;

       for(int k=1;k<=3;k++)
       {
          gmin[i-1]=gmin[i-1]-dt*FeLLF.get(i,1,k);
          ff[k-1]=-dt*(FeI.get(i,1,k)-FeLLF.get(i,1,k));
          if (ff[k-1]<0.0)
          {isgn[k-1]=1;}
          else{isgn[k-1]=0;}
          ffsum = ffsum + double(isgn[k-1])*ff[k-1];
       }
       //Comp. of Lambdas. We know that if we scale the flux F^i_k by \theta^i_k
       // Then we will have that \rho^{n+1}_i >0 if 0 \le \theta^i_k \le amin^i_k
       for (int k=0; k<3; k++)
       {
            if (isgn[k]==1)
                amin[i-1][k] = Min(-gmin[i-1]/(ffsum-eps),1.e0);
            else
                amin[i-1][k] = 1.e0;
       }

       if(q.get(i,1,1)<eps)
       {cout<<"cell averages negative "<<q.get(i,1,1)<<endl;exit(1);}
    }

    for (int i=NumPhysElems+1;i<=NumElems;i++)
    for(int k=0;k<3;k++)
    {amin[i-1][k]=1.0;}


    if(indicator>0)
    {//Some average values will be negative if we dont limit, so we must limit.
    // limiting on pressure
    for (int i=1; i<=NumPhysElems; i++)
    {
        double Area = Mesh.get_area_prim(i);
        for( int m=1; m <= meqn; m++ )  
            {qlf[m-1]=Area*q.get(i,m,1) - dt*(FeLLF.get(i,m,1)+FeLLF.get(i,m,2)+FeLLF.get(i,m,3));}

        rho    = qlf[0];
        u1     = qlf[1]/rho;
        u2     = qlf[2]/rho;
        u3     = qlf[3]/rho;
        energy = qlf[4];
        plow  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

        if (plow < 0.0)
            {cout << "Negative solution in Lax-Fridrichs flux "<< plow <<" "<<rho<<" "<<qlf[0]*qlf[4]-0.5*(qlf[1]*qlf[1]+qlf[2]*qlf[2]+qlf[3]*qlf[3])<<" "<<q.get(i,5,1)*q.get(i,1,1)-0.5*(q.get(i,2,1)*q.get(i,2,1)+q.get(i,3,1)*q.get(i,3,1)+q.get(i,4,1)*q.get(i,4,1))<< endl;exit(1);}


        for( int m=1; m <= meqn; m++ )
        {
            fhat_local[0][m-1] = FeI.get(i,m,1);
            fhat_local[1][m-1] = FeI.get(i,m,2);
            fhat_local[2][m-1] = FeI.get(i,m,3);
            flf_local[0][m-1]  =  FeLLF.get(i,m,1);
            flf_local[1][m-1]  =  FeLLF.get(i,m,2);
            flf_local[2][m-1]  =  FeLLF.get(i,m,3);
            qtmp[m-1] = q.get(i,m,1);
        }

        double rescale[2][2][2];
        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        for(int i3=0; i3<2; i3++)
        {    rescale[i1][i2][i3] = 1.0;}

        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        for(int i3=0; i3<2; i3++)
        { 
            aa[0] = double(i1)*amin[i-1][0];
            aa[1] = double(i2)*amin[i-1][1];
            aa[2] = double(i3)*amin[i-1][2];
            //compute the fluxes evaluated at the corners of the fesable region
            //where rho is positive
            for( int m=1; m <= meqn; m++ )
            {
                ff[0] = -dt*(aa[0]*fhat_local[0][m-1]+(1.0-aa[0])*flf_local[0][m-1]);
                ff[1] = -dt*(aa[1]*fhat_local[1][m-1]+(1.0-aa[1])*flf_local[1][m-1]);
                ff[2] = -dt*(aa[2]*fhat_local[2][m-1]+(1.0-aa[2])*flf_local[2][m-1]);
                qh[m-1] = Area*qtmp[m-1] + (ff[0]+ff[1]+ff[2]);
            }

            rho    = qh[0];
            u1     = qh[1]/rho;
            u2     = qh[2]/rho;
            u3     = qh[3]/rho;
            energy = qh[4];
            phigh  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

        //Solve for scalings of these flux values that allow positive density and pressure
        if (phigh < 0.0)
            {
                root1 = 0.0;
                root2 = 1.0;

                for( int m=1; m <= meqn; m++ )
                {   qtmp2[m-1] = 0.0;}

                // bisection
                for(int nn = 1; nn<10; nn++){
                    rate = (root1+root2)/2.0;
                    for( int m=1; m <= meqn; m++ )
                    {
                        qtmp2[m-1] = rate*qh[m-1]+(1.0-rate)*qlf[m-1];
                    }
                rho    = qtmp2[0];
                u1     = qtmp2[1]/rho;
                u2     = qtmp2[2]/rho;
                u3     = qtmp2[3]/rho;
                energy = qtmp2[4];
                phigh  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3))-eps;
                if (phigh < 0.0)
                    root2 = rate;
                else
                    root1 = rate;
                }
                if (phigh > 0.0)
                    rate = rate;
                else
                    rate = root1;

                rescale[i1][i2][i3] = rate;

/*
           for( int m=1; m <= meqn; m++ )
               qtmp2[m-1] = rate*qh[m-1]+(1.0-rate)*qlf[m-1];

           rho    = qtmp2[0];
           u1     = qtmp2[1]/rho;
           u2     = qtmp2[2]/rho;
           u3     = qtmp2[3]/rho;
           energy = qtmp2[4];
           phigh2  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

           if (phigh2 < 0.0){
               cout << " something is wrong (i,j)="<<i<<" "<<j<<" "
                    << phigh2<<" "<< phigh<<" plow="<<plow<< endl;
               cout << aa[0]<<" "<<aa[1]<<" "
                    <<" "<<aa[2]<<" "<<aa[3]<<endl;
               cout << rate<< " "<< i1<<i2<<i3<<i4 << endl;
           }
*/

            }
         }
        //find the minimum rescaling parameter for each edge
        double rescale2[3] = {1.0, 1.0, 1.0};
        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        {
            rescale2[0] = Min(rescale2[0],rescale[1][i1][i2]);
            rescale2[1] = Min(rescale2[1],rescale[i1][1][i2]);
            rescale2[2] = Min(rescale2[2],rescale[i1][i2][1]);
        }
            
        for(int i1=0; i1<3; i1++)
            {amin[i-1][i1] = rescale2[i1]*amin[i-1][i1];}
    }
    //find the minimum theta for each edge (each edge can appear at most two times in this loop)     
    for (int i=1;i<=NumEdges;i++)
    {
      thex[i-1]=1.0;
    }
    //Minimize over cells
    for (int i=1;i<=NumPhysElems;i++)
    {
       for(int k=1;k<=3;k++)
       {
         int ie=Mesh.get_tedge(i,k);//eI.get(i,k);
       if(ie>0)
       {
         //if(thex[ie-1]<1.0){cout<<ie<<" "<<thex[ie]<<endl;}
         thex[ie-1]=Min(thex[ie-1],amin[i-1][k-1]);
       }
       }

    }
    //Scale the fluxes as a linear combination of the high order flux
    //and the low order LLF flux
    for (int i=1;i<=NumPhysElems;i++)
    {
       for(int k=1;k<=3;k++)
       {
         int ie=Mesh.get_tedge(i,k);//eI.get(i,k);
       if(ie>0)
       {
         for( int m=1; m <= meqn; m++ )  
         {FeI.set(i,m,k,(FeI.get(i,m,k)-FeLLF.get(i,m,k))*thex[ie-1]+FeLLF.get(i,m,k));}
       }
       //else{cout<<i<<" bad things "<<endl;exit(1);}
       }

    }

}

}

