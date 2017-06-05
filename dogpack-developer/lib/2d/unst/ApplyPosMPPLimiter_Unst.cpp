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

void ApplyPosMPPLimiter_Unst( const double dt, const mesh& Mesh,const edge_data_Unst& EdgeData, 
        const dTensor1& smax,
        const dTensor3& aux, const dTensor3& q, 
        dTensor3& FeI )
{
    const int NumElems      = Mesh.get_NumElems();
    const int NumPhysElems  = Mesh.get_NumPhysElems();
    const int NumEdges      = Mesh.get_NumEdges();
    const int meqn          = q.getsize(2);
    const int kmax          = q.getsize(3);
    const int maux          = aux.getsize(2);
    const int space_order   = dogParams.get_space_order();


    const double eps = 1.0e-15;
    const double rho_min = eps;

    dTensor3 FeLLF(NumElems, meqn, 3 );  FeLLF.setall(0.);
    void ConstructL_LLF_Unst(const mesh& Mesh,
        const edge_data_Unst& EdgeData,
        const dTensor3& aux,     
        const dTensor3& q,       
        const dTensor1& smax, dTensor3& FeI);
    ConstructL_LLF_Unst(Mesh,EdgeData,aux,q,smax,FeLLF);
    
    double gmin[NumElems],amin[NumElems][3],thex[NumEdges],ff[3];
    int isgn[3];
   

    int indicator=0; 

    for (int i=1;i<=NumPhysElems;i++)
    {  

       double Area = Mesh.get_area_prim(i);
       double gtmp=Area*q.get(i,1,1)-dt*(FeI.get(i,1,1)+FeI.get(i,1,2)+FeI.get(i,1,3));
       if(gtmp<0.0){indicator=1;}

       gmin[i-1]=-rho_min+Area*q.get(i,1,1);
       double ffsum=0.0;

       for(int k=1;k<=3;k++)
       {
       gmin[i-1]=gmin[i-1]-dt*FeLLF.get(i,1,k);
       ff[k-1]=-dt*(FeI.get(i,1,k)-FeLLF.get(i,1,k));
       if (ff[k-1]<0.0)
       {   isgn[k-1]=1;}
       else{isgn[k-1]=0;}
       ffsum = ffsum + double(isgn[k-1])*ff[k-1];
       }
       //Comp. of Lambdas
       for (int k=0; k<3; k++)
       {
            if (isgn[k]==1)
                {amin[i-1][k] = Min(-gmin[i-1]/(ffsum-eps),1.e0);}
            else
                {amin[i-1][k] = 1.e0;}

       }

       if(q.get(i,1,1)<0.0)
       {cout<<"cell averages negative "<<q.get(i,1,1)<<endl;exit(1);}
    }
    cout<<indicator<<endl;
    if(indicator>0)
    {cout<<indicator<<endl;
    for (int i=1;i<=NumEdges;i++)
    {
      thex[i-1]=1.0;
    }
    //Minimize over cells
    for (int i=1;i<=NumElems;i++)
    {
       for(int k=1;k<=3;k++)
       {
         int ie=Mesh.get_tedge(i,k);
         thex[ie]=Min(thex[ie],amin[i-1][k-1]);
       }

    }

    for (int i=1;i<=NumElems;i++)
    {
       for(int k=1;k<=3;k++)
       {
         int ie=Mesh.get_tedge(i,k);
         FeI.set(i,1,k,(FeI.get(i,1,k)-FeLLF.get(i,1,k))*thex[ie]+FeLLF.get(i,1,k));
       }

    }
   }

}

