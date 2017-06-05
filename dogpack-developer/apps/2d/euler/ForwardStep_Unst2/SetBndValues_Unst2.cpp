#include "dogdefs.h"
#include <cmath>
#include "mesh.h"
#include "dog_math.h"
#include <iostream>
using namespace std;


void qinitfuncs(double x,double y,double &rho,double &press,double &u1,double &u2,double &u3,double &energy)
{

    const double gamma = 1.4;
    const double gm1   = 1.4 - 1.0;
    const double gp1   = 1.4 + 1.0;

    const double Mach  = 5.09;
    const double M2    = Mach*Mach;

const double vs = Mach*sqrt( gamma*press/rho );

        // Correct with the Mach number for incoming shock
        {
            press = press*(2.0*gamma*M2- gm1)/gp1;
            rho   = rho*gp1*M2 / (gm1*M2 + 2.0);
            u1    = vs*(1.0 - 1.4 / rho );
            energy = press/(gamma-1.0e0)
                + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);
        }
        if(press<0.0){press=1.0e-12;}
} 


// This is a user-supplied routine that sets the the boundary conditions
//
void SetBndValues_Unst(const mesh& Mesh, dTensor3* q, dTensor3* aux)
{   
    int meqn = q->getsize(2);
    int kmax = q->getsize(3);
    int maux = aux->getsize(2);
    int      NumElems = Mesh.get_NumElems();
    int  NumPhysElems = Mesh.get_NumPhysElems();
    int NumGhostElems = Mesh.get_NumGhostElems();
    int      NumNodes = Mesh.get_NumNodes();
    int  NumPhysNodes = Mesh.get_NumPhysNodes();
    int      NumEdges = Mesh.get_NumEdges();

    // ----------------------------------------
    // Loop over each ghost cell element and
    // place the correct information into
    // these elements
    // ----------------------------------------

        double  rho   =  1.0;
        double  u1    =  3.0;
        double  u2    =  0.0;
        double  u3    =  0.0;
        double  press =  1.0/1.4;
        double gamma  =  1.4;
        double energy = press/(gamma-1.0e0)
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);


    for (int i=1; i<=NumGhostElems; i++)
    {
        int j = Mesh.get_ghost_link(i);

       double x1=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,1),1);
       double x2=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,2),1);
       double x3=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,3),1);

       double y1=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,1),2);
       double y2=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,2),2);
       double y3=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,3),2);
 
       double xi=(x1+x2+x3)/3.0;
       double yi=(y1+y2+y3)/3.0;


       double xe,ye;
       int i1,i2,i3=-1;
       i3=0;
       for(int l=1;l<=3;l++)
       {
           double xo=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,l),1);
           double yo=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,l),2);
           for (int k=1;k<=3;k++)
           {
              
              double xp=Mesh.get_node(Mesh.get_tnode(j,k),1);
              double yp=Mesh.get_node(Mesh.get_tnode(j,k),2);
              if( fabs(xo-xp)<1.0e-14 && fabs(yo-yp)<1.0e-14)
              {
                  if(i3==0)
                  {i1=l;i3=1;}
                  else 
                  {i2=l;}
              }
           }
       }    

       if(i1==-1 || i2==-1)
       {
           printf("problem finding boundary edge");exit(1);
       }

           double xo1=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,i1),1);
           double yo1=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,i1),2);
           double xo2=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,i2),1);
           double yo2=Mesh.get_node(Mesh.get_tnode(i+NumPhysElems,i2),2);

       xe=0.5*(xo1+xo2);
       ye=0.5*(yo1+yo2);

       if(xi<1.0e-14)
       {
          q->set(i+NumPhysElems,1,1,rho);
          q->set(i+NumPhysElems,2,1,u1*rho);
          q->set(i+NumPhysElems,3,1,u2*rho);
          q->set(i+NumPhysElems,4,1,0.0);
          q->set(i+NumPhysElems,5,1,energy);
       }
       else if(xi>3.0-1.0e-14)
       {
          /*
          q->set(i+NumPhysElems,1,1,q->get(j,1,1));
          q->set(i+NumPhysElems,2,1,q->get(j,2,1));
          q->set(i+NumPhysElems,3,1,q->get(j,3,1));
          q->set(i+NumPhysElems,4,1,q->get(j,4,1));
          q->set(i+NumPhysElems,5,1,q->get(j,5,1));*/

        for (int m=1; m<=meqn; m++)
        {

            for (int k=1; k<=kmax; k++)
            {
                q->set(i+NumPhysElems,m,k, q->get(j,m,k) );
            }
        }

       }
       else if(yi>1.0-1.0e-14)
       {
         /*
          q->set(i+NumPhysElems,1,1,q->get(j,1,1));
          q->set(i+NumPhysElems,2,1,q->get(j,2,1));
          q->set(i+NumPhysElems,3,1,-q->get(j,3,1));
          q->set(i+NumPhysElems,4,1,q->get(j,4,1));
          q->set(i+NumPhysElems,5,1,q->get(j,5,1));*/
        for (int m=1; m<=meqn; m++)
        {
       
            for (int k=1; k<=kmax; k++)
            {
                q->set(i+NumPhysElems,m,k, q->get(j,m,k) );
            }
        } 

       }
       else if(xi<0.6+1.0e-14 && yi<1.0e-14)
       {
       /*   q->set(i+NumPhysElems,1,1,q->get(j,1,1));
          q->set(i+NumPhysElems,2,1,q->get(j,2,1));
          q->set(i+NumPhysElems,3,1,-q->get(j,3,1));
          q->set(i+NumPhysElems,4,1,q->get(j,4,1));
          q->set(i+NumPhysElems,5,1,q->get(j,5,1));
       */
        for (int m=1; m<=meqn; m++)
        {
       
            for (int k=1; k<=kmax; k++)
            {
                q->set(i+NumPhysElems,m,k, q->get(j,m,k) );
            }
        } 

       }
       else if(ye<0.2-1.0e-14 && fabs(xe-0.6)<1.0e-14)
       {
         /* q->set(i+NumPhysElems,1,1,q->get(j,1,1));
          q->set(i+NumPhysElems,2,1,-q->get(j,2,1));
          q->set(i+NumPhysElems,3,1,q->get(j,3,1));
          q->set(i+NumPhysElems,4,1,q->get(j,4,1));
          q->set(i+NumPhysElems,5,1,q->get(j,5,1));*/
        for (int m=1; m<=meqn; m++)
        {
       
            for (int k=1; k<=kmax; k++)
            {
                q->set(i+NumPhysElems,m,k, q->get(j,m,k) );
            }
        } 

       }
       else 
       {
         /* q->set(i+NumPhysElems,1,1,q->get(j,1,1));
          q->set(i+NumPhysElems,2,1,q->get(j,2,1));
          q->set(i+NumPhysElems,3,1,-q->get(j,3,1));
          q->set(i+NumPhysElems,4,1,q->get(j,4,1));
          q->set(i+NumPhysElems,5,1,q->get(j,5,1));*/
        for (int m=1; m<=meqn; m++)
        {
       
            for (int k=1; k<=kmax; k++)
            {
                q->set(i+NumPhysElems,m,k, q->get(j,m,k) );
            }
        } 

       }

        /*
        for (int m=1; m<=meqn; m++)
        {

            for (int k=2; k<=kmax; k++)
            {
                q->set(i+NumPhysElems,m,k, 0.0 );
            }
        }*/

        for (int m=1; m<=maux; m++)
        for (int k=1; k<=kmax; k++)
        {
            aux->set(i+NumPhysElems,m,k, aux->get(j,m,k) );
        }
    }

}
