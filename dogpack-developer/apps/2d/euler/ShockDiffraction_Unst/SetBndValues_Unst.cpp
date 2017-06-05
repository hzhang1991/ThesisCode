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
    for (int i=1; i<=NumGhostElems; i++)
    {
        int j = Mesh.get_ghost_link(i);

        for (int m=1; m<=meqn; m++)
        {
            q->set(i+NumPhysElems,m,1, q->get(j,m,1) );

        for (int k=2; k<=kmax; k++)
        {
            q->set(i+NumPhysElems,m,k, 0.0 );
        }
        }

        for (int m=1; m<=maux; m++)
        for (int k=1; k<=kmax; k++)
        {
            aux->set(i+NumPhysElems,m,k, aux->get(j,m,k) );
        }
    }

}
