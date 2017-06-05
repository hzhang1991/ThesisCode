#include "dogdefs.h"
#include "mesh.h"

// This routine takes a single Euler step on the solution.  The bulk of
// DoGPack uses the so-called low-storage methods, and therefore these "Euler
// steps" don't look like traditional Euler steps.
// In particular, the SSP solvers sometimes require convex combinations of two
// stage values q1, and q2, as well as a rhs that has been evaluated.

void AfterUpdateSoln_Unst(const mesh& Mesh, 
        dTensor3& aux, dTensor3& q,
        double dt, double beta);

void UpdateSoln_Unst(const double alpha1, const double alpha2, 
        const double beta, const double dt, const mesh& Mesh,
        dTensor3& aux, const dTensor3& qstar, 
        const dTensor3& Lstar, dTensor3& qnew)
{

    const int NumElems = qnew.getsize(1);
    const int meqn     = qnew.getsize(2);
    const int kmax     = qnew.getsize(3);
    const int NumPhysElems = Mesh.get_NumPhysElems();

    // Update solution
#pragma omp parallel for
    for (int i=1; i<=NumPhysElems; i++)
    for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
    {
        double tmp = alpha1*qstar.get(i,m,k) + alpha2*qnew.get(i,m,k)
                    + beta*dt*Lstar.get(i,m,k);

        qnew.set(i,m,k, tmp );
    }

    // Optional call to modify updated solution
    AfterUpdateSoln_Unst(Mesh,aux,qnew,dt,beta);    

}


// Update the solution using the constructed Lstar
void UpdateSoln_Unst(double g1,double g2, double g3, double delta, 
        double beta,double dt,const mesh& Mesh,
        dTensor3& aux,const dTensor3& qold, 
        const dTensor3& Lstar,
        dTensor3& q1, dTensor3& q2)
{

    const int NumElems = q1.getsize(1);
    const int meqn     = q1.getsize(2);
    const int kmax     = q1.getsize(3);
    const int NumPhysElems = Mesh.get_NumPhysElems();  

#pragma omp parallel for
    for (int i=1; i<=NumPhysElems; i++)
    for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
    {	  
        double s1 = q1.get(i,m,k);
        double s3 = qold.get(i,m,k);

        // update q2
        double s2 = q2.get(i,m,k) + delta*s1;
        q2.set(i,m,k, s2 );

        // update q
        double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.get(i,m,k);
        q1.set(i,m,k, tmp );
    }

    // Optional call to modify updated solution
    AfterUpdateSoln_Unst(Mesh,aux,q1,dt,beta);

}
