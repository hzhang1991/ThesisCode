#include "tensors.h"

void AfterUpdateSoln(const dTensorBC4& aux, dTensorBC4& q,
        double dt, double beta);

// ======= EAJMERGE: I moved this content to ~/lib/DogSolver.cpp
// on branch r1304, but when merging to head I restored this file
// to get applications working again.
void UpdateSoln(double alpha1, 
        double alpha2, 
        double beta, 
        double dt,
        const dTensorBC4& aux, 
        const dTensorBC4& qstar, const dTensorBC4& Lstar, 
        dTensorBC4& qnew)
{
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();    

    // Update solution
#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    for (int j=(2-mbc); j<=(my+mbc-1); j++)
    for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
    {
        double tmp = alpha1*qstar.get(i,j,m,k) + alpha2*qnew.get(i,j,m,k)
            + beta*dt*Lstar.get(i,j,m,k);

        qnew.set(i,j,m,k, tmp );
    }

    // Optional call to modify updated solution
    AfterUpdateSoln(aux,qnew,dt,beta);    
}

// Update the solution using the constructed Lstar
void UpdateSoln(double g1,
        double g2,
        double g3,
        double delta, 
        double beta,
        double dt,
        const dTensorBC4& aux,
        const dTensorBC4& qold,
        const dTensorBC4& Lstar,
        dTensorBC4& q1,
        dTensorBC4& q2)
{
    const int mx   = q1.getsize(1);
    const int my   = q1.getsize(2);
    const int meqn = q1.getsize(3);
    const int kmax = q1.getsize(4);
    const int mbc  = q1.getmbc();    

#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    for (int j=(2-mbc); j<=(my+mbc-1); j++)
    for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
    {

        double s1 = q1.get(i,j,m,k);
        double s3 = qold.get(i,j,m,k);

        // update q2
        double s2 = q2.get(i,j,m,k) + delta*s1;
        q2.set(i, j,m,k, s2 );

        // update q
        double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.get(i,j,m,k);
        q1.set(i, j, m, k, tmp );
    }

    // Optional call to modify updated solution
    // [I commented this out to get this to compile. -eaj]
    //AfterUpdateSoln(aux,q1,dt,beta);    

}
