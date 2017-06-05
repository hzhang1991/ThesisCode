#include "tensors.h"
#include "assert.h"

// Update the solution using the constructed Lstar
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensor2& node,const dTensorBC3& aux,
        const dTensorBC3& qstar, const dTensorBC3& Lstar,
        dTensorBC3& qnew)
{
    const int     mx = qnew.getsize(1);
    const int   meqn = qnew.getsize(2);
    const int   kmax = qnew.getsize(3);
    const int   maux = aux.getsize(2);
    const int mbc = qnew.getmbc();

#pragma omp parallel for
    for (int j=(2-mbc); j<=(mx+mbc-1); j++)
    for (int m=1; m<=meqn; m++)        
    for (int k=1; k<=kmax; k++)
    {
        double tmp = alpha1*qstar.get(j,m,k) + alpha2*qnew.get(j,m,k)
            + beta*dt*Lstar.get(j,m,k);
        qnew.set(j,m,k, tmp );
    }

    // Optional call to modify updated solution
    void AfterUpdateSoln(const dTensor2& node,
            const dTensorBC3& aux,
            dTensorBC3& q,
            double dt,
            double beta);
    AfterUpdateSoln(node,aux,qnew,dt,beta); 
}

// Update the solution using the constructed Lstar
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const dTensor2& node,const dTensorBC3& aux,
    const dTensorBC3& qold, const dTensorBC3& Lstar,
    dTensorBC3& q1, dTensorBC3& q2)
{
    const int     mx = q1.getsize(1);
    const int   meqn = q1.getsize(2);
    const int   kmax = q1.getsize(3);
    const int   maux = aux.getsize(2);
    const int mbc = q1.getmbc();

// TODO - do we want to replace these with vget? -DS 06/07/2014
#pragma omp parallel for
    for (int j=(2-mbc); j<=(mx+mbc-1); j++)
    for (int m=1; m<=meqn; m++)        
    for (int k=1; k<=kmax; k++)
    {

        double s1 = q1.get(j,m,k);
        double s3 = qold.get(j,m,k);

        // update q2
        double s2 = q2.get(j,m,k) + delta*s1;
        q2.set(j,m,k, s2 );

        // update q
        double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.get(j,m,k);
        q1.set(j,m,k, tmp );
    }

    // Optional call to modify updated solution
    void AfterUpdateSoln(const dTensor2& node,
            const dTensorBC3& aux,
            dTensorBC3& q,
            double dt,
            double beta);
    AfterUpdateSoln(node,aux,q1,dt,beta); 
}
