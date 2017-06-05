//#include "tensors.h"

#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "RKinfo.h"
#include "DogState.h"
using namespace std;


void StepSDC(const double& dt, const int method[], const dTensor2& node,
        dTensorBC1& smax, dTensorBC3& Lrhs,
        dTensorBC3& aux, dTensorBC3& qin, 
        dTensorBC3& qnew)
{

    void BeforeStep(double dt,const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    void AfterStep(double dt, const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    void ConstructL(const int[],const dTensor2&,dTensorBC3&,dTensorBC3&,
            dTensorBC3&,dTensorBC1& smax);
    void ApplyLimiter(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q, 
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    void RelaxLimiter(const dTensor2& node,dTensorBC3& aux,dTensorBC3& q);
    void EulerStepSDC(const double& dt, const dTensorBC3& aux, const dTensorBC3& qold, 
            const dTensorBC3& Lrhs, dTensorBC3& qnew);

    // needed for limiter
    void ProjectRightEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);
    void ProjectLeftEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);

    // Do all the stuff necessary for taking an euler step,then take an euler
    // step
    BeforeStep(dt, node, aux, qin);
    ConstructL(method,node,aux,qin,Lrhs,smax);
    EulerStepSDC(dt,aux,qin,Lrhs,qnew);
    if (dogParams.using_moment_limiter())
    { ApplyLimiter(node,aux,qnew,&ProjectRightEig,&ProjectLeftEig); }
    else if (dogParams.using_relax_limiter())
    {  RelaxLimiter(node,aux,qnew);  }
    AfterStep(dt,node,aux,qnew);

}

// Euler time stepping for the spectral deferred correction method
void EulerStepSDC(const double& dt, const dTensorBC3& aux, const dTensorBC3& qold, 
        const dTensorBC3& Lrhs, dTensorBC3& qnew)
{
    const int mbc  = qold.getmbc();
    const int imax = qold.getsize(1);
    const int mmax = qold.getsize(2);
    const int kmax = qold.getsize(3);

    // Euler time step loop
#pragma omp parallel for
    for (int i=(1-mbc); i<=(imax+mbc); i++)
        for (int m=1; m<=mmax; m++)
            for (int k=1; k<=kmax; k++)
            {
                double tmp = qold.get(i,m,k) + dt*Lrhs.get(i,m,k);
                qnew.set(i,m,k, tmp);
            }


}

void StepSDCRK2(const double& dt, const int method[], const dTensor2& node,
        dTensorBC1& smax, dTensorBC3& Lrhs, dTensorBC3& Lstar,
        dTensorBC3& aux, dTensorBC3& qin, dTensorBC3& qstar,
        dTensorBC3& qnew)
{

    void BeforeStep(double dt,const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    void AfterStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
    void ConstructL(const int[],const dTensor2&,dTensorBC3&,dTensorBC3&,
            dTensorBC3&,dTensorBC1& smax);
    void ApplyLimiter(const dTensor2&,const dTensorBC3&,dTensorBC3&,
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    void RelaxLimiter(const dTensor2& node,dTensorBC3& aux,dTensorBC3& q);
    void EulerStepSDC(const double& dt, const dTensorBC3& aux, const dTensorBC3& qold, 
            const dTensorBC3& Lrhs, dTensorBC3& qnew);

    // needed for limiter
    void ProjectRightEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);
    void ProjectLeftEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);

    void CopyQ(const dTensorBC3&,dTensorBC3&);

    const int mbc  = qin.getmbc();
    const int imax = qin.getsize(1);
    const int mmax = qin.getsize(2);
    const int kmax = qin.getsize(3);

    // Take a Single Euler Step to produce Qs
    BeforeStep(dt, node, aux, qin);
    ConstructL(method,node,aux,qin,Lrhs,smax);
    EulerStepSDC(dt,aux,qin,Lrhs,qstar);
    if (dogParams.using_moment_limiter())
    { ApplyLimiter(node,aux,qstar,&ProjectRightEig,&ProjectLeftEig); }
    else if (dogParams.using_relax_limiter())
    {  RelaxLimiter(node,aux,qstar);  }
    AfterStep(dt,node,aux,qstar);

    BeforeStep(dt, node, aux, qstar );
    ConstructL(method,node,aux,qstar,Lstar,smax);

    // Set the rhs to be the average of the two right hand sides
#pragma omp parallel for
    for (int i=(1-mbc); i<=(imax+mbc); i++)
        for (int m=1; m<=mmax; m++)
            for (int k=1; k<=kmax; k++)
            {
                double tmp = qin.get(i,m,k) + 0.5*dt*( Lrhs.get(i,m,k) + Lstar.get(i,m,k) );
                qnew.set(i,m,k,tmp);
            }
    if (dogParams.using_moment_limiter())
    { ApplyLimiter(node,aux,qnew,&ProjectRightEig,&ProjectLeftEig); }
    else if (dogParams.using_relax_limiter())
    {  RelaxLimiter(node,aux,qnew);  }
    AfterStep(dt,node,aux,qnew);

}

void StepSDCdeltaRK2(const double& dt, const int method[], const dTensor2& node,
        dTensorBC1& smax, 
        dTensorBC3& aux, 
        dTensorBC3& qstar,
        dTensorBC3& Lstar,
        dTensorBC3& L1, 
        dTensorBC3& L1new, 
        dTensorBC3& L2, 
        dTensorBC3& q1, 
        dTensorBC3& q2,
        int num, dTensorBC4& IL)
{

    // This function advances the error equation:
    //
    //   \dot{de} = -\dot{Q} + f(Q) + [ f(Q+de) - f(Q) ]
    //
    // And on output adds the correction term into Q.
    //
    // q1 is the old time value with the previous error already adde din.
    // q2 is the new value.  
    //
    // It is assumed that the previous error has already been added into q1

    void BeforeStep(double dt,const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    void AfterStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
    void ConstructL(const int[],const dTensor2&,dTensorBC3&,dTensorBC3&,
            dTensorBC3&,dTensorBC1& smax);
    void ApplyLimiter(const dTensor2&,const dTensorBC3&,dTensorBC3&,
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    void RelaxLimiter(const dTensor2& node,dTensorBC3& aux,dTensorBC3& q);

    // needed for limiter
    void ProjectRightEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);
    void ProjectLeftEig(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&);

    void CopyQ(const dTensorBC3&,dTensorBC3&);

    const int mbc  = q1.getmbc();
    const int imax = q1.getsize(1);
    const int mmax = q1.getsize(2);
    const int kmax = q1.getsize(3);

    // Take a Single Euler Step to produce Qs
    // BeforeStep(dt, node, aux, q1);
    // ConstructL(method,node,aux,q1,L1new,smax);


#pragma omp parallel for
    for (int i=(1-mbc); i<=(imax+mbc); i++)
        for (int m=1; m<=mmax; m++)
            for (int k=1; k<=kmax; k++)
            {
                double tmp = q1.get(i,m,k) + IL.get(i,m,k,num) + dt*( L1new.get(i,m,k) - L1.get(i,m,k) );
                qstar.set(i,m,k,tmp);
            }
    if (dogParams.using_moment_limiter())
    { ApplyLimiter(node,aux,qstar,&ProjectRightEig,&ProjectLeftEig); }
    else if (dogParams.using_relax_limiter())
    {  RelaxLimiter(node,aux,qstar);  }
    AfterStep(dt,node,aux,qstar);

    // Temporary storage
    ConstructL(method,node,aux,qstar,Lstar,smax);

    // Set the rhs to be the average of the two right hand sides
#pragma omp parallel for
    for (int i=(1-mbc); i<=(imax+mbc); i++)
        for (int m=1; m<=mmax; m++)
            for (int k=1; k<=kmax; k++)
            {
                double tmp = q1.get(i,m,k) + IL.get(i,m,k,num) + 0.5*dt*( (Lstar.get(i,m,k) - L2.get(i,m,k) ) 
                        + ( L1new.get(i,m,k) - L1.get(i,m,k) ) );
                q2.set(i,m,k,tmp);
            }
    if (dogParams.using_moment_limiter())
    { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
    else if (dogParams.using_relax_limiter())
    {  RelaxLimiter(node,aux,q2);  }
    AfterStep(dt,node,aux,q2);

}

