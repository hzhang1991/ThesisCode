#ifndef _DOGSOLVERK_UNST_H_
#define _DOGSOLVERK_UNST_H_

#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data_Unst.h"
#include "RKinfo.h"
#include "mesh.h"
#include "DogParams.h"
#include "DogStateUnst2.h"

// ------------------------------------------------------------
// Function definitions
void CopyQ_Unst(const dTensor3&,dTensor3&);
void ConSoln_Unst(const mesh& Mesh, 
        const dTensor3& aux, const dTensor3& q, 
        double t, string outputdir);

// Single RK time step functions:
void UpdateSoln_Unst(const double alpha1, const double alpha2, 
        const double beta, const double dt, const mesh& Mesh,
        const dTensor3& aux, const dTensor3& qstar, 
        const dTensor3& Lstar, dTensor3& qnew);
void UpdateSoln_Unst(double g1,double g2, double g3, double delta, 
        double beta,double dt,const mesh& Mesh,
        dTensor3& aux,const dTensor3& qold, 
        const dTensor3& Lstar,
        dTensor3& q1, dTensor3& q2);

void BeforeStep_Unst(const double,const mesh&,dTensor3&,dTensor3&);
void  AfterStep_Unst(const double,const mesh&,dTensor3&,dTensor3&);
void AfterFullTimeStep_Unst(const double dt, const mesh& Mesh,
        const dTensor3& auxold, const dTensor3& qold,
        const dTensor3& Lold, dTensor3& aux, dTensor3& q);
void AfterReject_Unst(const mesh& Mesh, const double dt, dTensor3& aux, dTensor3& q);
void UpdateSoln_Unst(const double alpha1, const double alpha2, 
        const double beta, const double dt, const mesh& Mesh,
        dTensor3& aux, const dTensor3& qstar, 
        const dTensor3& Lstar, dTensor3& qnew);

// RHS function
void ConstructL_Unst(const mesh& Mesh,
        const edge_data_Unst& EdgeData,
        dTensor3& aux, // SetBndValues modifies ghost cells
        dTensor3& q,   // SetBndValues modifies ghost cells
        dTensor3& Lstar,
        dTensor1& smax);

// Function for defining a CFL number:
double GetCFL_Unst(double dt, const mesh& Mesh,
                   const dTensor3& aux, const dTensor1& smax);
// ------------------------------------------------------------

#endif
