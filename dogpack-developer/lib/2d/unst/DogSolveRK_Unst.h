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
// Copy qin into qout.  ( see also: qout.copyfrom( qin ) )
void CopyQ_Unst(const dTensor3& qin, dTensor3& qout);
// ------------------------------------------------------------

// -------------------------------------------------------------------------- //
// Called once after each successful step.  This is useful for tracking scalar
// quantities of interest.  (e.g. mass conservation [which is automatic], or
// L2-norm of an electric field, etc.)
// -------------------------------------------------------------------------- //
void ConSoln_Unst(const mesh& Mesh, 
          const dTensor3& aux, const dTensor3& q, 
          double t, const char* outputdir);

// -------------------------------------------------------------------------- //
// Routines for taking a single (Euler) time step on the solution.  Most of the
// integrators written in DoGPack are the low-storage methods
// -------------------------------------------------------------------------- //
void UpdateSoln_Unst(const double alpha1, const double alpha2, const double beta, 
    const double dt, const mesh& Mesh,
    const dTensor3& aux, const dTensor3& qstar, const dTensor3& Lstar, 
    dTensor3& qnew);

void UpdateSoln_Unst(double g1, double g2, double g3, 
    double delta, double beta, double dt,
    const mesh& Mesh, dTensor3& aux,
    const dTensor3& qold, const dTensor3& Lstar,
    dTensor3& q1, dTensor3& q2);

void UpdateSoln_Unst(const double alpha1, const double alpha2, 
    const double beta, const double dt, 
    const mesh& Mesh, dTensor3& aux, 
    const dTensor3& qstar, const dTensor3& Lstar, 
    dTensor3& qnew);
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// Called once before each stage
// -------------------------------------------------------------------------- //
void BeforeStep_Unst( const double dt, const mesh& Mesh,
    dTensor3& aux, dTensor3& q);

// -------------------------------------------------------------------------- //
// Called once after each stage
// -------------------------------------------------------------------------- //
void  AfterStep_Unst( const double dt, const mesh& Mesh,
    dTensor3& aux, dTensor3& q);

// -------------------------------------------------------------------------- //
// Called once after completing a full time step.
// -------------------------------------------------------------------------- //
void AfterFullTimeStep_Unst(const double dt, const mesh& Mesh,
    const dTensor3& auxold, const dTensor3& qold, const dTensor3& Lold, 
    dTensor3& aux, dTensor3& q);

// -------------------------------------------------------------------------- //
// Called only if the time step was rejected
// -------------------------------------------------------------------------- //
void AfterReject_Unst(const mesh& Mesh, const double dt, 
    dTensor3& aux, dTensor3& q);

// -------------------------------------------------------------------------- //
// RHS function in the MOL formulation
// -------------------------------------------------------------------------- //
void ConstructL_Unst(const mesh& Mesh, const edge_data_Unst& EdgeData,
    dTensor3& aux, // SetBndValues modifies ghost cells
    dTensor3& q,   // SetBndValues modifies ghost cells
    dTensor3& Lstar, dTensor1& smax);

// -------------------------------------------------------------------------- //
// Function for defining a CFL number
// -------------------------------------------------------------------------- //
double GetCFL_Unst(double dt, const mesh& Mesh,
    const dTensor3& aux, const dTensor1& smax);
// -------------------------------------------------------------------------- //

#endif
