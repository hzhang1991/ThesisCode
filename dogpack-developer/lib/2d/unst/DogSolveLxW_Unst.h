#ifndef _DOGSOLVELXW_UNST_H_
#define _DOGSOLVELXW_UNST_H_

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
// Single step of the method.  See also UpdateSoln_Unst
// -------------------------------------------------------------------------- //
void StepLxW( double dt, const dTensor3& qold,
    const dTensor3& L, dTensor3& qnew );

// -------------------------------------------------------------------------- //
// Called once after each successful step.  This is useful for tracking scalar
// quantities of interest.  (e.g. mass conservation [which is automatic], or
// L2-norm of an electric field, etc.)
// -------------------------------------------------------------------------- //
void ConSoln_Unst(const mesh& Mesh, 
          const dTensor3& aux, const dTensor3& q, 
          double t, const char* outputdir);

// -------------------------------------------------------------------------- //
// Function for defining a CFL number
// -------------------------------------------------------------------------- //
double GetCFL_Unst(double dt, const mesh& Mesh,
    const dTensor3& aux, const dTensor1& smax);
// -------------------------------------------------------------------------- //

void LaxWendroff_Unst(double dt,
    const mesh& Mesh, const edge_data_Unst& EdgeData,
    dTensor3& aux,                  // SetBndValues modifies ghost cells
    dTensor3& q,                    // SetBndValues modifies ghost cells
    dTensor3& Lstar, dTensor1& smax);

#endif
