#ifndef _DOGSOLVE_HYBRID_H_
#define _DOGSOLVE_HYBRID_H_

#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data_Unst.h"
#include "RKinfo.h"
#include "mesh.h"
#include "DogParams.h"

#include "DogStateHybrid.h"
#include "DogStateUnst2.h"
#include "DogParamsCart2.h"

#include "QuadratureRules.h"

// Function definitions:
void DogSolveRK_Unst(const dTensor2& vel_vec,
    const mesh& Mesh, const edge_data_Unst& EdgeData,
        dTensor3& aux, dTensor3& qold, dTensor3& qnew, 
        const double tstart, const double tend, DogStateUnst2& dogStateUnst2);

// Silent Solvers for each slice
//
// `silent' simply means they don't output any of their data
//
double DogSolveRK_Unst_Parallel(
    const QuadratureRules& QuadFuncs,
    const mesh& Mesh, const edge_data_Unst& EdgeData,
    dTensorBC5& q, const double tstart, const double tend );
void DogSolveSL_Parallel( const QuadratureRules& QuadFuncs, const mesh& Mesh, 
    dTensorBC5& q, double tstart, double tend );

#endif
