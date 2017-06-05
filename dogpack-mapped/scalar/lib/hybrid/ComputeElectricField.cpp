#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "mesh.h"
#include "assert.h"

// This is the 'do nothing' example that can be swapped out later if desired.
// See VlasovPoisson/lib for the better routine.

// This routine simply glues together many of the routines that are already
// written in the Poisson solver library
//
// phi( 1:SubNumPhysNodes       ) is a scalar quantity.  
//
// E1 ( 1:NumElems, 1:kmax_unst ) is a vector quantity.
// E2 ( 1:NumElems, 1:kmax_unst ) is a vector quantity.
//
// See also: ConvertEfieldOntoDGbasis
//
void ComputeElectricField( const double t, const mesh& Mesh, const dTensorBC5& q,
    dTensor2& E1, dTensor2& E2)
{

    //
    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);
    const int     meqn = q.getsize(4);
    const int     kmax = q.getsize(5);

    const int space_order = dogParams.get_space_order();

    // unstructured parameters:
    const int kmax_unst    = E2.getsize(2);
    const int NumBndNodes  = Mesh.get_NumBndNodes();
    const int NumPhysNodes = Mesh.get_NumPhysNodes();
    const int SubFactor    = Mesh.get_SubFactor();

    
}
