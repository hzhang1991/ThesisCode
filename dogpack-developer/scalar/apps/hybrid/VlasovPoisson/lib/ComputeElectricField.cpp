#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "mesh.h"
#include <sstream>
#include <string>
#include "assert.h"
#include "SparseCholesky.h"

#include "VlasovParams.h"

//void ComputeEfield(const double t, const int space_order, const mesh& Mesh, const dTensor1& phi,
//    dTensor2& E1, dTensor2& E2);

void ComputeEfield(const int space_order, const mesh& Mesh, const dTensor1& phi,
    dTensor2& E1, dTensor2& E2);

// This routine simply glues together many of the routines that are already
// written in the Poisson solver library
//
// phi( 1:SubNumPhysNodes       ) is a scalar quantity.  
//
// E1 ( 1:NumElems, 1:kmax2d ) is a vector quantity.
// E2 ( 1:NumElems, 1:kmax2d ) is a vector quantity.
//
// See also: ConvertEfieldOntoDGbasis
void ComputeElectricField( const double t, const mesh& Mesh, const dTensorBC5& q,
    dTensor2& E1, dTensor2& E2)
{

    // Problem dimensions //
    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);
    const int     meqn = q.getsize(4);
    const int     kmax = q.getsize(5);

    const int space_order = dogParams.get_space_order();

    // unstructured parameters:
    const int kmax2d       = E2.getsize(2);
    const int NumBndNodes  = Mesh.get_NumBndNodes();
    const int NumPhysNodes = Mesh.get_NumPhysNodes();

    // Quick error check
    if( !Mesh.get_is_submesh() )
    {
        printf("ERROR: mesh needs to have subfactor set to %d\n", space_order);
        printf("Go to Unstructured mesh and remesh the problem\n");
        exit( 1 );
    }
    const int SubFactor    = Mesh.get_SubFactor();

    assert_eq( NumElems, Mesh.get_NumElems() );

    // -- Step 1: Compute rho -- //
    dTensor3 rho(NumElems, 1, kmax2d );

    // Compute moments performs: \int f(x,v) dv.
    //
    // In some cases, the Poisson solver needs to have background ions 
    // subtracted out, and the correct form for the right hand side is
    //
    //    rho = -\int f(x,v) dv + rho0,
    //
    // Where rho0 is the initial distribution.
    void ComputeDensity(const mesh& Mesh, const dTensorBC5& q, dTensor3& rho);
    ComputeDensity( Mesh, q, rho  );

    // Redefine rho as the correct rho:
    for( int n=1; n <= NumElems; n++ )
    {
        rho.set( n,1, 1, vlasovParams.rho0 - rho.get(n,1,1) );
        for( int k=2; k <= kmax2d; k++ )
        {
            rho.set( n, 1, k, -rho.get(n,1,k) );
        }
    }

// Compute the total charge in the domain:
const int  NumPhysElems = Mesh.get_NumPhysElems();
double rho_total = 0.0;
for( int n=1; n <= NumPhysElems; n++ )
{
    rho_total += rho.get(n,1,1) * Mesh.get_area_prim(n);
}
printf(" total rho = %2.15e\n", rho_total );

    // -- Step 2: Figure out how large phi needs to be -- //
    int SubNumPhysNodes = 0;
    int SubNumBndNodes  = 0;
    switch( dogParams.get_space_order() )
    {
        case 1:
            SubNumPhysNodes = NumPhysNodes;
            SubNumBndNodes  = NumBndNodes;
            break;

        case 2:
            SubNumPhysNodes = Mesh.get_SubNumPhysNodes();
            SubNumBndNodes  = Mesh.get_SubNumBndNodes();
            if(SubFactor!=2)
            {
                printf("\n");
                printf(" Error: for space_order = %i, need SubFactor = %i\n",space_order,2);
                printf("      SubFactor = %i\n",SubFactor);
                printf("\n");
                exit(1);
            }
            break;

        case 3:
            SubNumPhysNodes = Mesh.get_SubNumPhysNodes();
            SubNumBndNodes  = Mesh.get_SubNumBndNodes();
            if(SubFactor!=3)
            {
                printf("\n");
                printf(" Error: for space_order = %i, need SubFactor = %i\n",space_order,3);
                printf("      SubFactor = %i\n",SubFactor);
                printf("\n");
                exit(1);
            }
            break;

        default:
            printf("\n");
            printf(" ERROR in RunDogpack_unst.cpp: space_order value not supported.\n");
            printf("       space_order = %i\n",space_order);
            printf("\n");
            exit(1);
    }

    // local storage:
    dTensor1 rhs(SubNumPhysNodes);
    dTensor1 phi(SubNumPhysNodes);

    // Get Cholesky factorization matrix R
    //
    // TODO - this should be saved earlier in the code rather than reading
    // from file every time we with to run a Poisson solve!
    //
    SparseCholesky R(SubNumPhysNodes);
    string outputdir = dogParams.get_outputdir();
    R.init(outputdir);
    R.read(outputdir);

    // Create right-hand side vector
    void Rhs2D_unst(const int space_order,
            const mesh& Mesh, const dTensor3& rhs_dg,
            dTensor1& rhs);
    Rhs2D_unst(space_order, Mesh, rho, rhs);

    // Call Poisson solver  
    void PoissonSolver2D_unst(const int space_order,
            const mesh& Mesh,
            const SparseCholesky& R,
            const dTensor1& rhs,
            dTensor1& phi,
            dTensor2& E1,
            dTensor2& E2);
    PoissonSolver2D_unst(space_order, Mesh, R, rhs, phi, E1, E2);

}
