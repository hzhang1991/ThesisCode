#include <sstream>
#include <string>
#include "dog_math.h"           // for Max
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "mesh.h"
#include "assert.h"
#include "SparseCholesky.h"
#include "VlasovParams.h"
#include "ComputeElectricField.h"

// This routine has been borrowed from the library to include the applied
// electric field:
//
//     Ea = -omega0^2 {\bf x}.
//
// For this problem, we set omega0 = 1.
//

// Applied electrid field.  This is the only application that has one.  This
// electric field is of the form:
//
//    {\bf E} = -\omega0^2 {\bf x}
//
void AppliedElectricField(
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& E)
{

    const int numpts = xpts.getsize(1);
    const int meqn   =    E.getsize(2);

assert( meqn == 2 );

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        E.set(i, 1, 64.0*x );
        E.set(i, 2, 64.0*y );

    }

}

// This routine simply glues together many of the routines that are already
// written in the Poisson solver library
//
// phi( 1:SubNumPhysNodes    ) is a scalar quantity.  
//
// E1 ( 1:NumElems, 1:kmax2d ) is a vector quantity.
// E2 ( 1:NumElems, 1:kmax2d ) is a vector quantity.
//
// See also: ConvertEfieldOntoDGbasis
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
    const int kmax2d    = E2.getsize(2);
    const int NumBndNodes  = Mesh.get_NumBndNodes();
    const int NumPhysNodes = Mesh.get_NumPhysNodes();

    // Quick error check
    if( !Mesh.get_is_submesh() )
    {
        printf("ERROR: mesh needs to have subfactor set to %d\n", space_order);
        printf("Go to Unstructured mesh and remesh the problem\n");
        exit(-1);
    }
    const int SubFactor    = Mesh.get_SubFactor();

    assert_eq( NumElems, Mesh.get_NumElems() );

    // -- Step 1: Compute rho -- //
    dTensor3 rho(NumElems, 1, kmax2d );

    // Compute moments performs: \int f(x,v) dv.
    //
    // The Poisson solver needs to have background ions subtracted out, and
    // the correct form for the right hand side is actually,
    //
    //    rho = -\int f(x,v) dv + rho0,
    //
    // Where rho0 is the initial distribution.
    // 
    void ComputeMoments(const mesh& Mesh, const dTensorBC5& q, dTensor3& q_moments);
    ComputeMoments( Mesh, q, rho  );

    // Redefine rho as the correct rho:
    for( int n=1; n <= NumElems; n++ )
    {
        rho.set( n,1, 1, -rho.get(n,1,1) );
//      rho.set( n,1, 1, vlasovParams.rho0 - rho.get(n,1,1) );
        for( int k=2; k <= kmax2d; k++ )
        {
            rho.set( n, 1, k, -rho.get(n,1,k) );
        }
    }

    // -- Step 2: Figure out how large phi needs to be
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

    // NEW STUFF:
    //
    // Add in the extra applied electric field:
    dTensor3 Ea( NumElems, 2, kmax2d );
    
    // Used for the applied electric field
    L2Project_Unst( NULL, 1, NumElems, 
        space_order, space_order, space_order, space_order, 
        Mesh, &rho, &rho, &Ea, &AppliedElectricField );

    // Add in the result to the computed Electric field:
    for( int i=1; i <= NumElems; i++ )
    for( int k=1; k <= kmax2d; k++ )
    {
        E1.set( i, k, E1.get(i,k) + Ea.get(i,1,k) );
        E2.set( i, k, E2.get(i,k) + Ea.get(i,2,k) );
    }

// Compute an estimate of the largest electric field in the domain
//  double tmp = 0.;
//  for( int i=1; i <= NumPhysNodes; i++ )
//  {
//      tmp = Max( tmp, fabs( E1.get(i,1) ) );
//      tmp = Max( tmp, fabs( E2.get(i,1) ) );
//  }
//  printf("approximate max electric field = %2.5f\n", tmp );
//  exit(1);

}
