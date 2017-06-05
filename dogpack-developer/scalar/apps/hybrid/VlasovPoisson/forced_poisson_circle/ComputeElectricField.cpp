#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "mesh.h"
#include "assert.h"
#include "VlasovParams.h"

void ComputeEfield(const int space_order, const mesh& Mesh, const dTensor1& phi,
    dTensor2& E1, dTensor2& E2);

// In this application, we'll set E2 = 0, in SetLocalSpeeds.
//
// Exact solution: f = A(t) * (1-x^2-y^2) * exp( -vx^2 ) / sqrt(pi)
//
// Where, A(t) = ( 0.75 + 0.25*cos(2*pi*t) )
//
//
void ExactElectricField(
    const double t, 
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& exact_e)
{

    const int numpts = exact_e.getsize(1);
    const double A   = ( 0.75 + 0.25*cos(2.0*pi*t) );
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double r2 = pow(x,2) + pow(y,2);
        double r  = sqrt( r2 );

        exact_e.set(i,1, A*( -0.5*x + x*x*x/3. ) );
        exact_e.set(i,2, A*( -0.5*y + y*y*y/3. ) );

    }

}

// This routine simply glues together many of the routines that are already
// written in the Poisson solver library
//
// phi( 1:SubNumPhysNodes       ) is a scalar quantity.  
//
// E1 ( 1:NumElems, 1:kmax_unst ) is a vector quantity.
// E2 ( 1:NumElems, 1:kmax_unst ) is a vector quantity.
//
// See also: ConvertEfieldOntoDGbasis
void ComputeElectricField( const double t, const mesh& Mesh, const dTensorBC5& q,
    dTensor2& E1, dTensor2& E2)
{

    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);
    const int     meqn = q.getsize(4);
    const int     kmax = q.getsize(5);

    const int sorder = dogParams.get_space_order();

    // unstructured parameters:
    const int kmax2d       = E1.getsize(2);
    const int NumBndNodes  = Mesh.get_NumBndNodes();
    const int NumPhysNodes = Mesh.get_NumPhysNodes();

    // use the exact Electric field here:
    void L2Project_Unst(
        const double time,
        const dTensor2* vel_vec,
        const int istart, 
        const int iend, 
        const int QuadOrder, 
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,
        const mesh& Mesh, 
        const dTensor3* qin, 
        const dTensor3* auxin, 
        dTensor3* fout, 
        void (*Func)(const double t, const dTensor2* vel_vec, const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor2&));

    dTensor3 qtmp   (NumElems, 2, kmax2d );  qtmp.setall(0.);
    dTensor3 auxtmp (NumElems, 0, kmax2d );
    dTensor3 ExactE (NumElems, 2, kmax2d );
    L2Project_Unst( t, NULL, 1, NumElems, 
        sorder, sorder, sorder, sorder, Mesh, 
        &qtmp, &auxtmp, &ExactE, 
        &ExactElectricField );

    // Copy contents into E1 and E2:
    //
    // For now, we'll set E2 = 0.
    for( int n=1; n <= NumElems; n++ )
    for( int k=1; k <= kmax2d; k++ )
    {
        E1.set(n, k, ExactE.get(n,1,k) );
        E2.set(n, k, ExactE.get(n,2,k) );
    }

}

//  // Quick error check
//  if( !Mesh.get_is_submesh() )
//  {
//      printf("ERROR: mesh needs to have subfactor set to %d\n", space_order);
//      printf("Go to Unstructured mesh and remesh the problem\n");
//      exit(-1);
//  }
//  const int SubFactor    = Mesh.get_SubFactor();

//  assert_eq( NumElems, Mesh.get_NumElems() );

//  // -- Step 1: Compute rho -- //
//  dTensor3 rho(NumElems, 1, kmax_unst );
//  void ComputeDensity( const mesh& Mesh, const dTensorBC5& q, dTensor3& rho );
//  ComputeDensity( Mesh, q, rho );

//  // -- Step 2: Figure out how large phi needs to be
//  int SubNumPhysNodes = 0;
//  int SubNumBndNodes  = 0;
//  switch( dogParams.get_space_order() )
//  {
//      case 1:
//          SubNumPhysNodes = NumPhysNodes;
//          SubNumBndNodes  = NumBndNodes;
//          break;

//      case 2:
//          SubNumPhysNodes = Mesh.get_SubNumPhysNodes();
//          SubNumBndNodes  = Mesh.get_SubNumBndNodes();
//          if(SubFactor!=2)
//          {
//              printf("\n");
//              printf(" Error: for space_order = %i, need SubFactor = %i\n",space_order,2);
//              printf("      SubFactor = %i\n",SubFactor);
//              printf("\n");
//              exit(1);
//          }
//          break;

//      case 3:
//          SubNumPhysNodes = Mesh.get_SubNumPhysNodes();
//          SubNumBndNodes  = Mesh.get_SubNumBndNodes();
//          if(SubFactor!=3)
//          {
//              printf("\n");
//              printf(" Error: for space_order = %i, need SubFactor = %i\n",space_order,3);
//              printf("      SubFactor = %i\n",SubFactor);
//              printf("\n");
//              exit(1);
//          }
//          break;

//      default:
//          printf("\n");
//          printf(" ERROR in RunDogpack_unst.cpp: space_order value not supported.\n");
//          printf("       space_order = %i\n",space_order);
//          printf("\n");
//          exit(1);
//  }

//  // local storage:
//  dTensor1 rhs(SubNumPhysNodes);
//  dTensor1 phi(SubNumPhysNodes);

//  // Get Cholesky factorization matrix R
//  //
//  // TODO - this should be saved earlier in the code rather than reading
//  // from file every time we with to run a Poisson solve!
//  //
//  SparseCholesky R(SubNumPhysNodes);
//  string outputdir = dogParams.get_outputdir();
//  R.init(outputdir);
//  R.read(outputdir);

//  // Create right-hand side vector
//  void Rhs2D_unst(const int space_order,
//          const mesh& Mesh, const dTensor3& rhs_dg,
//          dTensor1& rhs);
//  Rhs2D_unst(space_order, Mesh, rho, rhs);

//  // Call Poisson solver  
//  void PoissonSolver2D_unst(const int space_order,
//          const mesh& Mesh,
//          const SparseCholesky& R,
//          const dTensor1& rhs,
//          dTensor1& phi,
//          dTensor2& E1,
//          dTensor2& E2);
//  PoissonSolver2D_unst(space_order, Mesh, R, rhs, phi, E1, E2);

//}
