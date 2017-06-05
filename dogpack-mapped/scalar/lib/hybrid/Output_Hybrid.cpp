#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "MonomialsToLegendre.h"
#include "mesh.h"
#include <sstream>
#include <string>
#include "assert.h"


void ComputeMoments( const mesh& Mesh, const dTensorBC5& q, dTensor3& q_moments)
{

    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);   assert_eq(NumElems,q_moments.getsize(1) );
    const int     meqn = q.getsize(4);   assert_eq(meqn, 1 );
    const int     kmax = q.getsize(5);
    const int kmax2d   = q_moments.getsize(3);

    const double v1_low = dogParamsCart2.get_xlow();
    const double v2_low = dogParamsCart2.get_ylow();

    const double dv1   = dogParamsCart2.get_dx();
    const double dv2   = dogParamsCart2.get_dy();

    const int num_moments = q_moments.getsize(2);

    // Integrate all of the moments
    int translate_k_rho[] = {1,2,3,6,7,8};
    int translate_k_dv1[] = { 4, 9,11 };
    int translate_k_dv2[] = { 5,10,12 };
#pragma omp parallel for
    for( int n=1; n <= NumElems; n++ )
    {

        // Compute particle Density 
        // (and any term that gets multiplied by a constant)
        for( int k2=1; k2 <= kmax2d; k2++ )
        {
            double rho    = 0.;   // particle density
            double rho_u1 = 0.;   // 1-momentum  \iint_{v1,v2} v1 * f
            double rho_u2 = 0.;   // 2-momentum  \iint_{v1,v2} v2 * f 
            for( int i=1; i <= mx; i++ )
            for( int j=1; j <= my; j++ )
            {
                double v1c = v1_low + (i-0.5)*dv1;
                double v2c = v2_low + (j-0.5)*dv2;

                rho    += q.get( i, j, n, 1, translate_k_rho[k2-1] );
                rho_u1 += v1c*q.get(i,j,n,1, translate_k_rho[k2-1] );
                rho_u2 += v2c*q.get(i,j,n,1, translate_k_rho[k2-1] );

            }
            q_moments.set(n, 1, k2, dogParamsCart2.get_prim_vol()*rho    );
            q_moments.set(n, 2, k2, dogParamsCart2.get_prim_vol()*rho_u1 );
            q_moments.set(n, 3, k2, dogParamsCart2.get_prim_vol()*rho_u2 );
        }

        // Compute particle momentum
        for( int k2=1; k2 <= dogParams.get_space_order(); k2++ )
        {

            double rho_u1 = 0.;   // 1-momentum  \iint_{v1,v2} v1 * f
            double rho_u2 = 0.;   // 2-momentum  \iint_{v1,v2} v2 * f 
            for( int i=1; i <= mx; i++ )
            for( int j=1; j <= my; j++ )
            {
                rho_u1 += dv1*q.get(i,j,n,1, translate_k_dv1[k2-1] );
                rho_u2 += dv2*q.get(i,j,n,1, translate_k_dv2[k2-1] );

            }
            q_moments.set(n, 2, k2, q_moments.get(n,2,k2) + dogParamsCart2.get_prim_vol()*rho_u1 );
            q_moments.set(n, 3, k2, q_moments.get(n,3,k2) + dogParamsCart2.get_prim_vol()*rho_u2 );
        }
    }
}

// Output solution on unstructured mesh
void Output_Hybrid(const mesh& Mesh, const dTensorBC5& q, 
        double t, int nframe, string outputdir)
{ 

    const int       mx = q.getsize(1);   assert_eq(mx, dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my, dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);
    const int     meqn = q.getsize(4);
    const int     kmax = q.getsize(5);
  
    // TODO - of course these need to be modified!
    const int NumMoments = 3;  

    int kmax2d;
    switch( dogParams.get_space_order() )
    {

        case 3:
            kmax2d = 6;
            break;

        case 2:
            kmax2d = 3;
            break;

        case 1:
            kmax2d = 1;
            break;
        default:
            unsupported_value_error( dogParams.get_space_order() );
    }

    // -- Compute moments we desire -- //
    dTensor3 q_moments(NumElems, NumMoments, kmax2d );
    ComputeMoments( Mesh, q, q_moments );

    // -- Write moments to output files -- //
    ostringstream basename;
    basename << setfill('0') << setw(4) << nframe;
    string fname1  = outputdir+"/q2dUnst"+basename.str();

    // Output Q values
    void WriteOutput_Unst(string fname, string varname, 
            const dTensor3& q, double t);
    WriteOutput_Unst(fname1, "q2dUnst", q_moments, t);

    // ----------------------------------------------------- //
    //
    // Print a single cartesian slice to file
    //
    // ----------------------------------------------------- //

    string fname_cart  = outputdir+"/q2dCart"+basename.str()+".dat";
    FILE* file = fopen(fname_cart.c_str(), "w");

    // cartesian moments that survive (see QuadratureRules.cpp)
    //
    // TODO - if we truly wanted slices, and not integrated control volumes, we
    // should be sampling the basis functions in place of this simple mapping.
    //
    // If we were to sum over all of configuration space, then this would
    // produce the correct answer, but here we're trying to take a slice ...
    //
    int ktrans[] = {1,4,5,13,14,15};

    // search for an element that's near the "center":
//  int nslice = 1;
//  while( 1 )
//  {
//      if( fabs( Mesh.get_node(nslice,1) ) < 0.1 && fabs(
//      Mesh.get_node(nslice,2) < 0.1 ) )
//      {
//          break;
//      }
//      nslice += 1;
//  }
    const int nslice = Mesh.get_NumPhysElems()/2;

    fprintf(file,"%24.16e\n",t);
    for (int k=1; k<=kmax2d; k++ )
    for (int j=1; j<=q.getsize(2); j++)
    for (int i=1; i<=q.getsize(1); i++)
    {
        double tmp = q.get(i, j, nslice, 1, ktrans[k-1]);
        fprintf(file, "%24.16e\n", tmp);
    }
    fclose(file);

    // Cartesian helper file 
    // (has different format than unstructured helper file)
    if( nframe == 0 )
    {

        string fn = outputdir+"/qhelp_cart.dat";
        // create qhelp_cart.dat if it doesn't exist:
        FILE* file = fopen(fn.c_str(),"w");

        fprintf(file,"%d\n", meqn);
        fprintf(file,"%d\n", 0 ); // maux
        fprintf(file,"%d\n", dogParams.get_nout() );
        fprintf(file,"%d\n", dogParams.get_space_order());
        fprintf(file,"%d\n", dogParamsCart2.get_mx() );
        fprintf(file,"%d\n", dogParamsCart2.get_my() );
        fprintf(file,"%f\n", dogParamsCart2.get_xlow() );
        fprintf(file,"%f\n", dogParamsCart2.get_xhigh() );
        fprintf(file,"%f\n", dogParamsCart2.get_ylow() );
        fprintf(file,"%f\n", dogParamsCart2.get_yhigh() );
        fprintf(file,"1\n"); // data format (ASCII==1, HDF5==5)
        fclose(file);

    }

    // TODO - print a "qhelp_unst.dat" file as well!

    void Output_Hybrid_Extra(const mesh& Mesh, const dTensorBC5& q, 
        double t, int nframe, string outputdir);
    Output_Hybrid_Extra(Mesh, q, t, nframe, outputdir);

}
