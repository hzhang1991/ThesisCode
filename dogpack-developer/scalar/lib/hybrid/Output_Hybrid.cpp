#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "MonomialsToLegendre.h"
#include "mesh.h"
#include <sstream>
#include <string>
#include "assert.h"

void ComputeDensity( const mesh& Mesh, const dTensorBC5& q, dTensor3& rho )
{

    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);
    const int     meqn = q.getsize(4);
    const int     kmax = q.getsize(5);

    const int kmax2d   = rho.getsize(3);

    const double dV = dogParamsCart2.get_prim_vol();

    // CG-solver only works with a single equation here:
//  assert_eq( rho.getsize(1), NumElems   );

    // mapping from the moments  that survive the integral, to those that DIE!
    int tranlate_k[] = {1,2,3,6,7,8};
    if( kmax2d == 6 || kmax2d == 3 || kmax2d == 1)
    { 
        const int me = 1;
#pragma omp parallel for
        for( int n=1; n <= NumElems; n++ )
        {
            for( int k2=1; k2 <= kmax2d; k2++ )
            {
 
                double tmp = 0.;
                for( int i=1; i <= mx; i++ )
                for( int j=1; j <= my; j++ )
                {
                    tmp += q.get(i,j,n,me, tranlate_k[k2-1] );
                }
                rho.set(n,me,k2, dV * tmp );
            }
        }
    }
    else
    {
        unsupported_value_error( kmax2d );
    }
}

/*
void ComputeDensity( const mesh& Mesh, const dTensorBC5& q, dTensor3& rho )
{

    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);   assert_eq(NumElems, rho.getsize(1) );
    const int     meqn = q.getsize(4);   assert_eq(meqn, 1 );
    const int     kmax = q.getsize(5);
    const int kmax2d   = rho.getsize(3);

    const double dV     = dogParamsCart2.get_prim_vol();

    const int num_moments = rho.getsize(2);
    assert_eq( num_moments, 1 );

    // Integrate all of the moments
    int translate_k_rho []  = { 1,2,3,6,7,8};

#pragma omp parallel for
    for( int n=1; n <= NumElems; n++ )
    {

        const int me = 1;

        // Density for triangle n
        dTensor1 irho(kmax2d);  irho.setall(0.);
        for( int k2=1; k2 <= kmax2d; k2++ )
        {
            for( int i=1; i <= mx; i++ )
            for( int j=1; j <= my; j++ )
            {
                irho.set( k2, irho.get(k2) + q.get(i,j,n,me,translate_k_rho[k2-1] ) );
            }
            rho.set( n, me, k2, dV*irho.get(k2) );
        }
    }

}
*/

void ComputeMoments( const mesh& Mesh, const dTensorBC5& q, dTensor3& q_moments)
{

    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);   assert_eq(NumElems, q_moments.getsize(1) );
    const int     meqn = q.getsize(4);   assert_eq(meqn, 1 );
    const int     kmax = q.getsize(5);
    const int kmax2d   = q_moments.getsize(3);

    const double v1_low = dogParamsCart2.get_xlow();
    const double v2_low = dogParamsCart2.get_ylow();

    const double dv1   = dogParamsCart2.get_dx();
    const double dv2   = dogParamsCart2.get_dy();

    const double dv1_sq = dv1*dv1;
    const double dv2_sq = dv2*dv2;
    const double dv12   = dv1*dv2;
    
    const double dV     = dogParamsCart2.get_prim_vol();
    const double one_third = 1./3.;

    const int num_moments = q_moments.getsize(2);
    q_moments.setall( 0. );

    // Poisson solver only requires the first moment
    if( num_moments == 1 )
    {
        ComputeDensity( Mesh, q, q_moments );
        return;
    }

    // Integrate all of the moments
//  int translate_k_rho []  = { 1,2,3,6,7,8};
//  int translate_k_mu  []  = { 4, 9,11    };
//  int translate_k_tau []  = { 5,10,12    };

#pragma omp parallel for
    for( int n=1; n <= NumElems; n++ )
    {

        const int me = 1;
//      double rho     = 0.;   // particle density
//      double mu      = 0.;   // \iint     mu * phi^{(k)}
//      double tau     = 0.;   // \iint    eta * phi^{(k)}
//      double mu2     = 0.;   // \iint   mu^2 * phi^{(k)}
//      double tau2    = 0.;   // \iint  tau^2 * phi^{(k)}
//      double mutau   = 0.;   // \iint mu*tau * phi^{(k)}

        // Locally integrated quantities:
        dTensor1 irho(kmax2d), imu(kmax2d), itau(kmax2d);
        dTensor1 imu2(kmax2d), itau2(kmax2d), imutau(kmax2d);

        for( int i=1; i <= mx; i++ )
        for( int j=1; j <= my; j++ )
        {

            // Reset the local quantities (for current cell)
            irho.setall  (0.);
            imu.setall   (0.);
            itau.setall  (0.);
            imu2.setall  (0.);
            itau2.setall (0.);
            imutau.setall(0.);

            // Constant part of the integral
            irho.set   (1, q.get(i,j,n,me,1) );

            switch( dogParams.get_space_order() )
            {
                case 3:
                irho.set(6, q.get(i,j,n,me,8) );
                irho.set(5, q.get(i,j,n,me,7) );
                irho.set(4, q.get(i,j,n,me,6) );
                irho.set(3, q.get(i,j,n,me,3) );
                irho.set(2, q.get(i,j,n,me,2) );

                imu.set  (3, q.get(i,j,n,me,11) / sq3 );
                imu.set  (2, q.get(i,j,n,me, 9) / sq3 );
                imu.set  (1, q.get(i,j,n,me, 4) / sq3 );

                itau.set (3, q.get(i,j,n,me,12) / sq3 );
                itau.set (2, q.get(i,j,n,me,10) / sq3 );
                itau.set (1, q.get(i,j,n,me, 5) / sq3 );

                imu2.set (6, one_third*irho.get(6));
                imu2.set (5, one_third*irho.get(5));
                imu2.set (4, one_third*irho.get(4));
                imu2.set (3, one_third*irho.get(3));
                imu2.set (2, one_third*irho.get(2));
                imu2.set (1, (4./15.)*q.get(i,j,n,me,14) + one_third*irho.get(1));

                itau2.set (6, one_third*irho.get(6));
                itau2.set (5, one_third*irho.get(5));
                itau2.set (4, one_third*irho.get(4));
                itau2.set (3, one_third*irho.get(3));
                itau2.set (2, one_third*irho.get(2));
                itau2.set (1, (4./15.)*q.get(i,j,n,me,15) + one_third*irho.get(1));

                imutau.set( 1, one_third*q.get(i,j,n,me,13) );
                break;

                case 2:
                irho.set(3, q.get(i,j,n,me,3) );
                irho.set(2, q.get(i,j,n,me,2) );

                imu.set  (1, q.get(i,j,n,me, 4) / sq3 );
                itau.set (1, q.get(i,j,n,me, 5) / sq3 );

                imu2.set (3, one_third*irho.get(3));
                imu2.set (2, one_third*irho.get(2));
                imu2.set (1, one_third*irho.get(1));

                itau2.set (3, one_third*irho.get(3));
                itau2.set (2, one_third*irho.get(2));
                itau2.set (1, one_third*irho.get(1));

            }

            // local coordinates (needed to know what to add in)
            const double v1c = v1_low + (i-0.5)*dv1;
            const double v2c = v2_low + (j-0.5)*dv2;

            for( int k2=1; k2 <= kmax2d; k2++ )
            {

                // Density
                double tmp = q_moments.get(n,1,k2) + dV*irho.get(k2);
                q_moments.set(n,1,k2, tmp );

                // 1-momentum  \iint_{v1,v2} v1 * f
                // v1    = dv1*mu/2 + v1c
                tmp = q_moments.get(n,2,k2) + dV*( v1c*irho.get(k2) + 0.5*dv1*imu.get(k2) );
                q_moments.set(n,2,k2, tmp );

                // 2-momentum  \iint_{v1,v2} v2 * f 
                // v2    = dv2*tau/2 + v2c
                tmp = q_moments.get(n,3,k2) + dV*( v2c*irho.get(k2) + 0.5*dv2*itau.get(k2) );
                q_moments.set(n,3,k2, tmp );

                // "Energy"

                // (1,1)-component of "energy"
                // v1**2 = dv1**2*mu**2/4 + dv1*mu*v1c + v1c**2
                tmp = q_moments.get(n,4,k2) + dV*( 0.25*dv1_sq*imu2.get(k2) + dv1*v1c*imu.get(k2) + v1c*v1c*irho.get(k2) );
                q_moments.set(n,4,k2, tmp );

                // (2,2)-component of "energy"
                // v2**2 = dv2**2*tau**2/4 + dv2*tau*v2c + v2c**2
                tmp = q_moments.get(n,5,k2) + dV*( 0.25*dv2_sq*itau2.get(k2) + dv2*v2c*itau.get(k2) + v2c*v2c*irho.get(k2) );
                q_moments.set(n,5,k2, tmp );

                // (1,2)-component of "energy"
                // v1*v2 = dv1*dv2*mu*tau/4 + dv1*mu*v2c/2 + dv2*tau*v1c/2 + v1c*v2c
                tmp = q_moments.get(n,6,k2) + dV*( 0.25*dv12*imutau.get(k2) + 0.5*( dv1*v2c*imu.get(k2) + dv2*v1c*itau.get(k2) ) + v1c*v2c*irho.get(k2) );
                q_moments.set(n,6,k2, tmp );

            }


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
 
    // TODO - this is hard-coded for the Vlasov-code
    const int NumMoments = 6;  

    int kmax2d_vec[] = {1,3,6};
    const int kmax2d = kmax2d_vec[ dogParams.get_space_order()-1 ];

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

        fprintf(file,"%d\n", meqn);                     // == 1 for VP
        fprintf(file,"%d\n", 0 );                       // maux
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
