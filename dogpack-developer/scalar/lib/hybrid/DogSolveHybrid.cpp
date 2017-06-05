#include "DogSolveHybrid.h"

// number of quadrature points used in each direction, and number of basis
// functions for each 2D problem
int num_quad_cart[] = {1,4,9,16,25};  // == M^2
int num_quad_unst[] = {1,3,6,12,16};  // slightly smaller than M^2
int kmax2d_vec[]    = {1,3,6,10,15};  // number of polynomials for each 2D problem

// Advance the solution q over time interval tstart to tend.
//
// All local information is allocated within this function.  The only part
// that gets shared are time values passed through DogStateHybrid.  This class
// should be modified to accept the state variable, q and aux in place of only
// containing time information as is currently the case.  (-DS)
void DogSolveHybrid( const mesh& Mesh, const edge_data_Unst& EdgeData, 
    dTensorBC5& q, const double tstart, const double tend )
{

    // Cartesian grid parameters:
    const int mx     = dogParamsCart2.get_mx();  assert_eq( mx, q.getsize(1) );
    const int my     = dogParamsCart2.get_my();  assert_eq( my, q.getsize(2) );

    // Number of triangles:
    const int NumElems = q.getsize(3);

    const int meqn = q.getsize(4);
    const int kmax = q.getsize(5);

// TODO - dogParams is used for the 'local' problems!
// We should add an extra section to get this parameter passed in correctly.
//  const double* cflv     = dogParams.get_cflv();
    double cflv[] = {0., 5., 50.5 };
    // double cflv[] = {0., 0.5, 0.05 };

    const int nv   = dogParams.get_nv();

    // define local variables
    int n_step = 0;
    double t   = tstart;
    double dt  = dogStateHybrid.get_initial_dt();

    const double CFL_max    = cflv[1];
    const double CFL_target = cflv[2];
    double cfl      = 0.0;
    double cfl_junk = 0.;
    double dtmin = dt;
    double dtmax = dt;

printf("DogSolveHybrid: Using des_cfl = %2.3f; ", CFL_target );
printf("this is currently hard coded in DogSolveHybrid.cpp\n");

    assert_eq( NumElems , Mesh.get_NumElems() ); // Number of total elements in mesh
    const int NumNodes = Mesh.get_NumNodes(); // Number of nodes in mesh
    const int NumEdges = Mesh.get_NumEdges(); // Number of edges in mesh 

    // we'll define the cfl number based on advection velocities:
    double vmax = Max( fabs(dogParamsCart2.get_xhigh()), fabs(dogParamsCart2.get_yhigh() ) );
    vmax = Max( vmax, fabs(dogParamsCart2.get_xlow()) );
    vmax = Max( vmax, fabs(dogParamsCart2.get_ylow()) );

    // Largest grid size in Cartesian space
    // double h    = Max( dogParamsCart2.get_dx(), dogParamsCart2.get_dy() );

    // TODO - does the mesh save its largest edge length?
    double max_vol = 0.;
    for( int n=1; n <= NumElems; n++ )
    { max_vol = Max( max_vol, 0.5*Mesh.get_area_prim(n) ); }
    double h = sqrt( max_vol );
   
    // Set up the quadrature rules that are used for each time step:
    QuadratureRules QuadFuncs;
    QuadFuncs.init( dogParams.get_space_order() );

    // initial time step ( based on [hard coded] desired CFL number )
    dt = Min( dogParams.get_max_dt(), CFL_target * h / vmax );

    // Main time stepping loop:
    while( t < tend )
    {

        // initialize time step
        n_step       = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            eprintf(" Error in DogSolveHybrid.cpp: "
                    " Exceeded allowed # of time steps \n"
                    "    n_step = %d\n"
                    "        nv = %d\n\n",
                    n_step,nv);
        }

        // TODO - choose uniform time steps in place of this? (-DS)

        // set current time
        double told = t;
        if (told+dt > tend)
        { dt = tend - told; }
        t = told + dt;

        // Current time and time step choice:
        dogStateHybrid.set_time(told);
        dogStateHybrid.set_dt(dt);

        // Take a full time step of size dt
        //
        // Hard-coded Strang splitting.  This way the time_order parameter
        // will actually refer to the RK time stepping time order
        //cfl_junk = DogSolveRK_Unst_Parallel( Mesh, EdgeData, q, told, told+0.5*dt );

        // -----------------------------------------------
        // Stage #1  	        
        // -----------------------------------------------
        cfl_junk = DogSolveRK_Unst_Parallel( QuadFuncs, Mesh, EdgeData, q, told, told+0.5*dt );

        // ------------------------------------------------
        // Stage #2
        // ------------------------------------------------
        DogSolveSL_Parallel( QuadFuncs, Mesh, q, told, t );

        // ------------------------------------------------
        // Stage #3
        // ------------------------------------------------
        cfl_junk = DogSolveRK_Unst_Parallel( QuadFuncs, Mesh, EdgeData, q, told+0.5*dt, t);

        // cfl number based on max velocity only:
        cfl = vmax * dt / h;

        // output time step information
        if (dogParams.get_verbosity()>0) 
        {
            printf("DogSolveHYBRID ... Step %5d"
                    "   CFL =%6.3f"
                    "   dt =%11.3e"
                    "   t =%11.3e\n",
                    n_step, cfl, dt, t);
        }

        // choose new time step
        if (cfl>0.0)
        {   
            dt    = Min(dogParams.get_max_dt(), dt*CFL_target/cfl);
            dtmin = Min(dt,dtmin);
            dtmax = Max(dt,dtmax);
        }
        else
        {
            dt = dogParams.get_max_dt();
        }

        // save the current time:
        dogStateHybrid.set_time(t);

    }

    // set initial time step for next call to DogSolveRK
    dogStateHybrid.set_initial_dt(dt);

}


// Function definitions:
double DogSolveRK_Unst_Parallel(
    const QuadratureRules& QuadFuncs,
    const mesh& Mesh, const edge_data_Unst& EdgeData,
    dTensorBC5& q,
    const double tstart, const double tend )
{

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    const int mx = dogParamsCart2.get_mx();    
    const int my = dogParamsCart2.get_my();
    const int NumElems  = Mesh.get_NumElems();
    assert_eq( mx,       q.getsize(1) );
    assert_eq( my,       q.getsize(2) );
    assert_eq( NumElems, q.getsize(3) );

    const int sorder = dogParams.get_space_order();

    double DogSolveRK_Unst_Quiet(
       const dTensor2* vel_vec,
       const mesh& Mesh, const edge_data_Unst& EdgeData,
       dTensor3& aux, dTensor3& qold, dTensor3& qnew, 
       const double tstart, const double tend, 
       DogStateUnst2& dogStateUnst2);

    const int kmax2d = kmax2d_vec     [sorder-1];
    const int meqn   = num_quad_cart  [sorder-1];

    const int maux = dogParams.get_maux();

    double cfl = 0.;

    // Get the quadrature points and weights for the cartesian slab:
    const dTensor2& spts_cart   = QuadFuncs.get_spts_cart();
    const int num_quad_pts_cart = spts_cart.getsize(1);

// This call is turned off because we're assuming that velocity space is much
// smaller compared to configuration (physical) space.
// 
// If you want to turn this back on, you should remove the pragma statements
// hidden inside lib/2d/unst/ConstructL_Unst.
// #pragma omp parallel for
    for( int i=1; i <= mx; i++ )
    {

        // time parameters used for this single slice:
        DogStateUnst2 dogStateUnst2;
        dogStateUnst2.init();
        dogStateUnst2.set_initial_dt(dogParams.get_initial_dt());

        // Storage for this thread:
        dTensor3 qold( NumElems, meqn, kmax2d );
        dTensor3 aux ( NumElems, maux, kmax2d );
        dTensor3 qnew( NumElems, meqn, kmax2d );

        dTensor2* vel_vec = new dTensor2(meqn,2);
        for( int j=1; j <= my; j++ )
        {

            // Evaluate the speeds, vx, and vy throughout the current cell at
            // each required quadrature point:
            const double xc = dogParamsCart2.get_xc(i);
            const double yc = dogParamsCart2.get_yc(j);
            for( int s=1; s <= num_quad_pts_cart; s++ )
            {  
                vel_vec->set(s,1, xc + 0.5*dx*spts_cart.get(s, 1 ) );
                vel_vec->set(s,2, yc + 0.5*dy*spts_cart.get(s, 2 ) );
            }

            // Copy data from q into qold:
            void ReadSlice( const QuadratureRules& QuadFuncs, const dTensorBC5& qin, const int i, const int j, dTensor3& qout);
            ReadSlice( QuadFuncs, q, i, j, qold );
            ReadSlice( QuadFuncs, q, i, j, qnew );

            // Write data into q:
            //
            // Warning: race condition if we try to select Max over all cfl
            // numbers ...  This is why it's now called 'junk_cfl'.
            double junk_cfl = DogSolveRK_Unst_Quiet(vel_vec, Mesh, EdgeData, aux, qold, qnew, tstart, tend, dogStateUnst2);
            void WriteSlice( const QuadratureRules& QuadFuncs, const dTensor3& qin, const int i, const int j, dTensorBC5& qout );
            WriteSlice( QuadFuncs, qnew, i, j, q );

        }
        delete vel_vec;

    }

    return 0.;

}

void DogSolveSL_Parallel( 
    const QuadratureRules& QuadFuncs,
    const mesh& Mesh, dTensorBC5& q, double tstart, double tend )
{

    const double dt     = tend-tstart;
    const int sorder    = dogParams.get_space_order();   

    // Cartesian grid parameters:
    const int mx = dogParamsCart2.get_mx();
    const int my = dogParamsCart2.get_my();
    const int mbc = q.getmbc();

    // unstructured grid parameters:
    const int NumElems = q.getsize(3);

    // number of 'equations' to be evolved per cell
    const int mpoints_unst = num_quad_unst [sorder-1];
    const int kmax2d       = kmax2d_vec    [sorder-1];

    // Function definitions:
    void StepAdvecCC(double dt, dTensorBC4& qold, const dTensor2& speeds, dTensorBC4& qnew );
    void ComputeElectricField( const double t, const mesh& Mesh, const dTensorBC5& q,
        dTensor2& E1, dTensor2& E2);
    void SetLocalSpeeds( const QuadratureRules& QuadFuncs, const int n, const dTensor2& E1, const dTensor2& E2, dTensor2& speeds );

    // Call a Poisson solve (if running the Vlasov codes):
    dTensor2 E1(NumElems, kmax2d );
    dTensor2 E2(NumElems, kmax2d );
    ComputeElectricField( 0.5*(tstart+tend), Mesh, q, E1, E2);

    // Print stuff to output (at the half time steps to avoid extra electric
    // field computations, if running the Vlasov codes):
    void PrintElectricField( 
        const double t, const mesh& Mesh, const dTensorBC5& q, 
        const dTensor2& E1, const dTensor2& E2 
        );
    PrintElectricField( 0.5*(tstart+tend), Mesh, q, E1, E2 );

#pragma omp parallel for
    for( int n=1; n <= NumElems; n++ )
    {

        // This is hot swappable if one wishes to set different speeds here:
        dTensor2 speeds ( mpoints_unst, 2 );
        SetLocalSpeeds  ( QuadFuncs, n, E1, E2, speeds );

        // Storage for this thread:
        dTensorBC4 qold( mx, my, mpoints_unst, kmax2d, mbc );
        dTensorBC4 qnew( mx, my, mpoints_unst, kmax2d, mbc );

        // Copy data from q into qold:
        void ReadSlice( const QuadratureRules& QuadFuncs, const dTensorBC5& qin, const int n, dTensorBC4& qout);
        ReadSlice( QuadFuncs, q, n, qold );

        // Single time step using Semi-Lagrangian scheme:
        StepAdvecCC(dt, qold, speeds, qnew );

        // Write data into q:
        //
        void WriteSlice( const QuadratureRules& QuadFuncs, const dTensorBC4& qin, const int n, dTensorBC5& qout );
        WriteSlice( QuadFuncs, qnew, n, q );

    }

}
