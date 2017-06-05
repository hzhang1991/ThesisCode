#include "RunDogpackHybrid.h"

int RunDogpackHybrid(string outputdir)
{

    // Get current time
    timeval start_time = get_utime();

    // Output title information
    printf("\n"
            "   ------------------------------------------------   \n"
            "   | DoGPack: The Discontinuous Galerkin Package  |   \n"
            "   | Developed by the research group of           |   \n"
            "   |            James A. Rossmanith               |   \n"
            "   |            Department of Mathematics         |   \n"
            "   |            Iowa State University             |   \n"
            "   ------------------------------------------------   \n\n");

    // Get parameters from parameters.ini
    dogParams.init();  

    // Time stepping information (for global solves)
    dogStateHybrid.init();
    dogStateHybrid.set_initial_dt(dogParams.get_initial_dt());

    // ----------------------------------------------- //
    // ------------- Unstructured Stuff -------------- //
    // ----------------------------------------------- //

    // Check to see if a mesh has been generated,
    //   if YES, then read in basic mesh parameters
    ifstream mesh_file("Unstructured_Mesh/mesh_output/mesh_params.dat", ios::in);
    if(mesh_file.is_open()!=1)
    {
        printf(" ERROR: file not found:"
                " 'Unstructured_Mesh/mesh_output/mesh_params.dat' \n"
                "   In order to run DoGPack in unstructured grid mode \n"
                "   you must first generate a mesh. You can do this by \n"
                "   following these steps: \n \n"
                "      (1) Type: $DOGPACK/scripts/create_unst2_dir \n"
                "      (2) Type: cd Unstructured_Mesh \n"
                "      (3) Modify the following files as desired:  \n"
                "            - input2D.data  \n"
                "            - SignedDistance.cpp  \n"
                "            - GridSpacing.cpp  \n"
                "            - MeshPreProcess.cpp  \n"
                "            - MeshPostProcess1.cpp  \n"
                "            - MeshPostProcess2.cpp  \n"
                "      (4) Run mesh generator by typing: mesh2d.exe  \n"
                "      (5) Visualize mesh using the 'plotmesh2.m' MATLAB script  \n \n");
        exit(1);
    }

    // Check to see that the correct inversion routines have been called from
    // Matlab.
//  ifstream matlab_file("matlab/R.dat", ios::in);
//  if(matlab_file.is_open()!=1)
//  {
//      printf( " ERROR: file not found:"
//              " 'matlab/R.dat' \n"
//              " Open Matlab and call the function CreateMatrix( Morder ) \n"
//              " Also, I need to write some more testing here! (-DS) \n" );
//      exit(1);
//  }
    // TODO - might as well read in the Sparse matrix routines here! (-DS)


    // Read-in in basic mesh parameters
    int NumElems,NumPhysElems,NumGhostElems,NumNodes;
    int NumPhysNodes,NumBndNodes,NumEdges,NumBndEdges;
    char buffer[256];
    mesh_file >> NumElems;
    mesh_file.getline(buffer,256);
    mesh_file >> NumPhysElems;
    mesh_file.getline(buffer,256);
    mesh_file >> NumGhostElems;
    mesh_file.getline(buffer,256);
    mesh_file >> NumNodes;
    mesh_file.getline(buffer,256);
    mesh_file >> NumPhysNodes;
    mesh_file.getline(buffer,256);
    mesh_file >> NumBndNodes;
    mesh_file.getline(buffer,256);
    mesh_file >> NumEdges;
    mesh_file.getline(buffer,256);
    mesh_file >> NumBndEdges;
    mesh_file.getline(buffer,256);
    mesh_file.close();

    // Initialize 2d unstructured parameters
    dogParamsUnst2.init(NumElems,
            NumPhysElems,
            NumGhostElems,
            NumNodes,
            NumPhysNodes,
            NumBndNodes,
            NumEdges,
            NumBndEdges,
            outputdir);  

    // Create and read-in entire unstructured mesh
    mesh Mesh(NumElems,NumPhysElems,NumNodes,NumPhysNodes,
            NumBndNodes,NumEdges,NumBndEdges);
    string mesh_dir = "Unstructured_Mesh/mesh_output";
    Mesh.InputMesh(mesh_dir);

    // create help file for plotting purposes
    string qhelp;
    qhelp=outputdir+"/qhelp.dat";
    dogParams.write_qhelp(qhelp.c_str());
    dogParamsUnst2.write_qhelp(qhelp.c_str());

    // Copy mesh into output directory
    RunMeshCopyScript(outputdir);

    // ----------------------------------------------- //
    // ------------- Structured Stuff ---------------- //
    // ----------------------------------------------- //
    dogParamsCart2.init();  // well, now wasn't that easy?
    const int mx  = dogParamsCart2.get_mx();
    const int my  = dogParamsCart2.get_my();
    const int mbc = dogParamsCart2.get_mbc();

    // Get application parameters    
    // This is currently a blank routine located in lib/2d/IniApp.cpp
    // InitApp_Unst(ini_doc,Mesh);

    // Dimension arrays
    const int meqn = dogParams.get_meqn();
    const int kmax = dogParams.get_kmax();

    const int space_order = dogParams.get_space_order();

    // initialize state of solver
    int    nstart = 0;
    double tstart = 0.;

    // -------------------------- 
    // Start new computation
    // -------------------------- 

    // Set initial data on computational grid
    void L2Project( const mesh& Mesh, dTensorBC5& q);
    dTensorBC5 q(mx, my, NumElems, 1, kmax, mbc);
    L2Project(Mesh,q);

    // Apply post processing to initial data
    //
    // For the Vlasov Routines, this also reads in the vlasovParams section
    AfterQinit_Unst(Mesh,q);

printf("*** We also believe ethat NDIMS = %d ***\n", NDIMS );

    // Output initial data to file
    dogStateHybrid.set_time(tstart);
    Output_Hybrid(Mesh, q, 0., 0, outputdir); 

    // Compute conservation and print to file
    //  ConSoln_Unst(Mesh,aux,qnew,0.0,outputdir);

    // Set edge information needed later by Riemann solver
    edge_data_Unst EdgeData(NumEdges);
    SetEdgeData_Unst(Mesh,space_order,space_order,EdgeData);

    // Main loop for time stepping
    const int nout = dogParams.get_nout();
    const double tfinal = dogParams.get_tfinal(); 
    double tend = tstart;
    const double dtout = (tfinal-tstart)/double(nout-nstart);
    for (int n=nstart+1; n<=nout; n++)
    {
        tstart = tend;	  
        tend = tstart + dtout;

        // Solve hyperbolic system from tstart to tend and output to file
        DogSolveHybrid (Mesh,EdgeData,q,tstart,tend);
        Output_Hybrid  (Mesh, q, tend, n, outputdir); 

        printf("DOGPACK: Frame %3d: at time t =%12.5e\n\n", n,tend);
    }

    // Get current time
    timeval end_time = get_utime();

    // Output elapsed time
    double diff_utime = timeval_diff(end_time, start_time);
    printf(" Total elapsed time in seconds = %11.5f\n\n", diff_utime);

    return 0;

}
