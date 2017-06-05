#include "dogdefs.h"
#include "IniDocument.h"
#include "DogParams.h"
#include "DogParamsUnst2.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include "DogStateUnst2.h"
#include "ext_time.h"               /* for get_utime and timeval_diff */
#include "RunDogpack_Unst.h"

int RunDogpack_Unst(string outputdir)
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

    // allocate state variables
    dogStateUnst2.init();
    dogStateUnst2.set_initial_dt(dogParams.get_initial_dt());

    // Check to see if a mesh has been generated,
    //   if YES, then read in basic mesh parameters
    ifstream mesh_file("Unstructured_Mesh/mesh_output/mesh_params.dat", ios::in);
    char buffer[256];
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

    // Read-in in basic mesh parameters
    int NumElems,NumPhysElems,NumGhostElems,NumNodes;
    int NumPhysNodes,NumBndNodes,NumEdges,NumBndEdges;
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

    // Get application parameters    
    InitApp_Unst(ini_doc,Mesh);

    // Dimension arrays
    const int meqn = dogParams.get_meqn();
    const int kmax = dogParams.get_kmax();
    const int space_order = dogParams.get_space_order();
    dTensor3 qnew(NumElems,meqn,kmax);
    dTensor3 qold(NumElems,meqn,kmax);
    dTensor1 smax(NumEdges);
    const int maux = dogParams.get_maux();
    dTensor3  aux(NumElems,maux,kmax);

    // initialize state of solver
    int nstart = 0;
    double tstart = 0.;
    if (dogParams.get_mrestart())
    {
        // -------------------------- 
        // Do a restart from old data
        // --------------------------
        nstart = dogParams.get_nrestart();

        // Set any auxiliary variables on computational grid
        if (maux>0)
        {  SetAux_Unst_restart(nstart,aux,outputdir);  }    

        // Set initial data on computational grid
        tstart = Qinit_Unst_restart(nstart,qnew,outputdir);     
    }
    else
    {
        // -------------------------- 
        // Start new computation
        // -------------------------- 

        // Set any auxiliary variables on computational grid
        // Set values using L2-projection
        if (maux>0)
        {  
            L2Project_Unst(1,NumElems,
                    space_order,space_order,space_order,space_order,		       
                    Mesh,&qnew,&aux,&aux,&AuxFuncWrapper);  
        }    

        // Set initial data on computational grid
        // Set values using L2-projection
        L2Project_Unst(1,NumElems,
                space_order,space_order,space_order,space_order,
                Mesh,&qnew,&aux,&qnew,&QinitFuncWrapper);

        // Apply post processing to initial data
        AfterQinit_Unst(Mesh,aux,qnew);

        // Output initial data to file
        dogStateUnst2.set_time(tstart);
        dogStateUnst2.set_dt(0.0);
        Output_Unst(Mesh,aux,qnew,tstart,nstart,outputdir);       
    }       

    // Compute conservation and print to file
    ConSoln_Unst(Mesh,aux,qnew,0.0,outputdir);

    // Set edge information needed later by Riemann solver
    edge_data_Unst EdgeData(NumEdges);
    SetEdgeData_Unst(Mesh,space_order,space_order,EdgeData);

    // Main loop for time stepping
    const int nout = dogParams.get_nout();
    const double tfinal = dogParams.get_tfinal(); 
    double tend   = tstart;
    const double dtout = (tfinal-tstart)/double(nout-nstart);
    for (int n=nstart+1; n<=nout; n++)
    {
        tstart = tend;	  
        tend = tstart + dtout;

        const string time_stepping_method = dogParams.get_time_stepping_method();
        // Solve hyperbolic system from tstart to tend
        if (time_stepping_method == "Runge-Kutta")
        {  
            // Runge-Kutta time-stepping  
            DogSolveRK_Unst(Mesh,EdgeData,aux,qold,qnew,
                    tstart,tend,outputdir);
        }
        else if (time_stepping_method == "SDC")
        {
            // Spectral deferred correction time-stepping
            printf("\n");
            printf(" SDC time-stepping has not yet been implemented ... \n");
            printf("\n");
            exit(1);
        }
        else if (time_stepping_method == "Lax-Wendroff")
        {
            // Lax-Wendroff time-stepping
            printf("\n");
            printf(" Lax-Wendroff time-stepping has not yet been implemented ... \n");
            printf("\n");
            exit(1);
        }
        else if (time_stepping_method == "User-Defined")
        {
            // User defined time-stepping scheme
            //	  DogSolveUser_Unst(Mesh,EdgeData,aux,qold,qnew,dogState,
            //                  tstart,tend,outputdir);
            printf("\n");
            printf(" DogSolveUser time-stepping has not yet been implemented ... \n");
            printf("\n");
            exit(1);
        }
        else
        {
            printf("\n");
            printf(" ERROR in RunDogpack_Unst.cpp: \n");
            printf("       Currently only Runge-Kutta time-stepping \n");
            printf("       has been implemented ... \n");
            printf("\n");
            exit(1);
        }

        // Output data to file
        Output_Unst(Mesh,aux,qnew,tend,n,outputdir); 

        // Done with solution from tstart to tend
        printf("DOGPACK: Frame %3d: at time t =%12.5e\n\n", n,tend);
    }

    // Get current time
    timeval end_time = get_utime();

    // Output elapsed time
    double diff_utime = timeval_diff(end_time, start_time);
    printf(" Total elapsed time in seconds = %11.5f\n\n", diff_utime);

    return 0;
}


void AuxFuncWrapper(const dTensor2& xpts,
        const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2,
        dTensor2& auxvals)
{
    void AuxFunc(const dTensor2& xpts, dTensor2& auxvals);
    AuxFunc(xpts,auxvals);
}

void QinitFuncWrapper(const dTensor2& xpts,
        const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2,
        dTensor2& qvals)
{
    void QinitFunc(const dTensor2& xpts, dTensor2& qvals);
    QinitFunc(xpts,qvals);
}
