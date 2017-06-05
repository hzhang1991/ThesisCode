#include "RunDogpack.h"

/*
 * Top level function to RunDogpack.  Briefly, this function calls the
 * following functions in the following order:
 *
 * 1.)  Iniitalize global structs dogParams and dogParamsCart1
 *
 * 2.) Call InitApp.  (additional application specific parameters)
 *
 * 3.) Write qhelp.dat to the output directory.  This is a total of two
 * functions calls: one on dogParams.write_qhelp and one on
 * dogParamsCart1.write_qhelp.
 *
 * 4.) Call GridSetup.  For a 1D grid, this sets up a single array called
 * node containing cell edges.
 * 
 * 5.) Call L2Project.  This Projects initial conditions onto basis functions.
 *
 * 6.) Call AfterQinit.  This is called once per simulation and can be used to
 * set up extra variables.
 *
 * 7.) Call Output - output initial conditions to the output directory.
 *
 * 8.) Call ConSoln - this is a call for saving 'conserved' quantities.  This
 * function is called once per time step.
 *
 * 9.) Run the main time stepping loop.  This consists of calling the
 * following two functions, once for each frame the user requested:
 *
 *     a.) Call DogSolve[TS-method], where TS-method is a valid time-stepping
 *     option.  (e.g. DogSolveRK, DogSolveSDC, DogSolveLxW, DogSolveUser).
 *
 *     b.) Call Output to print data to file
 *
 */
int RunDogpack(string outputdir)
{

    // Output title information
    cout << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << "   | DoGPack: The Discontinuous Galerkin Package  |   " << endl;
    cout << "   | Developed by the research group of           |   " << endl;
    cout << "   |            James A. Rossmanith               |   " << endl;
    cout << "   |            Department of Mathematics         |   " << endl;
    cout << "   |            Iowa State University             |   " << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << endl;

    // Get parameters
    dogParams.init();
    dogParamsCart1.init(ini_doc);
    cout << endl;

    // Get addtional parameters
    InitApp(ini_doc);
    cout << endl;

    // If we want to use the top-level solver, this routine needs to be written:
    //fetch_dogState().init();

    const string time_stepping_method = dogParams.get_time_stepping_method();
    const int&     nout     = dogParams.get_nout();
    const double&  tfinal   = dogParams.get_tfinal();
    double dtv[2+1];
    dtv[1] = dogParams.get_initial_dt();
    dtv[2] = dogParams.get_max_dt();
    const double*  cflv     = dogParams.get_cflv();
    const int      nv       = dogParams.get_nv();
    const int*     method   = dogParams.get_method();
    const int&     meqn     = dogParams.get_meqn();
    const int&     mdim     = dogParams.get_ndims();
    const int&     melems   = dogParamsCart1.get_melems();
    const int&     mbc      = dogParamsCart1.get_mbc();
    const double&  xlow     = dogParamsCart1.get_xlow();
    const double&  xhigh    = dogParamsCart1.get_xhigh();
    const double&  dx       = dogParamsCart1.get_dx();
    const int&     mrestart = dogParams.get_mrestart();
    const int mnodes = melems + 1;

    // Output meqn and nout for plotting purposes
    string qhelp;
    qhelp=outputdir+"/qhelp.dat";
    dogParams.write_qhelp(qhelp.c_str());
    dogParamsCart1.append_qhelp(qhelp.c_str());

    // Dimension arrays
    dTensor2      node(mnodes,mdim);
    dTensor1      prim_vol(melems);
    dTensorBC3    qnew(melems,meqn,method[1],mbc);
    dTensorBC3    qold(melems,meqn,method[1],mbc);
    dTensorBC1    smax(melems,mbc);
    dTensorBC3    aux(melems,iMax(method[6],1),method[1],mbc);

    // Construct 1D grid
    GridSetup(method[1],xlow,dx,node,prim_vol);

    // Set any auxiliary variables on computational grid
    // Set values and apply L2-projection
    if (method[6]>0)
    {  L2Project(0,1-mbc,melems+mbc,node,qnew,aux,aux,&AuxFunc);  }

    // Set initial data on computational grid
    // Set values and apply L2-projection
    L2Project(0,1-mbc,melems+mbc,node,qnew,aux,qnew,&QinitFunc);

    // Run AfterQinit to set any necessary variables
    AfterQinit(node,aux,qnew);

    // Output initial data to file
    // For each element, we output ``method[1]'' number of values
    Output(node,aux,qnew,0.0,0,outputdir);

    // Compute conservation and print to file
    ConSoln(method,node,aux,qnew,0.0,outputdir);

    // Main loop for time stepping
    double tstart = 0.0;
    double tend   = 0.0;
    double dtout = tfinal/double(nout);    
    for (int n=1; n<=nout; n++)
    {        
        tstart = tend;      
        tend = tstart + dtout;

        // Solve hyperbolic system from tstart to tend
        if (time_stepping_method == "Runge-Kutta")
        {  
            // Runge-Kutta time-stepping scheme
            DogSolveRK(node,prim_vol,aux,qold,qnew,smax,tstart,tend,
                    nv,method,dtv,cflv,outputdir);
        }
        else if (time_stepping_method == "SDC")
        {
            // Spectral deferred correction (SDC) time-stepping scheme
            DogSolveSDC(node,prim_vol,aux,qold,qnew,smax,tstart,tend,
                    nv,method,dtv,cflv,outputdir);
        }
        else if (time_stepping_method == "Lax-Wendroff")
        {
            // Lax-Wendroff time-stepping scheme
            DogSolveLxW(node,prim_vol,aux,qold,qnew,smax,tstart,tend,
                    nv,method,dtv,cflv,outputdir);
        }
        else if (time_stepping_method == "User-Defined")
        {
            // User-defined time-stepping scheme
            DogSolveUser(node,prim_vol,aux,qold,qnew,smax,tstart,tend,
                    nv,method,dtv,cflv,outputdir);
        }

        // Output data to file
        Output(node,aux,qnew,tend,n,outputdir);

        // Done with solution from tstart to tend
        cout << setprecision(5);
        cout << "DOGPACK: Frame " << setw(3) << n;
        cout << ": plot files done at time t =";
        cout << setw(12) << scientific << tend << endl;
        cout << endl;
    }

    return 0;
}

void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals)
{
    void QinitFunc(const dTensor1& xpts, dTensor2& qvals);
    QinitFunc(xpts,qvals);
}

void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals)
{
    void AuxFunc(const dTensor1& xpts, dTensor2& auxvals);
    AuxFunc(xpts,auxvals);
}
