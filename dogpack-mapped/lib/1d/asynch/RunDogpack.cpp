#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "IniDocument.h"
#include "DogParamsCart1.h"
#include "DogAsynch1d.h"
using namespace std;

int RunDogpack(string outputdir)
{
    // ------------------------------------------------------------
    // Function definitions
    void L2Project(int mopt, int istart, int iend, const DogAsynch1d& dogAsynch1d, dTensorBC3* Fout,
            void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));
    void Output(const DogAsynch1d& dogAsynch1d, double t, int nframe, string outputdir);
    void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
            const dTensor2& NOT_USED_2, dTensor2& qvals);
    void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
            const dTensor2& NOT_USED_2, dTensor2& auxvals);
    void ConSoln(const DogAsynch1d& dogAsynch1d, double t, string outputdir);
    void DogSolveLxW_synch(DogAsynch1d& dogAsynch1d, double tstart, double tend, string outputdir);
    void DogSolveLxW_asynch(DogAsynch1d& dogAsynch1d, double tstart, double tend, string outputdir);
    void InitApp(IniDocument& ini_doc);
    // ------------------------------------------------------------

    // Output title information
    cout << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << "   | DoGPack: The Discontinuous Galerkin Package  |   " << endl;
    cout << "   | Developed by the research group of           |   " << endl;
    cout << "   |            James A. Rossmanith               |   " << endl;
    cout << "   |            Department of Mathematics         |   " << endl;
    cout << "   |            University of Wisconsin - Madison |   " << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << endl;

    cout << "Asynch examples currently do not compile. " << endl;
    cout << "This was discovered 2/11/2014.  -DS "       << endl;
    cout << endl << "See $DOGPACK/lib/1d/async/RunDogpack" << endl;

    // Some needed local variables
//  int n;
//  double dtout,tstart,tend;

//  // Get parameters
//  dogParams.init(ini_doc);
//  dogParamsCart1.init(ini_doc);
//  cout << endl;

//  // Get addtional parameters
//  InitApp(ini_doc);

//  // Create dogAsynch1d and store in dogAsynch1d object
//  int mxin;
//  if (dogParamsCart1.get_read_grid())
//  {
//      ifstream meshfile("mesh.dat", ios::in);
//      meshfile >> mxin;
//      meshfile.close();
//  }
//  else
//  {
//      cout << endl;
//      mxin = dogParamsCart1.get_mx();
//  }

//  DogAsynch1d dogAsynch1d(mxin,dogParams.get_meqn(),dogParams.get_maux(),
//          dogParams.get_kmax(),dogParamsCart1.get_mbc());

//  if (dogParamsCart1.get_read_grid())
//  {    
//      ifstream meshfile("mesh.dat", ios::in);
//      meshfile >> mxin;
//      for (int i=1; i<=(mxin+1); i++)
//      {
//          double xtmp;
//          meshfile >> xtmp;
//          dogAsynch1d.node->set(i, xtmp );
//      }
//      meshfile.close();
//      dogParamsCart1.set_mx(mxin);
//      dogParamsCart1.set_xlims(dogAsynch1d.node->get(1),dogAsynch1d.node->get(mxin+1));
//      dogParamsCart1.reportParameters();
//      cout << endl;
//  }
//  else
//  {    
//      dogAsynch1d.InputMesh(dogParamsCart1.get_xlow(),dogParamsCart1.get_xhigh());    
//  }
//  dogAsynch1d.dx_compute();
//  dogAsynch1d.set_CFL_MAX(dogParams.get_cflv()[1]);
//  dogAsynch1d.set_CFL(dogParams.get_cflv()[2]);
//  dogAsynch1d.OutputMesh(outputdir);

//  // Some basic parameters    
//  const string time_stepping_method = dogParams.get_time_stepping_method();
//  const int&     nout     = dogParams.get_nout();
//  const double&  tfinal   = dogParams.get_tfinal();
//  double dtv[2+1];
//  dtv[1] = dogParams.get_initial_dt();
//  dtv[2] = dogParams.get_max_dt();
//  const double*  cflv     = dogParams.get_cflv();
//  int      nv       = dogParams.get_nv();
//  const int*     method   = dogParams.get_method();
//  const int&     meqn     = dogParams.get_meqn();
//  const int&     mdim     = dogParams.get_ndims();
//  const int&     mrestart = dogParams.get_mrestart();
//  int      nstart   = dogParams.get_nstart();    
//  const int&     mx       = dogParamsCart1.get_mx();
//  const int&     mbc      = dogParamsCart1.get_mbc();
//  const double&  xlow     = dogParamsCart1.get_xlow();
//  const double&  xhigh    = dogParamsCart1.get_xhigh();
//  const double&  dx       = dogParamsCart1.get_dx();
//  const int mnodes = mx + 1;
//  const int kmax = method[1];
//  const int maux = method[6];

//  // Output meqn and nout for plotting purposes
//  string qhelp;
//  qhelp=outputdir+"/qhelp.dat";
//  ofstream out_file(qhelp.c_str(), ios::out);

//  out_file << setprecision(16);
//  out_file << nout << endl << meqn << endl << maux << endl;
//  out_file << kmax << endl << mx << endl;
//  out_file << setw(24) << scientific << xlow << endl << setw(24)
//      << scientific << xhigh << endl << setw(24) 
//      << scientific << dx << endl;
//  out_file.close();

//  // Set any auxiliary variables on computational grid
//  // Set values and apply L2-projection
//  if (maux>0)
//  {  L2Project(0,1-mbc,mx+mbc,dogAsynch1d,dogAsynch1d.aux,&AuxFunc);  }

//  // Set initial data on computational grid
//  // Set values and apply L2-projection
//  L2Project(0,1-mbc,mx+mbc,dogAsynch1d,dogAsynch1d.q,&QinitFunc);

//  // Output initial data to file
//  // For each element, we output ``method[1]'' number of values
//  dogAsynch1d.OutputSoln(outputdir,0.0,0);

//  // Compute conservation and print to file
//  ConSoln(dogAsynch1d,0.0,outputdir);

//  // Main loop for time stepping
//  tstart = 0.0;
//  tend   = 0.0;
//  dtout = tfinal/double(nout);    
//  for (n=1; n<=nout; n++)
//  {        
//      tstart = tend;      
//      tend = tstart + dtout;

//      if (dtout>0.0e0)
//      {
//          // Solve hyperbolic system from tstart to tend
//          if (time_stepping_method == "Runge-Kutta")
//          {  
//              exit(1);
//          }
//          else if (time_stepping_method == "SDC")
//          {
//              exit(1);
//          }
//          else if (time_stepping_method == "LxW_synch")
//          {
//              // Lax-Wendroff time-stepping scheme          
//              DogSolveLxW_synch(dogAsynch1d,tstart,tend,outputdir);
//          }
//          else if (time_stepping_method == "LxW_asynch")
//          {
//              // Lax-Wendroff time-stepping scheme          
//              DogSolveLxW_asynch(dogAsynch1d,tstart,tend,outputdir);
//          }
//          else if (time_stepping_method == "User-Defined")
//          {
//              exit(1);
//          }
//      }

//      // Output data to file
//      dogAsynch1d.OutputSoln(outputdir,tend,n);

//      // Compute conservation and print to file
//      ConSoln(dogAsynch1d,tend,outputdir);

//      // Done with solution from tstart to tend
//      cout << setprecision(5);
//      cout << "DOGPACK: Frame " << setw(3) << n;
//      cout << ": plot files done at time t =";
//      cout << setw(12) << scientific << tend << endl;
//      cout << endl;
//  }

    return 0;
}
