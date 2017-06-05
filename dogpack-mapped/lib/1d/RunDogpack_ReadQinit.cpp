#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "IniDocument.h"
#include "DogParamsCart1.h"
using namespace std;

int RunDogpack(string outputdir)
{
  // ------------------------------------------------------------
  // Function definitions
  void GridSetup(int,double,double,dTensor2&,dTensor1&);
  void L2Project(int,int,int,const dTensor2&,const dTensorBC3&,
		 const dTensorBC3&,dTensorBC3&,
		 void (*Func)(const dTensor1&, const dTensor2&, 
			      const dTensor2&, dTensor2&));
  void Output(const dTensor2&,const dTensorBC3&,const dTensorBC3&,
	      double,int,string);
  void ReadQinit(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q);
  void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
	       const dTensor2& NOT_USED_2, dTensor2& auxvals);
  void AfterStep(double dt, const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
  void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
			 dTensorBC3& auxold, dTensorBC3& aux, 
			 dTensorBC3& qold, dTensorBC3& q);
  void ConSoln(const int method[], const dTensor2& node, const dTensorBC3& aux,
	       const dTensorBC3& q, double t, string outputdir);
  void DogSolveRK(const dTensor2&,const dTensor1&,
		  dTensorBC3&,dTensorBC3&,dTensorBC3&,
		  dTensorBC1&,double,double,int,const int[],
		  double[],const double[],string);
  void DogSolveSDC(const dTensor2&,const dTensor1&,
		   dTensorBC3&,dTensorBC3&,dTensorBC3&,
		   dTensorBC1&,double,double,int,const int[],
		   double[],const double[],string);
  void DogSolveLxW(const dTensor2&,const dTensor1&,
		   dTensorBC3&,dTensorBC3&,dTensorBC3&,
		   dTensorBC1&,double,double,int,const int[],
		   double[],const double[],string);
  void DogSolveUser(const dTensor2&,const dTensor1&,
		    dTensorBC3&,dTensorBC3&,dTensorBC3&,
		    dTensorBC1&,double,double,int,const int[],
		    double[],const double[],string);
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
  
  // Get parameters
  dogParams.init();
  dogParamsCart1.init(ini_doc);
  cout << endl;
  
  // Get addtional parameters
  InitApp(ini_doc);
  cout << endl;
  
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
  ofstream out_file(qhelp.c_str(), ios::out);
  
  out_file << setprecision(16);
  out_file << nout << endl << meqn << endl << method[6] << endl;
  out_file << method[1] << endl << melems << endl;
  out_file << setw(24) << scientific << xlow << endl << setw(24)
	   << scientific << xhigh << endl << setw(24) 
	   << scientific << dx << endl;
  out_file.close();
  
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
  // by reading it in from file
  ReadQinit(node,aux,qnew);
  
  // Run AfterStep to set any necessary variables
  AfterStep(0.0,node,aux,qnew);
  AfterFullTimeStep(0.0,node,prim_vol,aux,aux,qnew,qnew);
  
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
  
  return 1;
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
