#include "../defs.h"
#include "dog_math.h"
#include "GridPatch.h"

int RunDogpackHMM(int argc=0, char**argv=NULL)
{
    // ------------------------------------------------------------
    // Function definitions
    void ParseArguments(int argc,char**argv,string& outputdir);
    void GetParamsHMM(int&,int&,int[],double&,double[],double[],
		      int&,int&,double&,int&,int&,int&,int&,double&,
		      double&,int[],int&,int&,int&,double&,double&,string);
    void GridSetup(int,double,double,dTensor2&,dTensor1&,string,char[]);
    void SnapToMaster(GridPatch&,GridPatch&,iTensor1&,iTensor2&,iTensor2&);
    void L2Project(int,int,int,dTensor2,dTensorBC3,dTensorBC3,dTensorBC3&,
                   void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&));
    void Output(dTensor2,dTensorBC3,dTensorBC3,double,int,string,string);
    void ConSoln(dTensor2,int[],dTensorBC3,dTensorBC3,double,string);
    void DogSolveRK_HMM(double,double,int,double[],double[],double[],iTensor1,
			iTensor2,iTensor2,int[],dTensor2,dTensorBC3&,
			dTensorBC3&,dTensorBC1&,dTensor1,int[],dTensor2,
			dTensorBC3&,dTensorBC3&,dTensorBC1&,dTensor1,string);
    void DogSolveSDC(dTensor2,dTensor1,dTensorBC3&,dTensorBC3&,dTensorBC1&,
		     double,double,int,int[],double[],double[]);
    void QinitFunc(dTensor1,dTensor2,dTensor2,dTensor2&);
    void   AuxFunc(dTensor1,dTensor2,dTensor2,dTensor2&);
    void AfterStep(dTensor2,dTensorBC3&,dTensorBC3&);
    void QinitFunc_sub(dTensor1,dTensor2,dTensor2,dTensor2&);
    void   AuxFunc_sub(dTensor1,dTensor2,dTensor2,dTensor2&);
    void AfterStep_sub(dTensor2,dTensorBC3&,dTensorBC3&);
    // ------------------------------------------------------------

    // Output title information
    cout << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << "   | DoGPack: The Discontinuous Galerkin Package  |   " << endl;
    cout << "   | Written by James A. Rossmanith               |   " << endl;
    cout << "   |            Department of Mathematics         |   " << endl;
    cout << "   |            University of Wisconsin - Madison |   " << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << endl;

    // Some needed local variables
    int n,mx,mnodes;
    int mtmp1,mtmp2,mtmp3;
    int nout,nv,method[8],method_sub[8];
    double tfinal,dtv[3],dtv_sub[3],cflv[3];
    double dtout,tstart,tend,dx;
    int mbc,     meqn,     maux,     kmax;
    int mbc_sub, meqn_sub, maux_sub, kmax_sub;
    double xlow, xhigh, xlow_sub, xhigh_sub;
    string outputdir("output");

    // Parse arguments -- sets directory to which output will be sent,
    //                    the default is "output"
    ParseArguments(argc,argv,outputdir);
    
    // Get parameters from "dogpack.data"
    GetParamsHMM(nout,nv,method,tfinal,dtv,cflv,mx,mnodes,dx,
		 meqn,maux,kmax,mbc,xlow,xhigh,method_sub,
		 meqn_sub,maux_sub,kmax_sub,xlow_sub,xhigh_sub,
		 outputdir);
    dtv_sub[1] = dtv[1];
    dtv_sub[2] = dtv[2];

    // Figure out mbc_sub (number of ghost cells for SubGrid)
    if (method[2]<0)
      {  mbc_sub = mbc*abs(method[2]);  }
    else
      {  mbc_sub = mbc;  }
  
    // Output size parameters to screen
    cout << setprecision(4);
    cout << " MainGrid: " << endl;
    cout << "   Number of Equations:           " << meqn << endl;
    cout << "   Number of Spatial Dimensions:  " << 1 << endl;    
    cout << "   Order of Accuracy in Space:    " << method[1] << endl;
    cout << "   Order of Accuracy in Time:     " << method[2] << endl;
    if (method[3]==1)
    {  cout << "   Limiters:                      yes " << endl;  }
    else
    {  cout << "   Limiters:                      no  " << endl;  }
    cout << "   Number of Elements in Grid:    " << mx << endl;
    cout << "        Left Endpoint in Grid:    " << scientific <<  xlow << endl;
    cout << "       Right Endpoint in Grid:    " << scientific << xhigh << endl;
    cout << "        Number of Ghost Cells:    " << mbc << "  (on each boundary) " << endl;
    cout << endl;

    // Dimension GridPatch arrays
    GridPatch MainGrid(mx,xlow,xhigh,meqn,maux,kmax,mbc);
    GridPatch  SubGrid(1,xlow_sub,xhigh_sub,meqn_sub,maux_sub,kmax_sub,mbc_sub);
    iTensor1 map(1);
    iTensor2  leftBC(mbc_sub,2);
    iTensor2 rightBC(mbc_sub,2);
    
    // Snap SubGrid onto MainGrid
    SnapToMaster(MainGrid,SubGrid,map,leftBC,rightBC);    
				      
    // Arrays for the main grid
    dTensor2     node(mnodes,1);
    dTensor1 prim_vol(mx);
    dTensorBC3      q(mx,meqn,kmax,mbc);
    dTensorBC1   smax(mx,mbc);
    dTensorBC3    aux(mx,maux,kmax,mbc);    

    // Arrays for the sub grid
    int mx_sub       = SubGrid.getmx();
    int mnodes_sub   = mx_sub+1;
    double dx_sub    = SubGrid.getdx();
    xlow_sub  = SubGrid.getxlow();
    xhigh_sub = SubGrid.getxhigh();
    dTensor2     node_sub(mnodes_sub,1);
    dTensor1 prim_vol_sub(mx_sub);
    dTensorBC3      q_sub(mx_sub,meqn_sub,kmax_sub,mbc_sub);
    dTensorBC1   smax_sub(mx_sub,mbc_sub);
    dTensorBC3    aux_sub(mx_sub,maux_sub,kmax_sub,mbc_sub);

    // Output size parameters to screen
    cout << " SubGrid: " << endl;
    cout << "   Number of Equations:           " << meqn_sub << endl;
    cout << "   Number of Spatial Dimensions:  " << 1 << endl;
    cout << "   Order of Accuracy in Space:    " << method_sub[1] << endl;
    cout << "   Order of Accuracy in Time:     " << method_sub[2] << endl;
    if (method_sub[3]==1)
    {  cout << "   Limiters:                      yes " << endl;  }
    else
    {  cout << "   Limiters:                      no  " << endl;  }
    cout << "   Number of Elements in Grid:    " << mx_sub << endl;
    cout << "        Left Endpoint in Grid:    " << scientific <<  xlow_sub << endl;
    cout << "       Right Endpoint in Grid:    " << scientific << xhigh_sub << endl;
    cout << "        Number of Ghost Cells:    " << mbc_sub << "  (on each boundary) " << endl;
    cout << endl;

    // Output information to "qhelp.dat"
    string qhelp;
    qhelp=outputdir+"/qhelp.dat";
    ofstream out_file(qhelp.c_str(), ios::app);
    out_file << setprecision(16);
    out_file << meqn_sub << endl << nout << endl << method_sub[1] << endl;
    out_file << setw(24) << scientific << xlow_sub << endl << setw(24)
	     << scientific << xhigh_sub << endl << setw(24) 
	     << scientific << dx_sub << endl;
    out_file.close();


    // -------------------------
    // ***     MAIN GRID     ***
    // -------------------------

    // Construct main grid
    GridSetup(method[1],xlow,dx,node,prim_vol,outputdir,"grid.dat");
 
    // Set any auxiliary variables on computational grid
    // Set values and apply L2-projection
    if (method[6]>0)
    {  L2Project(0,1-mbc,mx+mbc,node,q,aux,aux,&AuxFunc);  }

    // Set initial data on computational grid
    // Set values and apply L2-projection
    L2Project(0,1-mbc,mx+mbc,node,q,aux,q,&QinitFunc);
    
    // Run AfterStep to set any necessary variables
    AfterStep(node,aux,q);
  
    // Output initial data to file
    // For each element, we output ``method[1]'' number of values
    Output(node,aux,q,0.0,0,outputdir,"q");

    // Compute conservation and print to file
    ConSoln(node,method,aux,q,0.0,outputdir);

    // ------------------------
    // ***     SUB GRID     ***
    // ------------------------

    // Construct main grid
    GridSetup(method_sub[1],xlow_sub,dx_sub,node_sub,prim_vol_sub,
	      outputdir,"grid_sub.dat");
 
    // Set any auxiliary variables on computational grid
    // Set values and apply L2-projection
    if (method_sub[6]>0)
    {  L2Project(0,1-mbc_sub,mx_sub+mbc_sub,node_sub,q_sub,aux_sub,aux_sub,&AuxFunc_sub);  }

    // Set initial data on computational grid
    // Set values and apply L2-projection
    L2Project(0,1-mbc_sub,mx_sub+mbc_sub,node_sub,q_sub,aux_sub,q_sub,&QinitFunc_sub);

    // Run AfterStep to set any necessary variables
    AfterStep_sub(node_sub,aux_sub,q_sub);
  
    // Output initial data to file
    // For each element, we output ``method_sub[1]'' number of values
    Output(node_sub,aux_sub,q_sub,0.0,0,outputdir,"qsub");


    // ----------------------------------
    // ***     TIME STEPPING LOOP     ***
    // ----------------------------------

    // Main loop for time stepping
    tstart = 0.0;
    tend   = 0.0;
    dtout = tfinal/double(nout);
    for (n=1; n<=nout; n++)
    {
        tstart = tend;	  
	tend = tstart + dtout;
	
	// Solve hyperbolic system from tstart to tend
	if (method[2]>1)
	{
	    // Spectral Deferred Correction time stepping
	    cout << " SDC Time-Stepping is currently not supported for HMM " << endl;
	    exit(1);
	}
	else
	{	  
	    // Runge-Kutta time stepping
 	    DogSolveRK_HMM(tstart,tend,nv,dtv,dtv_sub,cflv,map,leftBC,rightBC,
	  		   method,node,aux,q,smax,prim_vol,
	  		   method_sub,node_sub,aux_sub,q_sub,smax_sub,
	  		   prim_vol_sub,outputdir);
	}
    
	// Output data to file
	Output(node,aux,q,tend,n,outputdir,"q");
	Output(node_sub,aux_sub,q_sub,tend,n,outputdir,"qsub");
	
	// Done with solution from tstart to tend
	cout << setprecision(5);
	cout << "DOGPACK: Frame " << setw(3) << n;
	cout << ": plot files done at time t =";
	cout << setw(12) << scientific << tend << endl;
	cout << endl;
    }
    
    return 1;
   
}
