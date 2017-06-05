#include "../defs.h"
#include <sstream>
#include <string>

void GetParamsHMM(int& nout,int& nv,int method[],double& tfinal,
		  double dtv[],double cflv[],int& mx,int& mnodes,double& dx,
		  int& meqn,int& maux,int& kmax,int& mbc,double& xlow,double& xhigh,
		  int method_sub[],int& meqn_sub,int& maux_sub,int& kmax_sub,
		  double& xlow_sub,double& xhigh_sub, string outputdir)
{
    int itmp1,itmp2,itmp3;
    double dtmp1,dtmp2,dtmp3;

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
    ostringstream command;
    command
      << "if test -f startscript && test -x startscript;\n"
      << "then ./startscript " << outputdir << "\n"
      << "else ${DOGPACK}/scripts/startscript " << outputdir << "\n"
      << "fi";
    //cout << command.str() << endl;
    int exit_status = system(command.str().c_str());
    //// stop if error in script
    //if(exit_status!=0)
    //{
    //   exit(exit_status);
    //}

    string qhelp;
    qhelp=outputdir+"/qhelp.dat";
    ifstream read_file("dogpack.data", ios::in);
    ofstream out_file(qhelp.c_str(), ios::out);
    char buffer[256];

    if(read_file.is_open()!=1)
    {
        cout << " ERROR: file 'dogpack.data' not found " << endl << endl;
        exit(1);
    }

    // Read in parameters from redpack.data
    read_file >> nout;
    read_file.getline(buffer,256);
    read_file >> tfinal;
    read_file.getline(buffer,256);
    read_file >> dtv[1];
    read_file.getline(buffer,256);
    read_file >> dtv[2];
    read_file.getline(buffer,256);
    read_file >> cflv[1];
    read_file.getline(buffer,256);
    read_file >> cflv[2];
    read_file.getline(buffer,256);
    read_file >> nv;
    read_file.getline(buffer,256);
    read_file.getline(buffer,256);

    // Information for main grid
    read_file >> method[1];  
    read_file.getline(buffer,256);
    read_file >> method[2];
    read_file.getline(buffer,256);
    read_file >> method[3];
    read_file.getline(buffer,256);
    read_file >> method[4];
    read_file.getline(buffer,256);
    read_file >> method[5];
    read_file.getline(buffer,256);
    read_file >> method[6];
    read_file.getline(buffer,256);
    read_file >> method[7];
    read_file.getline(buffer,256);
    read_file >> meqn;
    read_file.getline(buffer,256);
    maux = method[6];
    kmax = method[1];    

    read_file >> mx;
    read_file.getline(buffer,256);
    mnodes = mx + 1;

    read_file >> mbc;
    read_file.getline(buffer,256);
  
    read_file >> xlow;
    read_file.getline(buffer,256);
    read_file >> xhigh;
    read_file.getline(buffer,256);
    read_file.getline(buffer,256);
  
    dx = (xhigh-xlow)/double(mx);

    // Information for sub grid
    method_sub[1] = method[1];
    method_sub[2] = method[2];
    read_file >> method_sub[3];
    read_file.getline(buffer,256);
    read_file >> method_sub[4];
    read_file.getline(buffer,256);
    read_file >> method_sub[5];
    read_file.getline(buffer,256);
    read_file >> method_sub[6];
    read_file.getline(buffer,256);
    read_file >> method_sub[7];
    read_file.getline(buffer,256);
    read_file >> meqn_sub;
    read_file.getline(buffer,256);
    maux_sub = method_sub[6];
    kmax_sub = method_sub[1];
    
    read_file >> xlow_sub;
    read_file.getline(buffer,256);
    read_file >> xhigh_sub;
    read_file.getline(buffer,256);

    read_file.close();    

    // Error check 1
    if (cflv[2] >= cflv[1])
    {
        cout << " ERROR: desired CFL must be strictly smaller than max CFL " 
	     << endl;
	cout << endl;
        exit(1);
    }

    // Error check 2
    if (method[5] > method[6])
    {
        cout << " ERROR: mcapa cannot be larger than maux " << endl;
	cout << "    mcapa = " << method[5] << endl;
	cout << "     maux = " << method[6] << endl;
	cout << endl;
        exit(1);
    }
  
    // Error check 3
    if (meqn < 1)
    {
        cout << " ERROR: meqn must be at least 1 " << endl;
	cout << "    meqn = " << meqn << endl;
	cout << endl;
        exit(1);
    }
  
    // Error check 4
    if (method[1]<1 || method[1]>5 )
    {
        cout << " ERROR: must have method[1] = 1, 2, 3, 4, or 5 " << endl;
        cout << endl;
        exit(1);
    }
  
    // Error check 5
    if ( method[2]==0 || method[2]<-3 || method[2]>5 )
    {
        cout << " ERROR: must have method[2] = -3, -2, -1, or 1, 2, 3, 4, 5 ";
	cout << endl;
        cout << endl;
        exit(1);
    }

    // Error check 6
    if (method_sub[5] > method_sub[6])
    {
        cout << " ERROR: mcapa_sub cannot be larger than maux_sub " << endl;
	cout << "    mcapa_sub = " << method_sub[5] << endl;
	cout << "     maux_sub = " << method_sub[6] << endl;
	cout << endl;
        exit(1);
    }
  
    // Error check 7
    if (meqn_sub < 1)
    {
        cout << " ERROR: meqn must be at least 1 " << endl;
	cout << "    meqn_sub = " << meqn_sub << endl;
	cout << endl;
        exit(1);
    }
  
    // Error check 8
    if (method_sub[1]<1 || method_sub[1]>5 )
    {
        cout << " ERROR: must have method_sub[1] = 1, 2, 3, 4, or 5 " << endl;
        cout << endl;
        exit(1);
    }
  
    // Error check 9
    if ( method_sub[2]==0 || method_sub[2]<-3 || method_sub[2]>5 )
    {
        cout << " ERROR: must have method_sub[2] = -3, -2, -1, or 1, 2, 3, 4, 5 ";
	cout << endl;
        cout << endl;
        exit(1);
    }
  
    // Output meqn and nout for plotting purposes (main grid)
    out_file << setprecision(16);
    out_file << meqn << endl << nout << endl << method[1] << endl;
    out_file << setw(24) << scientific << xlow << endl << setw(24)
	     << scientific << xhigh << endl << setw(24) 
	     << scientific << dx << endl;    
    out_file.close();
    
}
