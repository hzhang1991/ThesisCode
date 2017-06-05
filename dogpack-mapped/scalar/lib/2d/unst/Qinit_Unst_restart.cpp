#include "dogdefs.h"
#include "DogParams.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>

// Initial data from the file qXXXX_restart.dat (where XXXX=nstart)
double QinitRestartASCII_Unst(string fname, dTensor3& q)
{
    int i,j,m,k;
    int NumElems = q.getsize(1);
    int     meqn = q.getsize(2);
    int     kmax = q.getsize(3);
    double tmp,t;

    fname += ".dat";
    ifstream restart_file(fname.c_str(), ios::in);

    if(!restart_file.is_open())
    {
        cerr << "Qinit_Unst_restart: failed to open file " << fname << endl;
        exit(1);
    }

    // read in old data (written in Output.cpp)
    restart_file >> t;
    for (i=1; i<=NumElems; i++)
      {	    			
	for (m=1; m<=meqn; m++)
	  {
	    for (k=1; k<=kmax; k++)
	      {
		assert(restart_file >> tmp);
		q.set(i,m,k, tmp );
	      }
	  }
      }
    string word;
    restart_file >> word;
    assert(restart_file.eof());
    restart_file.close();
    return t;
}

double QinitRestart_Unst(int nstart, string varname, dTensor3& q, string outputdir)
{
    double retval;

    ostringstream basename;
    basename << setfill('0') << setw(4) << nstart;
    string fname  = outputdir+"/"+varname+basename.str()+"_restart";

    switch(dogParams.get_datafmt())
    {
    case ASCII:
      double QinitRestartASCII_Unst(string fname, dTensor3& q);
      retval = QinitRestartASCII_Unst(fname, q);
      break;
    default:
      cerr << "Qinit_Unst_restart.cpp datafmt: " << dogParams.get_datafmt()
           << " not supported" << endl;
      exit(1);
      break;
    }
    return retval;
}

double Qinit_Unst_restart(int nstart, dTensor3& q, string outputdir)
{
  return QinitRestart_Unst(nstart, "q", q, outputdir);
}

