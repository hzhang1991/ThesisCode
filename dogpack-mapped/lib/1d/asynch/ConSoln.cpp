#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "DogAsynch1d.h"
using namespace std;

void ConSoln(const DogAsynch1d& dogAsynch1d, double t, string outputdir)
{
  int     mx = dogAsynch1d.get_mx();
  int   meqn = dogAsynch1d.get_meqn();
  int   maux = dogAsynch1d.get_maux();
  string fname1 = outputdir+"/conservation.dat";
  ofstream write_file1,write_file2;
  dTensor1 qsum(meqn);
  dTensor1 res_sum(meqn);

  if (t==0) 
    {
      write_file1.open(fname1.c_str(), ofstream::out);
    }
  else
    {
      write_file1.open(fname1.c_str(), ofstream::app);
    }
  
  // -----------------
  // CONSERVATION
  // -----------------
  for (int m=1; m<=meqn; m++)
    {
      qsum.set(m,0.0);
      
      for (int i=1; i<=mx; i++)
	{
	  double dtmp = dogAsynch1d.dx->get(i);
	  double qtmp = dogAsynch1d.q->get(i,m,1);
	  
	  qsum.set(m, (qsum.get(m) + dtmp*qtmp) );
	}
    }
  
  write_file1 << setprecision(16);
  write_file1 << setw(24) << scientific << t << " ";
  for (int m=1; m<=meqn; m++)
    {
      if (abs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
      write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
  write_file1 << endl;
  
}
