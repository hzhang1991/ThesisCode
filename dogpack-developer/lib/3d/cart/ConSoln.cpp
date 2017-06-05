#include<iomanip>
#include<fstream>
#include<iostream>
#include<cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart3.h"
const char* get_outputdir();

using namespace std;

void ConSoln(const dTensorBC5& aux, const dTensorBC5& q, double t)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mz   = q.getsize(3);
  const int meqn = q.getsize(4);
  const int kmax = q.getsize(5);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(4);
  
  string fname1 = string(get_outputdir())+"/conservation.dat";
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
  const double dx = dogParamsCart3.get_dx();
  const double dy = dogParamsCart3.get_dy();
  const double dz = dogParamsCart3.get_dz();
  const double vol = dx*dy*dz;
  if (dogParams.get_mcapa()<1) // without capacity function
    {
      for (int m=1; m<=meqn; m++)
        {
          qsum.set(m,0.0);
          
          for (int i=1; i<=mx; i++)	    
            for (int j=1; j<=my; j++)
	      for (int k=1; k<=mz; k++)
		{
		  double qtmp = q.get(i,j,k,m,1);
		  
		  qsum.set(m, qsum.get(m) + vol*qtmp );
		}
        }
    }
  
  write_file1 << setprecision(16);
  write_file1 << setw(24) << scientific << t << " ";
  for (int m=1; m<=meqn; m++)
    {
      if (fabs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
      write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
  write_file1 << endl;
  
  write_file1.close();
  
}
