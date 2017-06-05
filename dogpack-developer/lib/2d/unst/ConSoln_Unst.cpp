#include "dogdefs.h"
#include "mesh.h"
#include "DogParams.h"
#include <fstream>

void ConSoln_Unst(const mesh& Mesh, 
		  const dTensor3& aux,
		  const dTensor3& q, 
		  double t, 
		  const char* outputdir)
{
  const int NumElems = q.getsize(1);
  const int     meqn = q.getsize(2);
  const int     kmax = q.getsize(3);
  const int     maux = aux.getsize(2);
  char fname1[1024];
  snprintf(fname1,1024,"%s/conservation.dat",outputdir);
  dTensor1 qsum(meqn);
  dTensor1 res_sum(meqn);
  const int NumPhysElems = Mesh.get_NumPhysElems();

  FILE* write_file1;
  if (t==0) 
    {  write_file1 = fopen(fname1,"w");  }
  else
    {  write_file1 = fopen(fname1,"a");  }

  // -----------------
  // CONSERVATION
  // -----------------
  if (dogParams.get_mcapa()<1) // without capacity function
    {
      for (int m=1; m<=meqn; m++)
        {
	  qsum.set(m,0.0);

	  for (int i=1; i<=NumPhysElems; i++)
            {             
	      double dtmp = Mesh.get_area_prim(i);
	      double qtmp = q.get(i,m,1);

	      qsum.set(m, (qsum.get(m) + dtmp*qtmp) );
            }
        }
    }

  fprintf(write_file1,"%24.16e  ",t);
  for (int m=1; m<=meqn; m++)
    {
      fprintf(write_file1,"%24.16e  ",qsum.get(m));
    }
  fprintf(write_file1,"\n");
  fclose(write_file1);
}
