#include "dogdefs.h"
#include "DogParams.h"
#include "MonomialsToLegendre.h"
#include "mesh.h"

// This simple wrapper function allows for variable names to be passed in to
// WriteOutput_Unst.
void Output_Unst(const mesh& Mesh, 
		 const dTensor3& aux, 
		 const dTensor3& q, 
		 double t, 
		 int nframe, 
		 const char* outputdir)
{ 
  void WriteOutput_Unst(char* fname, 
			const dTensor3* q, 
			double t);
  // Output Q values
  char fname1[1024];
  assert( snprintf(fname1, 1024, "%s/q%04d", outputdir, nframe) > 0 );
  WriteOutput_Unst(fname1, &q, t);

  if( dogParams.get_maux() > 0 ) 
    {  
      // Output aux values
      char fname2[1024];
      assert( snprintf(fname2, 1024, "%s/a%04d", outputdir, nframe) > 0 );
      WriteOutput_Unst(fname2, &aux, t);
    }

}

// Output solution on unstructured mesh
void WriteOutput_Unst(char* fname, 
		      const dTensor3* q,
		      double t)
{
  const int NumElems = q->getsize(1);
  const int     meqn = q->getsize(2);
  const int     kmax = q->getsize(3);

  char fname1[1024];
  assert( snprintf(fname1,1024,"%s.dat",fname) > 0 );

  FILE* write_file = fopen(fname1,"w");
  assert( write_file != NULL );

  fprintf(write_file,"%24.19e\n",t);
  for (int k=1; k<=kmax; k++)
    for (int m=1; m<=meqn; m++)      
      for (int i=1; i<=NumElems; i++)
	{
	  double tmp = q->get(i,m,k);
	  fprintf(write_file, "%24.19e\n", tmp);
	}
  fclose(write_file);
}
