#include "tensors.h"

//
// Aux data from the file aXXXX_restart.dat (where XXXX=nstart)
//
void SetAux_Unst_restart(int nstart,
			 dTensor3& aux,
			 const char* outputdir)
{
  double QinitRestart_Unst(int nstart,
			   char* varname,
			   dTensor3& q,
			   const char* outputdir);

  char* aname = new char[1];
  aname[0]='a';
  QinitRestart_Unst(nstart, aname, aux, outputdir);
  delete aname;
}
