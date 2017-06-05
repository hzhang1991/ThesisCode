#include "tensors.h"
#include <string>
using namespace std;

//
// Aux data from the file aXXXX_restart.dat (where XXXX=nstart)
//
void SetAux_Unst_restart(int nstart, dTensor3& aux, string outputdir)
{
  double QinitRestart_Unst(int nstart, string varname,
			   dTensor3& q, string outputdir);

  QinitRestart_Unst(nstart, "a", aux, outputdir);
}
