#include "dogdefs.h"
#include "L2Project.h"

void L2ProjectInitialCond(const int istart,
			  const int iend,
			  const int jstart,
			  const int jend,
			  const int kstart,
			  const int kend,
			  const int QuadOrder,
			  const int BasisOrder_qin,
			  const int BasisOrder_auxin,
			  const int BasisOrder_fout,    
			  const dTensorBC5* qin,
			  const dTensorBC5* auxin,
			  dTensorBC5* fout,
			  void (*Func)(const dTensor2&,
				       const dTensor2&,
				       const dTensor2&,
				       dTensor2&))
{
  L2Project(istart,iend,jstart,jend,kstart,kend,
	    QuadOrder,BasisOrder_qin,BasisOrder_auxin,
	    BasisOrder_fout,qin,auxin,fout,Func);
}
