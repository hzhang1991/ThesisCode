#include "tensors.h"
#include "diffusion.h"

// Optional call to function that adds extra piece into
// solution right-hand side (i.e., q_t = Lstar)
void LstarExtra(const dTensor2& node, 
		dTensorBC3& aux, 
		dTensorBC3& q, 
		dTensorBC3& Lstar)
{

}
