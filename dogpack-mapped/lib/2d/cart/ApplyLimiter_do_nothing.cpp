#include "Limiters.h"
#include "DogParams.h"
#include "debug.h"

void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
		  void (*ProjectRightEig)(int,const dTensor1&,
					  const dTensor1&,const dTensor2&,
					  dTensor2&),
		  void (*ProjectLeftEig)(int,const dTensor1&,
					 const dTensor1&,const dTensor2&,
					 dTensor2&))
{
}
