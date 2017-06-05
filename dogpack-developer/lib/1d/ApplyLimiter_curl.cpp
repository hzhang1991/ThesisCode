#include <cmath>
#include "dogdefs.h"

void ApplyLimiter(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q, 
		  void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
		  void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{
}
