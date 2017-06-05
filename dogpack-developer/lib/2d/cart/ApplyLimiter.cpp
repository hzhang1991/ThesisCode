#include "tensors.h"
#include "DogParams.h"
#include "Legendre2d.h"
#include "Limiters.h"

void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
    void (*ProjectRightEig)(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&),
    void (*ProjectLeftEig)(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&))
{
    // assumes that cell averages satisfy positivity
    ApplyLimiterKrivodonova(aux, q, ProjectRightEig, ProjectLeftEig);
}
