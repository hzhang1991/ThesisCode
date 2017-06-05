#include "dogdefs.h"

void ApplyLimiterKrivodonova(dTensorBC5& aux, 
                             dTensorBC5& q,
                             void (*ProjectRightEig)(int,
                                                     const dTensor1&,
                                                     const dTensor1&,
                                                     const dTensor2&,
                                                     dTensor2&),
                             void (*ProjectLeftEig)(int,
                                                    const dTensor1&,
                                                    const dTensor1&,
                                                    const dTensor2&,
                                                    dTensor2&));
