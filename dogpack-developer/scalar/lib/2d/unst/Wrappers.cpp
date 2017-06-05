#include "dogdefs.h"

void AuxFuncWrapper(
    const dTensor2* vel_vec,
    const dTensor2& xpts,
    const dTensor2& NOT_USED_1,
    const dTensor2& NOT_USED_2,
    dTensor2& auxvals)
{
    void AuxFunc(
        const dTensor2* vel_vec,
        const dTensor2& xpts, dTensor2& auxvals);
    AuxFunc(vel_vec, xpts,auxvals);
}

void QinitFuncWrapper(
    const dTensor2* vel_vec,
    const dTensor2& xpts,
    const dTensor2& NOT_USED_1,
    const dTensor2& NOT_USED_2,
    dTensor2& qvals)
{
    void QinitFunc(
        const dTensor2* vel_vec,
        const dTensor2& xpts, dTensor2& qvals);
    QinitFunc(vel_vec, xpts,qvals);
}
