#ifndef _LEGENDREPOLYS1D_H_
#define _LEGENDREPOLYS1D_H_

#include "dogdefs.h"
#include "dog_math.h"
#include "tensors.h"
#include "constants.h"
#include "Quadrature.h"

void evaluateLegendrePolys(const dTensor1& spts, dTensor2& phi );
void evaluateLegendrePolys( const double dx, const dTensor1& spts, dTensor2&
    phi, dTensor2& phi_x);

#endif
