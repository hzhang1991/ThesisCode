#ifndef _DOGSOLVE_TD_H__
#define _DOGSOLVE_TD_H__

#include <stdlib.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"

// -------------------------------------------------------------------------- //
// Function definitions
// -------------------------------------------------------------------------- //

void CopyQ(const dTensorBC4& qin,dTensorBC4& qout);

void LaxWendroff(double dt, 
    double alpha1,  double beta1,      
    dTensorBC4& aux, dTensorBC4& q,    
    dTensorBC4& Lstar, dTensorBC3& smax);

void LaxWendroffTD(double dt, 
    double alpha1,  double beta1,      // parameteres used for each derivative.
    dTensorBC4& aux1, dTensorBC4& q1,    // set bndy values modifies these
    double alpha2,  double beta2,      // parameteres used for each derivative.
    dTensorBC4& aux2, dTensorBC4& q2,    // set bndy values modifies these
    dTensorBC4& Lstar, dTensorBC3& smax);

void SetBndValues(dTensorBC4&,dTensorBC4&);

void StepLxW( double dt, const dTensorBC4& qold, 
    const dTensorBC4& L, dTensorBC4& qnew );
// -------------------------------------------------------------------------- //

#endif
