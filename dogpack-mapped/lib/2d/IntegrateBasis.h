#ifndef IntegrateBasis_h
#define IntegrateBasis_h
#include "tensors.h"

class Integrator
{

public:
void IntegrateBasis(
        const int istart,
        const int iend,
        const int jstart,
        const int jend,
        const int QuadOrder,
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,    
        const dTensorBC4* qin,
        const dTensorBC4* auxin,
        dTensorBC4* fout);

};


// -------------------------------------------------------------------------- //
//
// Single function call used to provide backwards compatability
//
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //

#endif // IntegrateBasis_h
