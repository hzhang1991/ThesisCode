///////////////////////////////////////////////////////////////////////////////
// 
// Header file for Semi-Lagrangian time stepping methods.  This should probably
// be moved to $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/, but then the linker
// needs to include these with the -I flag for each application that links to
// this.  -DS
//
///////////////////////////////////////////////////////////////////////////////

#ifndef _SLSTATE_H_
#define _SLSTATE_H_

#include "tensors.h"

typedef struct _sl_state_vars
{

    int stage_num;          // stage number currently being integrated

    double tn;              // this doesn't change from time tn to tnp1
    double dt;              // this is the total length of the desired time step
    double dtnow;           // length of time step about to be performed
    double t;               // "time" passed in

    // contains current time of each operator 
    // split_time.get( s, op_nun ) = time for stage "s" and operator op_num
    dTensor2* split_time;  
    dTensorBC4* aux1d;     //(max(mx,my), maux1d, kmax1d, # of derivatives, mbc, ndims-1);

    // pointer to q at time tn
    dTensorBC4* qold;      //(mx, my, meqn, kmax, mbc)
    dTensorBC4* qnew;      //(mx, my, meqn, kmax, mbc)

    dTensor1* dt_stages;   // pointer to time steps for each time step
    dTensor2* node1d;      // 1d nodal points (1:mpoints+1, 1)

} SL_state;
#endif
///////////////////////////////////////////////////////////////////////////////
