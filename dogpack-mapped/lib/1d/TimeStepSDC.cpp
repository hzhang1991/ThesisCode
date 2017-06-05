#include "constants.h"
#include "tensors.h"

void TimeStepSDC(int method2, 
        double t, 
        double dt, 
        dTensor1& dtvec, 
        dTensor1& tvec)
{
    // --------------------------------------------------------------
    // Select Gauss-Lobatto points in time and
    // construct the initial right-hand side for SDC
    switch(method2)
    {
        case 2: 
            dtvec.set(1, dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );  
            break;

        case 3:
            dtvec.set(1, 0.5*dt );
            dtvec.set(2, 0.5*dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );
            tvec.set(3, tvec.get(2) + dtvec.get(2) );
            break;

        case 4:
            dtvec.set(1, 0.25*dt );
            dtvec.set(2, 0.50*dt );
            dtvec.set(3, 0.25*dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );
            tvec.set(3, tvec.get(2) + dtvec.get(2) );
            tvec.set(4, tvec.get(3) + dtvec.get(3) );
            break;

        case 5:
            dtvec.set(1, (0.5 - 0.25*sq2) * dt );
            dtvec.set(2,       (0.25*sq2) * dt );
            dtvec.set(3,       (0.25*sq2) * dt );
            dtvec.set(4, (0.5 - 0.25*sq2) * dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );
            tvec.set(3, tvec.get(2) + dtvec.get(2) );
            tvec.set(4, tvec.get(3) + dtvec.get(3) );
            tvec.set(5, tvec.get(4) + dtvec.get(4) );
            break;

        case 6:
            dtvec.set(1, dt / 5.0);
            dtvec.set(2, dt / 5.0);
            dtvec.set(3, dt / 5.0);
            dtvec.set(4, dt / 5.0);
            dtvec.set(5, dt / 5.0);

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );
            tvec.set(3, tvec.get(2) + dtvec.get(2) );
            tvec.set(4, tvec.get(3) + dtvec.get(3) );
            tvec.set(5, tvec.get(4) + dtvec.get(4) );
            tvec.set(6, tvec.get(5) + dtvec.get(5) );
            break;

    }
    // --------------------------------------------------------------
}
