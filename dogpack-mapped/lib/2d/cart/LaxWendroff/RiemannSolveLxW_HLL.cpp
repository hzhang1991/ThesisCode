#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "dogdefs.h"

// This is the HLL(E) method.  The HLLE method defines left and right
// speeds by
//
// s1 := minimum speed of waves for left state and intermediate state and
// s2 := maximum speed of waves for right state and intermediate state,
//       where intermediate state is theoretically the Roe average of
//       the left and right states and in practice is often the simple
//       average.
//
// (The original HLL method defines s1 and s2 as extreme characteristic
// speeds of Riemann solutions.)
//
// In case s1 < 0 < s2 approximate the Riemann solution by two
// shockwaves of speeds s1 and s2 separating a single constant
// intermediate state Q* from the left state Ql and the right state Qr.
// Define Fl := F(Ql) and Fr := F(Qr).
// Conservation (Rankine-Hugoniot) at the shock fronts
// gives a system of two equations which can be solved
// for the two unknowns Q* and F*:
//
// left  shock: (Fr - F*) = (Qr - Q*) s2
// right shock: (F* - Fl) = (Q* - Ql) s1
//
// Solving gives:
//
//   Q* = (Fr - Fl + s1 Ql - s2 Qr)/(s1-s2)  and
//   F* = (Fl + Fr - s1 Ql - s2 Qr + (s1+s2)Q*)/2.
//

double RiemannSolveLxW(const dTensor1& nvec,
    const dTensor1& xedge,
    const dTensor1& Ql,   const dTensor1& Qr,
    const dTensor1& Auxl, const dTensor1& Auxr,
    const dTensor1& ffl,  const dTensor1& ffr,
    dTensor1& Fl, dTensor1& Fr)
{

    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();

    // HLLE parameters
    double s1    = 0.;
    double s2    = 0.;
    double Qstar = 0.;

    // Calculate minimum and maximum HLLE speeds
    void SetWaveSpd(const dTensor1& nvec,
                    const dTensor1& xedge,
                    const dTensor1& Ql,
                    const dTensor1& Qr,
                    const dTensor1& Auxl,
                    const dTensor1& Auxr,
                    double& s1,
                    double& s2);
    SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

    int mcase = -1;

    // Calculate Fluxes
    if (fabs(s1)<=1.0e-12 && fabs(s2)<=1.0e-12)
    { mcase = 1; }
    else
    {
        if (s1*s2 > 0.0)
        {
            if (s1>0)
            { mcase = 1; }
            else
            { mcase = 2; }
        }
        else
        { mcase = 3; }
    }

    switch ( mcase )
    {
        case 1:  // both s1 and s2 are positive (or both zero)

            for (int m=1; m<=meqn; m++)
            { 
                Fl.set(m, ffl.get(m) );
                Fr.set(m, ffl.get(m) );
            }	
            break;

        case 2:  // both s1 and s2 are negative

            for (int m=1; m<=meqn; m++)
            { 
                Fl.set(m, ffr.get(m) );
                Fr.set(m, ffr.get(m) );
            }	
            break;

        case 3:  // s1 is negative and s2 is positive

            for (int m=1; m<=meqn; m++)
            {
                Qstar = (ffr.get(m)-ffl.get(m) 
                        + s1*Ql.get(m) - s2*Qr.get(m))/(s1-s2);
                Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                            + (s1+s2)*Qstar - s1*Ql.get(m) - s2*Qr.get(m)) );
                Fr.set(m, Fl.get(m) );
            }	
            break;
        default:
            printf("you did not find a valid case\n");
            exit(1);
    }

    double smax_edge = Max(fabs(s1), fabs(s2));
    return smax_edge;
}
