#include <cmath>
#include "tensors.h"
#include "dog_math.h"

// Local Lax-Friedrichs solver used for Lax Wendroff (and two-derivative)
// time stepping.
//
// The default library function evaluates the flux function at cell edges.
// The Lax-Wendroff method ( as well as the multiderivative schemes ) rely on
// a variation of the actual flux function, and so there is no need to
// recompute those values again here.
//
// To use the HLL(E) solver, one needs to replace the linker in order to link to
// this.  Add the following to your Makefile:
//
//     RiemannSolveLxW  = $(LIB2D)/LaxWendroff/RiemannSolveLxW_HLL
//
// if you do not wish to use the LLF solver.
//
double RiemannSolveLxW(const dTensor1& nvec,const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr, const dTensor1& Auxl,
        const dTensor1& Auxr, const dTensor1& ffl, const dTensor1& ffr,
        dTensor1& Fl, dTensor1& Fr)
{

    const int meqn = Ql.getsize();

    // Calculate minimum and maximum HLLE speeds
    double s1,s2;
    void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
                    const dTensor1& Ql, const dTensor1& Qr,
                    const dTensor1& Auxl, const dTensor1& Auxr,
                    double& s1, double& s2);
    SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

    double smax_edge = Max(fabs(s1),fabs(s2));

    // Calculate Fluxes (LLF Riemman solver)
    for(int m=1; m<=meqn; m++)
    {
        Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                    + smax_edge*(Ql.get(m) - Qr.get(m)) ) );
        Fr.set(m, Fl.get(m) );
    }

    return smax_edge;
}
