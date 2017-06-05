#include <cmath>
#include "tensors.h"
#include "dog_math.h"

// Upwind method used for LxW flux function.
// We already know the LxW flux function, so no need to compute that here...
double RiemannSolveLxWupwind(const dTensor1& xedge,
			     const dTensor1& Ql, const dTensor1& Qr, const dTensor1& Auxl,
			     const dTensor1& Auxr, const dTensor1& ffl, const dTensor1& ffr,
			     dTensor1& Fl, dTensor1& Fr,
			     void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,
						const dTensor1&,const dTensor1&,double&,double&))
{
    int m;
    double smax_edge = 0.0e0;
    int meqn = Ql.getsize();
    double s1,s2;

    // Calculate minimum and maximum HLLE speeds
    SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1,s2);
    smax_edge = Max(fabs(s1),fabs(s2));

    // Calculate Fluxes (LLF Riemman solver)
//    for (m=1; m<=meqn; m++)
//    {
//        Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
//                         + smax_edge*(Ql.get(m) - Qr.get(m)) ) );
//        Fr.set(m, Fl.get(m) );
//    }

    // Upwind method for meq == 1
    // TODO there's a bug here with comparing s1 > 0 versus s1 >= 0 for meqn == 1
    double qstar = 0.0;
    for(m=1; m<=meqn; m++)
    {
        //if( s1 >= 0 && s2 > 0)
        if( s1 >= 0 ) // && s2 > 0)
        {   // positive speeds
            Fl.set(m, ffl.get(m) );
            Fr.set(m, Fl.get(m) );
        //} else if ( s2 <= 0 && s1 < 0)
        } else if ( s1 < 0 ) // && s1 < 0)
        {   // all negative speeds
            Fl.set(m, ffr.get(m) );
            Fr.set(m, Fl.get(m) );
        } 
//        else // this case shouldn't occur for meqn == 1
//        {    // positive and negative speeds (HLLE Reimman solver)
//             cout << "     bad case here!! " << endl;
//             qstar = (Qr.get(m)-Ql.get(m)+s1*Ql.get(m)-s2*Qr.get(m))/(s1-s2);
//             Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) + (s1+s2)*qstar 
//                               - s1*Ql.get(m) + s2*Qr.get(m) ) );
//                 + smax_edge*(Ql.get(m) - Qr.get(m)) ) );
//             Fr.set(m, Fl.get(m) );
//        }
    }

    return smax_edge;
}
