#include <cmath>
#include "fcmp.h"
#include "constants.h"
//
// test whether two doubles are equal
//
bool test_equal(double a, double b)
{
    return !fcmp(a,b,EPSILON);
}

// Truncating such very small numbers to zero was done in
// clawpack at some point for some reason (perhaps because a
// 3-digit exponent adds another character to the total width or
// perhaps because an old version of matlab fscanf ignored the
// third digit?).  I think that we should not use this function
// when writing out the state to allow for the possibility of
// restarting with exactly the same data down to the last bit.
// 
double matlab_fix_double(double val)
{
    if(fabs(val) < 1.0e-99) return 0.;
    return val;
}

#if 0
// test whether two doubles are equal
//
bool test_equal(double a, double b)
{
    double size = abs(a)+abs(b);
    double diff = abs(a-b);
    if(size < 1)
    { if(diff<EPSILON) return true;
    }
    else
    {
        double ratio = abs(b-a)/size;
        if(ratio<EPSILON) return true;
    }
    return false;
}
#endif
