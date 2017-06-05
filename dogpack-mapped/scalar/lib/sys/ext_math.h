#ifndef ext_math_h
#define ext_math_h
#include <cmath>
#include <float.h>

// override any previous definitions?
// (If you do not trust the performance or reliability of
// vendor-supplied versions in math.h.)
//#undef isnan
//#undef isinf
//#undef isfinite
//
#ifndef isnan
// How do I test for -funsafe-math-optimizations?
// #define isnan(x) false
# define isnan(x) \
    (sizeof (x) == sizeof (long double) ? isnan_ld (x) \
     : sizeof (x) == sizeof (double) ? isnan_d (x) \
     : isnan_f (x))
// These definitions work only if -funsafe-math-optimizations is not set;
// otherwise cannot do this inline (would have to implement functions in a
// an object module compiled without -funsafe-math-optimizations,
// tending to defeat the desired code acceleration).
//
// These definitions work only if -funsafe-math-optimizations is not set:
static inline int isnan_f  (float       x) { return (x != x); }
static inline int isnan_d  (double      x) { return (x != x); }
static inline int isnan_ld (long double x) { return (x != x); }
#if 0
// These definitions work only if -funsafe-math-optimizations is set
// (but they still do not work if -O1 is set).
static inline int isnan_f  (volatile float       x) { return (x < x); }
static inline int isnan_d  (volatile double      x) { return (x < x); }
static inline int isnan_ld (volatile long double x) { return (x < x); }
// So these definitions would work unless -O1 is set.
static inline int isnan_f(volatile float  x) { return (!(x<x)^(x==x)); }
static inline int isnan_d(volatile double x) { return (!(x<x)^(x==x)); }
static inline int isnan_ld(volatile long double x) { return (!(x<x)^(x==x)); }
#endif
#endif // isnan

#ifndef isinf
# define isinf(x) \
    (sizeof (x) == sizeof (long double) ? isinf_ld (x) \
     : sizeof (x) == sizeof (double) ? isinf_d (x) \
     : isinf_f (x))
//
static inline int isinf_f  (float       x) { return (x>FLT_MAX || x<-FLT_MAX); }
static inline int isinf_d  (double      x) { return (x>DBL_MAX || x<-DBL_MAX); }
static inline int isinf_ld (long double x) { return (x>LDBL_MAX || x<-LDBL_MAX); }
//
// the following versions do not work if gcc flag
// -funsafe-math-optimizations is set (e.g. if -ffast-math is set)
//
//static inline int isinf_f  (float       x) { return isnan (x - x); }
//static inline int isinf_d  (double      x) { return isnan (x - x); }
//static inline int isinf_ld (long double x) { return isnan (x - x); }
#endif // isinf

#ifndef isfinite
# define isfinite(x) \
    (sizeof (x) == sizeof (long double) ? isfinite_ld (x) \
     : sizeof (x) == sizeof (double) ? isfinite_d (x) \
     : isfinite_f (x))
//
static inline int isfinite_f  (float       x) { return (x<=FLT_MAX && x>=-FLT_MAX); }
static inline int isfinite_d  (double      x) { return (x<=DBL_MAX && x>=-DBL_MAX); }
static inline int isfinite_ld (long double x) { return (x<=LDBL_MAX && x>=-LDBL_MAX); }
#endif // isfinite

#endif // ext_math_h
