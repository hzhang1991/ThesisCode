#ifndef DOG_MATH_H
#define DOG_MATH_H
#include <cmath>
// This header extends math.h

// function declarations
// (only declaring the functions from dog_math.cpp that actually get used)

//////////////////////
// inline functions //
//////////////////////

// sgnum function.  undefined behavior for x near zero
inline double sgn( double x )
{ return x<0?-1.:1.; }

inline double Max(double a, double b)
{ return a>b?a:b; }

inline double Min(double a, double b)
{ return a<b?a:b; }

inline int iMax(int a, int b)
{ return a>b ? a : b; }

inline int iMin(int a, int b)
{ return a<b ? a : b; }

inline int iMod(int n, int m)
{
    int modval = n%m;
    return n%m < 0 ? (modval)+m : modval;
}

inline double minmod(double a, double b)
{
    if(std::signbit(a))
    {
      if(std::signbit(b)) return Max(a,b);
      return 0.;
    }
    if(std::signbit(b)) return 0.;
    return Min(a,b);
}
inline double minmod(double a, double b, double c)
{
  if(std::signbit(a))
  {
    if(std::signbit(b) && std::signbit(c)) return Max(Max(a,b),c);
    return 0;
  }
  if(std::signbit(b) || std::signbit(c)) return 0;
  return Min(Min(a,b),c);
}

// poor man's Factorial function //
// Note: if n < 0, this function returns 1 //
inline int Factorial(int n)
{
    int result = 1;
    while( n > 1 )
    { result *= n--; }
    return result;
}
#endif // DOG_MATH_H
