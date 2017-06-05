#ifndef Interval_h
#define Interval_h
#include "assert.h"
#include "debug.h"

// This class might cause executable bloat
// but inline might be needed for performance.
//
class Interval
{

 private:
  double a;
  double b;
 private:
  static inline double Max(double a, double b)
  { return a>b?a:b; }

  static inline double Min(double a, double b)
  { return a<b?a:b; }

 public:
  Interval(){}
  Interval(double a_in, double b_in):a(a_in),b(b_in)
  {
    // all methods assume and should maintain this condition
    assert_le(a,b);
    //assert(a<=b);
  }
  double max() const {return b;}
  double min() const {return a;}
  void set_max(double in) {b=in;}
  void set_min(double in) {a=in;}
  // could define e.g. operator+= instead of these
  Interval& operator+=(Interval const& x)
  {
    a+=x.a; b+=x.b;
    return *this;
  }
  Interval& operator+=(double x)
  {
    a+=x; b+=x;
    return *this;
  }
  void add(Interval const& x) { a+=x.a; b+=x.b; }
  void subtract(Interval const& x) { a-=x.b; b-=x.a; }
  void subtract(double x) { a-=x; b-=x; }
  // multiplication
  Interval square() const
  {
    Interval out;
    if(a>=0.)
    {
      out.a = a*a;
      out.b = b*b;
    }
    else if(b<=0.)
    {
      out.a = b*b;
      out.b = a*a;
    }
    else
    {
      out.a = 0.;
      double temp = Max(-a,b);
      out.b = temp*temp;
    }
    return out;
  }
  Interval operator*(double c) const
  {
    Interval z;
    if(c<0.)
    {
      z.a = b*c;
      z.b = a*c;
    }
    else
    {
      z.a = a*c;
      z.b = b*c;
    }
    return z;
  }
  Interval operator*(const Interval& c) const
  {
    const double& x1=a;
    const double& x2=b;
    const double& c1=c.a;
    const double& c2=c.b;
    Interval z;
    double& z1=z.a;
    double& z2=z.b;
    if(x1>=0.) // pos * something
    {
      if(c1>=0.) // pos * pos
      {
        z1 = x1*c1; // small
        z2 = x2*c2; // big
      }
      else if(c2<=0.) // pos * neg
      {
        z1=x2*c1; // big
        z2=x1*c2; // small
      }
      else // pos * 0 (x2 is big, controls)
      {
        z1=x2*c1;
        z2=x2*c2;
      }
    }
    else if(x2<=0.) // neg * something
    {
      if(c2<=0.) // neg * neg
      {
        z1=x2*c2; // small
        z2=x1*c1; // big
      }
      else if(c1>=0.) // neg * pos
      {
        z1=x1*c2; // big
        z2=x2*c1; // small
      }
      else // neg * 0 (x1 is big, controls & flips)
      {
        z1=x1*c2;
        z2=x1*c1;
      }
    }
    else // 0 * something
    {
      if(c1>=0.) // 0 * pos (c2 controls)
      {
        z1=x1*c2;
        z2=x2*c2;
      }
      else if(c2<=0.) // 0 * neg (c1 controls, flips)
      {
        z1=x2*c1;
        z2=x1*c1;
      }
      // 0-centered intervals persist, so it is reasonable
      // to expect them to be common, but not if intervals
      // are small (which occurs as mesh is refined).
      // if 0-centered intervals are common
      // we could accelerate Interval arithmetic
      // by adding signum to the class.
      // (nested switch instead of if,else if, else.;
      // so only two conditional branchings in all cases
      // instead of as many as the four required to get here)
      // But maintaining signum adds expense e.g. to addition.
      else // 0 * 0
      {
        // so in this case we still
        // get 6 branchings and 4 multiplications
        z2 = Max(x2*c2,x1*c1);
        z1 = Min(x1*c2,x2*c1);
      }
    }
    return z;
  }
  Interval& operator*=(Interval const& x)
  {
    *this = (*this) * x;
    return *this;
  }
  // less efficient version of multiplication;
  // used to test correctness of more efficient version.
  void multiply(Interval const& x)
  {
    // involves four multiplications and 6 branchings
    // in all cases.  Can reduce this except for
    // 0-centered intervals.
    //
    const double& c=x.a;
    const double& d=x.b;
    double ac = a*c;
    double ad = a*d;
    double bc = b*c;
    double bd = b*d;
    a = ac;
    a = Min(a,ad);
    a = Min(a,bc);
    a = Min(a,bd);
    b = ac;
    b = Max(b,ad);
    b = Max(b,bc);
    b = Max(b,bd);
  }
  Interval reciprocal() const
  {
    Interval temp;
    if(a>0.||b<0.)
    {
      temp.a = 1./b;
      temp.b = 1./a;
    }
    else
    {
      eprintf("division by interval containing zero: [%f, %f", a,b);
    }
    return temp;
  }
  Interval& operator/=(Interval const& x)
  {
    *this = (*this) * x.reciprocal();
    return *this;
  }

  // cannot divide by zero
  void divide(Interval const& x)
  {
    (*this)*=x.reciprocal();
    // multiply(x.reciprocal());
  }

  //
  bool operator>=(double in) const
  {
    if(a>=in) return true;
    return false;
  }

  //
  void print() const
  {
    printf("[%+f,%+f]",a,b);
  }
};

class IntervalArray
{

  Interval*vec;
  int size;

 private:
  void operator=(const IntervalArray&); // disable assignment
  IntervalArray(const IntervalArray&); // disable copy constructor

 public:
  IntervalArray(int _s1):size(_s1)
  {
    vec = new Interval[size];
  }
  ~IntervalArray()
  {
    delete [] vec;
  }
  const Interval& get(int n1)const
  {
    assert_ge(n1,0); assert_lt(n1,size);
    return vec[n1];
  }
  Interval& fetch(int n1)
  {
    assert_ge(n1,0); assert_lt(n1,size);
    return vec[n1];
  }
  void set(int n1, const Interval& p)
  {
    assert_ge(n1,0); assert_lt(n1,size);
    vec[n1] = p;
  }

};

inline Interval operator+(const Interval& lhs, double rhs)
{
  Interval temp = lhs;
  temp+=rhs;
  return temp;
}

inline Interval operator+(const Interval& lhs, const Interval& rhs)
{
  Interval temp = lhs;
  //temp+=rhs;
  temp.add(rhs);
  return temp;
}

inline Interval operator-(const Interval& lhs, double rhs)
{
  Interval temp = lhs;
  temp.subtract(rhs);
  return temp;
}

inline Interval operator-(const Interval& lhs, const Interval& rhs)
{
  Interval temp = lhs;
  temp.subtract(rhs);
  return temp;
}

inline Interval operator*(double lhs, const Interval& rhs)
{
  return rhs*lhs;
}

inline Interval operator/(const Interval& lhs, const Interval& rhs)
{
  return lhs*rhs.reciprocal();
  // Interval temp = lhs;
  // temp.divide(rhs);
  // return temp;
}

#endif
