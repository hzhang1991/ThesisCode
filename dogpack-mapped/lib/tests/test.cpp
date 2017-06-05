// Generic test program for users to use to test snippets of code for DogPack
//
// To compile and run:
//
//   shell_prompt> make test
//   shell_prompt> ./test
//
#include <dogdefs.h>
#include <float.h>
#include <math.h>
//#include "ext_math.h"

#if 0
struct C
{
  ~C(){dprint("destructor called");}
};

void f()
{
  dprint("throwing an exception");
  throw 1;
}

void f2()
{
  C c;
  f();
}

void g()
{
  try{ f2(); }
  catch(int j){
    dprintf("Caught int %d",j);
  }
}
#endif

#if 0
#include <stdlib.h> // for rand, srand
#include <time.h> // to seed
double rand_double()
{
  double out;
  // produces a number in [0, 1)
  double scale_inv = 1./(RAND_MAX+1.);
  //out = double(rand())/(RAND_MAX+1.);
  out = double(rand())*scale_inv;
  return out;
}

void seed_rand()
{
  srand((unsigned)time(NULL)); 
}

void test_set_B_orth_system()
{
  seed_rand();
  const double B1 = rand_double();
  const double B2 = rand_double();
  const double B3 = rand_double();
  const double Bmag = sqrt(B1*B1+B2*B2+B3*B3);
  double a[4][4];
  set_B_orth_system(a,B1,B2,B3,Bmag);
}
#endif

#if 0
void test_inf_and_nan()
{
  volatile double z = 0.;
  volatile double o = 1.;
  volatile double t = 2.;
  volatile double n = z/z;
  volatile double p = 1./z;
  volatile double m = -1./z;
  dprint(n)
  dprint(n == n)
  dprint(n < n)
  dprint(isnan(z))
  dprint(isnan(t))
  dprint(isnan(n))
  dprint(isinf(-o/z))
  dprint(isinf(o/z))
  dprint(isinf(o/o))
  dprint(isinf(o))
  dprint(isinf(z))
  dprint(isnan(-o/z))
  dprint(isnan(o/z))
  dprint(isnan(o/o))
  dprint(isnan(o))
  dprint(isnan(z))
  dprint(isfinite(o))
  dprint(isfinite(z))
  dprint(isfinite(n))
  dprint(isfinite(p))
  dprint(isfinite(m))
  
  dprint(isinf(-1./0.));
  dprint(isinf(1./0.));
  dprint(isinf(1.));
  dprint(isinf(0./0.));
  dprint(0./0.);
  dprint(0./0. > 3.);
  dprint(1./0. <= 2*DBL_MAX);
  dprint(0./0. <= 2*DBL_MAX);
  dprint(n <= 2*DBL_MAX);
  dprint(isnan(0./0.));
  dprint(isnan(1./0.));
  dprint(isnan(-1./0.));
  dprint(isnan(1.));
  //dprint(isunordered(n,o))
}
#endif

// probably this is compiler-dependent,
// since inline is supposed to be just a hint
inline void* myalloca(size_t size)
{
    return alloca(size);
}

// there does not seem to be any mechanism to allocat
// the data of a class from the stack that hides the
// details from the user.
void test_alloc_new()
{
  int i[6];
  //int* j = new int[6];
  // this does not do what I want
  int* j1 = (int*) myalloca(sizeof(int)*6);
  int* j2 = (int*) alloca(sizeof(int)*6);
  volatile int s = 6;
  int* k = new int[s];
  int m[6];
  printf("i =%p\n",i);
  printf("j1=%p\n",j1);
  printf("j2=%p\n",j2);
  printf("k =%p\n",k);
  printf("m =%p\n",m);
  delete [] k;
}

int main()
{
  test_alloc_new();
  //test_inf_and_nan();
  //test_member_reference();
  //g();
  //test_set_B_orth_system();
}
