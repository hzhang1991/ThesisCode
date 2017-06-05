#include <assert.h>
#include <math.h>
#include <ext_time.h>
#include <float.h>
#include <limits>
#include <complex>

#define make_test(args...) \
{ \
  double sum_time = 0.; \
  double sum_sqrs = 0.; \
  int num_samples = 6; \
  int num_iterations = 10000; \
  for(int i=1;i<=num_samples;i++) \
  { \
    volatile bool B; \
    volatile double b; \
    volatile double c=.3; \
    double one_over_max = 1./num_iterations; \
    timeval start_time = get_utime(); \
    volatile double a; \
    for(a=one_over_max;a<1.;a+=one_over_max) \
    { \
      /* do it a bunch of times so that the for loop
         overhead is relatively negligible. */ \
      args; \
      args; \
      args; \
      args; \
      args; \
      args; \
      args; \
      args; \
      args; \
      args; \
    } \
    timeval end_time = get_utime(); \
    double diff_utime = timeval_diff(end_time, start_time); \
    sum_time += diff_utime; \
    sum_sqrs += pow(diff_utime,2); \
    /*printf(" utime = %.5f", diff_utime);*/ \
    /*printf(" for %s\n",#args);*/ \
  } \
  double avg_time = sum_time/num_samples; \
  double var_time = (sum_sqrs - sum_time*avg_time)/(num_samples-1); \
  double std_time = sqrt(var_time); \
  double percent = 100*std_time/avg_time; \
  printf(" utime = %7.3f +/- %05.2f%%", 1000*avg_time, percent); \
  printf(" for %s\n",#args); \
}

void test_double_math()
{
  printf("\n");
  printf("testing powers\n");
  printf("\n");
  make_test(b = pow(a,0.6));
  make_test(b = exp(log(a)*0.6));
  make_test(b = pow(a,0.5));
  make_test(b = sqrt(a));
  make_test(b = cbrt(a));
  printf("\n");
  make_test(b = exp(log(a)*5.));
  make_test(b = pow(a,5)); // IEEE: takes 10 times as long
  make_test(b = a*a*a*a*a); // IEEE forces 5 multiplications
  make_test(b = (a*a)*(a*a)*a); // this is done with 4 multiplications
  make_test(b = pow(pow(a,2),2)*a);
  printf("\n");
  make_test(b = exp(log(a)*3.));
  make_test(b = pow(a,3));
  make_test(b = a*a*a);
  printf("\n");
  make_test(b = exp(log(a)*2.));
  make_test(b = pow(a,2));
  make_test(b = a*a);
  make_test(b = a*c);
  make_test(b = a*a);
  make_test(b = c*a);
  make_test(b = 2.718281828*a);
  make_test(b = 2*a);
  make_test(b = 0.5*a);
  make_test(b = a+a);
  make_test(b = a+c);
  make_test(b = a+a);
  make_test(b = a);
  printf("\n");
  printf("testing division\n");
  printf("\n");
  // division takes 7 times as long as mult
  // unless compiling with -funsafe-math-optimizations
  make_test(b = a/c);
  make_test(b = a/3.);
  make_test(b = a*(1./3.));
  printf("\n");
  printf("timing elementary functions\n");
  printf("\n");
  make_test(b=log(a));
  make_test(b=exp(a));
  make_test(b=cos(a));
  make_test(b=sin(a));
  make_test(b=cosh(a));
  make_test(b=sinh(a));
  make_test(b=tan(a));
  make_test(b=acos(a));
  make_test(b=asin(a));
  make_test(b=atan(a));
  make_test(b=sqrt(a));
  printf("\n");
  // these checks are five to ten times more expensive
  make_test(B=isinf(a));
  make_test(B=isnan(a));
  make_test(B=std::isinf(a));
  make_test(B=std::isnan(a));
  make_test(B=(a>=-std::numeric_limits<double>::infinity() && a<=std::numeric_limits<double>::infinity()));
  // these are much cheaper than using standard library routines
  make_test(B=!(a==a));
  // supported after gcc 4.2, I think.
  //#pragma GCC diagnostic push
  //#pragma GCC diagnostic ignored "-Wno-div-by-zero"
  //make_test(B=(a>=-1./0. && a<=1./0.));
  make_test(B=(a>=-2*DBL_MAX && a<=2*DBL_MAX));
  //#pragma GCC diagnostic pop
  make_test(B=(a>=-DBL_MAX && a<=DBL_MAX));
  printf("\n");
}
  
int main()
{
  test_double_math();
}
