#ifndef time_commands_h
#define time_commands_h

#include <ext_time.h>

// one can use this to report execution time of a command

#define time_test(args...) \
{ \
  double sum_time = 0.; \
  double sum_sqrs = 0.; \
  int num_samples = 6; \
  int num_iterations = 10000; \
  for(int i=1;i<=num_samples;i++) \
  { \
    volatile double b; \
    volatile double c=.3; \
    double one_over_max = 1./num_iterations; \
    timeval start_time = get_utime(); \
    volatile double a; \
    for(a=one_over_max;a<1.;a+=one_over_max) \
    { \
      /* do it a bunch of times so that the for loop
         overhead is relatively negligible. */ \
      {args;} \
      {args;} \
      {args;} \
      {args;} \
      {args;} \
      {args;} \
      {args;} \
      {args;} \
      {args;} \
      {args;} \
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

#endif // time_commands_h
