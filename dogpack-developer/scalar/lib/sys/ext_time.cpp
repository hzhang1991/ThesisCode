#include "ext_time.h"

timeval get_utime()
{
    // struct rusage r_usage; getrusage(RUSAGE_SELF,&r_usage);
    // return r_usage.ru_utime;
    timeval ret;
    gettimeofday(&ret,NULL);
    return ret;
};

timeval timeval_subtract (const timeval &x, const timeval &y_in)
{
    /* Perform the carry for the later subtraction by updating y. */
    timeval y = y_in;
    if (x.tv_usec < y.tv_usec) {
        int nsec = (y.tv_usec - x.tv_usec) / 1000000 + 1;
        y.tv_usec -= 1000000 * nsec;
        y.tv_sec += nsec;
    }
    if (x.tv_usec - y.tv_usec > 1000000) {
        int nsec = (x.tv_usec - y.tv_usec) / 1000000;
        y.tv_usec += 1000000 * nsec;
        y.tv_sec -= nsec;
    }

    /* Compute the time remaining to wait.
       tv_usec is certainly positive. */
    timeval result;
    result.tv_sec  = x.tv_sec  - y.tv_sec;
    result.tv_usec = x.tv_usec - y.tv_usec;
    return result;
}

double timeval2double(const timeval &in)
{
    double retval = in.tv_sec + in.tv_usec/1000000.;
    return retval;
}

double timeval_diff(const timeval &x, const timeval &y)
{
    return timeval2double(timeval_subtract(x, y));
}

