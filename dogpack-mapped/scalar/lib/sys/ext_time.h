#ifndef EXT_TIME_H
#define EXT_TIME_H
// conceptually this header extends time.h
#include <sys/time.h>     /* for getrusage() */
#include <sys/resource.h> /* for getrusage() */
#include <stddef.h>       /* for NULL */

timeval get_utime();
double timeval_diff(const timeval &x, const timeval &y);
#endif // EXT_TIME_H
