
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifndef CHECK_BOUNDS
#define CHECK_BOUNDS
#endif
#ifndef NDIMS
#define NDIMS 0
#endif
#include "tensors.h"
#include "assert.h"
#include "debug.h"

//#define myprintf(e, args...) printf(e, ##args)
#define myprintf(e, args...) ((void)0)

double timeval_subtract (struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  struct timeval result;
  result.tv_sec = x->tv_sec - y->tv_sec;
  result.tv_usec = x->tv_usec - y->tv_usec;

  double retval = result.tv_sec + result.tv_usec/1000000.;
  return retval;
}

#define toc(task) report_time_fileLine(__FILE__,__LINE__,task);

static timeval start_utime;

timeval get_utime()
{
  struct rusage r_usage; getrusage(RUSAGE_SELF,&r_usage);
  return r_usage.ru_utime;
};

void tic()
{
  // struct rusage r_usage; getrusage(RUSAGE_SELF,&r_usage);
  // start_utime = r_usage.ru_utime;
  start_utime = get_utime();
}

void report_time_fileLine(const char *filename, int line, const char* task)
{
  // struct rusage r_usage; getrusage(RUSAGE_SELF,&r_usage);
  // timeval utime = r_usage.ru_utime;
  timeval utime = get_utime();
  double change_utime = timeval_subtract(&utime,&start_utime);
  printf("  time elapsed: %12.6f for %-18s at %s:%d\n",
    change_utime, task, filename, line);
}

int b1=2;
int b2=3;
int b3=-1;
int b4=-2;
int b5=1;
int b6=5;

#define mult 10
int s1=mult+1;
int s2=mult+2;
int s3=mult+3;
int s4=mult+4;
int s5=mult+5;
int s6=mult+5;

int e1=b1+s1-1;
int e2=b2+s2-1;
int e3=b3+s3-1;
int e4=b4+s4-1;
int e5=b5+s5-1;
int e6=b6+s6-1;

// === section: multidimensional tensor base classes ===

void test_dTensor6d()
{
    printf("=== testing dTensor6d ===\n");
    dTensor6d t(s1, s2, s3, s4, s5, s6,
                b1, b2, b3, b4, b5, b6);
    dTensor6d u(s1, s2, s3, s4, s5, s6,
                b1, b2, b3, b4, b5, b6);
    t.setall(0.);
    u.setall(0.);

    // reversing order increases time to access idx by 5 times
    // (.055 vs. .220)
    tic();
    int idx = 0;
    for(int i6=b6;i6<=e6;i6++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5,i6);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx reversed");

    tic();
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i6=b6;i6<=e6;i6++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5,i6);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx");

    // this initialization requires 1.12 s
    // (4.25 if order is reversed)
    tic();
    double a=0.;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i6=b6;i6<=e6;i6++)
    {
      t.set(i1,i2,i3,i4,i5,i6, a++);
    }
    myprintf("  a=%f\n",a);
    myprintf("  t.numel()=%d\n",t.numel());
    toc("t.set");

    // this initialization requires 0.82 s
    tic();
    double a2=0.;
    for(int i=0; i<u.numel(); i++) u.vset(i,a2++);
    myprintf("  a2=%f\n",a2);
    toc("u.vset");

    tic();
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i6=b6;i6<=e6;i6++)
    {
      assert(t.get(i1,i2,i3,i4,i5,i6)==u.get(i1,i2,i3,i4,i5,i6));
    }
    toc("t.get==u.get");

    tic();
    //for(int k=0; k<u.numel(); k++) assert_eq(t.vget(k),u.vget(k));
    myprintf("  a2=%f\n",a2);
    toc("t.vget==u.vget");

    // this initialization requires 1.12 s
    // (4.25 if order is reversed)
    tic();
    a=0.;
    for(int i6=b6;i6<=e6;i6++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      // for some bizarre reason this is faster
      // when I increment a than when I don't.
      t.set(i1,i2,i3,i4,i5,i6, a++);
    }
    myprintf("  a=%f\n",a);
    toc("t.set reversed");

    // requires .05 s, .66 s if order is reversed
    tic();
    double vb;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i6=b6;i6<=e6;i6++)
    {
      vb=t.get(i1,i2,i3,i4,i5,i6);
    }
    toc("t.get");

    // why does reversing the order seem to make no difference?
    // requires .05 s, .66 s if order is reversed
    tic();
    double va;
    for(int i6=b6;i6<=e6;i6++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      va=t.get(i1,i2,i3,i4,i5,i6);
    }
    toc("t.get reversed");

    // requires 0.04 s
    tic();
    double vc;
    for(int k=0;k<t.numel();k++)
    {
      vc+=t.vget(k);
    }
    toc("t.vget");

    // for some reason a single loop with the
    // same number of iterations requires half
    // as much time as nested loops (.036 vs. .076)
    tic();
    int ja=0;
    int jb=0;
    for(int k=0;k<t.numel();k++)
    {
      ja++;
      jb+=ja;
    }
    myprintf("  jb=%d\n",jb);
    toc("simple loop");

    tic();
    int jc=0;
    int jd=0;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    for(int i6=b6;i6<=e6;i6++)
    {
      jc++;
      jd+=jc;
    }
    myprintf("  jd=%d\n",jd);
    toc("full loop");
}

void test_dTensor5d()
{
    printf("=== testing dTensor5d ===\n");
    dTensor5d t(s1, s2, s3, s4, s5, 
                b1, b2, b3, b4, b5);
    dTensor5d u(s1, s2, s3, s4, s5, 
                b1, b2, b3, b4, b5);
    t.setall(0.);
    u.setall(0.);

    // reversing order increases time to access idx by 5 times
    // (.055 vs. .220)
    tic();
    int idx = 0;
    for(int i5=b5;i5<=e5;i5++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx reversed");

    tic();
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx");

    // this initialization requires 1.12 s
    // (4.25 if order is reversed)
    tic();
    double a=0.;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    {
      t.set(i1,i2,i3,i4,i5, a++);
      // u.set(i1,i2,i3,i4,i5, t.get(i1,i2,i3,i4,i5));
    }
    myprintf("  a=%f\n",a);
    myprintf("  t.numel()=%d\n",t.numel());
    toc("t.set");

    // this initialization requires 0.82 s
    tic();
    double a2=0.;
    for(int i=0; i<u.numel(); i++) u.vset(i,a2++);
    myprintf("  a2=%f\n",a2);
    toc("u.vset");

    tic();
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    {
      assert(t.get(i1,i2,i3,i4,i5)==u.get(i1,i2,i3,i4,i5));
    }
    toc("t.get==u.get");

    tic();
    //for(int k=0; k<u.numel(); k++) assert_eq(t.vget(k),u.vget(k));
    myprintf("  a2=%f\n",a2);
    toc("t.vget==u.vget");

    // this initialization requires 1.12 s
    // (4.25 if order is reversed)
    tic();
    a=0.;
    for(int i5=b5;i5<=e5;i5++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      // for some bizarre reason this is faster
      // when I increment a than when I don't.
      t.set(i1,i2,i3,i4,i5, a++);
    }
    myprintf("  a=%f\n",a);
    toc("t.set reversed");

    // requires .05 s, .66 s if order is reversed
    tic();
    double vb;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    {
      vb=t.get(i1,i2,i3,i4,i5);
    }
    toc("t.get");

    // why does reversing the order seem to make no difference?
    // requires .05 s, .66 s if order is reversed
    tic();
    double va;
    for(int i5=b5;i5<=e5;i5++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      va=t.get(i1,i2,i3,i4,i5);
    }
    toc("t.get reversed");

    // requires 0.04 s
    tic();
    double vc;
    for(int k=0;k<t.numel();k++)
    {
      vc+=t.vget(k);
    }
    toc("t.vget");

    // for some reason a single loop with the
    // same number of iterations requires half
    // as much time as nested loops (.036 vs. .076)
    tic();
    int ja=0;
    int jb=0;
    for(int k=0;k<t.numel();k++)
    {
      ja++;
      jb+=ja;
    }
    myprintf("  jb=%d\n",jb);
    toc("simple loop");

    tic();
    int jc=0;
    int jd=0;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    for(int i5=b5;i5<=e5;i5++)
    {
      jc++;
      jd+=jc;
    }
    myprintf("  jd=%d\n",jd);
    toc("full loop");
}

void test_dTensor4d()
{
    printf("=== testing dTensor4d ===\n");
    dTensor4d t(s1, s2, s3, s4, 
                b1, b2, b3, b4);
    dTensor4d u(s1, s2, s3, s4, 
                b1, b2, b3, b4);

    tic();
    int ja=0;
    int jb=0;
    for(int k=0;k<t.numel();k++)
    {
      ja++;
      jb+=ja;
    }
    myprintf("  jb=%d\n",jb);
    toc("simple loop");

    tic();
    int jc=0;
    int jd=0;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    {
      jc++;
      jd+=jc;
    }
    myprintf("  jd=%d\n",jd);
    toc("full loop");


    tic();
    int idx = 0;
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      idx = t.getidx(i1,i2,i3,i4);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx reversed");

    tic();
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    {
      idx = t.getidx(i1,i2,i3,i4);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx");

    tic();
    double a=0.;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    {
      t.set(i1,i2,i3,i4, a++);
    }
    myprintf("  a=%f\n",a);
    myprintf("  t.numel()=%d\n",t.numel());
    toc("t.set");

    // this initialization requires 0.82 s
    tic();
    double a2=0.;
    for(int i=0; i<u.numel(); i++) u.vset(i,a2++);
    myprintf("  a2=%f\n",a2);
    toc("u.vset");

    tic();
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    {
      assert(t.get(i1,i2,i3,i4)==u.get(i1,i2,i3,i4));
    }
    toc("t.get==u.get");

    tic();
    //for(int k=0; k<u.numel(); k++) assert(t.vget(k)==u.vget(k));
    myprintf("  a2=%f\n",a2);
    toc("t.vget==u.vget");

    // this initialization requires 1.12 s
    // (4.25 if order is reversed)
    tic();
    a=0.;
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      // for some bizarre reason this is faster
      // when I increment a than when I don't.
      t.set(i1,i2,i3,i4, a++);
    }
    myprintf("  a=%f\n",a);
    toc("t.set reversed");

    // requires .05 s, .66 s if order is reversed
    tic();
    double vb;
    for(int i1=b1;i1<=e1;i1++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i4=b4;i4<=e4;i4++)
    {
      vb=t.get(i1,i2,i3,i4);
    }
    toc("t.get");

    // why does reversing the order seem to make no difference?
    // requires .05 s, .66 s if order is reversed
    tic();
    double va;
    for(int i4=b4;i4<=e4;i4++)
    for(int i3=b3;i3<=e3;i3++)
    for(int i2=b2;i2<=e2;i2++)
    for(int i1=b1;i1<=e1;i1++)
    {
      va=t.get(i1,i2,i3,i4);
    }
    toc("t.get reversed");

    // requires 0.04 s
    tic();
    double vc;
    for(int k=0;k<t.numel();k++)
    {
      vc+=t.vget(k);
    }
    toc("t.vget");
}

// === section: multidimensional 1-based tensor classes ===

void test_dTensor5()
{
    printf("=== testing dTensor5 ===\n");
    dTensor5 t(s1, s2, s3, s4, s5);
    dTensor5 u(s1, s2, s3, s4, s5);

    // reversing order increases time to access idx by 5 times
    // (.055 vs. .220)
    tic();
    int idx = 0;
    for(int i5=1;i5<=s5;i5++)
    for(int i4=1;i4<=s4;i4++)
    for(int i3=1;i3<=s3;i3++)
    for(int i2=1;i2<=s2;i2++)
    for(int i1=1;i1<=s1;i1++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx reversed");

    tic();
    for(int i1=1;i1<=s1;i1++)
    for(int i2=1;i2<=s2;i2++)
    for(int i3=1;i3<=s3;i3++)
    for(int i4=1;i4<=s4;i4++)
    for(int i5=1;i5<=s5;i5++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx");

    // this initialization requires 1.12 s
    // (4.25 if order is reversed)
    tic();
    double a=0.;
    for(int i1=1;i1<=s1;i1++)
    for(int i2=1;i2<=s2;i2++)
    for(int i3=1;i3<=s3;i3++)
    for(int i4=1;i4<=s4;i4++)
    for(int i5=1;i5<=s5;i5++)
    {
      t.set(i1,i2,i3,i4,i5, a++);
    }
    myprintf("  a=%f\n",a);
    myprintf("  t.numel()=%d\n",t.numel());
    toc("t.set");

    // this initialization requires 0.82 s
    tic();
    double a2=0.;
    for(int i=0; i<u.numel(); i++) u.vset(i,a2++);
    myprintf("  a2=%f\n",a2);
    toc("u.vset");

    tic();
    for(int i1=1;i1<=s1;i1++)
    for(int i2=1;i2<=s2;i2++)
    for(int i3=1;i3<=s3;i3++)
    for(int i4=1;i4<=s4;i4++)
    for(int i5=1;i5<=s5;i5++)
    {
      assert(t.get(i1,i2,i3,i4,i5)==u.get(i1,i2,i3,i4,i5));
    }
    toc("t.get==u.get");

    tic();
    //for(int k=0; k<u.numel(); k++) assert(t.vget(k)==u.vget(k));
    myprintf("  a2=%f\n",a2);
    toc("t.vget==u.vget");

    // this initialization requires 1.12 s
    // (4.25 if order is reversed)
    tic();
    a=0.;
    for(int i5=1;i5<=s5;i5++)
    for(int i4=1;i4<=s4;i4++)
    for(int i3=1;i3<=s3;i3++)
    for(int i2=1;i2<=s2;i2++)
    for(int i1=1;i1<=s1;i1++)
    {
      // for some bizarre reason this is faster
      // when I increment a than when I don't.
      t.set(i1,i2,i3,i4,i5, a++);
    }
    myprintf("  a=%f\n",a);
    toc("t.set reversed");

    // requires .05 s, .66 s if order is reversed
    tic();
    double vb;
    for(int i1=1;i1<=s1;i1++)
    for(int i2=1;i2<=s2;i2++)
    for(int i3=1;i3<=s3;i3++)
    for(int i4=1;i4<=s4;i4++)
    for(int i5=1;i5<=s5;i5++)
    {
      vb=t.get(i1,i2,i3,i4,i5);
    }
    toc("t.get");

    // why does reversing the order seem to make no difference?
    // requires .05 s, .66 s if order is reversed
    tic();
    double va;
    for(int i5=1;i5<=s5;i5++)
    for(int i4=1;i4<=s4;i4++)
    for(int i3=1;i3<=s3;i3++)
    for(int i2=1;i2<=s2;i2++)
    for(int i1=1;i1<=s1;i1++)
    {
      va=t.get(i1,i2,i3,i4,i5);
    }
    toc("t.get reversed");

    // requires 0.04 s
    tic();
    double vc;
    for(int k=0;k<t.numel();k++)
    {
      vc+=t.vget(k);
    }
    toc("t.vget");
}

// === section: multidimensional boundary condition (BC) tensor classes ===

void test_dTensorBC5()
{
    printf("=== testing dTensorBC5 ===\n");
    int mbc=2;
    int S1=s1-2*mbc;
    int S2=s2-2*mbc;
    int S3=s3-2*mbc;
    int S4=s4;
    int S5=s5;

    dTensorBC5 t(S1, S2, S3, S4, S5, mbc, 3);
    dTensorBC5 u(S1, S2, S3, S4, S5, mbc, 3);

    tic();
    int idx = 0;
    for(int i5=1;i5<=S5;i5++)
    for(int i4=1;i4<=S4;i4++)
    for(int i3=1-mbc;i3<=S3+mbc;i3++)
    for(int i2=1-mbc;i2<=S2+mbc;i2++)
    for(int i1=1-mbc;i1<=S1+mbc;i1++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx reversed");

    tic();
    for(int i1=1-mbc;i1<=S1+mbc;i1++)
    for(int i2=1-mbc;i2<=S2+mbc;i2++)
    for(int i3=1-mbc;i3<=S3+mbc;i3++)
    for(int i4=1;i4<=S4;i4++)
    for(int i5=1;i5<=S5;i5++)
    {
      idx = t.getidx(i1,i2,i3,i4,i5);
    }
    myprintf("  idx=%d\n",idx);
    toc("t.getidx");

    tic();
    double a=0.;
    for(int i1=1-mbc;i1<=S1+mbc;i1++)
    for(int i2=1-mbc;i2<=S2+mbc;i2++)
    for(int i3=1-mbc;i3<=S3+mbc;i3++)
    for(int i4=1;i4<=S4;i4++)
    for(int i5=1;i5<=S5;i5++)
    {
      t.set(i1,i2,i3,i4,i5, a++);
    }
    myprintf("  a=%f\n",a);
    myprintf("  t.numel()=%d\n",t.numel());
    toc("t.set");

    tic();
    double a2=0.;
    for(int i=0; i<u.numel(); i++) u.vset(i,a2++);
    myprintf("  a2=%f\n",a2);
    toc("u.vset");

    tic();
    for(int i1=1-mbc;i1<=S1+mbc;i1++)
    for(int i2=1-mbc;i2<=S2+mbc;i2++)
    for(int i3=1-mbc;i3<=S3+mbc;i3++)
    for(int i4=1;i4<=S4;i4++)
    for(int i5=1;i5<=S5;i5++)
    {
      assert(t.get(i1,i2,i3,i4,i5)==u.get(i1,i2,i3,i4,i5));
    }
    toc("t.get==u.get");

    tic();
    //for(int k=0; k<u.numel(); k++) assert(t.vget(k)==u.vget(k));
    myprintf("  a2=%f\n",a2);
    toc("t.vget==u.vget");

    tic();
    a=0.;
    for(int i5=1;i5<=S5;i5++)
    for(int i4=1;i4<=S4;i4++)
    for(int i3=1-mbc;i3<=S3+mbc;i3++)
    for(int i2=1-mbc;i2<=S2+mbc;i2++)
    for(int i1=1-mbc;i1<=S1+mbc;i1++)
    {
      t.set(i1,i2,i3,i4,i5, a++);
    }
    myprintf("  a=%f\n",a);
    toc("t.set reversed");

    // requires .05 s, .66 s if order is reversed
    tic();
    double vb;
    for(int i1=1-mbc;i1<=S1+mbc;i1++)
    for(int i2=1-mbc;i2<=S2+mbc;i2++)
    for(int i3=1-mbc;i3<=S3+mbc;i3++)
    for(int i4=1;i4<=S4;i4++)
    for(int i5=1;i5<=S5;i5++)
    {
      vb=t.get(i1,i2,i3,i4,i5);
    }
    toc("t.get");

    // why does reversing the order seem to make no difference?
    tic();
    double va;
    for(int i5=1;i5<=S5;i5++)
    for(int i4=1;i4<=S4;i4++)
    for(int i3=1-mbc;i3<=S3+mbc;i3++)
    for(int i2=1-mbc;i2<=S2+mbc;i2++)
    for(int i1=1-mbc;i1<=S1+mbc;i1++)
    {
      va=t.get(i1,i2,i3,i4,i5);
    }
    toc("t.get reversed");

    // requires 0.04 s
    tic();
    double vc;
    for(int k=0;k<t.numel();k++)
    {
      vc+=t.vget(k);
    }
    toc("t.vget");
}

void test_dTensor1()
{
  int s1=4;
  dTensor1 t1(s1);
  for(int i=1; i<t1.getsize(); i++) t1.set(i,i);
  for(int i=1; i<t1.getsize(); i++) printf("t1.get(i)=%f\n",t1.get(i));
}

void test_uninitialized()
{
    dTensorBase t(20);
    dprint(t.vget(3));
}

int main()
{

    test_dTensor6d();

    test_dTensor5d();
    test_dTensor5();
    test_dTensorBC5();

    test_dTensor4d();

    test_uninitialized();

}
