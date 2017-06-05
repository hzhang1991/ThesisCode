#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifndef CHECK_BOUNDS
#define CHECK_BOUNDS
#endif
#ifndef NDIMS
#define NDIMS 4
#endif
#include "tensors.h"
#include "assert.h"
#include "debug.h"

//#define myprintf(e, args...) printf(e, ##args)
#define myprintf(e, args...) ((void)0)

int e1   = 3;
int e2   = 2;
int e3   = 1;
int e4   = 3;
int e5   = 1;
int e6   = 3;
int mbc  = 1;

void test_dTensorBC6()
{

    printf("ndims = %d\n", NDIMS);
    printf("=== testing dTensorBC6 ===\n");
    dTensorBC6   t(e1,e2,e3,e4,e5,e6,mbc);
    dTensor6     u(e1,e2,e3,e4,e5,e6);
    t.setall(0.);
    u.setall(0.);

    int jc=0;
    for(int i1=1-mbc;  i1<=e1+mbc; i1++)
    for(int i2=1-mbc;  i2<=e2+mbc; i2++)
    for(int i3=1-mbc;  i3<=e3+mbc; i3++)
    for(int i4=1-mbc;  i4<=e4+mbc; i4++)
    for(int i5=1;      i5<=e5;     i5++)
    for(int i6=1;      i6<=e6;     i6++)
    {
        jc++;
        t.set(i1,i2,i3,i4,i5,i6, jc);
    }

    jc=0;
    for(int i1=1;  i1<=e1; i1++)
    for(int i2=1;  i2<=e2; i2++)
    for(int i3=1;  i3<=e3; i3++)
    for(int i4=1;  i4<=e4; i4++)
    for(int i5=1;  i5<=e5; i5++)
    for(int i6=1;  i6<=e6; i6++)
    {
        jc++;
        u.set(i1,i2,i3,i4,i5,i6, jc);
    }


    const int numelbc = t.numel();
    printf("numelbc = %d\n", numelbc);
    for( int k=0; k < numelbc; k++ )
    {
        printf("t = %f\n", t.vget(k) );
    }

    const int numelnbc = u.numel();
    printf("numelnbc = %d\n", numelnbc);
    for( int k=0; k < numelnbc; k++ )
    {
        printf("u = %f\n", u.vget(k) );
    }

}

int main()
{

    test_dTensorBC6();

}
