// turn on CHECK_BOUNDS before including tensors.h
// so that the non-inline versions of get and set are declared.
// This way they are implemented in this file so that
// elsewhere in the code one can choose whether to use
// the in-line or bounds-checking versions of get and set
#ifndef CHECK_BOUNDS
#define CHECK_BOUNDS
#endif
#include "tensors.h"
#include "debug.h"
#include "assert.h"
#include<new>
#define bounds_check() expensive_bounds_check();
//#define bounds_check() cheap_bounds_check();
#define PARALLEL_SIZE 512
#include <limits> // for numeric_limits

// class iTensor2

iTensor2::iTensor2(int numRows, int numColumns)
{
    assert_ge(numRows,0);
    assert_ge(numColumns,0);

    size = numRows*numColumns;
    vec = new int[size];
    rows = numRows;
    columns = numColumns;
}

iTensor2::iTensor2(const iTensor2& anotheriTensor2)
{
    rows = anotheriTensor2.rows;
    columns = anotheriTensor2.columns;
    size = anotheriTensor2.size;
    vec = new int[size];

    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor2.vec[i]; }
    else
#pragma omp parallel for
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor2.vec[i]; }
}

iTensor2::~iTensor2()
{
    delete[] vec;
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_gt(n1,0); assert_le(n1,rows); \
assert_gt(n2,0); assert_le(n2,columns); 
#define cheap_bounds_check() \
    assert_printf(k<size && k>=0, "n1=%d, n2=%d",n1,n2);

const int& iTensor2::get(int n1, int n2) const
{
    int k = (n1-1)*columns + (n2-1);
    bounds_check();
    return vec[k];
}

void iTensor2::set(int n1, int n2, int value)
{
    int k = (n1-1)*columns + (n2-1);
    bounds_check();
    vec[k] = value;
}

// class iTensor3

iTensor3::iTensor3(int n1, int n2, int n3)
{
    assert_printf(n1>0 && n2>0 && n3>0, "n1=%d, n2=%d, n3=%d",n1,n2,n3);
    size = n1*n2*n3;
    vec = new int[size];
    numElem1 = n1;
    numElem2 = n2;
    numElem3 = n3;
}

iTensor3::iTensor3(const iTensor3& anotheriTensor3)
{
    numElem1 = anotheriTensor3.numElem1;
    numElem2 = anotheriTensor3.numElem2;
    numElem3 = anotheriTensor3.numElem3;
    size = anotheriTensor3.size;
    vec = new int[size];

    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor3.vec[i]; }
    else
#pragma omp parallel for
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor3.vec[i]; }
}

iTensor3::~iTensor3()
{
    delete[] vec;
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_gt(n1,0); assert_le(n1,numElem1); \
assert_gt(n2,0); assert_le(n2,numElem2); \
assert_gt(n3,0); assert_le(n3,numElem3); 
#define cheap_bounds_check() \
    assert_printf(k<size && k>=0, "n1=%d, n2=%d, n3=%d",n1,n2,n3);

const int& iTensor3::get(int n1, int n2, int n3) const
{
    int k = ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1);
    bounds_check();
    return vec[k];
}

void iTensor3::set(int n1, int n2, int n3, int value)
{
    int k = ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1);
    bounds_check();
    vec[k] = value;
}

// === major section: itensors (arrays of int) ===

// class iTensorBase

void iTensorBase::init()
{
    assert_printf(size>=0, "size=%d", size);
    vec = new int[size];
}

iTensorBase::iTensorBase(const iTensorBase& in) :
    size(in.size)
{
    vec = new int[size];
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++) vec[i] = in.vec[i];
    else
#pragma omp parallel for
        for (int i=0; i<size; i++) vec[i] = in.vec[i];
}

iTensorBase::~iTensorBase()
{
    delete[] vec;
}

const int& iTensorBase::vget(int k) const
{
    assert_ge(k,0); assert_lt(k,size);
    return vec[k];
}

int& iTensorBase::vfetch(int k)
{
    assert_ge(k,0); assert_lt(k,size);
    return vec[k];
}

void iTensorBase::vset(int k, int value)
{
    assert_ge(k,0); assert_lt(k,size);
    vec[k]=value;
}

void iTensorBase::setall(int value)
{
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++) vec[i] = value;
    else
#pragma omp parallel for
        for (int i=0; i<size; i++) vec[i] = value;
}

// === major section: dtensors (arrays of double) ===

// class dTensorBase

void dTensorBase::init()
{
    assert_printf(size>=0, "size=%d", size);
    vec = new double[size];
#ifdef CHECK_INIT
    setall(std::numeric_limits<double>::quiet_NaN());
    // setall(0./0); // initialize all entries to nan
#endif
}

dTensorBase::dTensorBase(const dTensorBase& in, CopyMode::Enum copyMode) :
    size(in.size)
{
    switch(copyMode)
    {
        // should make vec counted_ptr first
        //case CopyMode::SHALLOW:
        //  vec(in.vec);
        //  break;
        case CopyMode::DIMS:
            vec = new double[size];
            break;
        case CopyMode::DEEP:
            vec = new double[size];
            copyfrom(in);
            break;
        default:
            unsupported_value_error(copyMode);
    }
}

void dTensorBase::copyfrom(const dTensorBase& in)
{
    if(this!=&in)
    {
        assert_eq(size,in.size);
        if(size<PARALLEL_SIZE)
            for (int i=0; i<size; i++) vec[i] = in.vec[i];
        else
#pragma omp parallel for
            for (int i=0; i<size; i++) vec[i] = in.vec[i];
    }
}

dTensorBase::~dTensorBase()
{
    delete[] vec;
}

const double& dTensorBase::vget(int k) const
{
    assert_ge(k,0); assert_lt(k,size);
#ifdef CHECK_INIT
    if(!(vec[k]==vec[k])) eprintf("vec[%d]=%24.16e",k,vec[k]);
    // double ret = vec[k];
    // assert_eq(ret,ret);
#endif

    return vec[k];
}

double& dTensorBase::vfetch(int k)
{
    assert_ge(k,0); assert_lt(k,size);
#ifdef CHECK_INIT
    if(!(vec[k]==vec[k])) eprintf("vec[%d]=%24.16e",k,vec[k]);
#endif
    return vec[k];
}

void dTensorBase::vset(int k, double value)
{
    assert_ge(k,0); assert_lt(k,size);
#ifdef CHECK_INIT
    if(!(value==value)) eprintf("for k=%d value=%24.16e",k,value);
#endif
    vec[k]=value;
}

void dTensorBase::setall(double value)
{
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++) vec[i] = value;
    else
#pragma omp parallel for
        for (int i=0; i<size; i++) vec[i] = value;
}

void dTensorBase::check()
{
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++){
            assert_printf(vec[i]==vec[i],"vec[%d]=%24.16e",i,vec[i]);
        }
    else
#pragma omp parallel for
        for (int i=0; i<size; i++){
            assert_printf(vec[i]==vec[i],"vec[%d]=%24.16e",i,vec[i]);
        }
}

// --- section: multidimensional tensor base classes ---

// class dTensor6d

// check parameters and update derived parameters.
void dTensor6d::init()
{
    assert_printf(s1>=0 && s2>=0 && s3>=0 && s4>=0 && s5>=0,
            "b1=%d, b2=%d, b3=%d, b4=%d, b5=%d, b6=%d,"
            "s1=%d, s2=%d, s3=%d, s4=%d, s5=%d, s6=%d",
            b1,b2,b3,b4,b5,b6,
            s1,s2,s3,s4,s5,s6);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    e3 = b3 + s3 - 1;
    e4 = b4 + s4 - 1;
    e5 = b5 + s5 - 1;
    e6 = b6 + s6 - 1;
    size = s1*s2*s3*s4*s5*s6;
    dTensorBase::init();
}

dTensor6d::dTensor6d(
        int s1i, int s2i, int s3i, int s4i, int s5i, int s6i,
        int b1i, int b2i, int b3i, int b4i, int b5i, int b6i ) :
        b1(b1i), b2(b2i), b3(b3i), b4(b4i), b5(b5i), b6(b6i),
        s1(s1i), s2(s2i), s3(s3i), s4(s4i), s5(s5i), s6(s6i)
{ init(); }

dTensor6d::dTensor6d(const dTensor6d& in, CopyMode::Enum copyMode) 
    : dTensorBase(in, copyMode),
    b1(in.b1), b2(in.b2), b3(in.b3), b4(in.b4), b5(in.b5), b6(in.b6),
    e1(in.e1), e2(in.e2), e3(in.e3), e4(in.e4), e5(in.e5), e6(in.e6),
    s1(in.s1), s2(in.s2), s3(in.s3), s4(in.s4), s5(in.s5), s6(in.s6) { }

    // alias for assert_eq
    //#define ae(a,b) assert_op(==,a,b);
#define ae assert_eq
void dTensor6d::copyfrom(const dTensor6d& in)
{
    ae(b1,in.b1); ae(b2,in.b2); ae(b3,in.b3); ae(b4,in.b4); ae(b5,in.b5); ae(b6,in.b6);
    ae(e1,in.e1); ae(e2,in.e2); ae(e3,in.e3); ae(e4,in.e4); ae(e5,in.e5); ae(e6,in.e6);
    ae(s1,in.s1); ae(s2,in.s2); ae(s3,in.s3); ae(s4,in.s4); ae(s5,in.s5); ae(s6,in.s6);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2); \
assert_ge(n3,b3); assert_le(n3,e3); \
assert_ge(n4,b4); assert_le(n4,e4); \
assert_ge(n5,b5); assert_le(n5,e5); \
assert_ge(n6,b6); assert_le(n6,e6); 
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d, n3=%d, n4=%d, n5=%d, n6=%d", n1, n2, n3, n4, n5, n6);

const double& dTensor6d::get(int n1,int n2,int n3,int n4,int n5,int n6) const
{
    int k = getidx(n1,n2,n3,n4,n5,n6);
    bounds_check();
    return vec[k];
}

double& dTensor6d::fetch(int n1,int n2,int n3,int n4,int n5,int n6)
{
    int k = getidx(n1,n2,n3,n4,n5,n6);
    bounds_check();
    return vec[k];
}


void dTensor6d::set(int n1,int n2,int n3,int n4,int n5,int n6,double value)
{
    int k = getidx(n1,n2,n3,n4,n5,n6);
    bounds_check();
    vec[k] = value;
}


// class dTensor5d

// check parameters and update derived parameters.
void dTensor5d::init()
{
    assert_printf(s1>=0 && s2>=0 && s3>=0 && s4>=0 && s5>=0,
            "b1=%d, b2=%d, b3=%d, b4=%d, b5=%d"
            "s1=%d, s2=%d, s3=%d, s4=%d, s5=%d",
            b1,b2,b3,b4,b5,
            s1,s2,s3,s4,s5);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    e3 = b3 + s3 - 1;
    e4 = b4 + s4 - 1;
    e5 = b5 + s5 - 1;
    size = s1*s2*s3*s4*s5;
    dTensorBase::init();
}

dTensor5d::dTensor5d(
        int s1i, int s2i, int s3i, int s4i, int s5i,
        int b1i, int b2i, int b3i, int b4i, int b5i ) :
    s1(s1i), s2(s2i), s3(s3i), s4(s4i), s5(s5i), 
    b1(b1i), b2(b2i), b3(b3i), b4(b4i), b5(b5i)
{ init(); }

dTensor5d::dTensor5d(const dTensor5d& in, CopyMode::Enum copyMode) 
    : dTensorBase(in, copyMode),
    s1(in.s1), s2(in.s2), s3(in.s3), s4(in.s4), s5(in.s5),
    b1(in.b1), b2(in.b2), b3(in.b3), b4(in.b4), b5(in.b5),
    e1(in.e1), e2(in.e2), e3(in.e3), e4(in.e4), e5(in.e5)
     { }

    // alias for assert_eq
    //#define ae(a,b) assert_op(==,a,b);
#define ae assert_eq
void dTensor5d::copyfrom(const dTensor5d& in)
{
    ae(b1,in.b1); ae(b2,in.b2); ae(b3,in.b3); ae(b4,in.b4); ae(b5,in.b5);
    ae(e1,in.e1); ae(e2,in.e2); ae(e3,in.e3); ae(e4,in.e4); ae(e5,in.e5);
    ae(s1,in.s1); ae(s2,in.s2); ae(s3,in.s3); ae(s4,in.s4); ae(s5,in.s5);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2); \
assert_ge(n3,b3); assert_le(n3,e3); \
assert_ge(n4,b4); assert_le(n4,e4); \
assert_ge(n5,b5); assert_le(n5,e5); 
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d, n3=%d, n4=%d, n5=%d", n1, n2, n3, n4, n5);

const double& dTensor5d::get(int n1,int n2,int n3,int n4,int n5) const
{
    int k = getidx(n1,n2,n3,n4,n5);
    bounds_check();
    return vec[k];
}

double& dTensor5d::fetch(int n1,int n2,int n3,int n4,int n5)
{
    int k = getidx(n1,n2,n3,n4,n5);
    bounds_check();
    return vec[k];
}


void dTensor5d::set(int n1,int n2,int n3,int n4,int n5,double value)
{
    int k = getidx(n1,n2,n3,n4,n5);
    bounds_check();
    vec[k] = value;
}

// class dTensor4d

// check parameters and update derived parameters.
void dTensor4d::init()
{
    assert_printf(s1>=0 && s2>=0 && s3>=0 && s4>=0,
            "b1=%d, b2=%d, b3=%d, b4=%d"
            "s1=%d, s2=%d, s3=%d, s4=%d",
            b1,b2,b3,b4,
            s1,s2,s3,s4);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    e3 = b3 + s3 - 1;
    e4 = b4 + s4 - 1;
    size = s1*s2*s3*s4;
    dTensorBase::init();
}

dTensor4d::dTensor4d(
        int s1i, int s2i, int s3i, int s4i,
        int b1i, int b2i, int b3i, int b4i ) :
    s1(s1i), s2(s2i), s3(s3i), s4(s4i), 
    b1(b1i), b2(b2i), b3(b3i), b4(b4i)
{ init(); }

dTensor4d::dTensor4d(const dTensor4d& in, CopyMode::Enum copyMode)
    : dTensorBase(in, copyMode),
    s1(in.s1), s2(in.s2), s3(in.s3), s4(in.s4),
    b1(in.b1), b2(in.b2), b3(in.b3), b4(in.b4),
    e1(in.e1), e2(in.e2), e3(in.e3), e4(in.e4)
     { }

void dTensor4d::copyfrom(const dTensor4d& in)
{
    ae(b1,in.b1); ae(b2,in.b2); ae(b3,in.b3); ae(b4,in.b4);
    ae(e1,in.e1); ae(e2,in.e2); ae(e3,in.e3); ae(e4,in.e4);
    ae(s1,in.s1); ae(s2,in.s2); ae(s3,in.s3); ae(s4,in.s4);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2); \
assert_ge(n3,b3); assert_le(n3,e3); \
assert_ge(n4,b4); assert_le(n4,e4);
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d, n3=%d, n4=%d", n1, n2, n3, n4);

const double& dTensor4d::get(int n1,int n2,int n3,int n4) const
{
    int k = getidx(n1,n2,n3,n4);
    bounds_check();
    return vec[k];
}

double& dTensor4d::fetch(int n1,int n2,int n3,int n4)
{
    int k = getidx(n1,n2,n3,n4);
    bounds_check();
    return vec[k];
}

void dTensor4d::set(int n1,int n2,int n3,int n4,double value)
{
    int k = getidx(n1,n2,n3,n4);
    bounds_check();
    vec[k] = value;
}

// class dTensor3d

// check parameters and update derived parameters.
void dTensor3d::init()
{
    assert_printf(s1>=0 && s2>=0 && s3>=0,
            "b1=%d, b2=%d, b3=%d"
            "s1=%d, s2=%d, s3=%d",
            b1,b2,b3,
            s1,s2,s3);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    e3 = b3 + s3 - 1;
    size = s1*s2*s3;
    dTensorBase::init();
}

dTensor3d::dTensor3d(
        int s1i, int s2i, int s3i,
        int b1i, int b2i, int b3i ) :
    s1(s1i), s2(s2i), s3(s3i), 
    b1(b1i), b2(b2i), b3(b3i)
{ init(); }

dTensor3d::dTensor3d(const dTensor3d& in, CopyMode::Enum copyMode)
    : dTensorBase(in, copyMode),
    b1(in.b1), b2(in.b2), b3(in.b3),
    e1(in.e1), e2(in.e2), e3(in.e3),
    s1(in.s1), s2(in.s2), s3(in.s3) { }

dTensor3d::dTensor3d(const dTensor3d& in) : dTensorBase(in),
    s1(in.s1), s2(in.s2), s3(in.s3),
    b1(in.b1), b2(in.b2), b3(in.b3),
    e1(in.e1), e2(in.e2), e3(in.e3) { }

void dTensor3d::copyfrom(const dTensor3d& in)
{
    ae(b1,in.b1); ae(b2,in.b2); ae(b3,in.b3);
    ae(e1,in.e1); ae(e2,in.e2); ae(e3,in.e3);
    ae(s1,in.s1); ae(s2,in.s2); ae(s3,in.s3);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2); \
assert_ge(n3,b3); assert_le(n3,e3);
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d, n3=%d", n1, n2, n3);

const double& dTensor3d::get(int n1,int n2,int n3) const
{
    int k = getidx(n1,n2,n3);
    bounds_check();
    return vec[k];
}

void dTensor3d::set(int n1,int n2,int n3,double value)
{
    int k = getidx(n1,n2,n3);
    bounds_check();
    vec[k] = value;
}

// class dTensor2d

// check parameters and update derived parameters.
void dTensor2d::init()
{
    assert_printf(s1>=0 && s2>=0,
            "b1=%d, b2=%d"
            "s1=%d, s2=%d",
            b1,b2,
            s1,s2);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    size = s1*s2;
    dTensorBase::init();
}

dTensor2d::dTensor2d(
        int s1i, int s2i,
        int b1i, int b2i ) :
    s1(s1i), s2(s2i), 
    b1(b1i), b2(b2i)
{ init(); }

dTensor2d::dTensor2d(const dTensor2d& in) : dTensorBase(in),
    s1(in.s1), s2(in.s2),
    b1(in.b1), b2(in.b2),
    e1(in.e1), e2(in.e2)
     { }

void dTensor2d::copyfrom(const dTensor2d& in)
{
    ae(b1,in.b1); ae(b2,in.b2);
    ae(e1,in.e1); ae(e2,in.e2);
    ae(s1,in.s1); ae(s2,in.s2);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2);
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d", n1, n2);

const double& dTensor2d::get(int n1,int n2) const
{
    int k = getidx(n1,n2);
    bounds_check();
    return vec[k];
}

double& dTensor2d::fetch(int n1,int n2)
{
    int k = getidx(n1,n2);
    bounds_check();
    return vec[k];
}

void dTensor2d::set(int n1,int n2,double value)
{
    int k = getidx(n1,n2);
    bounds_check();
    vec[k] = value;
}

void dTensor1d::copyfrom(const dTensor1d& in)
{
    ae(b1,in.b1);
    dTensorBase::copyfrom(in);
}

// --- section: multidimensional 1-based tensor classes ---

// class dTensorBC6

void dTensorBC6::copyfrom(const dTensorBC6& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor6d::copyfrom(in);
}

dTensorBC6::dTensorBC6(const dTensorBC6& in, CopyMode::Enum copyMode)
: dTensor6d(in,copyMode),
  S1(in.S1), S2(in.S2), S3(in.S3), S4(in.S4), S5(in.S5), S6(in.S6),
  mbc(in.mbc), ndims(in.ndims)
{ }

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC6::dTensorBC6(int S1i, int S2i, int S3i, int S4i, int S5i, int S6i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), S3(S3i), S4(S4i), S5(S5i), S6(S6i),
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i; s3=S3i; s4=S4i; s5=S5i;  s6=S6i;
    b1=1;   b2=1;   b3=1;   b4=1;   b5=1;    b6=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 6: b6=1-mbc; s6=S6+2*mbc;
        case 5: b5=1-mbc; s5=S5+2*mbc;
        case 4: b4=1-mbc; s4=S4+2*mbc;
        case 3: b3=1-mbc; s3=S3+2*mbc;
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}


// class dTensorBC5

void dTensorBC5::copyfrom(const dTensorBC5& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor5d::copyfrom(in);
}

dTensorBC5::dTensorBC5(const dTensorBC5& in, CopyMode::Enum copyMode)
: dTensor5d(in,copyMode),
  S1(in.S1), S2(in.S2), S3(in.S3), S4(in.S4), S5(in.S5),
  mbc(in.mbc), ndims(in.ndims)
{ }

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC5::dTensorBC5(int S1i, int S2i, int S3i, int S4i, int S5i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), S3(S3i), S4(S4i), S5(S5i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i; s3=S3i; s4=S4i; s5=S5i; 
    b1=1;   b2=1;   b3=1;   b4=1;   b5=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 5: b5=1-mbc; s5=S5+2*mbc;
        case 4: b4=1-mbc; s4=S4+2*mbc;
        case 3: b3=1-mbc; s3=S3+2*mbc;
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC4

void dTensorBC4::copyfrom(const dTensorBC4& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor4d::copyfrom(in);
}

dTensorBC4::dTensorBC4(const dTensorBC4& in, CopyMode::Enum copyMode)
: dTensor4d(in,copyMode),
  S1(in.S1), S2(in.S2), S3(in.S3), S4(in.S4), 
  mbc(in.mbc), ndims(in.ndims)
{ }

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC4::dTensorBC4(int S1i, int S2i, int S3i, int S4i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), S3(S3i), S4(S4i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i; s3=S3i; s4=S4i;
    b1=1;   b2=1;   b3=1;   b4=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 4: b4=1-mbc; s4=S4+2*mbc;
        case 3: b3=1-mbc; s3=S3+2*mbc;
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC3

void dTensorBC3::copyfrom(const dTensorBC3& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor3d::copyfrom(in);
}

dTensorBC3::dTensorBC3(const dTensorBC3& in, CopyMode::Enum copyMode)
  : dTensor3d(in,copyMode),
  S1(in.S1), S2(in.S2), S3(in.S3),
  mbc(in.mbc), ndims(in.ndims)
{ }

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC3::dTensorBC3(int S1i, int S2i, int S3i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), S3(S3i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i; s3=S3i;
    b1=1;   b2=1;   b3=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 3: b3=1-mbc; s3=S3+2*mbc;
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC2

void dTensorBC2::copyfrom(const dTensorBC2& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor2d::copyfrom(in);
}

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC2::dTensorBC2(int S1i, int S2i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i;
    b1=1;   b2=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC1

void dTensorBC1::copyfrom(const dTensorBC1& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor1d::copyfrom(in);
}

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC1::dTensorBC1(int S1i,
        int mbcin, int ndimsin) :
    S1(S1i),
    mbc(mbcin), ndims(ndimsin)
{
    size=S1i;
    b1=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 1: b1=1-mbc; size=S1+2*mbc;
        case 0: ;
    }
    init();
}

