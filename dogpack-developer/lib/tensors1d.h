#ifndef tensors1d_h
#define tensors1d_h

#ifndef NDIMS
#include <dimdefs.h> // for NDIMS (number of spatial dimensions)
#endif
#if (!defined(NDIMS) || !(NDIMS>=0))
#error "NDIMS must be defined"
#endif

// What's the difference between DEEP and DIMS copy? (-DS)
namespace CopyMode
{
    enum Enum
    {
        DEEP = 1,
        DIMS = 2,
    };
}

class iTensorBase
{
    protected:
        void init();
        iTensorBase(){};
        iTensorBase(const iTensorBase& in);
    public:
        iTensorBase(int size_in):size(size_in){init();}
        ~iTensorBase();
        void setall(int);
        const int numel() const { return size; }
#ifdef CHECK_BOUNDS
        const int& vget(int k) const;
        int& vfetch(int k);
        void vset(int k, int value);
#else
        const int& vget(int k) const
        {
            return vec[k];
        }
        int& vfetch(int k) {
            return vec[k];
        }
        void vset(int k, int value){
            vec[k]=value;
        }
#endif
    protected:
        int* vec;
        int size;
};

class iTensor1d : public iTensorBase
{
    // data
    protected:
        int b1;
        // methods
    protected:
        iTensor1d(){};
    public:
        int getidx(int n1) const
        {
            int k = n1-b1;
            return k;
        }
        int getsize() const { return size; }
    public:
        // constructor takes size and initial index in each dimension
        iTensor1d( int s1i, int b1i ) : b1(b1i) { size=s1i; init(); }
        iTensor1d(const iTensor1d& in) : iTensorBase(in), b1(in.b1) {}

        const int& get(int n1) const
        { return vget(n1-b1); }
        void set(int n1, int value)
        { return vset(n1-b1, value); }
};

class iTensor1 : public iTensor1d
{
    public:
        iTensor1(int s1i) :
            iTensor1d(s1i,1) { }
};

// === major section: dtensors (arrays of double) ===

class dTensorBase
{
    private: // disabled
        dTensorBase& operator=(const dTensorBase& in);
    protected:
        void init();
        dTensorBase(){};
        dTensorBase(const dTensorBase& in, CopyMode::Enum copyMode=CopyMode::DEEP);
        void copyfrom(const dTensorBase& in);
    public:
        dTensorBase(int size_in):size(size_in){init();}
        ~dTensorBase();
        void check();
        void setall(double);
        const int numel() const { return size; }
#ifdef CHECK_BOUNDS
        const double& vget(int k) const;
        double& vfetch(int k);
        void vset(int k, double value);
#else
        const double& vget(int k) const { return vec[k]; }
        double& vfetch(int k) {
            return vec[k];
#ifdef CHECK_INIT
            if(!(vec[k]==vec[k])) eprintf("vec[%d]=%24.16e",k,vec[k]);
#endif
        }
        void vset(int k, double value){vec[k]=value;}
#endif
    protected:
        double* vec;
        int size;
};

class dTensor1d : public dTensorBase
{
    // data
    protected:
        int b1;
        // methods
    private: // disabled
        dTensor1d& operator=(const dTensor1d& in);
    protected:
        dTensor1d(){};
    public:
        int getidx(int n1) const
        {
            int k = n1-b1;
            return k;
        }
        int getsize() const { return size; }
    public:
        // constructor takes size and initial index in each dimension
        dTensor1d( int s1i, int b1i ) : b1(b1i) { size=s1i; init(); }
        dTensor1d(const dTensor1d& in) : dTensorBase(in), b1(in.b1) {}
        void copyfrom(const dTensor1d& in);

        const double& get(int n1) const
        { return vget(n1-b1); }
        double& fetch(int n1)
        { return vfetch(n1-b1); }
        void set(int n1, double value)
        { return vset(n1-b1, value); }
};

class dTensor1 : public dTensor1d
{
    public:
        dTensor1(int s1i) :
            dTensor1d(s1i,1) { }
        void copyfrom(const dTensor1& in){ dTensor1d::copyfrom(in); }
    private: // disabled
        dTensor1& operator=(const dTensor1& in){ copyfrom(in); return *this; }
};

class dTensorBC1 : public dTensor1d
{
    private: // disabled
        dTensorBC1& operator=(const dTensorBC1& in);
    public:
        dTensorBC1(int s1i,
                int mbc, int ndims=1);
        void copyfrom(const dTensorBC1& in);
        int getmbc() const {return mbc;}
    private:
        int S1;
        int mbc;
        int ndims;
};

// methods for tensors
int count_fields(const char* str_in, char field_sep);
bool str_into_tensor(const char* str, iTensorBase& t, char field_sep);
bool str_into_tensor(const char* str, dTensorBase& t, char field_sep);
#endif
