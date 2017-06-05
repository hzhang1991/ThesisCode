// Including this header file allows access to all the 1D paramters described in
// the [grid] section of any 1D parameters.ini file.
//
// The singleton dogParamsCart1 allows access to these parameters.  For example,
// after including DogParamsCart1.h, one may accesses mx, the number of grid cells
// via:
//
//      dogParamsCart1.get_mx()
//
// There are currently four parameters in the [grid] section for 1D
// applications:
//
//          mx (==melems), mbc, xlow, xhigh and dx.
//
//
// See also: DogParams.h, DogParamsCart2.h, IniDocument.h (global parser)

#ifndef _DOGPARAMSCART1_H_
#define _DOGPARAMSCART1_H_

class IniDocument;   // The parser used in DoGPack for reading parameters.ini

struct DogParamsCart1
{

    // -- Methods -- //
    public:

        bool get_is_initialized(){return is_initialized;}
        DogParamsCart1(){is_initialized=false;}

        void init(IniDocument& ini_doc);
        void append_qhelp(const char* filename);

        // User feedback for the welcome screen
        void reportParameters();

        // Accessors
        //
        // Note: melems == mx for backwards compatability
        const int   & get_mx()      const{ return  mx;   }
        const int   & get_melems()  const{ return  mx;   }
        const int   & get_mbc()     const{ return  mbc;  }
        const double& get_xlow()    const{ return  xlow; }
        const double& get_xhigh()   const{ return  xhigh;}
        const double& get_dx()      const{ return  dx;   }

        // Setters
        void set_mx(int mx_in){ mx = mx_in; }
        void set_xlims(double,double);

    private:
        void checkParameters();
        void setDerivedParameters();

    // -- Fields -- //
    private:
        bool is_initialized;    // flag defining whether or not the section has been read
        int mx;                 // Number of elements (==melems)
        int mbc;                // number of ghost cells on each side
        double xlow;            // Lower end of domain
        double xhigh;           // Right hand endpoint of domain
        double dx;              // Derived from other parameters

};

extern DogParamsCart1 dogParamsCart1;

#endif
