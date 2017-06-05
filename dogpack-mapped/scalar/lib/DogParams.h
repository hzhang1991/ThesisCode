#ifndef _DOGPARAMS_H_
#define _DOGPARAMS_H_

// Including this header file allows access to all the paramters described in
// the [dogParams] section of any parameters.ini input file.
//
// The singleton dogParams allows access to these parameters.  For example,
// after including DogParams.h, one may accesses ndims, the number of
// dimensions in the problem through the following call:
//
//      dogParams.get_ndims()
//

#include <stdio.h> // for FILE
#define MAX_NUM_SUBSYSTEMS 10

using namespace std;

// What's this namespace used for? (-DS)
namespace SplittingType
{
    enum Enum{ none=0, sf=1, fs=2, sfs=3, fsf=4 };
}

enum DataFmt {ASCII=1,HDF5=5};

struct DogParams
{

    private:
        bool is_initialized;
        int ndims;
        char* mesh_type;
        int nout;
        double tfinal;
        double  dtv[2+1];
        double cflv[2+1];
        int nv;
        char* time_stepping_method;
        char* limiter_method;
        int method[7+1];
        SplittingType::Enum splitting;
        bool flux_term;
        int meqn;
        int nrestart;
        int maintained_restart;
        int report_frame_idx;
        DataFmt datafmt;
        int nout_per_restart;
        int how_many_plot_resolutions; int* nout_per_plot;
        // number of time steps per call to DogSolve*
        // at which ConSoln should print output (negative=always)
        // (see e.g. apps/plasma/2d/twofluid/g10/ConSoln.cpp for example of use)
        int num_subintervals;
        //
        // these have been deprecated
        double time; // should use DogState::time
        double dt; // should use DogSolver::dt
        //
        int use_divfree;
        int how_many_vectors_divfree; int* which_compnt_divfree;
        int ic_quad_order; // james: this is new thing
        //
        // derived parameters
        //
        int kmax;
        int kmax_divfree;
        double frame_interval; // time interval per call to DogSolve*
        double frame_subinterval; // time interval per ConSoln output (0=always)
        int how_many_generic_components; int* generic_components;

    //
    // methods
    //
    private:
        void setDerivedParameters();
        void set_kmax();
        void checkParameters();
        void checkParameters1();
        void checkParameters2();
        void checkParameters4();
        void reportParameters();
        void resetParameters();

    // Static methods moved from DogSolver:
    // (TODO - declaring these as static doesn't let this compile!  is this
    // because DogParams is a struct and not a class? (-DS)
    public:
//      static char* outputdir;
//      static void parse_arguments(int argc, char**argv);
        void parse_arguments(int argc, char**argv);
        const char* get_outputdir(){ return outputdir; }

    private:
//      static void set_outputdir(const char*);
        void set_outputdir(const char*);
        char* outputdir;

    public:

        // Destructor
        ~DogParams();

        // Default Constructor using initialization lists 
        // (e.g. http://www.cprogramming.com/tutorial/initialization-lists-c++.html )
        DogParams():
            mesh_type(NULL),
            time_stepping_method(NULL),
            limiter_method(NULL),
            is_initialized(false), nout_per_plot(NULL),
            use_divfree(0), how_many_vectors_divfree(0), which_compnt_divfree(NULL),
            how_many_generic_components(0), generic_components(NULL),
            kmax(0),kmax_divfree(0), //withPyClawPlotting(0),
            flux_term(true), splitting(SplittingType::none)
        {}

        void init();
        void write_qhelp(const char* filename);
        void write_output_parameters_ini(const char* filename);

        bool get_is_initialized()               const{ return is_initialized;}
        int   get_ndims()                       const{ return ndims;}
        const char* get_mesh_type()             const{ return mesh_type;}
        const char* get_time_stepping_method()  const{ return time_stepping_method;}
        const char* get_limiter_method()        const{ return limiter_method;}
        int   get_nout()                        const{ return nout;}
        double get_tfinal()                     const{ return tfinal;}
        const double* get_cflv()                const{ return cflv;}
        int   get_nv()                          const{ return nv;}

        // get_dtv() is deprecated
        const double* get_dtv() const{ return  dtv;}
        // use these instead
        double get_initial_dt() const{ return dtv[1]; }
        double get_max_dt()     const{ return dtv[2]; }
        // get_method is deprecated
        const int* get_method() const{ return  method;}

        // use these instead
        int    get_space_order()      const{ return method[1];}
        int    get_time_order()       const{ return method[2];}
        int    get_use_limiter()      const{ return method[3];}
        int    get_verbosity()        const{ return method[4];}
        int    get_mcapa()            const{ return method[5];}
        int    get_maux()             const{ return method[6];}
        int    get_source_term()      const{ return method[7];}

        SplittingType::Enum get_splitting() const{ return splitting;}
        bool   get_flux_term()        const{ return flux_term;}
        int    get_meqn()             const{ return  meqn;}
        int    get_mrestart()         const{ if(nrestart >= 0) return 1; return 0;}
        int    get_nrestart()         const{ return  nrestart;}
        int    get_nstart_for_time(double t)const{ return int(t/tfinal*nout+1e-14);}
        int    get_maintained_restart() const{ return maintained_restart;}
        int    get_report_frame_idx() const{ return  report_frame_idx;}
        DataFmt get_datafmt()         const{ return  datafmt;}
        int get_num_subintervals()    const{ return  num_subintervals;}
        //
        double get_frame_interval()   const{ return frame_interval;}
        double get_frame_subinterval()const{ return frame_subinterval;}
        int get_nout_per_restart() const{ return nout_per_restart;}
        int get_how_many_plot_resolutions() const{ return how_many_plot_resolutions;}
        const int* get_nout_per_plot() const{ return nout_per_plot;}
        int get_kmax() const{ return kmax;}
        int get_kmax_divfree() const {return kmax_divfree;}
        int get_use_divfree() const { return use_divfree; }
        int get_how_many_vectors_divfree()
            const { return how_many_vectors_divfree; }
        const int* get_which_compnt_divfree()
            const{ return which_compnt_divfree; }
        int get_how_many_generic_components()
            const{ return how_many_generic_components; }
        int get_generic_component(int i)
            const{ return generic_components[i]; }

        // public methods for setting should not exist in this class
        //
        // deprecated; this is in the DogState class.
        void set_time(double in) { time = in; }
        // deprecated; this is in the DogSolver class.
        void set_dt(double in)   { dt = in; }

        // deprecated
        double get_time() const { return time; }
        double get_dt()   const { return dt; }

        // convenience functions
        bool using_moment_limiter();
        bool using_viscosity_limiter();
        bool using_relax_limiter();
        bool using_positive_limiter();

        // initial condition helper
        int get_ic_quad_order() const{ return ic_quad_order; }
};

extern DogParams dogParams;
#endif
