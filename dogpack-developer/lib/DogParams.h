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
// The bulk work of defining these parameters is carried out by the init()
// function. Parameters are passed around as a global variable.

#include <stdio.h> // for FILE
#define MAX_NUM_SUBSYSTEMS 10

// What's this namespace used for? (-DS)
namespace SplittingType
{
    enum Enum{ none=0, sf=1, fs=2, sfs=3, fsf=4 };
}

enum DataFmt {ASCII=1,HDF5=5};

struct DogParams
{

    public:

        // This function does the bulk work of scanning the section labeled
        // [dogParams], and saving the information located there inside this
        // singleton struct.
        //
        // This needs to be called at least once, and future calls are
        // ignored, so it's OK to call this extra times.
        //
        void init();
        bool get_is_initialized()               const{ return is_initialized;}

        // files included in the output directory.  These are used for helping
        // the visualization routines and saving what parameters were used for
        // creating the simulation.
        void write_qhelp(const char* filename);
        void write_output_parameters_ini(const char* filename);

        int   get_ndims()                       const{ return ndims;}
        const char* get_mesh_type()             const{ return mesh_type;}
        const char* get_time_stepping_method()  const{ return time_stepping_method;}
        const char* get_limiter_method()        const{ return limiter_method;}
        int   get_nout()                        const{ return nout;}
        double get_tfinal()                     const{ return tfinal;}
        const double* get_cflv()                const{ return cflv;}
        int   get_nv()                          const{ return nv;}


        // These calls are depracated, but are included to retain backwards
        // compatability:
        const int* get_method() const{ return  method;}
        const double* get_dtv() const{ return  dtv;}

        // These accessors are safer, given that the can't modify the arrays
        // method[] and dtv[]
        double get_initial_dt()       const{ return dtv[1]; }
        double get_max_dt()           const{ return dtv[2]; }
        int    get_space_order()      const{ return method[1];}
        int    get_time_order()       const{ return method[2];}
        int    get_use_limiter()      const{ return method[3];}
        int    get_verbosity()        const{ return method[4];}
        int    get_mcapa()            const{ return method[5];}
        int    get_maux()             const{ return method[6];}
        int    get_source_term()      const{ return method[7];}

        // ?? 
        SplittingType::Enum get_splitting() const{ return splitting;}

        // ??
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
        int get_nout_per_restart()    const{ return nout_per_restart;}
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

        // TODO remove these depracated calls, and make sure the rest of the
        // code still functions:
        // deprecated; this is in the DogState class.
        void set_time(double in) { time = in; }

        // deprecated; this is in the DogSolver class.
        void   set_dt(double in) { dt = in; }
        double get_time() const  { return time; }
        double get_dt()   const  { return dt; }

        // convenience functions (derived from setting, e.g. 
        // moment_limiter = 1 in the [dogParams] section of parameters.ini. )
        bool using_moment_limiter();
        bool using_viscosity_limiter();
        bool using_relax_limiter();
        bool using_positive_limiter();

        // pyclaw stuff
        int  get_withPyClawPlotting() const{ return withPyClawPlotting; }
        void set_withPyClawPlotting(const int in) { withPyClawPlotting=in; }

        // quadrature order used for projecting the initial conditions
        int get_ic_quad_order() const{ return ic_quad_order; }

        // Destructor
        ~DogParams();

        // Default Constructor using initialization lists 
        // (e.g. http://www.cprogramming.com/tutorial/initialization-lists-c++.html )
        //
        // 0 and NULL are equivalent here.
        DogParams():
            is_initialized(false),
            mesh_type(NULL),
            time_stepping_method(NULL),
            limiter_method(NULL),         
            splitting(SplittingType::none),
            flux_term(true),
	    nout_per_plot(NULL),
            use_divfree(0), 
	    how_many_vectors_divfree(0), 
            which_compnt_divfree(NULL),
            withPyClawPlotting(0),            
            kmax(0),
            kmax_divfree(0),      
            how_many_generic_components(0), 
            generic_components(NULL)            
    {}

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

        // these have been deprecated
        double time;                        // should use DogState::time
        double dt;                          // should use DogSolver::dt

        int use_divfree;
        int how_many_vectors_divfree; int* which_compnt_divfree;
        int ic_quad_order; // james: this is new thing

        // flag indicating whether or not to print contents to the output
        // directory that conform with PyClaw's plotting routines:
        int withPyClawPlotting;

        // derived parameters
        //
        int kmax;
        int kmax_divfree;
        double frame_interval;              // time interval per call to DogSolve*
        double frame_subinterval;           // time interval per ConSoln output (0=always)
        int how_many_generic_components; 
        int* generic_components;

    //
    // methods
    //
    private:
        void setDerivedParameters();
        void set_kmax();
        void set_kmax_divfree();
        void checkParameters();
        void checkParameters1();
        void checkParameters2();
        void checkParameters3();
        void checkParameters4();
        void reportParameters();
        void resetParameters();

};

extern DogParams dogParams;

#endif
