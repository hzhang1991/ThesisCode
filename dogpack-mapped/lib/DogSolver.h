#ifndef _DOG_SOLVER_H__
#define _DOG_SOLVER_H__

#include "DogState.h"
#include "debug.h"

// Used to create output directory and copy a few files (Makefiles,
// parameters.ini, ... ) into the output folder.
void RunStartScript(int ndims);

class PositivityState;
class dTensorBase;
class dTensor1;
class dTensor2;

// -------------------------------------------------------------------------- //
//
// The main DoGPack solver class from which all solvers are derived
//
// -------------------------------------------------------------------------- //
class DogSolver
{

    // saving and initialization
    public:

        // ------------------------------------------------------------------ //
        // This is the first function that gets called from anywhere in the
        // code.  Its purpose is to perform the following tasks:
        //
        //    1.) Parse user supplied arguments.  These can be used to identify
        //        an output directory that is different from the default output
        //        directory.
        //
        //        TODO: APPARENTLY THIS CAN ALSO READ IN DEBUG OPTIONS! (-DS)
        //
        //    2.) Run the start script, lib/RunStartScript.cpp, which
        //        calls $(DOGPACK)/scripts/startscript.  The purpose of that
        //        script is to create the output directory for storing data from
        //        the simulation, as well as copy enough files to that directory
        //        which can be used to reproduce the simulation.
        //
        //    3.) Call RunDogpack().
        //
        // ------------------------------------------------------------------ //
        int main(int argc, char* argv[]);

        // ------------------------------------------------------------------ //
        // RunDogpack is the final routine that DogSolver::main(int,char**)
        // calls.  This routine is the top level time manager for running the
        // simulation, and plays as referee during the simulation.  In brief, 
        // its purpose is to perform the following tasks:
        //
        //    1.) Print a welcome message to standard output.
        //
        //    2.) Ask dogParams to read and save the [dogParams] section of the
        //    "parameters.ini" file.  This global variable is accessible by
        //    including "DogParams.h" from anywhere in the code.
        //
        //    3.) Add extra information to the output directory through a single
        //    call to initOutputDirectories().
        //
        //    4.) Ask DogState to initialize themselves, or restart from old
        //    data.
        //
        //    5.) Run a single loop for each call to Output().  In each
        //    instance, a single call to DogSolve(tstart, tend) is made to
        //    advance the solution over a large interval.
        //
        //    6.) Print a timing message describing how long it took to run the
        //    code, and finally return.
        //
        // ------------------------------------------------------------------ //
        int RunDogpack();

        // ------------------------------------------------------------------ //
        // parser to read input arguments.  Current valid options are 
        // "-o outputdir" and "-d X", for an output directory and debugging
        // level.
        static void parse_arguments(int argc,char**argv);

        // ------------------------------------------------------------------ //
        // Main solver.  Right now, this is essentially hard coding a
        // Runge-Kutta method.  SDC methods are called by running DogSolveSDC in
        // place of DogSolve.
        // ------------------------------------------------------------------ //
        void DogSolve(const double tstart, const double tend);

        // I'm not exactly sure what this call does.  It looks like it tells the
        // state variable who its solver is.  (-DS)
        virtual void init();

        // TODO - these should be eliminated by incorporating them into DogSolve
        // virtual void DogSolveRK   (double tstart, double tend)=0;
        virtual void DogSolveLxW  (double tstart, double tend)=0;
        virtual void DogSolveSDC  (double tstart, double tend)=0;
        virtual void DogSolveUser (double tstart, double tend)=0;

        double get_dt() const { return dt;}
        const PositivityState& get_positivityState() const
        { return *positivityState; }
        PositivityState& fetch_positivityState() 
        { return *positivityState; }

        // constructor (0 is equivalent to NULL for pointers)
        DogSolver():state(0),state_old(0),dt(0), next_report_time(0),next_frame_time(0), positivityState(0)
        {
            set_solver(this);
            set_time_hack(0.);
        }

        // Destructor
        virtual ~DogSolver();

        // methods
        static void set_solver(DogSolver* in);
        static void set_time_hack(double in) { time_hack = in;}
        static double get_time_hack() { return time_hack; }
        static const DogSolver& get_solver() { return *solver; }
        static DogSolver& fetch_solver() { return *solver; }
        static const char* get_outputdir() {return outputdir;}

        // saving and initialization
        void write_time_state(int nframe, const char* outputdir) const;
        void read_time_state (int nframe, const char* outputdir);


    protected:
        double set_dt(double in) { dt=in; }
        const DogState& get_state() const { return *state; }
        const DogState& get_state_old() const { return *state_old; }
        DogState& fetch_state() { return *state; }
        DogState& fetch_state_old() { return *state_old; }
        void set_state(DogState* in){state=in;}
        void set_state_old(DogState* in){state_old=in;}

    protected:
        virtual void truncate_logfiles(double time)const;
        static  void truncate_logfile(double time, const char* logfile_basename);
        virtual void initParams();
        virtual void initOutputDirectory(int idx, const char* framedir);
        virtual void saveState();
        virtual void revertToSavedState();

    protected:
        double get_time() const { return state->get_time(); }

    private:

        // -- Data -- //
        double dt;
        double next_report_time;
        double next_frame_time;
        PositivityState* positivityState;
        DogState *state;
        DogState *state_old;

        void set_next_frame_time(double in){next_frame_time=in;}


        virtual void reset_smax()=0;

        // function that can be used to print extra information to output
        // directories.  For example, the 2D-Cartesian code uses this to save
        // plotting help files, in a single file called qhelp.dat.
        virtual void initOutputDirectories()=0;

    private:
        // "callbacks" overridden only in library
        // (make these private in e.g. DogSolverCart2)
        //
        virtual double GetCFL(double dt)const=0;
        //
        void reportAfterStep(const DogState& state); // This function (usually) calls ConSoln

    private:

        // These three functions are 'pure virtual function's, meaning classes
        // that inherit from this class must implement these functions.  In
        // addition, this means that this class cannot be instantiated.
        virtual void write_restart(int n) const=0;
        virtual void Output(int n)const=0;
        virtual void Restart(int nstart)=0;

    // RK methods - I truly don't think these methods belong here, otherwise
    // this class becomes too cumbersome.  It looks like this class is supposed
    // to be the guy who plays referee for the whole simulation, and shouldn't
    // need or have to know the details about which time stepping option is
    // being applied.  Just my two cents (-DS).
    private:
        virtual void AdvanceStageRK(int mstage, double dt,
                const DogState& in_Qstar, DogState& in_Qnew)=0;
        virtual void AdvanceStageRK2(int mstage, double dt,
                const DogState& in_Qold, DogState& in_Q1, DogState& in_Q2)=0;
        virtual void AdvanceStageRK4( double dt,
                const DogState& in_Qold, DogState& in_Q1, DogState& in_Q2)=0;
        int advanceTimeStageRK( int mstage, const double dt,
                DogState& Qinit, const DogState& Qstar, DogState& Qnew);
        int advanceTimeStageRK2( const int mstage, double dt,
                const DogState& Qold, DogState& Q1, DogState& Q2);
        int advanceTimeStageRK4(double dt,
                const DogState& Qold, DogState& Q1, DogState& Q2);
        int advanceTimeStepRK(double dt);

        // There are a few applications in 2d/twofluid that split the source
        // term off of the main solver.
        virtual int advanceSplitTimeStep(double dt);
        //
        int advanceFullTimeStep(double dt);

    private: // disabled (currently singleton)
        DogSolver(const DogSolver& in);
        DogSolver& operator=(const DogSolver& in);

        // static data
        //
    private:
        static DogSolver* solver;
        static double time_hack;
        static char* outputdir;

    private:
        static void set_outputdir(const char*);

}; // end class DogSolver

// convenience global accessor
const char* get_outputdir(); // { return DogSolver::get_outputdir(); }

class dTensor1;
class dTensor2;

// -------------------------------------------------------------------------- //
// RungeKutta solver
//
// I'm not convinced this is an appropriate spot to place this.
// Is this module supposed to describe the time integrators?  If so, then
// maybe.  We really want to make sure we retain the option of switching (or
// adding) different time integrator options to the code.  (-DS)
//
// -------------------------------------------------------------------------- //
class RKsolver
{

    public:

        // Constructor.  All pointers are by default initialized to NULL, and
        // private ints are initialized to 0.
        //
        // These are later set with a single call to SetRKinfo.cpp
        RKsolver():
            mstage(0),
            num_stages(0),
            alpha1 (NULL),
            alpha2 (NULL),
            beta   (NULL),
            gamma  (NULL),
            delta  (NULL)
        {}

        // Destructor:
        ~RKsolver();

        // TODO : WHY ISN'T AUX PASSED IN HERE TOO? (-DS)
        void advanceTimeStage( const int mstage, const double dt,
                const double told, const dTensorBase& qold, const dTensorBase& L,
                double &tnew, dTensorBase& qnew);
        void advanceTimeStage2(const int mstage, const double dt,
                const double told, const dTensorBase& qold,
                double& t1, dTensorBase& q1, const dTensorBase& L,
                double& t2, dTensorBase& q2);

        // Only used in the 4th-order Runge-Kutta method.
        void advanceTimeStage4( const double dt,
                const double told, const dTensorBase& qold,
                double& t1, dTensorBase& q1,
                double& t2, dTensorBase& q2);

        // This routine sets the time order.  In doing so, it derives the
        // following quantities which are computed and saved:
        //
        //      (*) stage number, mstage.
        //      (*) number of stages, num_stages.
        //      (*) Values for alpha1, alpha2 and beta in methods of orders 
        //          1-4.
        //      (*) Values for delta, gamma and beta in 5th order scheme.
        void set(int time_order);

    private:

        // stage number:
        int mstage;

        // Number of stages for a given order of accuracy.  This is derived
        // and saved when calling ::set(int time_order).
        int num_stages;

        // These parameters are used for all RK methods of orders 1-4.  See
        // Ketcheson's paper.
        dTensor1* alpha1;
        dTensor1* alpha2;
        dTensor1* beta;

        // These two are needed for 5th order Stepping
        dTensor2* gamma;
        dTensor1* delta;

};

// -------------------------------------------------------------------------- //
//
// -------------------------------------------------------------------------- //
class DogSolverTB : public DogSolver
{

    protected:
        void set_smax(dTensorBase* in){smax = in;}
        void set_L(dTensorBase* in){L = in;}

    // initialization
    protected:

        // Classes who inherit from this top class overwrite this initializer:
        virtual void init();

        // Constructor:
        DogSolverTB() : L(0),smax(0) {}

    private:

        // -- Methods -- // 

        // Accessors
        const dTensorBase& get_L()const{return *L;}

        // Methods
        virtual void reset_smax();
        virtual void AdvanceStageRK(int mstage, double dt,
                const DogState& in_Qstar, DogState& in_Qnew);
        virtual void AdvanceStageRK2(int mstage, double dt,
                const DogState& in_Qold, DogState& in_Q1, DogState& in_Q2);

        // Used in Runge-Kutta solver.  This advances the "fourth"-stage in
        // the low-storage fourth-order scheme.
        virtual void AdvanceStageRK4( double dt,
                const DogState& in_Qold, DogState& in_Q1, DogState& in_Q2);
    
        // -- Fields -- //
        dTensorBase* L;
        dTensorBase* smax;
        RKsolver rk;

};
// -------------------------------------------------------------------------- //

#endif
