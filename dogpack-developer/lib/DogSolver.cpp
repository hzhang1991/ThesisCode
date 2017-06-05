// -------------------------------------------------------------------------- //
// DogSolver.cpp.    TODO - add some documentation describing the purpose of 
//                          this module.
//
// This module was originally put together to define a MOL formulation that
// could be used everywhere in the code.  However, this solver was never pushed
// past the 2D-Cartesian part of the code.  For example, 2D unstructured code
// still calls DogSolveRK_Unst, the 1D code still calls main_global which calls
// RunDogPack, etc.
//
// This "solver" currently contains a handful of parsing functions, various
// routines required to implement an RK-solver, calls to positivity state 
// and other items.
// -------------------------------------------------------------------------- //

// These header files are not included in DogSolver.h in order to speed up
// compilation time.
#include <stdio.h> // for printf, sscanf, fopen
#include <cstring> // for strdup
#include <string>
#include <unistd.h> // for ftruncate
#include <stdlib.h> // for free, malloc
#include <sys/stat.h> // for mkdir
#include <sys/errno.h> // for errno, EEXIST
#include <fcntl.h> // for open
#include <dirent.h> // for readdir
#include "getdelim.h"
#include "dog_io.h" // for copyFile()
#include "assert.h"
#include "debug.h"
#include "dog_math.h"
#include "linked_ptr.h"
#include "DogState.h"
#include "DogSolver.h"
#include "IniDocument.h"
#include "DogParams.h"
#include "Positivity.h"
#include "tensors.h" // needed by RKsolver and DogSolverTB
#include "ext_time.h" /* for get_utime and timeval_diff */
#include <iostream>
using namespace std;

using std::string;

// -------------------------------------------------------------------------- //
// Single function to read ndims, the number of dimensions from the [dogParams]
// section of parameters.ini.
static int get_ndims()
{
    int ndims;
    const char* s_ndims = ini_doc["dogParams"]["ndims"].c_str();
    if(!sscanf(s_ndims, "%d" ,&ndims))
        eprintf("invalid_value: s_ndims = %s", s_ndims);
    return ndims;
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// This is the `global main', that every piece of the code calls.
// See notes in DogSolver.h.
int DogSolver::main(int argc, char* argv[])
{
    ini_doc.initFromFile("parameters.ini");
    DogSolver::parse_arguments(argc,argv);
    RunStartScript(get_ndims());
    return RunDogpack();
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// This routine runs the code.  It's responsible for printing a welcome message
// to standard output, and plays the referee role for running the simulation.
//
// See the notes in DogSolver.h
//
// -------------------------------------------------------------------------- //
int DogSolver::RunDogpack()
{
    // Get current time
    timeval start_time = get_utime();

    // ------------------------------------------------------------

    // Print a welcome message.
    printf("\n"
            "   ------------------------------------------------   \n"
            "   | DoGPack: The Discontinuous Galerkin Package  |   \n"
            "   | Developed by the research group of           |   \n"
            "   |            James A. Rossmanith               |   \n"
            "   |            Department of Mathematics         |   \n"
            "   |            Iowa State University             |   \n"
            "   ------------------------------------------------   \n\n");

    // Ask the global variable DogParams dogParams to read and save 
    // the [dogParams] section of parameters.ini.
    //
    // This variable is accessible by adding #includde "DogParams.h" in any 
    // file in the code.
    initParams();

    // initialize solver based on ini parameters
    init();

    // populate output directory with more data.  This routine is a virtual
    // routine, and can be overwritten by inheriting classes.  The 2D-Cartesian
    // code uses this call to save qhelp.dat.
    initOutputDirectories();

    // initialize state of solver
    int    nstart = 0;
    double tstart = 0.;

    if (dogParams.get_mrestart())
    {
        // --------------------------
        // Do a restart from old data
        // --------------------------
        Restart(dogParams.get_nrestart());
        tstart = get_state().get_time();
        nstart = dogParams.get_nstart_for_time(tstart);
        dprintf1("restarting at tstart=%f, frame=%d, nrestart=%d",
                tstart,nstart,dogParams.get_nrestart());
        truncate_logfiles(tstart);
    }
    else
    {
        // --------------------------
        // Start new computation
        // --------------------------

        // Set initial time
        fetch_state().set_time(tstart);
        set_time_hack(tstart);

        fetch_state().InitState();

        // Apply post processing (e.g. smoothing) to initial data
        fetch_state().AfterInitState();

        // Output initial frame to file
        Output(nstart);

        // Compute conservation and print to file
        reportAfterStep(get_state());
    }

    // loop for Output(). In this loop:
    //   tstart is maintained as time at beginning of current loop iteration
    //   tend is maintained as anticipated time at end of current loop iteration
    //   n is maintained as current loop iteration number
    const int nout      = dogParams.get_nout();
    const double tfinal = dogParams.get_tfinal();
    double tend         = tstart;
    const double dtout  = tfinal/double(nout);
    for (int n=nstart+1; n<=nout; n++)
    {

        tstart = tend;
        tend = dtout*n;
        assert_ge(tend, tstart);

        const string time_stepping_method = dogParams.get_time_stepping_method();
        // Solve hyperbolic system from tstart to tend
        // (This select case should really be done at the time-stepping level)
        if (time_stepping_method == "Runge-Kutta")
        {

            // Main call:
            DogSolve(tstart,tend);
            // DogSolveRK(tstart,tend); // deprecated

        }
        else if (time_stepping_method == "SDC")
        {
            // Spectral deferred correction (SDC) time-stepping scheme
            DogSolveSDC(tstart,tend); // should be incorporated into DogSolve()
        }
        else if (time_stepping_method == "Lax-Wendroff")
        {
            DogSolveLxW(tstart,tend); // should be incorporated into DogSolve()
        }
        else if (time_stepping_method == "User-Defined")
        {
            // User defined time-stepping scheme
            DogSolveUser(tstart,tend); // should be incorporated into DogSolve()
        }
        else
        {
            unsupported_value_error(time_stepping_method.c_str());
        }

        // Output data to file
        fetch_state().set_time(tend);

        Output(n);

        // Done with solution from tstart to tend
        printf("DOGPACK: Frame %3d: at time t =%12.5e\n\n", n,tend);
    }

    // Get current time
    timeval end_time = get_utime();

    // Output elapsed time
    double diff_utime = timeval_diff(end_time, start_time);
    printf(" Total elapsed time in seconds = %11.5f\n\n", diff_utime);

    return 0;
}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// I don't declare this (e.g. as a static method) in DogSolver in
// order to avoid #include <string> (#include <iosfwd> would be
// insufficient since returning a reference is not an option due
// to dynamic memory allocation). -eaj
//
string get_framedir(int movie_idx)
{

    assert_lt(movie_idx,10);
    assert_gt(movie_idx,0);
    char movie_id[2];
    snprintf(movie_id,2,"%d",movie_idx);

    // initOutputDirectories() creates this directory
    // on startup and populates it with sufficient information
    // to plot data (parameters.ini, out_parameters.ini)
    return string(get_outputdir()) + "/movie" + movie_id;
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Section: DogSolver methods
// -------------------------------------------------------------------------- //

char* DogSolver::outputdir   = 0;
DogSolver* DogSolver::solver = NULL;
double DogSolver::time_hack  = 0.0;

// convenience accessor
const char* get_outputdir() { return DogSolver::get_outputdir(); }

void DogSolver::set_outputdir (const char* arg)
{
    free(outputdir);
    outputdir = strdup(arg);
}

void DogSolver::set_solver(DogSolver* in)
{
    // require this class to be singleton
    assert(!solver);
    solver = in;
}

// subsection: saving, initialization, and cleanup

DogSolver::~DogSolver()
{
    if(solver) solver=0;
    assert(state); {delete state; state=0;}
    assert(state_old); {delete state_old; state_old=0;}
    delete positivityState;
}

void DogSolver::revertToSavedState()
{
    state->copyfrom(*state_old);
}

void DogSolver::saveState()
{
    state_old->copyfrom(*state);
}

// -------------------------------------------------------------------------- //
// This function parses user supplied arguments when running the code.
//   e.g., this function parses the section "-o output_some_other_directory"
//   when the user calls "./dog.exe -o output_some_other_directory"
// 
// TODO - apparently DEBUG options can be specified from the command line.
// Because there's no documentation for this, I don't know what all options are
// available (-DS).
// -------------------------------------------------------------------------- //
void DogSolver::parse_arguments(int argc,char**argv)
{

    // flag used for telling if arguments were sucessfully parsed.
    int show_usage=0;

    // default output directory:
    char* outputdir=(char*)"output";

    // parse arguments
    if(argc>1)
    {
        for(int arg_idx=1; arg_idx<argc; arg_idx++)
        {
            if(argv[arg_idx][0]=='-')
            {

                // "-o" specifies which output directory to use.
                if(argv[arg_idx][1]=='o' && ++arg_idx < argc)
                { outputdir=argv[arg_idx]; }

                // "-d" specifies which debug level to use.   What are the legal
                // options that can be used here? (-DS)
                else if(argv[arg_idx][1]=='d' && ++arg_idx < argc)
                {

                    int debug_level;
                    int success = sscanf(argv[arg_idx],"%d",&debug_level);

                    if(success!=1)
                    { goto show_usage; }

                    DebugLevel::set(debug_level);
                    dprintf("set debug level to %d",DebugLevel::get());
                }
                else{ goto show_usage; }
            }
            else { goto show_usage; }
        }
    }
    set_outputdir(outputdir);
    return;

show_usage:
    {
        // TODO - add in the debug options as part of this helper statement!
        printf(
                "usage: %s [-o outputdir]\n"
                "  -o outputdir : use outputdir as output directory [default: output]\n",
                argv[0]);
        exit(1);
    }
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
void DogSolver::init()
{

    dt=dogParams.get_initial_dt();

    // initialize positivity indicator
    assert(!positivityState);
    positivityState = new PositivityState;

    // finish initializing states
    //
    // What do these calls do? (-DS)  Tell the state variables who its solver
    // is?
    //
    assert(&fetch_state());
    fetch_state().set_solver(*this);
    if(!&get_state_old())
    { set_state_old(get_state().clone(CopyMode::DIMS)); }

}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
void DogSolver::initParams()
{ dogParams.init(); }
// -------------------------------------------------------------------------- //


//////// Q: Can we move the parsing functions elsewhere? (-DS) /////////////////

// -------------------------------------------------------------------------- //
// this could be renamed and made static
void DogSolver::initOutputDirectory(int idx, const char* framedir)
{
    int error = mkdir( framedir, 0777 );
    // complain if could not create directory
    // (and if directory does not already exist)
    if(error<0 && errno!=EEXIST)
    {   eprintf("could not create directory %s", framedir); }

    // copy parameters.ini into framedir
    string outFile = string(framedir)+"/parameters.ini";
    copyFile("parameters.ini",outFile.c_str());
}

static ssize_t my_getline(char** lineptr, size_t* n, FILE* stream)
{
    ssize_t ret = getdelim (lineptr, n, '\n', stream);
    return ret;
}

void DogSolver::truncate_logfile(double time, const char* logfile_basename)
{
    dprintf3("attempting to truncate file %s",logfile_basename);

    string filename = string(get_outputdir())+"/"+logfile_basename;
    const char* logfile = filename.c_str();
    // open the file

    size_t bytes_read;
    size_t nbytes = 1024;
    char *line;

    long len;
    char *buf;
    int fd = open(logfile,O_RDWR); // fd = file descriptor
    if(fd<=0)
    {
        dprintf3("Evidently there is no logfile to truncate named: %s\n",logfile);
        return;
    }
    FILE *fp; // file pointer
    fp = fdopen(fd, "r+");
    assert(fp);

    line = (char *) malloc (nbytes + 1);

    // find the first line on which the first number exceeds time
    // and truncate the file
    size_t total_bytes_read = 0;
    double first_double=-1.;
    while(bytes_read = my_getline (&line, &nbytes, fp), bytes_read!=-1)
    {
        // if the line begin with a number that exceeds time truncate
        {
            const char* ptr=line;
            double old_first_double = first_double;
            int retval = sscanf(ptr,"%lf",&first_double);
            if(retval==EOF)
            {
                dprintf3("no need to truncate file %s", logfile);
                break;
            }

            // detect if the file does not actually conform to the format
            // of a log file and print an error (could change this to
            // print a warning and skip this file).
            {
                if(retval!=1)
                {
                    eprintf("not truncating file %s"
                            "\n   because the following line does not begin with a number:"
                            "\n   line=%s",
                            logfile_basename, line);
                    break;
                }
                // in a log file the first number of each line
                // (the time) should be increasing
                assert_gt(first_double,old_first_double);
            }

            if(first_double>time)
            {
                // truncate file
                dprintf3("truncating logfile %s beginning with line\n%s", logfile,line);

                // make sure this truncation didn't causea n error
                if( ftruncate(fd,total_bytes_read) == -1 )
                {
                    printf("Error while calling ftruncate");
                    exit(1);
                }
                break;
            }
        }
        total_bytes_read += bytes_read;
    }

    // clean up
    // close(fd); // equivalent to the following
    fclose(fp);
    free(line);
}

static bool is_numeric_character(char c)
{
    if('0' <= c && c <= '9') return true;
    return false;
}

void DogSolver::truncate_logfiles(double time)const
{
    // The hack in this file could be removed if we used
    // a suffix .log which is only used for log files
    //
    // truncate all log files. (A log file has a .log or .dat
    // extension, is not an array file, and each line begins with a
    // number which is supposed to be the time).
    DIR *dp = opendir (get_outputdir());

    if (dp != NULL)
    {
        struct dirent *ep;     
        while (ep == readdir (dp))
        {
            //puts (ep->d_name);
            const char * filename = ep->d_name;
            int length = strlen(filename);
            if(length<4) continue;

            bool is_log = !strncmp(".log",&filename[length-4],4);
            bool is_dat = !strncmp(".dat",&filename[length-4],4);
            if(!is_log && !is_dat) continue;

            // hack to determine if a .dat file is really a log file
            if(is_dat)
            {
                bool is_array = is_dat && (filename[0]=='a' || filename[0]=='q')
                    && is_numeric_character(filename[1])
                    && is_numeric_character(filename[2])
                    && is_numeric_character(filename[3])
                    && is_numeric_character(filename[4]);
                bool is_qhelp = !strcmp("qhelp.dat",filename);
                if(!is_array && !is_qhelp) is_log = true;
            }
            if(is_log)
            {
                dprintf("truncating log file %s",filename);
                truncate_logfile(time,filename);
            }
        }
        (void) closedir (dp);
    }
    else
        eprintf("Couldn't open directory %s", get_outputdir());
}

// -------------------------------------------------------------------------- //
// reading and writing routines
// -------------------------------------------------------------------------- //

string get_varframe_filename(const char* dir, const char* varname, int nframe)
{
    assert_ge(nframe,0);
    assert_le(nframe,9999);
    char frame_num[5];
    snprintf(frame_num,5,"%04d",nframe);
    return string(dir) + "/" + varname + frame_num;
}

static FILE* open_state_file(const char* outputdir,int nframe,const char* mode)
{
    string filename = get_varframe_filename(outputdir,"t",nframe)+".ini";
    FILE* file = fopen(filename.c_str(),mode);
    if(!file) eprintf("Could not open file %s for %s",
            filename.c_str(), (mode[0]=='w'?"writing":"reading"));
    return file;
}

// -------------------------------------------------------------------------- //
// This should save information sufficient to reconstruct the
// state of the solver.
// -------------------------------------------------------------------------- //
void DogSolver::write_time_state(int nframe, const char* outputdir) const
{
    FILE* file = open_state_file(outputdir,nframe,"w");
    fprintf(file,"time = %24.16e\n", get_time());
    fprintf(file,"dt = %24.16e\n", get_dt());
    fprintf(file,"next_report_time = %24.16e\n", next_report_time);
    positivityState->write_state(file);
    //fprintf(file,"next_frame_time = %24.16e\n", next_frame_time);
    // for an exact restart we should also have positivityState
    // save itself.
    fclose(file);
}

// -------------------------------------------------------------------------- //
// This is a crude way to read the time state back in; the proper
// way would be to use the ini parser (so that we can insert new
// data and maintain backward compatibility).
// -------------------------------------------------------------------------- //
void DogSolver::read_time_state(int nframe, const char* outputdir)
{
    FILE* file = open_state_file(outputdir,nframe,"r");
    double time;

    // make sure all these things were read correctly:
    assert( fscanf(file,"time = %lf\n", &time) );
    assert( fscanf(file,"dt = %lf\n", &dt)     );
    assert( fscanf(file,"next_report_time = %lf\n", &next_report_time) );

    positivityState->read_state(file);

    // We do not actually need to do this since it gets overwritten
    // at the beginning of DogSolve
    //fscanf(file,"next_frame_time = %lf\n", &next_frame_time);
    set_time_hack(time);
    fclose(file);
    assert_le(time, next_report_time);

}

// -------------------------------------------------------------------------- //
//
// SubSection: DogSolver RK methods
//
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Advance solution forward one full time frame.
//
// This is the old DogSolveRK method
//
// -------------------------------------------------------------------------- //
void DogSolver::DogSolve(const double tstart, const double tend)
{
    set_next_frame_time(tend);

    assert_eq(tstart, get_state().get_time());

    // maximum allowable time steps for the simulation:
    const int nv = dogParams.get_nv();

    // define local variables
    int n_step = 0;
    double t  = tstart;
    double dt = get_dt();

    // CFL parameters
    const double* cflv = dogParams.get_cflv();
    const double CFL_max    = cflv[1];
    const double CFL_target = cflv[2];

    // Q: where does storage for the stage values get allocated? (-DS).
    // In the old DogSolveRK method (e.g. see 1D code), it was allocated here,
    // and then recycled throughout the time stepping.
    //
    while (t<tend)
    {
        // initialize time step
        bool time_step_completed = false;
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            eprintf(" Error in DogSolve.cpp: "
                    " Exceeded allowed # of time steps \n"
                    "    n_step = %d\n"
                    "        nv = %d\n\n",
                    n_step,nv);
        }

        // save the current state
        //
        // copy qnew into qold
        saveState();

        // keep trying until we get a dt that does not violate CFL condition
        while (!time_step_completed)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }

            // what the time will be at the conclusion of this step
            t = told + dt;

            if(told != get_state_old().get_time())
            {
                dprintf("told1=%24.16e",told);
                dprintf("told2=%24.16e",get_state_old().get_time());
            }
            assert_almost_eq(told, get_state_old().get_time());
            assert_almost_eq(told, get_state().get_time());
            set_dt(dt);

            set_time_hack(told);

            // Set initial maximum wave speeds to zero
            reset_smax();

            // -------------------------------------------------------------------------- //
            // *** main_call ***        
            //
            // In the DogSolveUser case, this is the spot that throws an
            // error.  In principle, this is the only spot that should be
            // modified, but how do we know what storage to allocate? (-DS)
            //        
            int error = advanceFullTimeStep(dt); 
            if(!get_state().check_valid_state()) eprintf("invalid state");
            // -------------------------------------------------------------------------- //

            // compute cfl number
            const double flux_dt
                = dogParams.get_splitting()==SplittingType::fsf ? 0.5*dt : dt;
            double cfl = GetCFL(flux_dt); // accesses smax

            // output time step information
            if (dogParams.get_verbosity()>0) 
            {
                printf("DogSolve2D ... Step %5d"
                        "   CFL =%6.3f"
                        "   dt =%11.3e"
                        "   t =%11.3e\n",
                        n_step,cfl,dt,t);
            }

            // Currently we do not provide the user a mechanism
            // to supply a cfl number for the source term or
            // the split-off operator, so the new time step is
            // determined entirely by the CFL number based on the
            // wave speeds from the Riemann solves.
            //
            // choose new time step
            if (cfl>0.0)
            {   
                const double cfl_target
                    = CFL_target*get_positivityState().get_cflFactor();
                double new_dt = dt*cfl_target/cfl;
                dt = Min(dogParams.get_max_dt(),new_dt);
                dprintf3("new time step: %f",dt);
            }
            else if(error)
            {
                dt = get_positivityState().get_dt();
            }
            else
            {
                dt = dogParams.get_max_dt();
            }

            // see whether to accept or reject this step
            double cfl_max = CFL_max*get_positivityState().get_cflFactor();
            bool satisfied_cfl = (cfl <= cfl_max);
            if (!satisfied_cfl)
            {
                if(dogParams.get_verbosity()>0)
                    printf("DogSolve2D rejecting step..."
                            "CFL number too large\n");
            }
            if(error)
            {
                if (dogParams.get_verbosity()>0)
                {
                    printf("DogSolve2D rejecting step..."
                            "positivity violation\n");
                }
                // take the newly suggested dt.
                //double dt_old = dt;
                //dt = get_positivityState().get_dt();
                //dprintf1("shrunk dt from %6.3f to %6.3f",dt_old,dt);
            }

            if(!error && satisfied_cfl)
            { 
                time_step_completed = true;
                // we do this to prevent the time from
                // drifting by machine epsilon because of
                // inexactness in the solver.
                fetch_state().set_time(t);
                set_time_hack(t);
            }

            if(!time_step_completed) // reject time step
            {   
                // revert to old state
                //
                revertToSavedState();        
                fetch_state().AfterReject(dt);
                t = get_state().get_time();
            }
        }

        // do any extra work
        fetch_state().AfterFullTimeStep(dt);
        if(!get_state().check_valid_state()) eprintf("invalid state");

        //write_restart(9990);
        // compute conservation and print to file
        //
        reportAfterStep(get_state());
    }

    // set initial time step for next call to DogSolveRK
    set_dt(dt);
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// poor-man's "exception" propagation
//
// Why aren't we simply throwing exceptions? (-DS)
//
#define ret_if_err(args...)                \
{ int retcode = args; if(retcode) return retcode; }
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
//
// This call serves as an intermediate step to determine if the user
// implemented a splitting option.
//
// Unless the user implements a time splitting option, this routine
// simply calls advanceTimeStepRK().
//
// -------------------------------------------------------------------------- //
int DogSolver::advanceFullTimeStep(double dt)
{

    using namespace SplittingType;
    switch(dogParams.get_splitting())
    {
        default:
            unsupported_value_error(dogParams.get_splitting())

        case none:

                // TODO: without splitting turned on, we don't want to only
                // select advanceTimeStepRK().
                ret_if_err(advanceTimeStepRK(dt));        
                break;

        case sf: // 1st-order Godunov splitting
                ret_if_err(advanceSplitTimeStep(dt));
                ret_if_err(advanceTimeStepRK(dt));
                break;

        case fs: // 1st-order Godunov splitting
                ret_if_err(advanceTimeStepRK(dt));
                ret_if_err(advanceSplitTimeStep(dt));
                break;

    }
    return 0;
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Take a full Runge-Kutta time step of length dt
//
// returns 0 on success
// called e.g. from DogSolve()
// -------------------------------------------------------------------------- //
int DogSolver::advanceTimeStepRK(double dt)
{
    //dprintf("solving  flux term for dt = %e",dt);
    DogState& Qnew = fetch_state();
    DogState& Qold = fetch_state_old();

    if(dogParams.get_time_order()>=4)
    {
        // if split term was called then Qold!=Qnew, so we would need
        // to allocate another register; but if using high-order in
        // time then probably we are not using splitting anyway.
        // It is conceivable that one might use high-order in time
        // for stability. To enable splitting in this case we would
        // use auto_ptr and allocate Qold and initialize to Qnew if
        // necessary.
        assert_eq(dogParams.get_splitting(),SplittingType::none);
    }

    // in first-order case there is no need to allocate a secondary q register
    if ( dogParams.get_time_order()==1 )
    {      
        ret_if_err(advanceTimeStageRK(1,dt,Qnew,Qnew,Qnew));    
        return 0;
    }

    // ==================================================================== //
    // allocate ancillary register
    linked_ptr<DogState> Qstar_ptr(Qnew.clone(CopyMode::DEEP));
    DogState& Qstar = *Qstar_ptr;
    // ==================================================================== //

    // In most cases (orders 2,3, and 4): 
    //
    // advanceTimeStageRK( dt, q1, q2, q ) produces:
    //     qstar = alpha1 q1 + alpha2 q2 + beta dt L( q1 ).
    //
    switch ( dogParams.get_time_order() )
    {
        case 2: // Second order in time (2-stages)
            {

                // This second order method is supposed to set the following:
                //
                //     qstar = qnew + dt L( qnew )
                //     qnew  = (0.5)*( qnew + qstar ) + 0.5*dt*L( qstar )
                //

                //qstar.setall(0.); // evidently unnecessary
                ret_if_err(advanceTimeStageRK(1, dt, Qnew,  Qnew,  Qstar ));
                ret_if_err(advanceTimeStageRK(2, dt, Qstar, Qstar, Qnew  ));
            }
            break;

        case 3: // Third order, 3-stage TVD method of Osher and Shu.
            // 
            // See: "Total variation diminising Runge-Kutta schemes", Gottlieb and Shu, 1999.
            // http://www.ams.org/journals/mcom/1998-67-221/S0025-5718-98-00913-2/
            //

            //qstar.setall(0.); // evidently unnecessary

            ret_if_err(advanceTimeStageRK( 1, dt, Qnew,  Qnew,  Qstar));
            ret_if_err(advanceTimeStageRK( 2, dt, Qstar, Qnew,  Qstar));
            ret_if_err(advanceTimeStageRK( 3, dt, Qstar, Qstar, Qnew));
            break;

        case 4: // Fourth order in time (10-stages)
            // Low-storage implementation of the optimal SSP(10,4) method.  See:
            // "Highly efficient strong stability-preserving Runge-Kutta methods with 
            // low-storage implementations", Ketcheson, (2008).
            //
            // In particular, we follow the low-storage implementation
            // suggested in Pseudocode 3.  There's a single "stage" that
            // reassigns each register with a linear combination of the
            // two registers being used.
            {
                // rename the registers
                DogState& Q1 = Qnew; // the current state
                DogState& Q2 = Qstar; // the ancillary register

                // Stage: 1,2,3,4, and 5 (these all look like Euler steps)
                for (int s=1; s<=5; s++)
                {
                    ret_if_err(advanceTimeStageRK(s,dt,Q1,Q1,Q1));
                }

                // This step swaps values between the two registers, and
                // doesn't need to evaluate the right hand side.  One can
                // think of this as a "reorganizing" step.
                ret_if_err(advanceTimeStageRK4(dt,Qold,Q1,Q2));

                // Stage: 6,7,8, and 9 (these all look like Euler steps)
                for (int s=6; s<=9; s++)
                {
                    ret_if_err(advanceTimeStageRK(s,dt,Q1,Q1,Q1));
                }

                // Stage: 10
                ret_if_err(advanceTimeStageRK(10,dt,Q1,Q2,Q1));
            }
            break;

        case 5: // Fifth order in time (8-stages)  See same paper from 4th
            // order scheme.
            {
                DogState& Q1 = Qnew;
                DogState& Q2 = Qstar; // Q2.q_setall(0.); (no longer needed)

                for (int s=1; s<=8; s++)
                {
                    ret_if_err(advanceTimeStageRK2(s,dt,Qold,Q1,Q2));
                }
            }
            break;

        default: unsupported_value_error(dogParams.get_time_order());
    }
    return 0;
}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// Take a single Eulerian step of length dt.  These are used for the
// intermediate stage values in a full Runge-Kutta method.
//
// The Runge-Kutta methods implemented in DoGPack are all low storage.
//
// TODO: where is Qstar allocated? (-DS)
//
// returns 0 if successful
//
// called e.g. from DogSolve -> advanceFullTimeStep
//                           -> advanceTimeStepRK   -> advanceTimeStageRK
// -------------------------------------------------------------------------- //
int DogSolver::advanceTimeStageRK(int mstage,
        const double dt,
        DogState& Qinit, // used to set L and smax
        const DogState& Qstar,
        DogState& Qnew)
{
    // Preliminary calls:
    Qinit.BeforeStage(dt);
    set_time_hack(Qinit.get_time()); // needed by ConstructL
    //Qinit.SetBndValues();

    // Set L and smax:
    Qnew.SetBndValues();
    // Check for limiters (moment limiter is applied after each stage):
    if(dogParams.using_moment_limiter())
    { Qnew.ApplyLimiter(); }

    Qinit.ConstructL();

    // Take an Eulerian step:
    // (Why do we need the DogStateTB class? -DS)
    //
    // main_call: update soln using L
    AdvanceStageRK(mstage,dt,Qstar,Qnew);
    // in case user needs it
    set_time_hack(Qnew.get_time()); 

    // Check positivity of cell average:
    fetch_positivityState().clearFlags();
    Qnew.AfterAdvanceStage();

    // Extra optional calls used after each stage:
    Qnew.AfterStage(dt);

    return 0;

}


// Time step advance that does not require evaluation of L
// and does not require limiting (because it generates a convex
// combination of already-limited states)
// 
int DogSolver::advanceTimeStageRK4(double dt,
        const DogState& Qold, // state at start of full time step
        DogState& Q1, // i.e. Qnew, used to set L; defines "current time"
        DogState& Q2) // i.e. Qstar, the ancillary register
{
    assert_eq(dogParams.get_time_order(),4);
    while(true)
    {
        set_time_hack(Q1.get_time()); // (in case BeforeStage needs it)
        Q1.BeforeStage(dt);

        // The original code in cart/DogSolveRK.cpp
        // was equivalent to making the following call with
        // none of the other code in this method.
        AdvanceStageRK4(dt,Qold,Q1,Q2);
        set_time_hack(Q1.get_time()); // (in case AfterStage needs it)

        // issue: AdvanceStageRK4 just defines Q1 to be a convex
        // combination of Qold and Q1, so probably no limiting is
        // necessary, but BeforeStage and AfterStage conceivably
        // could cause limiting to be needed.
        //
        //positivityState->clearFlags();
        //Q1.AfterAdvanceStage(); // e.g. check positivity
        //if(positivityState->get_repeatStepFlag()) return 1;
        //if(positivityState->get_repeatStageFlag()) continue;
        //
        //if (dogParams.using_moment_limiter())
        //  ApplyLimiter(Q1);
        Q1.AfterStage(dt);
        return 0;
    }
}

// it is expected that qnew and q1 point to the same data
int DogSolver::advanceTimeStageRK2(const int mstage,
        double dt,
        const DogState& Qold, // state at start of full time step
        DogState& Q1, // i.e. Qnew, used to set L; defines "current time"
        DogState& Q2) // i.e. Qstar, the ancillary register
{
    while(true)
    {
        set_time_hack(Q1.get_time()); // needed by ConstructL
        Q1.BeforeStage(dt);
        Q1.SetBndValues(); // set boundary conditions
        Q1.ConstructL(); // set L and smax
        AdvanceStageRK2(mstage,dt,Qold,Q1,Q2);
        set_time_hack(Q1.get_time()); // (in case user needs it)
        //
        positivityState->clearFlags();
        Q1.AfterAdvanceStage(); // e.g. check positivity
        if(positivityState->get_repeatStepFlag()) return 1;
        if(positivityState->get_repeatStageFlag()) continue;
        //
        if (dogParams.using_moment_limiter())
            Q1.ApplyLimiter();
        Q1.AfterStage(dt);
        return 0;
    }
}



// -------------------------------------------------------------------------- //
// default implementations of virtual methods
//
// (It is conceivable that the user might want to make use of
// Qstar, but almost always the user will leave the defaults as
// is except to override DogState::advanceSplitTimeStep().)
//
// There were three methods here that I removed, because they were never
// called from anywhere in the code:
//
//     advanceBeforeTimeStep, advanceAfterTimeStep, and
//     advanceBetweenTimeSteps.
//
// There are a few applications in plasma/2d/twofluid
// that use this method, so I'll leave it be:
int DogSolver::advanceSplitTimeStep(double dt)
{ return fetch_state().advanceSplitTimeStep(dt); }
// -------------------------------------------------------------------------- //



// -------------------------------------------------------------------------- //
// This is not currently used for anything.
static int steps_that_have_been_accepted=0;
void increment_steps_that_have_been_accepted()
{ steps_that_have_been_accepted++; }
int get_steps_that_have_been_accepted()
{ return steps_that_have_been_accepted; }

static string get_request_flag_filename()
{
    int frame_idx = dogParams.get_report_frame_idx();
    string flagfilename = get_varframe_filename(
            get_outputdir(), "frame_request_flag.", frame_idx);
    return flagfilename;
}

static bool frame_requested()
{
    int frame_idx = dogParams.get_report_frame_idx();
    if(frame_idx < 0) return false;
    string filename = get_request_flag_filename();
    if(flagfile_check(filename.c_str()))
    {
        return true;
    }
    return false;
}

void DogSolver::reportAfterStep(const DogState& state)
{
    if(frame_requested())
        // matlab wants a plot at the current time
    {
        int frame_idx = dogParams.get_report_frame_idx();
        write_restart(frame_idx);
        // remove request flag (signaling matlab to respond)
        string filename = get_request_flag_filename();
        flagfile_remove(filename.c_str());
    }

    // determine whether we have reached the next subinterval
    const double EPSILON = 1e-16;
    increment_steps_that_have_been_accepted();
    // restrict reporting time-steps to
    // DogParams::num_subintervals many per frame
    //
    // have we attained the time of a new frame?
    if(state.get_time()-next_frame_time >=-EPSILON)
    {
        next_report_time = next_frame_time + dogParams.get_frame_subinterval();
        // this is not really necessary since it will get reset
        // anyway after this call
        next_frame_time += dogParams.get_frame_interval();
    }
    else if(state.get_time() - next_report_time >= -EPSILON)
    {
        next_report_time += dogParams.get_frame_subinterval();
    }
    else
    {
        // do not generate a report
        return;
    }

    state.ReportAfterStep();
}

// -------------------------------------------------------------------------- //
// Section: RKsolver methods
// -------------------------------------------------------------------------- //

// copied from lib/2d/cart/DogSolveRK.cpp
//
// Both lib/1d/DogSolveRK and lib/2d/unst/DogSolveRK_Unst do not use these
// calls.  
//
// See: lib/2d/unst/RunDogpack_Unst.cpp to see the main body of the code.
void RKsolver::set(int time_order)
{
    RKsolver& rk = *this;

    rk.mstage = 0;

    switch( time_order )
    {
        case 5:
            rk.num_stages = 8;
            break;
        case 4:
            rk.num_stages = 10;
            break;
        default:
            rk.num_stages = time_order;
    }

    rk.alpha1 = new dTensor1(rk.num_stages);
    rk.alpha2 = new dTensor1(rk.num_stages);
    rk.beta   = new dTensor1(rk.num_stages);

    // 5th order stuff 
    rk.delta  = new dTensor1(rk.num_stages);
    rk.gamma  = new dTensor2(3, rk.num_stages);

    switch(time_order)
    {
        default:      
            unsupported_value_error(time_order);

        case 1: // first-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            break;

        case 2: // second-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            rk.alpha1->set(2, 0.5 );
            rk.alpha2->set(2, 0.5 );
            rk.beta->set(  2, 0.5 );

            break;

        case 3: // third-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            rk.alpha1->set(2, 0.75 );
            rk.alpha2->set(2, 0.25 );
            rk.beta->set(  2, 0.25 );

            rk.alpha1->set(3, 2.0e0/3.0e0 );
            rk.alpha2->set(3, 1.0e0/3.0e0 );
            rk.beta->set(  3, 2.0e0/3.0e0 );

            break;

        case 4: // fourth-order

            for (int i=1; i<=9; i++)
            {
                rk.alpha1->set(i, 1.0 );
                rk.alpha2->set(i, 0.0 );
                rk.beta->set(  i, 1.0/6.0 );
            }

            rk.alpha1->set(10, 1.0 );
            rk.alpha2->set(10, 3.0/5.0 );
            rk.beta->set(  10, 1.0/10.0 );

            break;

        case 5: // 5th-order

            rk.delta->set(1, 1.0e0);
            rk.delta->set(2, 1.528486658778845e00);
            rk.delta->set(3, 4.720094096662784e-02);
            rk.delta->set(4, 8.801244253465348e-01);
            rk.delta->set(5, 1.019066090228480e+00);
            rk.delta->set(6, 1.049772291176110e+01);
            rk.delta->set(7, -4.254616508506826e+00);
            rk.delta->set(8, 0.0);

            rk.gamma->set(1, 1, 0.0);
            rk.gamma->set(1, 2, -1.552288007713033e+01);
            rk.gamma->set(1, 3,  4.127286635722417e-01);
            rk.gamma->set(1, 4, -1.011819196331377e+00);
            rk.gamma->set(1, 5, -2.765748383780848e-01);
            rk.gamma->set(1, 6,  5.075770311217778e-02);
            rk.gamma->set(1, 7,  6.999810478513669e+00);
            rk.gamma->set(1, 8, -1.114908881433104e+01);

            rk.gamma->set(2, 1,  1.0);
            rk.gamma->set(2, 2,  6.534691420958578e+00);
            rk.gamma->set(2, 3,  2.280056542904473e-01);
            rk.gamma->set(2, 4,  1.308684311397668e+00);
            rk.gamma->set(2, 5,  4.769419552531064e-01);
            rk.gamma->set(2, 6, -6.368809762042849e-03);
            rk.gamma->set(2, 7,  9.339446057238532e-02);
            rk.gamma->set(2, 8,  9.556626047962331e-01);

            rk.gamma->set(3, 1, 0.0);
            rk.gamma->set(3, 2, 0.0);
            rk.gamma->set(3, 3, 0.0);
            rk.gamma->set(3, 4, -2.510747784045939e+00);
            rk.gamma->set(3, 5, -8.576822794622042e-01);
            rk.gamma->set(3, 6,  1.044599944472272e+00);
            rk.gamma->set(3, 7, -7.000810861049136e+00);
            rk.gamma->set(3, 8,  1.906311811144179e+00);

            rk.beta->set(1,  8.653258038183180e-02);
            rk.beta->set(2,  9.544677980851571e-01);
            rk.beta->set(3,  2.651941386774408e-01);
            rk.beta->set(4,  2.736914413910379e-01);
            rk.beta->set(5,  5.999778649323600e-01);
            rk.beta->set(6,  4.433177471748104e-03);
            rk.beta->set(7,  5.309971130968292e-03);
            rk.beta->set(8,  5.830861806762871e-01);
            break;
    }
}

RKsolver::~RKsolver()
{
    delete alpha1;
    delete alpha2;
    delete beta;
    delete gamma;
    delete delta;
}

void RKsolver::advanceTimeStage(const int mstage,
        const double dt,
        const double told, const dTensorBase& qold,
        const dTensorBase& L,
        double &tnew, dTensorBase& qnew)
{
    const double alpha1  = this->alpha1->get(mstage);
    const double alpha2  = this->alpha2->get(mstage);
    const double beta_dt = this->beta->get(mstage)*dt;

    // Update solution
    //
    const int numel = qnew.numel();
    assert_eq(qold.numel(),numel);
    assert_eq(L.numel(),numel);

    if(alpha2==0.) // to accelerate performance
    {
        assert_eq(alpha1,1.);

#pragma omp parallel for
        for(int v=0;v<numel;v++)
        {
            double tmp = qold.vget(v) + beta_dt*L.vget(v);
            qnew.vset(v, tmp );
        }
    }
    else
    {

#pragma omp parallel for
        for(int v=0;v<numel;v++)
        {
            double tmp = alpha1*qold.vget(v) + 
                alpha2*qnew.vget(v) + beta_dt*L.vget(v);
            qnew.vset(v, tmp );
        }
    }
    // advance time (for which L=1)
    tnew = alpha1*told + alpha2*tnew + beta_dt;
}

// -------------------------------------------------------------------------- //
// This routine is used for the fourth-order SSP(10,4) time integrator only.
//
// This is the one exception to the normal "Euler" steps, where the two
// registers are reassigned to each other by taking linear cominations from 
// each other.  This does not require any evaluations of a right hand side
// function.
//
// See 'Pseudocode 3.' in "Highly efficient strong stability-preserving
// Runge-Kutta methods with low-storage implementations", (2008).
//
// -------------------------------------------------------------------------- //
void RKsolver::advanceTimeStage4(const double dt,
        const double told, const dTensorBase& qold,
        double& t1, dTensorBase& q1,
        double& t2, dTensorBase& q2)
{
    const int numel = qold.numel();
    assert_eq(numel, q1.numel());
    assert_eq(numel, q2.numel());

    t2 = (told + 9.0*t1)*(1./25.);
    t1 = 15.0*t2 - 5.0*t1; // i.e., t1 = .6*told+.4*t1

#pragma omp parallel for
    for (int v=0; v<numel;v++)
    {
        const double s1 = q1.vget(v);
        const double tmp = (qold.vget(v) + 9.0*s1)*(1./25.0);
        q2.vset(v, tmp ); // proportional to a convex combination

        // the following is q1.vset(v, .6*qold.vget()+.4*s1), a convex combination
        q1.vset(v, 15.0*tmp - 5.0*s1 );
    }
}

// -------------------------------------------------------------------------- //
//
// Update the solution using the constructed L
// (used in fifth-order method only)
//
// -------------------------------------------------------------------------- //
void RKsolver::advanceTimeStage2(const int mstage,
        const double dt,
        const double told, const dTensorBase& qold,
        double& t1, dTensorBase& q1,
        const dTensorBase& L, // was set using q1
        double& t2, dTensorBase& q2)
{

    const double g1 = gamma->get(1,mstage);
    const double g2 = gamma->get(2,mstage);
    const double g3 = gamma->get(3,mstage);
    const double delta = this->delta->get(mstage);
    const double beta = this->beta->get(mstage);


    const int numel = L.numel();
    assert_eq(numel, qold.numel());
    assert_eq(numel, q1.numel());
    assert_eq(numel, q2.numel());

    // to save expense first test for simple cases
    if(mstage==1)
    {
        // because we deal separately with mstage=1
        // it is unnecessary to initialize t2 and q2 to zero.
        //
        assert_eq(g1,0.);
        assert_eq(g2,1.);
        assert_eq(g3,0.);
        assert_eq(delta,1.);

        // advance time (L=1)
        {
            t2 = t1;
            t1 = t2 + beta*dt;
        }
        // advance state
#pragma omp parallel for
        for(int v=0;v<numel;v++)
        {
            // update q2
            const double s2 = q1.vget(v);
            q2.vset(v, s2);

            // update q
            const double tmp = s2 + beta*dt*L.vget(v);
            q1.vset(v, tmp);
        }
    }
    else if(g3==0.)
    {
        // advance time (L=1)
        {
            t2 = t2 + delta*t1;
            t1 = g1*t1 + g2*t2 + beta*dt;
        }
        // advance state
#pragma omp parallel for
        for(int v=0;v<numel;v++)
        {
            // update q2
            const double s1 = q1.vget(v);
            const double s2 = q2.vget(v) + delta*s1;
            q2.vset(v, s2);

            // update q
            const double tmp = g1*s1 + g2*s2 + beta*dt*L.vget(v);
            q1.vset(v, tmp);
        }
    }
    else // do the full method
    {
        // advance time (L=1)
        {
            t2 = t2 + delta*t1;
            t1 = g1*t1 + g2*t2 + g3*told + beta*dt;
        }
        // advance state
#pragma omp parallel for
        for(int v=0;v<numel;v++)
        {
            // update q2
            const double s1 = q1.vget(v);
            const double s2 = q2.vget(v) + delta*s1;
            q2.vset(v, s2);

            // update q
            const double s3 = qold.vget(v);
            const double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*L.vget(v);
            q1.vset(v, tmp);
        }
    }
}

// -------------------------------------------------------------------------- //
// Section: DogSolverTB methods
// -------------------------------------------------------------------------- //

// What does the acronym "TB" stand for? (-DS)

// This was called UpdateSoln in DogSolveRK.cpp;
// I changed the name to reflect the fact that
// it is specific to the Runge-Kutta method. -eaj
//
void DogSolverTB::AdvanceStageRK(int mstage,
        double dt,
        const DogState& in_Qstar,
        DogState& in_Qnew)
{
    const DogStateTB& Qstar = (const DogStateTB&) in_Qstar;
    DogStateTB& Qnew        =       (DogStateTB&) in_Qnew;

    double tstar = Qstar.get_time();
    double tnew  = Qnew.get_time();
    rk.advanceTimeStage(mstage,dt,
            tstar,Qstar.get_q(),
            get_L(),
            tnew,Qnew.fetch_q());

    Qnew.set_time(tnew);
}


// Update the solution using the constructed L
void DogSolverTB::AdvanceStageRK2(int mstage,
        double dt,
        const DogState& in_Qold,   // state at start of full time step
        DogState& in_Q1,           // always used to set L; defines "current time"
        DogState& in_Q2)           // the ancillary register
{
    const DogStateTB& Qold = (const DogStateTB&) in_Qold;
    DogStateTB& Q1 = (DogStateTB&) in_Q1;
    DogStateTB& Q2 = (DogStateTB&) in_Q2;

    //
    rk.advanceTimeStage2(mstage,dt,
            Qold.get_time(),Qold.get_q(),
            Q1.fetch_time(),Q1.fetch_q(),
            get_L(),
            Q2.fetch_time(),Q2.fetch_q());
}

void DogSolverTB::AdvanceStageRK4(double dt,
        const DogState& in_Qold,    // state at start of full time step
        DogState& in_Q1,            // always used to set L; defines "current time"
        DogState& in_Q2)            // the ancillary register
{
    const DogStateTB& Qold = (const DogStateTB&) in_Qold;
    DogStateTB& Q1 = (DogStateTB&) in_Q1;
    DogStateTB& Q2 = (DogStateTB&) in_Q2;

    // This method is used for the fourth-order SSP(10,4) time integrator
    // only.
    rk.advanceTimeStage4(dt,
            Qold.get_time(),Qold.get_q(),
            Q1.fetch_time(),Q1.fetch_q(),
            Q2.fetch_time(),Q2.fetch_q());
}

void DogSolverTB::reset_smax()
{
    smax->setall(0.);
}

void DogSolverTB::init()
{
    DogSolver::init();
    rk.set(dogParams.get_time_order());
}
