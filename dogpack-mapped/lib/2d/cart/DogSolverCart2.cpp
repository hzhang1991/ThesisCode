#include "dogdefs.h"
#include "DogParams.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"
#include "DogParamsCart2.h"
#include "dog_math.h" // for Min
#include <string>

// idx == 0 specifies main output directory;
// idx >= 1 specifies movie subdirectory
// (called from PlasmaSolverCart2::initOutputDirectory)
void DogSolverCart2::initOutputDirectory(int idx, const char* framedir)
{
    DogSolver::initOutputDirectory(idx, framedir);

    int mx = dogParamsCart2.get_mx();
    int my = dogParamsCart2.get_my();
    if(idx!=0)
    {
        mx = dogParamsCart2.get_plot_mx(idx);
        my = dogParamsCart2.get_plot_my(idx);
    }
    dogParamsCart2.write_plotting_help_files(mx, my, framedir);
}

void DogSolverCart2::initOutputDirectories()
{
    // ensure that the top-level output directory exists
    // and is populated with plotting information
    initOutputDirectory(0, get_outputdir());

    // ensure that output subdirectories exists and are
    // populated with plotting information
    for(int idx=1;idx<=dogParams.get_how_many_plot_resolutions(); idx++)
    {
        string get_framedir(int movie_idx);
        string framedir = get_framedir(idx);
        initOutputDirectory(idx, framedir.c_str());
    }
}

void DogSolverCart2::saveState()
{
    fetch_state_old().copyfrom(get_state());
}

void DogSolverCart2::revertToSavedState()
{
    fetch_state().copyfrom(get_state_old());
}

void DogSolverCart2::init()
{
    if(!&get_state())
    {
        set_state(new DogStateCart2);
        fetch_state().init();		
    }

    DogSolverTB::init();

    // parameters specific to the 2D-Cartesian problem:
    // These come from the [grid] section of a parameters.ini file.
    const int mx   = dogParamsCart2.get_mx();
    const int my   = dogParamsCart2.get_my();
    const int meqn = dogParams.get_meqn();
    const int kmax = dogParams.get_kmax();
    const int mbc  = dogParamsCart2.get_mbc();

    // maximum speeds.  Again specific to 2D-Cartesian problem.
    smax   = new dTensorBC3(mx,my,2,mbc);
    DogSolverTB::set_smax(smax);

    // What's this L used for? Is this a single stage for the rhs of an RK time
    // stepper? (-DS)
    //
    // From lib/DogSolver.h, the header for this call is: 
    //     void set_L(dTensorBase* in){L = in;}
    // In the class DogSolverTB, there one of the three private field members
    // is dTensorBase* L.
    //
    DogSolverTB::set_L( L );
    L = new dTensorBC4(mx,my,meqn,kmax,mbc);
    L->setall(0.);
    DogSolverTB::set_smax(smax);
    DogSolverTB::set_L(L);
}

void DogSolverCart2::initParams()
{
    DogSolver::initParams();
    dogParamsCart2.init();
}

// this automatically calls ~DogSolver, which
// will take care of calling delete on fetch_state()
// and fetch_state_old()
//
DogSolverCart2::~DogSolverCart2()
{
    delete smax;
    delete L;
}

void DogSolverCart2::write_restart(int num_restart_frame) const
{
    write_time_state(num_restart_frame, get_outputdir());
    get_state().write_frame(num_restart_frame,get_outputdir());
}

// called from DogSolver
void DogSolverCart2::Output(int noutput) const
{
    // write restart file if appropriate
    if(noutput%dogParams.get_nout_per_restart()==0)
    {
        const int num_restart_frame = noutput/dogParams.get_nout_per_restart();
        write_restart(num_restart_frame);
    }
    if(dogParams.get_maintained_restart()>=0)
    {    
        write_restart(dogParams.get_maintained_restart());
    }

    get_state().write_output(noutput);

    // to provide backward compatability (-DS)
    void Output_Extra(const DogSolverCart2& solver, int n);
    Output_Extra( *this, noutput );

}

// called by global method that user can replace
void DogSolverCart2::Restart(int nrestart)
{
    read_time_state(nrestart, get_outputdir());
    fetch_state().read_frame(nrestart);
    fetch_state().SetBndValues();
}

// In earlier versions of this method q was not const and
// SetBndValues was called on q. User is now responsible to call
// SetBndValues before calling this method
//
//void DogSolverCart2::ConstructL(const DogState& state_in)
//{
//    const DogStateCart2& state = (DogStateCart2&) state_in;
//    dTensorBC3& smax = fetch_smax();
//    dTensorBC4& L = fetch_L();
//    //
//    const dTensorBC4& aux = state.get_aux();
//    const dTensorBC4& q = state.get_q();
//
//    ::ConstructL(aux,q,L,smax);
//}
