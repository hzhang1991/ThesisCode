#include "dogdefs.h"
#include "DogParams.h"
#include "DogStateCart4.h"
#include "DogSolverCart4.h"
#include "DogParamsCart4.h"
#include "dog_math.h" // for Min
#include "OutSliceCart4.h"
#include <string>

// idx == 0 specifies main output directory;
// idx >= 1 specifies movie subdirectory
void DogSolverCart4::initOutputDirectory(int idx, const char* framedir)
{

    DogSolver::initOutputDirectory(idx, framedir);

    int mx = dogParamsCart4.get_mx();
    int my = dogParamsCart4.get_my();
    int mz = dogParamsCart4.get_mz();
    int mw = dogParamsCart4.get_mw();

    if(idx!=0)
    {
        mx = dogParamsCart4.get_plot_mx(idx);
        my = dogParamsCart4.get_plot_my(idx);
        mz = dogParamsCart4.get_plot_mz(idx);
        mw = dogParamsCart4.get_plot_mw(idx);
    }

    dogParamsCart4.write_plotting_help_files(mx, my, mz, mw, framedir);
    outSliceCart4.write_plotting_help_files(framedir);
}

void DogSolverCart4::initOutputDirectories()
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

void DogSolverCart4::saveState()
{
    fetch_state_old().copyfrom(get_state());
}

void DogSolverCart4::revertToSavedState()
{
    fetch_state().copyfrom(get_state_old());
}

void DogSolverCart4::init()
{
    if(!&get_state())
    {
        set_state(new DogStateCart4);
        fetch_state().init();
    }

    DogSolverTB::init();

    // parameters specific to the 4D-Cartesian problem:
    // These come from the [grid] section of a parameters.ini file.
    const int mx   = dogParamsCart4.get_mx();
    const int my   = dogParamsCart4.get_my();
    const int mz   = dogParamsCart4.get_mz();
    const int mw   = dogParamsCart4.get_mw();
    const int meqn = dogParams.get_meqn();
    const int kmax = dogParams.get_kmax();
    const int mbc  = dogParamsCart4.get_mbc();

    // maximum speeds.  Again specific to 3D-Cartesian problem.
    smax   = new dTensorBC5(mx, my, mz, mw, 4, mbc, 4);
    DogSolverTB::set_smax(smax);

    // What's this L used for? Is this a single stage for the rhs of an RK time
    // stepper? (-DS)
    DogSolverTB::set_L(L);
    L = new dTensorBC6(mx, my, mz, mw, meqn, kmax, mbc, 4); L->setall(0.);
    DogSolverTB::set_smax(smax);
    DogSolverTB::set_L(L);
}

void DogSolverCart4::initParams()
{
    DogSolver::initParams();
    dogParamsCart4.init();    
    outSliceCart4.init();
}

// this automatically calls ~DogSolver, which
// will take care of calling delete on fetch_state()
// and fetch_state_old()
//
DogSolverCart4::~DogSolverCart4()
{
    delete smax;
    delete L;
}

void DogSolverCart4::write_restart(int num_restart_frame) const
{
    write_time_state(num_restart_frame, get_outputdir());
    get_state().write_frame(num_restart_frame,get_outputdir());
}

// called from DogSolver
void DogSolverCart4::Output(int noutput) const
{
    if (outSliceCart4.get_output4D()==1)
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
        void Output_Extra(const DogSolverCart4& solver, int n);
        Output_Extra( *this, noutput );
    }

    if (outSliceCart4.get_numslices()>0)
    {
        // output slices
        outSliceCart4.write_output(noutput,
                get_state().get_time(),
                get_state().get_q(),
                get_state().get_aux());
    }
}

// called by global method that user can replace
void DogSolverCart4::Restart(int nrestart)
{
    read_time_state(nrestart, get_outputdir());
    fetch_state().read_frame(nrestart);
    fetch_state().SetBndValues();
}
