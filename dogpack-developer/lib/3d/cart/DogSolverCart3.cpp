#include "dogdefs.h"
#include "DogParams.h"
#include "DogStateCart3.h"
#include "DogSolverCart3.h"
#include "DogParamsCart3.h"
#include "dog_math.h" // for Min
#include "OutSliceCart3.h"
#include <string>

using namespace std;

// idx == 0 specifies main output directory;
// idx >= 1 specifies movie subdirectory
void DogSolverCart3::initOutputDirectory(int idx, const char* framedir)
{
  DogSolver::initOutputDirectory(idx, framedir);
  
  int mx = dogParamsCart3.get_mx();
  int my = dogParamsCart3.get_my();
  int mz = dogParamsCart3.get_mz();
  
  if(idx!=0)
    {
      mx = dogParamsCart3.get_plot_mx(idx);
      my = dogParamsCart3.get_plot_my(idx);
      mz = dogParamsCart3.get_plot_mz(idx);
    }
  
  dogParamsCart3.write_plotting_help_files(mx, my, mz, framedir);
  outSliceCart3.write_plotting_help_files(framedir);
}

void DogSolverCart3::initOutputDirectories()
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

void DogSolverCart3::saveState()
{
  fetch_state_old().copyfrom(get_state());
}

void DogSolverCart3::revertToSavedState()
{
  fetch_state().copyfrom(get_state_old());
}

void DogSolverCart3::init()
{
  if(!&get_state())
    {
      set_state(new DogStateCart3);
      fetch_state().init();
    }
  
  DogSolverTB::init();
  
  // parameters specific to the 3D-Cartesian problem:
  // These come from the [grid] section of a parameters.ini file.
  const int mx   = dogParamsCart3.get_mx();
  const int my   = dogParamsCart3.get_my();
  const int mz   = dogParamsCart3.get_mz();
  const int meqn = dogParams.get_meqn();
  const int kmax = dogParams.get_kmax();
  const int mbc  = dogParamsCart3.get_mbc();

  // maximum speeds.  Again specific to 3D-Cartesian problem.
  smax   = new dTensorBC4(mx,my,mz,3,mbc,3);
  DogSolverTB::set_smax(smax);
  
  // What's this L used for? Is this a single stage for the rhs of an RK time
  // stepper? (-DS)
  DogSolverTB::set_L(L);
  L = new dTensorBC5(mx,my,mz,meqn,kmax,mbc,3); L->setall(0.);
  DogSolverTB::set_smax(smax);
  DogSolverTB::set_L(L);
}

void DogSolverCart3::initParams()
{
  DogSolver::initParams();
  dogParamsCart3.init();    
  outSliceCart3.init();
}

// this automatically calls ~DogSolver, which
// will take care of calling delete on fetch_state()
// and fetch_state_old()
//
DogSolverCart3::~DogSolverCart3()
{
  delete smax;
  delete L;
}

void DogSolverCart3::write_restart(int num_restart_frame) const
{
  write_time_state(num_restart_frame, get_outputdir());
  get_state().write_frame(num_restart_frame,get_outputdir());
}

// called from DogSolver
void DogSolverCart3::Output(int noutput) const
{
  if (outSliceCart3.get_output3D()==1)
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
      void Output_Extra(const DogSolverCart3& solver, int n);
      Output_Extra( *this, noutput );
    }
  
  if (outSliceCart3.get_numslices()>0)
    {
      // output slices
      outSliceCart3.write_output(noutput,
				 get_state().get_time(),
				 get_state().get_q(),
				 get_state().get_aux());
    }
}

// called by global method that user can replace
void DogSolverCart3::Restart(int nrestart)
{
  read_time_state(nrestart, get_outputdir());
  fetch_state().read_frame(nrestart);
  fetch_state().SetBndValues();
}
