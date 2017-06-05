#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogAsynch1d.h"
using namespace std;

void DogSolveLxW_synch(DogAsynch1d& dogAsynch1d, double tstart, double tend, string outputdir)
{  
  // Output space/time evolution info
  ostringstream fname;
  fname << outputdir << "/" << "spacetime.dat";
  ofstream spt_file;
  if (fabs(tstart)<=1.0e-15) 
    {
      spt_file.open(fname.str().c_str(), ios::out );
    }
  else
    {
      spt_file.open(fname.str().c_str(), ios::app );
    }

  // starting time
  double dt_step = Min(dogAsynch1d.get_min_dt_step(),tend-tstart);
  double tstart_step = tstart;
  double tend_step = tstart + dt_step;

  // ----------------------------------------------------
  // time-step entire mesh from  t = tstart  to  t = tend
  // ----------------------------------------------------
  while (dt_step>0.0)
    {
      // Initialize
      dogAsynch1d.initialize_for_time_stepping(tstart_step,tend_step);

      // --------------------------------------------------------------
      // time-step entire mesh from  t = tstart_step  to  t = tend_step
      // --------------------------------------------------------------
      while (dogAsynch1d.get_num_to_be_updated()>0)
	{	
	  // get next updatable element
	  int i = dogAsynch1d.get_next_updatable_element();
	  
	  // update current element
	  cout << " i = " << setw(4) << i << "    ";
	  cout << " time = " << setw(10) << scientific << setprecision(5) << dogAsynch1d.time->get(i);
	  dogAsynch1d.LxW_update(i);
	  cout << "    time = " << setw(10) << scientific << setprecision(5) << dogAsynch1d.time->get(i);
	  cout << "    CFL = " << setw(10) << scientific << setprecision(5) 
	       << 2.0*dogAsynch1d.dt->get(i)/dogAsynch1d.dx->get(i) << endl;
	  
	  spt_file << setw(4) << i << "      ";
	  spt_file << setw(12) << scientific << dogAsynch1d.time_old->get(i) << "      ";
	  spt_file << setw(12) << scientific << dogAsynch1d.time->get(i) << "      ";
	  spt_file << setw(12) << scientific << dogAsynch1d.node->get(i) << "      ";
	  spt_file << setw(12) << scientific << dogAsynch1d.node->get(i+1) << endl;
	}

      dt_step = Min(dogAsynch1d.get_min_dt_step(),tend-tend_step);
      tstart_step = tend_step;
      tend_step = tstart_step + dt_step;
    }

}
