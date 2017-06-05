#include "DogParams.h"
#include "DogParamsCart3.h"
#include "DogStateCart3.h"
#include "FaceData.h"
#include "L2Project.h"
#include "dogdefs.h"
#include <string>
#include "debug.h"
#include <float.h> // for DBL_MAX

const char* get_outputdir();   // from DogSolver.h

// === read/write routines ===

static void WriteStateASCII(string fname, const dTensorBC5& q, double t)
{  
  fname+=".dat";
  //
  // output values
  //
  FILE* file = fopen(fname.c_str(),"w");
  fprintf(file,"%24.16e\n",t);
  
  for (int ell=1; ell<=q.getsize(5); ell++)
    for (int m=1; m<=q.getsize(4); m++)
      for (int k=1; k<=q.getsize(3); k++)
	for (int j=1; j<=q.getsize(2); j++)
	  for (int i=1; i<=q.getsize(1); i++)
	    {
	      double tmp = q.get(i,j,k,m,ell);
	      fprintf(file,"%24.16e\n",tmp);
	    }
  fclose(file);
}

string get_varframe_filename(const char* dir, const char* varname, int nframe);

void WriteStateArray(const char* framedir, 
		     const char* varname,
		     const dTensorBC5& q, 
		     double t, 
		     int nframe)
{
    string fname = get_varframe_filename(framedir,varname,nframe);
    switch(dogParams.get_datafmt())
    {
        case ASCII:
            WriteStateASCII(fname, q, t);
            break;
        case HDF5:
        default:
            unsupported_value_error(dogParams.get_datafmt());
            break;
    }
}

DogStateCart3::DogStateCart3(const DogStateCart3&in, CopyMode::Enum copyMode) :
    DogStateTB(in)
{
  
    q = in.q->clone(copyMode);
    DogStateTB::set_q(q);

    aux = in.aux->clone(copyMode);
    DogStateTB::set_aux(aux);

}

void DogStateCart3::write_frame(int nframe, const char* framedir) const
{
    //DogState::write(nframe,framedir);
    const double t = get_time();
    WriteStateArray(framedir, "q", *q, t, nframe);
    if(aux->numel() > 0)
    {
        dprintf3("printing aux array");
        WriteStateArray(framedir, "a", *aux, t, nframe);
    }
}

// Initial data from the file qXXXX_restart.dat (where XXXX=nstart)
static double ReadStateArrayASCII(string fname, dTensorBC5& q)
{
    fname += ".dat";
    // ifstream restart_file(fname.c_str(), ios::in);
    FILE* restart_file = fopen(fname.c_str(), "r");
    //if(!restart_file.is_open())
    if(!restart_file) eprintf("could not open for reading: %s\n", fname.c_str());

    double t;
    int numconv = fscanf(restart_file,"%lf",&t);
    assert(numconv==1);
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int mz   = q.getsize(3);
    const int meqn = q.getsize(4);
    const int kmax = q.getsize(5);
    for (int ell=1; ell<=kmax; ell++)
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=mz; k++)   
	  for (int j=1; j<=my; j++)   
	    for (int i=1; i<=mx; i++)
	      {
		double datum;
		int numconv = fscanf(restart_file,"%lf",&datum);
		if(numconv!=1) eprintf("scan failed at "
				       "i=%d, j=%d, m=%d, k=%d in file %s", i,j,m,k,fname.c_str());
		q.set(i,j,k,m,ell, datum);
	      }
    
    double tmp;
    numconv = fscanf(restart_file,"%lf",&tmp);
    assert_eq(numconv,EOF);
    fclose(restart_file);
    return t;
}
double ReadStateArray(const char* dir, 
		      const char* varname,
		      dTensorBC5& q, 
		      int nstart)
{
    double retval;

    string fname = get_varframe_filename(dir,varname,nstart);

    switch(dogParams.get_datafmt())
    {
        case ASCII:
            retval = ReadStateArrayASCII(fname, q);
            break;
        case HDF5:            
        default:
            unsupported_value_error(dogParams.get_datafmt());
            break;
    }
    return retval;
}

void DogStateCart3::read_frame(int nframe)
{
    //DogState::read(nframe, get_outputdir());
    set_time(ReadStateArray(get_outputdir(), "q", *q, nframe));
    if (dogParams.get_maux()>0)
    {
        const double t = ReadStateArray(get_outputdir(), "a", *aux, nframe);
        assert_eq(t, get_time());
    }
}

// Write output to files
void DogStateCart3::write_output(int noutput) const
{
    // write plot frame
    for(int idx = 1; idx<=dogParams.get_how_many_plot_resolutions(); idx++)
    {
        if(noutput%dogParams.get_nout_per_plot()[idx]==0)
        {
            const int num_plot = noutput/dogParams.get_nout_per_plot()[idx];
            // also need to pass idx=plot_type_idx
            const int plot_mx = dogParamsCart3.get_plot_mx(idx);
            const int plot_my = dogParamsCart3.get_plot_my(idx);
	    const int plot_mz = dogParamsCart3.get_plot_mz(idx);
            string get_framedir(int movie_idx);
            string framedir = get_framedir(idx);
            DogStateCart3 coarseState;
            coarseState.init_coarser_state(*this,plot_mx,plot_my,plot_mz);
            coarseState.write_frame(num_plot,framedir.c_str());
            //DogSolver::write_time_info(num_plot,framedir.c_str());
        }
    }
}

// === init/cleanup routines ===

DogStateCart3::~DogStateCart3()
{
    delete q;
    delete aux;
}

void DogStateCart3::init()
{
    const int mx = dogParamsCart3.get_mx();
    const int my = dogParamsCart3.get_my();
    const int mz = dogParamsCart3.get_mz();
    init(mx,my,mz);
}

void DogStateCart3::init(int mx, int my, int mz)
{
    const int meqn  = dogParams.get_meqn();
    const int maux  = dogParams.get_maux();
    const int kmax  = dogParams.get_kmax();
    const int mbc   = dogParamsCart3.get_mbc();

    // Dimension and initialize arrays
    //
    assert(q==0); // should not call init() twice
    assert(aux==0); // should not call init() twice
    q = new dTensorBC5(mx,my,mz,meqn,kmax,mbc,3);
    DogStateTB::set_q(q);
    q->setall(0.);
    aux = new dTensorBC5(mx,my,mz,maux,kmax,mbc,3);
    DogStateTB::set_aux(aux);
    aux->setall(0.);
}

// copyfrom checks that dimensions agree
void DogStateCart3::copyfrom(const DogStateCart3& in)
{
  q->copyfrom(*in.q);
  aux->copyfrom(*in.aux);
  DogState::copyfrom(in);
}

void DogStateCart3::init_coarser_state(const DogStateCart3& in, int mx, int my, int mz)
{
  init(mx,my,mz);
  //ProjectFineToCoarse(*(in.q),*q);
  //ProjectFineToCoarse(*(in.aux),*aux);
  DogState::copyfrom(in);
}

// Reverting to top-level function call for debug_check_condition(), which
// always returns true.  It looks like this was used to double check the
// positivity of the solution as a debugging tool, not as a tool for selecting
// time step size.
//bool DogStateCart3::debug_check_condition() const
//{

    //double check_positivity_of_cell_averages(const dTensorBC4& q);
    //double minval = check_positivity_of_cell_averages(get_q());
    //if(minval<0.)
    //{
    //  return false;
    //}
//  return true;
//}

bool DogStateCart3::check_valid_state() const
{
    const dTensorBC5& aux = *DogStateCart3::aux;
    const dTensorBC5& q = *DogStateCart3::q;
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int mz   = q.getsize(3);
    const int meqn = q.getsize(4);
    const int kmax = q.getsize(5);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(4);

    if(!debug_check_condition())
    {
        dprintf1("debug condition violated");
        return false;
    }

    bool okay=true;
    int bad_i=0, bad_j=0, bad_m=0, bad_k=0, bad_ell=0;

// TODO: it would be nice to only use a single pragma loop here  (-DS)
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
	  for (int k=1; k<=mz; k++)
            for (int m=1; m<=meqn; m++)
	      for (int ell=1; ell<=kmax; ell++)
                {
		  const double val = q.get(i,j,k,m,ell);
		  //assert_isfinite(val);
		  bool is_finite = val >= -DBL_MAX && val <= DBL_MAX;
		  if(!is_finite)
		    // a block preceeded by a critical directive
		    // of a given name can only be executed by a
		    // single thread at a time
		    // (omitting "(NaN)" would use default name)
#pragma omp critical (not_finite)
                    {
                        okay=false;
                        Wprintf("q(%d,%d,%d,%d) is not finite", i,j,k,m,ell);
                        bad_i = i;
                        bad_j = j;
			bad_k = k;
                        bad_m = m;
                        bad_ell = ell;
                    }
                }

    bool aux_okay = true;
    int aux_i=0, aux_j=0, aux_k=0, aux_m=0, aux_ell=0;
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
	  for (int k=1; k<=mz; k++)
            for (int m=1; m<=maux; m++)
	      for (int ell=1; ell<=kmax; ell++)
                {
		  const double val = aux.get(i,j,k,m,ell);
		  //assert_isfinite(val);
		  bool is_finite = val >= -DBL_MAX && val <= DBL_MAX;
		  if(!is_finite)
#pragma omp critical (not_finite)
                    {
		      aux_okay=false;
		      Wprintf("aux(%d,%d,%d,%d) is not finite", i,j,k,m,ell);
		      aux_i = i;
		      aux_j = j;
		      aux_k = k;
		      aux_m = m;
		      aux_ell = ell;
                    }
                }
    if(!okay || !aux_okay)
      {
        write_frame(9999,get_outputdir());
        if(!okay)
            errmsg_printf("at time %e q(%d,%d,%d,%d) is not finite", get_time(),
			  bad_i,bad_j,bad_k,bad_m,bad_ell);
        if(!aux_okay)
            errmsg_printf("at time %e aux(%d,%d,%d,%d) is not finite", get_time(),
			  aux_i,aux_j,aux_k,aux_m,aux_ell);
    }
    return true;
}

#include "DogSolverCart3.h"
//
// user "callbacks"
// (AppStateCart3 should override these deprecated implementations)
//
void DogStateCart3::SetBndValues()
{

  // Hackish way to get time passed into the boundary routine.
  // What's the correct way to do this?
  dogParams.set_time( get_time() );
  
  // map onto deprecated global method
  void SetBndValues(dTensorBC5&, dTensorBC5&);
  SetBndValues(fetch_q(),fetch_aux());
}
void DogStateCart3::BeforeStage(double dt)
{
    // map onto deprecated global method
    //void ::BeforeStep(double, dTensorBC4&, dTensorBC4&, DogSolverCart3&);
    ::BeforeStep(dt,fetch_aux(),fetch_q(),fetch_solver());
}

// User is responsible to call SetBndValues before calling this method
//
void DogStateCart3::ConstructL()
{
    dTensorBC4& smax = fetch_solver().fetch_smax();
    dTensorBC5& L    = fetch_solver().fetch_L();
    ::ConstructL(fetch_aux(), fetch_q(), L, smax);
}

void DogStateCart3::ApplyLimiter()
{
    // map onto deprecated global method
    void ProjectRightEig(int,
			 const dTensor1&, 
			 const dTensor1&,
			 const dTensor2&,
			 dTensor2&);
    void ProjectLeftEig(int,
			const dTensor1&, 
			const dTensor1&,
			const dTensor2&,
			dTensor2&);
    ::ApplyLimiter(fetch_aux(),fetch_q(),ProjectRightEig,ProjectLeftEig);
}

void DogStateCart3::AfterStage(double dt)
{
    // map onto deprecated global method
    ::AfterStep(dt,fetch_aux(),fetch_q(),fetch_solver());
}

void DogStateCart3::AfterReject(double dt)
{
    // map onto deprecated global method
    ::AfterReject(dt,fetch_aux(),fetch_q(),fetch_solver());
}

void DogStateCart3::AfterFullTimeStep(double dt)
{
    // map onto deprecated global method
    ::AfterFullTimeStep(fetch_solver());
}

void DogStateCart3::ReportAfterStep() const
{
    // map onto deprecated global method
    ::ConSoln(get_aux(),get_q(),get_time());
}

void DogStateCart3::AfterInitState()
{
    // map onto deprecated global method
    ::AfterQinit(fetch_solver());
}

// Wrapper functions so we can make cleaner calls to user supplied functions
static void AuxFunc(const dTensor2& xpts,
		    const dTensor2& NOT_USED_1,
		    const dTensor2& NOT_USED_2,
		    dTensor2& auxvals)
{
  void AuxFunc(const dTensor2& xpts,dTensor2& auxvals);
  AuxFunc(xpts,auxvals);
}

static void QinitFunc(const dTensor2& xpts,
		      const dTensor2& NOT_USED_1,
		      const dTensor2& NOT_USED_2,
		      dTensor2& qvals)
{
  void QinitFunc(const dTensor2& xpts,dTensor2& qvals);
  QinitFunc(xpts,qvals);
}

// Routine to project initial conditions onto state vector
void DogStateCart3::InitState()
{
  // Fetch gives write access
  dTensorBC5& qnew = fetch_q();
  dTensorBC5& aux  = fetch_aux();  
  
  // Cartesian parameters
  const int mx  = dogParamsCart3.get_mx();
  const int my  = dogParamsCart3.get_my();
  const int mz  = dogParamsCart3.get_mz();
  const int mbc = dogParamsCart3.get_mbc();
  
  // parameters used for all solvers
  const int ic_quad_order = dogParams.get_ic_quad_order();
  const int kmax          = dogParams.get_kmax();
  const int space_order   = dogParams.get_space_order();
  
  // Set any auxiliary variables on computational grid
  // Set values using L2-projection
  void AuxFunc(const dTensor2& xpts,
	       const dTensor2& NOT_USED_1,
	       const dTensor2& NOT_USED_2,
	       dTensor2& auxvals);
  if (dogParams.get_maux()>0)
    {
      L2Project(1-mbc,mx+mbc,
		1-mbc,my+mbc,
		1-mbc,mz+mbc,
		ic_quad_order,
		space_order,
		space_order,
		space_order,
		&qnew,
		&aux,
		&aux,
		&AuxFunc);  
    }    
    
  // Set initial data on computational grid using L2-projection
  void QinitFunc(const dTensor2& xpts,
		 const dTensor2& NOT_USED_1,
		 const dTensor2& NOT_USED_2,
		 dTensor2& qvals);
  L2ProjectInitialCond(1-mbc,mx+mbc,
		       1-mbc,my+mbc,
		       1-mbc,mz+mbc,
		       ic_quad_order,
		       space_order,
		       space_order,
		       space_order,
		       &qnew,
		       &aux,
		       &qnew,
		       &QinitFunc);
}
