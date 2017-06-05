#include "DogParams.h"
#include "DogParamsCart4.h"
#include "DogStateCart4.h"
#include "FaceData.h"
#include "L2Project.h"
#include "dogdefs.h"
#include <string>
#include "debug.h"
#include <float.h> // for DBL_MAX

const char* get_outputdir();   // from DogSolver.h

// === read/write routines ===

static void WriteStateASCII(string fname, const dTensorBC6& q, double t)
{  
    fname+=".dat";
    //
    // output values
    //
    FILE* file = fopen(fname.c_str(),"w");
    fprintf(file,"%24.16e\n",t);

    for (int ell=1; ell<=q.getsize(6); ell++)
    for (int m=1; m<=q.getsize(5); m++)
    for (int l=1; l<=q.getsize(4); l++)
    for (int k=1; k<=q.getsize(3); k++)
    for (int j=1; j<=q.getsize(2); j++)
    for (int i=1; i<=q.getsize(1); i++)
    {
        double tmp = q.get(i,j,k,l,m,ell);
        fprintf(file,"%24.16e\n",tmp);
    }
    fclose(file);
}

string get_varframe_filename(const char* dir, const char* varname, int nframe);

void WriteStateArray(const char* framedir, 
        const char* varname,
        const dTensorBC6& q, 
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

DogStateCart4::DogStateCart4(const DogStateCart4&in, CopyMode::Enum copyMode) :
    DogStateTB(in)
{

    q = in.q->clone(copyMode);
    DogStateTB::set_q(q);

    aux = in.aux->clone(copyMode);
    DogStateTB::set_aux(aux);

}

void DogStateCart4::write_frame(int nframe, const char* framedir) const
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
static double ReadStateArrayASCII(string fname, dTensorBC6& q)
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
    const int mw   = q.getsize(4);
    const int meqn = q.getsize(5);
    const int kmax = q.getsize(6);
    for (int ell=1; ell<=kmax; ell++)
    for (int m=1; m<=meqn; m++)
    for (int l=1; l<=mw; l++)   
    for (int k=1; k<=mz; k++)   
    for (int j=1; j<=my; j++)   
    for (int i=1; i<=mx; i++)
    {
        double datum;
        int numconv = fscanf(restart_file,"%lf",&datum);
        if(numconv!=1) eprintf("scan failed at "
                "i=%d, j=%d, k=%d, l=%d, m=%d, k=%d in file %s", i,j,k,l,m, fname.c_str());
        q.set(i,j,k,l,m,ell, datum);
    }

    double tmp;
    numconv = fscanf(restart_file,"%lf",&tmp);
    assert_eq(numconv,EOF);
    fclose(restart_file);
    return t;
}
double ReadStateArray(const char* dir, 
        const char* varname,
        dTensorBC6& q, 
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

void DogStateCart4::read_frame(int nframe)
{
    set_time(ReadStateArray(get_outputdir(), "q", *q, nframe));
    if (dogParams.get_maux()>0)
    {
        const double t = ReadStateArray(get_outputdir(), "a", *aux, nframe);
        assert_eq(t, get_time());
    }
}

// Write output to files
void DogStateCart4::write_output(int noutput) const
{
    // write plot frame
    for(int idx = 1; idx<=dogParams.get_how_many_plot_resolutions(); idx++)
    {
        if(noutput%dogParams.get_nout_per_plot()[idx]==0)
        {
            const int num_plot = noutput/dogParams.get_nout_per_plot()[idx];

            // also need to pass idx=plot_type_idx
            const int plot_mx = dogParamsCart4.get_plot_mx(idx);
            const int plot_my = dogParamsCart4.get_plot_my(idx);
            const int plot_mz = dogParamsCart4.get_plot_mz(idx);
            const int plot_mw = dogParamsCart4.get_plot_mw(idx);

            string get_framedir(int movie_idx);
            string framedir = get_framedir(idx);
            DogStateCart4 coarseState;
            coarseState.init_coarser_state(*this,plot_mx,plot_my,plot_mz,plot_mw);
            coarseState.write_frame(num_plot,framedir.c_str());
        }
    }
}

// === init/cleanup routines ===

DogStateCart4::~DogStateCart4()
{
    delete q;
    delete aux;
}

void DogStateCart4::init()
{
    const int mx = dogParamsCart4.get_mx();
    const int my = dogParamsCart4.get_my();
    const int mz = dogParamsCart4.get_mz();
    const int mw = dogParamsCart4.get_mw();
    init(mx,my,mz,mw);
}

void DogStateCart4::init(int mx, int my, int mz, int mw)
{
    const int meqn  = dogParams.get_meqn();
    const int maux  = dogParams.get_maux();
    const int kmax  = dogParams.get_kmax();
    const int mbc   = dogParamsCart4.get_mbc();

    // Dimension and initialize arrays
    //
    assert(q==0); // should not call init() twice
    assert(aux==0); // should not call init() twice
    q = new dTensorBC6(mx,my,mz,mw,meqn,kmax,mbc,4);
    DogStateTB::set_q(q);
    q->setall(0.);
    aux = new dTensorBC6(mx,my,mz,mw,maux,kmax,mbc,4);
    DogStateTB::set_aux(aux);
    aux->setall(0.);
}

// copyfrom checks that dimensions agree
void DogStateCart4::copyfrom(const DogStateCart4& in)
{
    q->copyfrom(*in.q);
    aux->copyfrom(*in.aux);
    DogState::copyfrom(in);
}

void DogStateCart4::init_coarser_state(const DogStateCart4& in, 
    int mx, int my, int mz, int mw)
{
    init(mx,my,mz,mw);
    DogState::copyfrom(in);
}

bool DogStateCart4::check_valid_state() const
{

    const dTensorBC6& aux = *DogStateCart4::aux;
    const dTensorBC6& q   = *DogStateCart4::q;
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int mz   = q.getsize(3);
    const int mw   = q.getsize(4);
    const int meqn = q.getsize(5);
    const int kmax = q.getsize(6);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(5);

    if(!debug_check_condition())
    {
        dprintf1("debug condition violated");
        return false;
    }

    bool okay=true;
    int bad_i=0, bad_j=0, bad_m=0, bad_k=0, bad_l=0, bad_ell=0;

    // TODO: it would be nice to only use a single pragma loop here  (-DS)
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
    for (int k=1; k<=mz; k++)
    for (int l=1; l<=mw; l++)
    for (int m=1; m<=meqn; m++)
    for (int ell=1; ell<=kmax; ell++)
    {
        const double val = q.get(i,j,k,l,m,ell);
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
            Wprintf("q(%d,%d,%d,%d,%d) is not finite", i,j,k,l,m,ell);
            bad_i = i;
            bad_j = j;
            bad_k = k;
            bad_l = l;
            bad_m = m;
            bad_ell = ell;
        }
    }

    bool aux_okay = true;
    int aux_i=0, aux_j=0, aux_k=0, aux_l=0, aux_m=0, aux_ell=0;
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
    for (int k=1; k<=mz; k++)
    for (int l=1; l<=mw; l++)
    for (int m=1; m<=maux; m++)
    for (int ell=1; ell<=kmax; ell++)
    {
        const double val = aux.get(i,j,k,l,m,ell);
        //assert_isfinite(val);
        bool is_finite = val >= -DBL_MAX && val <= DBL_MAX;
        if(!is_finite)
#pragma omp critical (not_finite)
        {
            aux_okay=false;
            Wprintf("aux(%d,%d,%d,%d,%d) is not finite", i,j,k,l,m,ell);
            aux_i = i;
            aux_j = j;
            aux_k = k;
            aux_l = l;
            aux_m = m;
            aux_ell = ell;
        }
    }
    if(!okay || !aux_okay)
    {
        write_frame(9999,get_outputdir());
        if(!okay)
            errmsg_printf("at time %e q(%d,%d,%d,%d,%d) is not finite", get_time(),
                    bad_i,bad_j,bad_k,bad_l,bad_m,bad_ell);
        if(!aux_okay)
            errmsg_printf("at time %e aux(%d,%d,%d,%d,%d) is not finite", get_time(),
                    aux_i,aux_j,aux_k,aux_l,aux_m,aux_ell);
    }
    return true;
}

#include "DogSolverCart4.h"
//
// user "callbacks"
// (AppStateCart4 should override these deprecated implementations)
//
void DogStateCart4::SetBndValues()
{

    // Hackish way to get time passed into the boundary routine.
    // What's the correct way to do this?
    dogParams.set_time( get_time() );

    // map onto deprecated global method
    void SetBndValues(dTensorBC6&, dTensorBC6&);
    SetBndValues(fetch_q(),fetch_aux());
}

void DogStateCart4::BeforeStage(double dt)
{
    ::BeforeStep(dt, fetch_solver());
}

// User is responsible to call SetBndValues before calling this method
//
void DogStateCart4::ConstructL()
{
    dTensorBC5& smax = fetch_solver().fetch_smax();
    dTensorBC6& L    = fetch_solver().fetch_L();
    ::ConstructL(fetch_aux(), fetch_q(), L, smax);
}

void DogStateCart4::ApplyLimiter()
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
    ::ApplyLimiter(fetch_aux(), fetch_q(), ProjectRightEig, ProjectLeftEig);
}

void DogStateCart4::AfterStage(double dt)
{
    // map onto deprecated global method
    ::AfterStep(dt, fetch_solver());
}

void DogStateCart4::AfterReject(double dt)
{
    // map onto deprecated global method
    ::AfterReject(dt, fetch_solver());
}

void DogStateCart4::AfterFullTimeStep(double dt)
{
    // map onto deprecated global method
    ::AfterFullTimeStep(fetch_solver());
}

void DogStateCart4::ReportAfterStep() const
{
    // map onto deprecated global method
//  ::ConSoln(get_aux(),get_q(),get_time());
//  void ConSoln(const DogSolverCart4& solver);
//  ConSoln(get_solver());
// TODO - get something to work here
}

void DogStateCart4::AfterInitState()
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
//
// TODO - write this routine!
void DogStateCart4::InitState()
{
    // Fetch gives write access
    dTensorBC6& qnew = fetch_q();
    dTensorBC6& aux  = fetch_aux();  

    // Cartesian parameters
    const int mx  = dogParamsCart4.get_mx();
    const int my  = dogParamsCart4.get_my();
    const int mz  = dogParamsCart4.get_mz();
    const int mw  = dogParamsCart4.get_mw();
    const int mbc = dogParamsCart4.get_mbc();

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
//  if (dogParams.get_maux()>0)
//  {
//      L2Project(
//          1-mbc,mx+mbc,
//          1-mbc,my+mbc,
//          1-mbc,mz+mbc,
//          1-mbc,mw+mbc,
//          ic_quad_order,
//          space_order,
//          space_order,
//          space_order,
//          &qnew,
//          &aux,
//          &aux,
//          &AuxFunc);  
//  }    

    // Set initial data on computational grid using L2-projection
    void QinitFunc(const dTensor2& xpts,
            const dTensor2& NOT_USED_1,
            const dTensor2& NOT_USED_2,
            dTensor2& qvals);
// TODO
//  L2ProjectInitialCond(
//      1-mbc,mx+mbc,
//      1-mbc,my+mbc,
//      1-mbc,mz+mbc,
//      1-mbc,mw+mbc,
//      ic_quad_order,
//      space_order,
//      space_order,
//      space_order,
//      &qnew,
//      &aux,
//      &qnew,
//      &QinitFunc);
}
