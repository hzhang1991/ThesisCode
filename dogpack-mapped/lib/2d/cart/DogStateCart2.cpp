#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogStateCart2.h"
#include "edge_data.h"
#include "dogdefs.h"
#include "L2Project.h"
#include <string>
#include "debug.h"
#include <float.h> // for DBL_MAX

const char* get_outputdir();   // from DogSolver.h

// === mesh coarsening/refinement projection routines ===

// compute average-integral of product of fine basis with coarse basis
//
// This needs to be changed to work with second-order basis
//
static void compute_avg_fb_cb(dTensor4& fc, double& dx, double& dy)
{
    fc.setall(0.); // so we only have to fill in the nonzero inner products
    const int mkx = fc.getsize(2);
    const int mky = fc.getsize(3);

    // compute width of fine cells in canonical coordinates of coarse cells.
    dx = 2./mkx;
    dy = 2./mky;
    for(int kx=1; kx<=mkx;kx++)
        for(int ky=1; ky<=mky;ky++)
        {
            const double xk = -1+kx*dx-dx/2.;
            const double yk = -1+ky*dy-dy/2.;
            switch(dogParams.get_space_order())
            {
                case 3:
                    if(mky!=2) // exactly zero in case mkx==2
                        fc.set(1,kx,ky,6, (sq5/8.)*(12.*yk*yk+dy*dy-4.));
                    if(mkx!=2) // exactly zero in case mkx==2
                        fc.set(1,kx,ky,5, (sq5/8.)*(12.*xk*xk+dx*dx-4.));
                    fc.set(1,kx,ky,4, 3.*xk*yk);
                case 2:
                    fc.set(1,kx,ky,3, sq3*yk);
                    fc.set(1,kx,ky,2, sq3*xk);
                case 1:
                    fc.set(1,kx,ky,1, 1.);
                    break;
                default: unsupported_value_error(dogParams.get_space_order());
            }
            //
            switch(dogParams.get_space_order())
            {
                case 3:
                    fc.set(6,kx,ky,6, dy*dy/4.);
                    fc.set(5,kx,ky,5, dx*dx/4.);
                    fc.set(4,kx,ky,4, dx*dy/4.);
                    fc.set(3,kx,ky,4, (sq3/2.)*xk*dy);
                    fc.set(3,kx,ky,6, (sq15/2.)*yk*dy);
                    fc.set(2,kx,ky,5, (sq15/2.)*xk*dx);
                    fc.set(2,kx,ky,4, (sq3/2.)*yk*dx);
                case 2:
                    fc.set(3,kx,ky,3, dy/2.);
                    fc.set(2,kx,ky,2, dx/2.);
                case 1:
                    break;
                default: unsupported_value_error(dogParams.get_space_order());
                         //
            }
        }
}

// This shows agreement with the integrals
// computed numerically in dogpack_notes/refine2D.m
//
static void test_avg_fb_cb()
{
    const int mbas = dogParams.get_kmax();
    const int mkx = 4;
    const int mky = 4;
    dTensor4 avg_fb_cb(mbas,mkx,mky,mbas);
    double dx,dy;
    compute_avg_fb_cb(avg_fb_cb,dx,dy);

    dprintf("=== displaying inner products of fine basis with coarse basis ===");
    for(int kx=1; kx<=mkx;kx++)
        for(int ky=1; ky<=mky;ky++)
        {
            printf("for kx=%d, ky=%d\n",kx,ky);
            for(int fi=1;fi<=mbas;fi++)
            {
                for(int ci=1;ci<=mbas;ci++)
                {
                    printf(" %12.4g", avg_fb_cb.get(fi,kx,ky,ci));
                }
                printf("\n");
            }
        }
    exit(0);
}

// crb_fb = transpose of avg_fb_cb times dx*dy/4.
static void compute_crb_fb(dTensor4& crb_fb)
{
    const int mbas = dogParams.get_kmax();
    assert_eq(mbas, crb_fb.getsize(1));
    assert_eq(mbas, crb_fb.getsize(4));
    const int mkx = crb_fb.getsize(2);
    const int mky = crb_fb.getsize(3);
    dTensor4 avg_fb_cb(mbas,mkx,mky,mbas);
    double dx,dy;
    compute_avg_fb_cb(avg_fb_cb,dx,dy);
    const double dxdy_4 = dx*dy/4.;

    for(int ci=1; ci<=mbas; ci++)
        for(int kx=1; kx<=mkx; kx++)
            for(int ky=1; ky<=mky; ky++)
                for(int fi=1; fi<=mbas; fi++)
                {
                    crb_fb.set(ci,kx,ky,fi, avg_fb_cb.get(fi,kx,ky,ci)*dxdy_4);
                }
}

static void ProjectFineToCoarse(const dTensorBC4& qf,dTensorBC4& qc)
{
    const int mx   = qf.getsize(1);
    const int my   = qf.getsize(2);
    const int meqs = qf.getsize(3);
    const int mbas = qf.getsize(4);
    assert_eq(meqs,qc.getsize(3));
    assert_eq(qf.getsize(4),qc.getsize(4));

    const int plot_mx   = qc.getsize(1);
    const int plot_my   = qc.getsize(2);

    const int plot_rx = mx/plot_mx;
    const int plot_ry = my/plot_my;
    assert_eq(mx%plot_mx,0);
    assert_eq(my%plot_my,0);

    //qc.setall(0.); // probably unnecessary

    dTensor4 crb_fb(mbas,plot_rx,plot_ry,mbas);
    compute_crb_fb(crb_fb);

    // compute innerprod of coarse reciprocal basis with fine basis
    //
#if( NDIMS==2 )
#pragma omp parallel for
#endif
    // iterate over coarse cells
    for(int mx=1;mx<=plot_mx;mx++)
        for(int my=1;my<=plot_my;my++)
            for(int me=1;me<=meqs;me++)
                for(int ci=1;ci<=mbas;ci++)
                {
                    double acc=0.;
                    // sum over all fine cells within coarse cell ...
                    int nx = (mx-1)*plot_rx+1;
                    for(int kx=1;kx<=plot_rx;kx++,nx++)
                    {
                        int ny = (my-1)*plot_ry+1;
                        for(int ky=1;ky<=plot_rx;ky++,ny++)
                        {
                            //assert_eq(nx, kx+(mx-1)*plot_rx);
                            //assert_eq(ny, ky+(my-1)*plot_ry);
                            // ... and over all basis elements
                            for(int fi=1;fi<=ci;fi++) /* ci rather than mbas saves computation */
                            {
                                {
                                    // accelerate? (6 out of 21 entries are zero,
                                    // 8 if plot_mx and plot_my equal 2)
                                    acc += crb_fb.get(ci,kx,ky,fi)*qf.get(nx,ny,me,fi);
                                }
                            }
                        }
                    }
                    qc.set(mx,my,me,ci, acc);
                }
}

// === read/write routines ===

static void WriteStateASCII(string fname, const dTensorBC4& q, double t)
{  
    fname+=".dat";
    //
    // output values
    //
    FILE* file = fopen(fname.c_str(),"w");
    fprintf(file,"%24.16e\n",t);

    for (int k=1; k<=q.getsize(4); k++)
    for (int m=1; m<=q.getsize(3); m++)
    for (int j=1; j<=q.getsize(2); j++)
    for (int i=1; i<=q.getsize(1); i++)
    {
        double tmp = q.get(i,j,m,k);
        fprintf(file,"%24.16e\n",tmp);
    }
    fclose(file);

}

string get_varframe_filename(const char* dir, const char* varname, int nframe);

void WriteStateArray(const char* framedir, const char* varname,
        const dTensorBC4& q, double t, int nframe)
{
    string fname = get_varframe_filename(framedir,varname,nframe);
    switch(dogParams.get_datafmt())
    {
        case ASCII:
            WriteStateASCII(fname, q, t);
            break;
        case HDF5:
            void WriteStateHDF5(string fname, string varname,
                    const dTensorBC4& q, double t);
            WriteStateHDF5(fname, varname, q, t);
            break;
        default:
            unsupported_value_error(dogParams.get_datafmt());
            break;
    }
}

DogStateCart2::DogStateCart2(const DogStateCart2&in, CopyMode::Enum copyMode) :
    DogStateTB(in)
{
  
    q = in.q->clone(copyMode);
    DogStateTB::set_q(q);

    aux = in.aux->clone(copyMode);
    DogStateTB::set_aux(aux);

}

void DogStateCart2::write_frame(int nframe, const char* framedir) const
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
static double ReadStateArrayASCII(string fname, dTensorBC4& q)
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
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    for (int k=1; k<=kmax; k++)
        for (int m=1; m<=meqn; m++)
            for (int j=1; j<=my; j++)   
                for (int i=1; i<=mx; i++)
                {
                    double datum;
                    int numconv = fscanf(restart_file,"%lf",&datum);
                    if(numconv!=1) eprintf("scan failed at "
                            "i=%d, j=%d, m=%d, k=%d in file %s", i,j,m,k,fname.c_str());
                    q.set(i,j,m,k, datum);
                }

    double tmp;
    numconv = fscanf(restart_file,"%lf",&tmp);
    assert_eq(numconv,EOF);
    fclose(restart_file);
    return t;
}
double ReadStateArray(const char* dir, const char* varname,
        dTensorBC4& q, int nstart)
{
    double retval;

    string fname = get_varframe_filename(dir,varname,nstart);

    switch(dogParams.get_datafmt())
    {
        case ASCII:
            retval = ReadStateArrayASCII(fname, q);
            break;
        case HDF5:
            double ReadStateArrayHDF5(string fname, string varname, dTensorBC4& q);
            retval = ReadStateArrayHDF5(fname, varname, q);
            break;
        default:
            unsupported_value_error(dogParams.get_datafmt());
            break;
    }
    return retval;
}

void DogStateCart2::read_frame(int nframe)
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
void DogStateCart2::write_output(int noutput) const
{
    // write plot frame
    for(int idx = 1; idx<=dogParams.get_how_many_plot_resolutions(); idx++)
    {
        if(noutput%dogParams.get_nout_per_plot()[idx]==0)
        {
            const int num_plot = noutput/dogParams.get_nout_per_plot()[idx];
            // also need to pass idx=plot_type_idx
            const int plot_mx = dogParamsCart2.get_plot_mx(idx);
            const int plot_my = dogParamsCart2.get_plot_my(idx);
            string get_framedir(int movie_idx);
            string framedir = get_framedir(idx);
            DogStateCart2 coarseState;
            coarseState.init_coarser_state(*this,plot_mx,plot_my);
            coarseState.write_frame(num_plot,framedir.c_str());
            //DogSolver::write_time_info(num_plot,framedir.c_str());
        }
    }
}

// === init/cleanup routines ===

DogStateCart2::~DogStateCart2()
{
    delete q;
    delete aux;
}

void DogStateCart2::init()
{
    const int mx = dogParamsCart2.get_mx();
    const int my = dogParamsCart2.get_my();
    init(mx,my);
}

void DogStateCart2::init(int mx, int my)
{
    const int meqn  = dogParams.get_meqn();
    const int maux  = dogParams.get_maux();
    const int kmax  = dogParams.get_kmax();
    const int mbc   = dogParamsCart2.get_mbc();

    // Dimension and initialize arrays
    //
    assert(q==0); // should not call init() twice
    assert(aux==0); // should not call init() twice
    q = new dTensorBC4(mx,my,meqn,kmax,mbc);
    DogStateTB::set_q(q);
    q->setall(0.);
    aux = new dTensorBC4(mx,my,maux,kmax,mbc);
    DogStateTB::set_aux(aux);
    aux->setall(0.);
}

// copyfrom checks that dimensions agree
void DogStateCart2::copyfrom(const DogStateCart2& in)
{
    q->copyfrom(*in.q);
    aux->copyfrom(*in.aux);
    DogState::copyfrom(in);
}

void DogStateCart2::init_coarser_state(const DogStateCart2& in, int mx, int my)
{
    init(mx,my);
    ProjectFineToCoarse(*(in.q),*q);
    ProjectFineToCoarse(*(in.aux),*aux);
    DogState::copyfrom(in);
}

// Reverting to top-level function call for debug_check_condition(), which
// always returns true.  It looks like this was used to double check the
// positivity of the solution as a debugging tool, not as a tool for selecting
// time step size.
//bool DogStateCart2::debug_check_condition() const
//{

    //double check_positivity_of_cell_averages(const dTensorBC4& q);
    //double minval = check_positivity_of_cell_averages(get_q());
    //if(minval<0.)
    //{
    //  return false;
    //}
//  return true;
//}

bool DogStateCart2::check_valid_state() const
{
    const dTensorBC4& aux = *DogStateCart2::aux;
    const dTensorBC4& q = *DogStateCart2::q;
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();

    //
    // Something wrong here -James (3/13/2013)
    //
    if(!debug_check_condition())
    {
        dprintf1("debug condition violated");
        return false;
    }

    bool okay=true;
    int bad_i=0, bad_j=0, bad_m=0, bad_k=0;

// TODO: it would be nice to only use a single pragma loop here  (-DS)
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
            for (int m=1; m<=meqn; m++)
                for (int k=1; k<=kmax; k++)
                {
                    const double val = q.get(i,j,m,k);
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
                        Wprintf("q(%d,%d,%d,%d) is not finite \n", i,j,m,k);
                        printf("kmax=%d %d\n",k,kmax);
                        bad_i = i;
                        bad_j = j;
                        bad_m = m;
                        bad_k = k;
                    }
                }

    bool aux_okay = true;
    int aux_i=0, aux_j=0, aux_m=0, aux_k=0;
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
            for (int m=1; m<=aux.getsize(3); m++)
                for (int k=1; k<=kmax; k++)
                {
                    const double val = aux.get(i,j,m,k);
                    //assert_isfinite(val);
                    bool is_finite = val >= -DBL_MAX && val <= DBL_MAX;
                    if(!is_finite)
#pragma omp critical (not_finite)
                    {
                        aux_okay=false;
                        Wprintf("aux(%d,%d,%d,%d) is not finite", i,j,m,k);
                        aux_i = i;
                        aux_j = j;
                        aux_m = m;
                        aux_k = k;
                    }
                }
    if(!okay || !aux_okay)
    {
        write_frame(9999,get_outputdir());
        if(!okay)
            errmsg_printf("at time %e q(%d,%d,%d,%d) is not finite", get_time(),
                    bad_i,bad_j,bad_m,bad_k);
        if(!aux_okay)
            errmsg_printf("at time %e aux(%d,%d,%d,%d) is not finite", get_time(),
                    aux_i,aux_j,aux_m,aux_k);
    }
    return true;
}

#include "DogSolverCart2.h"
//
// user "callbacks"
// (AppStateCart2 should override these deprecated implementations)
//
void DogStateCart2::SetBndValues()
{

// Hackish way to get time passed into the boundary routine.
// What's the correct way to do this?
dogParams.set_time( get_time() );

    // map onto deprecated global method
    void SetBndValues(dTensorBC4&, dTensorBC4&);
    SetBndValues(fetch_q(),fetch_aux());
}
void DogStateCart2::BeforeStage(double dt)
{
    // map onto deprecated global method
    //void ::BeforeStep(double, dTensorBC4&, dTensorBC4&, DogSolverCart2&);
    ::BeforeStep(dt,fetch_aux(),fetch_q(),fetch_solver());
}

// User is responsible to call SetBndValues before calling this method
//
void DogStateCart2::ConstructL()
{
    dTensorBC3& smax = fetch_solver().fetch_smax();
    dTensorBC4& L    = fetch_solver().fetch_L();

    // TODO: I have no idea where this thing is declared ... -DS
    // Why can't I change its definition to accept non const arguments for aux
    // and q?
//  ::ConstructL(get_aux(), get_q(), L, smax);
    ::ConstructL(fetch_aux(), fetch_q(), L, smax);

}

void DogStateCart2::ApplyLimiter()
{
    // map onto deprecated global method
    void ProjectRightEig(int,const dTensor1&, const dTensor1&,const dTensor2&,
            dTensor2&);
    void ProjectLeftEig(int,const dTensor1&, const dTensor1&,const dTensor2&,
            dTensor2&);
    ::ApplyLimiter(fetch_aux(),fetch_q(),ProjectRightEig,ProjectLeftEig);
}

void DogStateCart2::AfterStage(double dt)
{
    // map onto deprecated global method
    ::AfterStep(dt,fetch_aux(),fetch_q(),fetch_solver());
}

void DogStateCart2::AfterReject(double dt)
{
    // map onto deprecated global method
    ::AfterReject(dt,fetch_aux(),fetch_q(),fetch_solver());
}

void DogStateCart2::AfterFullTimeStep(double dt)
{
    // map onto deprecated global method
    ::AfterFullTimeStep(fetch_solver());
}

void DogStateCart2::ReportAfterStep() const
{
    // map onto deprecated global method
    ::ConSoln(get_aux(),get_q(),get_time());
}

void DogStateCart2::AfterInitState()
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
void DogStateCart2::InitState()
{
    // Fetch gives write access
    dTensorBC4& qnew = fetch_q();
    dTensorBC4& aux  = fetch_aux();  

    // Cartesian parameters
    const int mx  = dogParamsCart2.get_mx();
    const int my  = dogParamsCart2.get_my();
    const int mbc = dogParamsCart2.get_mbc();

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
        L2Project(1-mbc,mx+mbc,1-mbc,my+mbc,
                ic_quad_order,space_order,
                space_order,space_order,
                &qnew,&aux,&aux,&AuxFunc);  
    }    


    // Set initial data on computational grid using L2-projection
    void QinitFunc(const dTensor2& xpts,
            const dTensor2& NOT_USED_1,
            const dTensor2& NOT_USED_2,
            dTensor2& qvals);
    L2Project(1-mbc,mx+mbc,1-mbc,my+mbc,
            ic_quad_order,space_order,
            space_order,space_order,
            &qnew,&aux,&qnew,&QinitFunc);
}
