#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart4.h"
#include "L2ProjectInline4d.h"
#include "Legendre4d.h"

// -------------------------------------------------------------
// Routine for computing the L2-projection of an input function
// onto an orthonormal Legendre basis
// -------------------------------------------------------------
static void L2ProjectAdd(const bool add,
             const int istart,
             const int iend,
             const int jstart,
             const int jend,
             const int kstart,
             const int kend,
             const int lstart,
             const int lend,
             const int QuadOrder,
             const int BasisOrder_qin,
             const int BasisOrder_auxin,
             const int BasisOrder_fout,    
             const dTensorBC6* qin,
             const dTensorBC6* auxin,
             dTensorBC6* fout,
             void (*Func)(),
             bool has_data, 
             void* data)
{

/*
  if(!has_data) assert(!data);
  
  // dx, dy, dz
  const double dx = dogParamsCart4.get_dx();
  const double dy = dogParamsCart4.get_dy();
  const double dz = dogParamsCart4.get_dz();
  
  // mbc
  const int mbc = qin->getmbc();
  assert_eq(mbc,auxin->getmbc());
  assert_eq(mbc,fout->getmbc());

  // qin variable
  const int       mx = qin->getsize(1);
  const int       my = qin->getsize(2);
  const int       mz = qin->getsize(3);
  const int     meqn = qin->getsize(4);
  const int kmax_qin = qin->getsize(5);
  int ktmp = (BasisOrder_qin*(BasisOrder_qin+2)*(BasisOrder_qin+1))/6;
  assert_eq(kmax_qin,ktmp);
  
  // auxin variable
  assert_eq(mx,auxin->getsize(1));
  assert_eq(my,auxin->getsize(2));
  assert_eq(mz,auxin->getsize(3));
  const int       maux = auxin->getsize(4);
  const int kmax_auxin = auxin->getsize(5);
  ktmp = (BasisOrder_auxin*(BasisOrder_auxin+2)*(BasisOrder_auxin+1))/6;
  assert_eq(kmax_auxin,ktmp);

  // fout variables
  //  TODO - why assume this has the same size?  What if we want to only
  //  project onto a single line, e.g. for padding boundary cell data? (-DS)
  assert_eq(mx,fout->getsize(1));  
  assert_eq(my,fout->getsize(2));
  assert_eq(mz,fout->getsize(3));
  const int mcomps_out = fout->getsize(4);
  const int  kmax_fout = fout->getsize(5);
  assert_eq(kmax_fout,(BasisOrder_fout*(BasisOrder_fout+2)*(BasisOrder_fout+1))/6);
  
  // starting and ending indeces
  assert_ge(istart,1-mbc);
  assert_le(iend,mx+mbc);
  assert_ge(jstart,1-mbc);
  assert_le(jend,my+mbc);
  assert_ge(kstart,1-mbc);
  assert_le(kend,mz+mbc);

  // number of quadrature points
  assert_ge(QuadOrder,1);
  assert_le(QuadOrder,20);
  int QuadOrder_MOD = QuadOrder;
  switch(QuadOrder)
    {
    case 7:
      QuadOrder_MOD = 8;
      break;
    case 9:
      QuadOrder_MOD = 10;
      break;
    case 11:
      QuadOrder_MOD = 12;
      break;
    case 13:
      QuadOrder_MOD = 14;
      break;
    case 15:
      QuadOrder_MOD = 16;
      break;
    case 17:
      QuadOrder_MOD = 18;
      break;
    case 19:
      QuadOrder_MOD = 20;
      break;
    }
  const int mpoints = QuadOrder_MOD*QuadOrder_MOD*QuadOrder_MOD;
  
  // set quadrature point and weight information
  void SetQuadWgtsPts(const int, dTensor1&, dTensor2&);
  dTensor1 wgt(mpoints);
  dTensor2 spts(mpoints, 3);
  SetQuadWgtsPts(QuadOrder_MOD, wgt, spts);
  
  // Loop over each quadrature point to construct Legendre polys
  const int kmax = iMax(iMax(kmax_qin,kmax_auxin),kmax_fout);
  dTensor2 phi(mpoints, kmax);
  void SetLegendrePolys(const int, const int, const dTensor2&, dTensor2&);
  SetLegendrePolys(mpoints, kmax, spts, phi);
  
  // For efficiency compute weight*phi and then take transpose
  dTensor2 wgt_phi_transpose(kmax,mpoints);
  for(int mp=1;mp<=mpoints;mp++)
    for(int k=1;k<=kmax;k++)
      {
        wgt_phi_transpose.set(k,mp, wgt.get(mp)*phi.get(mp,k) );
      }

  // ------------------------------------------------------------- //
  // Loop over every grid cell indexed by user supplied parameters //
  // described by istart...iend, jstart...jend                     // 
  // ------------------------------------------------------------- //
#pragma omp parallel for
  for (int i=istart; i<=iend; i++)
    {
      
      // Local storage for q, aux and xpts (all passed into user supplied
      // function)
      dTensor2   qvals(mpoints, meqn);
      dTensor2 auxvals(mpoints, maux);
      dTensor2    xpts(mpoints, 3);

      // Flux function and its Jacobian:
      dTensor2   fvals(mpoints, mcomps_out);

      //find center of current cell
      const double xc = dogParamsCart4.get_xc(i);
      
      for (int j=jstart; j<=jend; j++)
    {
      //find center of current cell
      const double yc = dogParamsCart4.get_yc(j);

      for (int k=kstart; k<=kend; k++)
        {
          //find center of current cell
          const double zc = dogParamsCart4.get_zc(k);
      
          // Compute q, aux and fvals at each Gaussian quadrature point
          // for this current cell indexed by (i,j,k)
          // Save results into dTensor2 qvals, auxvals and fvals.
          L2ProjectInline4d::set_vals_at_each_Gaussian_quadrature_point(i, j, k,
                                        mpoints, meqn, maux, 
                                        kmax_qin, kmax_auxin, 
                                        xc, yc, zc, dx, dy, dz, 
                                        spts, phi, qin, auxin, 
                                        xpts, qvals, auxvals);
          // Evaluate Func at Gaussian quadrature point
          if(has_data)
        {
          ((void (*)(const dTensor2&, const dTensor2&, 
                 const dTensor2&, dTensor2&,  void* data))
           Func)(xpts, qvals, auxvals, fvals, data);
        }
          else
        {
          ((void (*)(const dTensor2&, const dTensor2&, 
                 const dTensor2&, dTensor2&))
           Func)(xpts, qvals, auxvals, fvals);
        }
          
          // Evaluate integral on current cell (project onto Legendre basis)
          // using Gaussian Quadrature for the integration
          L2ProjectInline4d::integrate_on_current_cell(add, i, j, k,
                               mcomps_out, kmax_fout, mpoints, 
                               wgt_phi_transpose, 
                               fvals, fout);
        }
    }
    }
*/

}

// L2Project
void L2Project(
    const int istart, const int iend,
    const int jstart, const int jend,
    const int kstart, const int kend,
    const int lstart, const int lend,
    const int QuadOrder,
    const int BasisOrder_qin,
    const int BasisOrder_auxin,
    const int BasisOrder_fout,    
    const dTensorBC6* qin,
    const dTensorBC6* auxin,
    dTensorBC6* fout,
    void (*Func)(const dTensor2&,const dTensor2&,
        const dTensor2&,dTensor2&))
{
    L2ProjectAdd(false,istart,iend,jstart,jend,kstart,kend,lstart,lend,
            QuadOrder,BasisOrder_qin,BasisOrder_auxin,
            BasisOrder_fout,qin,auxin,fout,
            ((void(*)()) Func),false,0);
}

void L2Project(
    const int istart, const int iend, 
    const int jstart, const int jend,
    const int kstart, const int kend,
    const int lstart, const int lend,
    const dTensorBC6& q,
    const dTensorBC6& aux, 
    dTensorBC6& Fout,
    void (*Func)(const dTensor2& xpts,
        const dTensor2& qvals,
        const dTensor2& auxvals,
        dTensor2& source))
{
    const int space_order = dogParams.get_space_order();
    L2ProjectAdd(false,
         istart, iend,
         jstart, jend,
         kstart, kend,
         lstart, lend,
         space_order,
         space_order,
         space_order,
         space_order,
         &q,
         &aux,
         &Fout,
         ((void(*)()) Func),false,0);
}

// L2Project version to support a callback function that
// receives extra data registered by the user.
//
// This mechanism could be used to pass DogSolver::solver if
// we someday decide that it should no longer be a singleton;
// for now it is accessible via DogSolver::get_solver();
//
// This mechanism can be used to pass time increment
// information. Note that the time increment is available
// via DogSolver::get_solver().get_dt(), and the time of the
// current stage is set in DogSolver::advanceTimeStageRK() and
// is available via a call to DogSolver::get_time_hack().
//
void L2Project_extra(
    const int istart, const int iend,
    const int jstart, const int jend,
    const int kstart, const int kend,
    const int lstart, const int lend,
    const int QuadOrder,
    const int BasisOrder_qin,
    const int BasisOrder_auxin,
    const int BasisOrder_fout,    
    const dTensorBC6* qin,
    const dTensorBC6* auxin,
    dTensorBC6* fout,
    void (*Func)(const dTensor2&,const dTensor2&,
        const dTensor2&,dTensor2&, void* data),
    void* data)
{
    L2ProjectAdd(false,
            istart, iend,
            jstart, jend,
            kstart, kend,
            lstart, lend,
            QuadOrder,
            BasisOrder_qin,
            BasisOrder_auxin,
            BasisOrder_fout,
            qin,
            auxin,
            fout,
            ((void(*)()) Func),true,data);
}

void L2Project_extra(
    const int istart, const int iend, 
    const int jstart, const int jend, 
    const int kstart, const int kend, 
    const int lstart, const int lend,
    const dTensorBC6& q,
    const dTensorBC6& aux, 
    dTensorBC6& Fout,
    void (*Func)(const dTensor2&,const dTensor2&,
         const dTensor2&,dTensor2&, void* data),
    void* data)
{

/*
  const int space_order = dogParams.get_space_order();
  L2ProjectAdd(false,
           istart,
           iend,
           jstart,
           jend,
           kstart,
           kend,
           space_order,
           space_order,
           space_order,
           space_order,
           &q,
           &aux,
           &Fout,
           ((void(*)()) Func),true,data);
*/

}
