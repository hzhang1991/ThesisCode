#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart3.h"
#include "Legendre3d.h"
#include "L2ProjectInline3d.h"

// -------------------------------------------------------------------------- //
// All-purpose routine for computing the L2-projection
// of various functions onto Gradient (e.g divergence) of the Legendre basis.
//
// The routine L2ProjectGradAdd assumes that there is information stored in 
// fout that the user doesn't want overwritten.  In this case, start by 
// projecting Func onto the basis functions, and then "add" the moments onto
// fout.
//
// This is the case in "Part III: compute intra-element contributions" of
// ConstructL, where we already have inter-element constributions saved.
//
// The form of Func is inconsistent with what is currently in L2Project.  The
// reason being that this routine was only ever used to project FluxFunc onto
// the basis functions.  The idea is the following:
//
// A gradient of phi produces three components: grad( phi ) = < phi_x, phi_y, phi_z >.
//
// Therefore, to "project" this onto the basis functions, we need to have
// access to something with three components: f = < f1, f2, f3 >.
//
// The fourth argument to Func sets fvals( 1:numpts, 1:MEQN, 1:NDIMS ), which
// supplies us with the necessary arguments for f.
//
// At the end of the day, we end up computing the following quantity:
//
//    fout( istart:iend, jstart:jend, kstart:kend, 1:length(f), k ) = 
//
//              1 / dV \iiint grad( phi^{(k)} ) \cdot f dx dy dz
//
// -------------------------------------------------------------------------- //
void L2ProjectGradAdd(const int istart, 
        const int iend, 
        const int jstart, 
        const int jend,
        const int kstart, 
        const int kend,
        const int QuadOrder, 
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,
        const dTensorBC5* qin,
        const dTensorBC5* auxin, 
        dTensorBC5* fout,
        void Func(const dTensor2&, const dTensor2&,
            const dTensor2&, dTensor3&))
{
    // dx, dy, and dz
    const double dx   = dogParamsCart3.get_dx();
    const double dy   = dogParamsCart3.get_dy();
    const double dz   = dogParamsCart3.get_dz();
    const double xlow = dogParamsCart3.get_xlow();
    const double ylow = dogParamsCart3.get_ylow();
    const double zlow = dogParamsCart3.get_zlow();

    // mbc
    const int mbc = qin->getmbc();
    assert_eq(mbc, auxin->getmbc());
    assert_eq(mbc, fout->getmbc());

    // qin variable
    const int       mx = qin->getsize(1);
    const int       my = qin->getsize(2);
    const int       mz = qin->getsize(3);
    const int     meqn = qin->getsize(4);
    const int kmax_qin = qin->getsize(5);
    assert_eq(kmax_qin, (BasisOrder_qin*(BasisOrder_qin+1)*(BasisOrder_qin+2))/6);

    // auxin variable
    assert_eq(mx, auxin->getsize(1));
    assert_eq(my, auxin->getsize(2));
    assert_eq(mz, auxin->getsize(3));
    const int       maux = auxin->getsize(4);
    const int kmax_auxin = auxin->getsize(5);
    assert_eq(kmax_auxin, (BasisOrder_auxin*(BasisOrder_auxin+1)*(BasisOrder_auxin+2))/6);

    // fout variables
    assert_eq(mx, fout->getsize(1));
    assert_eq(my, fout->getsize(2));
    assert_eq(mz, fout->getsize(3));
    const int mcomps_out = fout->getsize(4);
    const int  kmax_fout = fout->getsize(5);
    assert_eq(kmax_fout, (BasisOrder_fout*(BasisOrder_fout+1)*(BasisOrder_fout+2))/6);

    // starting and ending indeces
    assert_ge(istart, 1-mbc);
    assert_le(iend, mx+mbc);
    assert_ge(jstart, 1-mbc);
    assert_le(jend, my+mbc);
    assert_ge(kstart, 1-mbc);
    assert_le(kend, mz+mbc);

    // number of quadrature points
    assert_ge(QuadOrder, 1);
    assert_le(QuadOrder, 6);
    const int mpoints = iMax((QuadOrder-1)*(QuadOrder-1)*(QuadOrder-1), 1);

    // trivial case
    if ( BasisOrder_fout==1 )
    {
        fout->setall(0.);
        return;
    }

    // set quadrature point and weight information
    void SetQuadWgtsPts(const int, dTensor1&, dTensor2&);
    dTensor1 wgt(mpoints);
    dTensor2 spts(mpoints, 3);
    SetQuadWgtsPts(QuadOrder-1, wgt, spts);

    // Loop over each quadrature point to construct Legendre polys
    const int kmax = iMax(iMax(kmax_qin, kmax_auxin), kmax_fout);
    dTensor2   phi(mpoints, kmax);
    dTensor2 phi_x(mpoints, kmax_fout);
    dTensor2 phi_y(mpoints, kmax_fout);
    dTensor2 phi_z(mpoints, kmax_fout);

    void SetLegendrePolys(const int, 
            const int, 
            const dTensor2&, 
            dTensor2&);
    SetLegendrePolys(mpoints, kmax, spts, phi);

    void SetLegendrePolysGrad(
            const double dx, const double dy, const double dz,
            const int mpoints, const int kmax_fout,
            const dTensor2& spts, dTensor2& phi_x, dTensor2& phi_y, dTensor2& phi_z);
    SetLegendrePolysGrad(dx, dy, dz, mpoints, kmax_fout, spts, phi_x, phi_y, phi_z);

    dTensor2 wgt_phi_x_tr(kmax_fout, mpoints);
    dTensor2 wgt_phi_y_tr(kmax_fout, mpoints);
    dTensor2 wgt_phi_z_tr(kmax_fout, mpoints);
    iTensor1 nonzero_flag(kmax_fout);  
    for(int k=1; k<=kmax_fout; k++)
    {
        iTensor1 nonzero_tmp(3);
        nonzero_tmp.set( 1, 0 );
        nonzero_tmp.set( 2, 0 );
        nonzero_tmp.set( 3, 0 );
        nonzero_flag.set(k, 0 );

        for(int mp=1; mp<=mpoints; mp++)
        {
            double xval = wgt.get(mp)*phi_x.get(mp, k);
            double yval = wgt.get(mp)*phi_y.get(mp, k);
            double zval = wgt.get(mp)*phi_z.get(mp, k);
            wgt_phi_x_tr.set(k, mp, xval );
            wgt_phi_y_tr.set(k, mp, yval );
            wgt_phi_z_tr.set(k, mp, zval );
            if(xval) nonzero_tmp.set(1, 1 );
            if(yval) nonzero_tmp.set(2, 1 );
            if(zval) nonzero_tmp.set(3, 1 );      
        }

        switch(100*nonzero_tmp.get(1) + 10*nonzero_tmp.get(2) + nonzero_tmp.get(3))
        {
            case 0:   // 000
                nonzero_flag.set(k, 1 );
                break;
            case 100: // 100
                nonzero_flag.set(k, 2 );
                break;
            case 10:  // 010
                nonzero_flag.set(k, 3 );
                break;
            case 1:   // 001
                nonzero_flag.set(k, 4 );
                break;
            case 110: // 110
                nonzero_flag.set(k, 5 );
                break;
            case 101: // 101
                nonzero_flag.set(k, 6 );
                break;
            case 11:  // 011
                nonzero_flag.set(k, 7 );
                break;
            case 111: // 111
                nonzero_flag.set(k, 8 );
                break;      
        }
    }

    // ----------------------------------
    // Loop over all elements of interest
    // ----------------------------------  
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {
        // Local storage for q, aux and xpts (all passed into user supplied
        // function)
        dTensor2 xpts(mpoints, 3);
        dTensor2 qvals(mpoints, meqn);
        dTensor2 auxvals(mpoints, maux);

        // Point values of output function
        dTensor3 fvals(mpoints, mcomps_out, 3);

        //find center of current cell
        const double xc = xlow + (double(i)-0.5)*dx;

        for (int j=jstart; j<=jend; j++)
        {
            //find center of current cell
            const double yc = ylow + (double(j)-0.5)*dy;

            for (int k=kstart; k<=kend; k++)
            {
                //find center of current cell
                const double zc = zlow + (double(k)-0.5)*dz;

                // Compute q, aux and fvals at each Gaussian quadrature point
                // for this current cell indexed by (i,j,k)
                // Save results into dTensor2 qvals, auxvals and fvals.
                L2ProjectInline3d::set_vals_at_each_Gaussian_quadrature_point(i, j, k,
                        mpoints, 
                        meqn, 
                        maux, 
                        kmax_qin, 
                        kmax_auxin, 
                        xc, yc, zc, dx, dy, dz, 
                        spts, phi, qin, auxin, 
                        xpts, qvals, auxvals);

                // Call user-supplied function to set fvals
                Func(xpts, qvals, auxvals, fvals);

                // Evaluate integral on current cell (project onto gradient Legendre basis)
                // using Gaussian Quadrature for the integration          
                L2ProjectInline3d::integrate_on_current_cell_grad(i, j, k,
                        mcomps_out, 
                        kmax_fout, 
                        mpoints, 
                        nonzero_flag,
                        wgt_phi_x_tr,
                        wgt_phi_y_tr,
                        wgt_phi_z_tr,
                        fvals,
                        fout);

            }
        }
    }


}// end of function L2ProjectGrad

void L2ProjectGradAdd(int istart, 
        int iend, 
        int jstart, 
        int jend,
        int kstart,
        int kend,
        const dTensorBC5& qin,
        const dTensorBC5& auxin, 
        dTensorBC5& Fout,
        void Func(const dTensor2&, const dTensor2&,
            const dTensor2&, dTensor3&))
{
    const int space_order = dogParams.get_space_order();
    L2ProjectGradAdd(istart, 
            iend, 
            jstart, 
            jend, 
            kstart, 
            kend,
            space_order, 
            space_order, 
            space_order, 
            space_order,
            &qin, 
            &auxin, 
            &Fout, 
            Func);
}
