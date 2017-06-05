#ifndef _L2PROJECTINLINE3D_H_
#define _L2PROJECTINLINE3D_H_

#include "tensors.h"

// declare inner methods in-line to avoid function call overhead
//
namespace L2ProjectInline3d
{

    inline void set_vals_at_each_Gaussian_quadrature_point(
        const int i, const int j, const int k,
        const int mpoints,
        const int meqn,
        const int maux,
        const int kmax_qin,
        const int kmax_auxin,
        const double xc, const double yc, const double zc,
        const double dx, const double dy, const double dz,
        const dTensor2& spts,
        const dTensor2& phi,
        const dTensorBC5* qin,
        const dTensorBC5* auxin,
        dTensor2& xpts, dTensor2& qvals, dTensor2& auxvals)
    {
        for (int m=1; m<=(mpoints); m++)
        {
            // grid point (xi,eta) in the unit square
            const double& xi   = spts.get(m,1);
            const double& eta  = spts.get(m,2);
            const double& zeta = spts.get(m,3);

            // location of each gaussian quadrature point in the square with
            // width dx and height dy. this will be passed to SourceTermFunc
            xpts.set( m,1, xc + 0.5*dx*xi   );
            xpts.set( m,2, yc + 0.5*dy*eta  );
            xpts.set( m,3, zc + 0.5*dz*zeta );

            // Solution values (q) at each grid point
            for (int me=1; me<=meqn; me++)
            {
                double val=qin->get(i,j,k,me,1);
                for (int n=2; n<=kmax_qin; n++)
                {
                    val += phi.get(m,n) * qin->get(i,j,k,me,n);
                }
                qvals.set(m,me, val );
            }

            // Auxiliary values (aux) at each grid point
            for (int ma=1; ma<=maux; ma++)
            {
                double val=auxin->get(i,j,k,ma,1);
                for (int n=2; n<=kmax_auxin; n++)
                {
                    val += phi.get(m,n) * auxin->get(i,j,k,ma,n);
                }
                auxvals.set(m,ma, val);
            }
        }
    }

    // This method performs a single projection of a function onto the basis
    // functions.
    //
    // Parameters:
    // ----------
    //
    //    bool add : flag to determine if we want to add the solution to
    //    whatever currently exists in the large tensor, fout.  If true, this
    //    adds the solution to fout, if false, it overwrites whatever is
    //    stored in fout.
    //
    //    (i,j,k)    : index of cell to be integrated.
    //
    //    mcomps_out : number of components to be written to.  This block
    //    writes to fout( i, j, k, 1:mcomps_out, 1:kmax_fout ).
    //
    //    kmax_fout  : number of polynomials for fout
    //
    //    mpoints    : number of quadrature points used for the integration.
    //
    //    wgt_phi_transpose  : product of quadrature weights and basis
    //                         functions evaluate at each point.  The purpose
    //                         of the `transpose' is to speed up the code.
    //                         wgt_phi_transpose( 1:kmax_fout, 1:mpoints ),
    //                         having the most aggressive index on the outer
    //                         loop is the fastest.
    //
    //    fvals      : pointwise values of the function being integrated.
    //                 These are of the form, fvals( 1:mpoints, 1:? )
    //                      
    // Returns:
    // --------
    //
    //    fout       : tensor being written to.  If add = true, it adds the
    //                 result to whatever is currently there, otherwise it
    //                 overwrites the result.
    //
    inline void integrate_on_current_cell(
        const bool add,
        const int i, const int j, const int k,
        const int mcomps_out,
        const int kmax_fout, 
        const int mpoints,
        const dTensor2& wgt_phi_transpose,
        const dTensor2& fvals,
        dTensorBC5* fout)
    {
        for (int meq=1; meq<=mcomps_out; meq++)
        for (int n=1; n<=kmax_fout; n++)
        {
            double tmp = 0.;
            int fvals_idx = fvals.getidx(1,meq);
            for (int mp=1; mp<=mpoints; mp++)
            {
                // -------------------------------------------------- //
                // This is the most executed block in the whole code.
                // See notes in lib/2d/cart/L2ProjectInline.h for why 
                // we use phi_transpose
                // -------------------------------------------------- //
                tmp += wgt_phi_transpose.get(n,mp)*fvals.vget(fvals_idx);
                fvals_idx += mcomps_out;
            }
            if(add)
            {  fout->fetch(i,j,k,meq,n) += tmp*0.125;  }
            else
            {  fout->set(i,j,k,meq,n, tmp*0.125);  }
        } 
    }

    // Version for projecting onto the gradient of the basis functions.
    inline void integrate_on_current_cell_grad(
        const int i, const int j, const int k,
        const int mcomps_out,
        const int kmax_fout, 
        const int mpoints,
        const iTensor1& nonzero_flag,
        const dTensor2& wgt_phi_x_tr, const dTensor2& wgt_phi_y_tr, const dTensor2& wgt_phi_z_tr,
        const dTensor3& fvals, dTensorBC5* fout)
    {
        for (int mi=1; mi<=mcomps_out; mi++)
        {
            // This is the second-most executed block in the code after
            // the loop in L2Project.
            // can start with 2 since phi_*(:, 1)=0
            for (int k2=2; k2<=kmax_fout; k2++)
            {
                double tmp = 0.;
                // switch statement to avoid unnecessary computations
                switch( nonzero_flag.get(k2) )
                {
                    case 1: // 000
                        eprintf("this case should only arise when k2=%d equals 1", k2);

                    case 2: // 100
                        for (int mp=1; mp<=mpoints; mp++)
                        {
                            tmp += fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp);
                        }
                        break;

                    case 3: // 010
                        for (int mp=1; mp<=mpoints; mp++)
                        {
                            tmp += fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp);
                        }
                        break;

                    case 4: // 001
                        for (int mp=1; mp<=mpoints; mp++)
                        {
                            tmp += fvals.get(mp, mi, 3)*wgt_phi_z_tr.get(k2, mp);
                        }
                        break;

                    case 5: // 110
                        for (int mp=1; mp<=mpoints; mp++)
                        {
                            tmp += ( fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp) +
                                    fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp) );
                        }
                        break;

                    case 6: // 101
                        for (int mp=1; mp<=mpoints; mp++)
                        {
                            tmp += ( fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp) +
                                    fvals.get(mp, mi, 3)*wgt_phi_z_tr.get(k2, mp) );
                        }
                        break;

                    case 7: // 011
                        for (int mp=1; mp<=mpoints; mp++)
                        {
                            tmp += ( fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp) +
                                    fvals.get(mp, mi, 3)*wgt_phi_z_tr.get(k2, mp) );
                        }
                        break;

                    case 8: // 111
                        for (int mp=1; mp<=mpoints; mp++)
                        {
                            tmp += ( fvals.get(mp, mi, 1)*wgt_phi_x_tr.get(k2, mp) +
                                    fvals.get(mp, mi, 2)*wgt_phi_y_tr.get(k2, mp) +
                                    fvals.get(mp, mi, 3)*wgt_phi_z_tr.get(k2, mp) );
                        }
                        break;
                    default:
                        // should never get here
                        eprintf("error in L2ProjectInline3d::integrate_on_current_cell_grad, nonzero_flag.get(%i) = %i\n",
                                k2,nonzero_flag.get(k2));
                }
                fout->fetch(i, j, k, mi, k2) += tmp*0.125;
            }
        }
    }

}

#endif
