#ifndef _L2PROJECTINLINE4D_H_
#define _L2PROJECTINLINE4D_H_

#include "tensors.h"

// declare inner methods in-line to avoid function call overhead
//
namespace L2ProjectInline4d
{

    // This routine samples the basis functions at a list of points provided by
    // spts.
    //
    // Input:
    // ------
    //
    //      i,j,k,l - index into the grid describing what element to sample
    //      mpoints - number of points used for sampling (should be const. with
    //                spts)
    //      [others] 
    //
    // Returns:
    // --------
    //
    //     xpts( 1:mpts, 1:4 ) - in 4D, we have x,y,z,w for the second
    //                           coefficient
    //
    //    qvals( 1:mpts, 1:meqn ) - the solution at each point.
    //  auxvals( 1:mpts, 1:maux ) - the auxiliary function at each point.
    //
    inline void set_vals_at_each_Gaussian_quadrature_point(
        const int i, const int j, const int k, const int l,
        const int mpoints,
        const int meqn,
        const int maux,
        const int kmax_qin,
        const int kmax_auxin,
        const double xc, const double yc, const double zc, const double wc,
        const double dx, const double dy, const double dz, const double dw,
        const dTensor2& spts,
        const dTensor2& phi,
        const dTensorBC6* qin,
        const dTensorBC6* auxin,
        dTensor2& xpts,
        dTensor2& qvals,
        dTensor2& auxvals)
    {
        for (int m=1; m<=(mpoints); m++)
        {
            // grid point (xi,eta,zeta,tau) in the unit square
            const double& xi   = spts.get(m,1);
            const double& eta  = spts.get(m,2);
            const double& zeta = spts.get(m,3);
            const double&  tau = spts.get(m,4);

            // location of each gaussian quadrature point in the square with
            // width dx and height dy. this will be passed to SourceTermFunc
            xpts.set( m,1, xc + 0.5*dx*xi   );
            xpts.set( m,2, yc + 0.5*dy*eta  );
            xpts.set( m,3, zc + 0.5*dz*zeta );
            xpts.set( m,4, zc + 0.5*dw*tau  );

            // Solution values (q) at each grid point
            for (int me=1; me<=meqn; me++)
            {
                double val = qin->get(i,j,k,l,me,1);
                for (int n=2; n<=kmax_qin; n++)
                {
                    val += phi.get(m,n) * qin->get(i,j,k,l,me,n);
                }
                qvals.set(m,me, val );
            }

            // Auxiliary values (aux) at each grid point
            for (int ma=1; ma<=maux; ma++)
            {
                double val=auxin->get(i,j,k,l,ma,1);
                for (int n=2; n<=kmax_auxin; n++)
                {
                    val += phi.get(m,n) * auxin->get(i,j,k,l,ma,n);
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
    inline void integrate_on_current_cell(const bool add,
            const int i, 
            const int j,
            const int k,
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
                    // This is the most executed block in the whole code.
                    //
                    // We reduce the number of multiplications executed here
                    // in each pass from four to one by defining wgt_phi_transpose
                    // and accessing the fvals index directly using vget.
                    // Using wgt_phi_transpose seems to shave about 6% off the execution time.
                    // Using vget shaves at most 1% (if anything). Times with -O4:
                    // nothing: 10.64012 9.82373 9.76175 9.79337 9.78146
                    // wgt_phi_transpose: : 9.27372 9.37958 9.25974 9.25850 9.24290
                    // and vget: 9.19254 9.26525 9.19923 9.15430 9.26599 9.28980 9.16268 9.16859
                    //
                    //tmp += wgt.get(mp)*phi.get(mp,k)*fvals.get(mp,meq);
                    //tmp += wgt_phi_transpose.get(k,mp)*fvals.get(mp,meq);
                    tmp += wgt_phi_transpose.get(n,mp)*fvals.vget(fvals_idx);
                    fvals_idx += mcomps_out;
                }
                if(add)
                {  fout->fetch(i,j,k,meq,n) += tmp*0.125;  }
                else
                {  fout->set(i,j,k,meq,n, tmp*0.125);  }
            } 
    }

    inline void integrate_on_current_cell_grad(const int i, 
            const int j,
            const int k,
            const int mcomps_out,
            const int kmax_fout, 
            const int mpoints,
            const iTensor1& nonzero_flag,
            const dTensor2& wgt_phi_x_tr,
            const dTensor2& wgt_phi_y_tr,
            const dTensor2& wgt_phi_z_tr,
            const dTensor3& fvals,
            dTensorBC5* fout)
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
                        eprintf("error in L2ProjectInline4d::integrate_on_current_cell_grad, nonzero_flag.get(%i) = %i\n",
                                k2,nonzero_flag.get(k2));
                }
                fout->fetch(i, j, k, mi, k2) += tmp*0.125;
            }
        }
    }

}

#endif
