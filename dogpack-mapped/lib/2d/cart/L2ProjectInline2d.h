#ifndef _L2PROJECTINLINE2D_H_
#define _L2PROJECTINLINE2D_H_

    void mapc2p(double& xc,double& yc);
    double mapx(double xi,double eta,double xp1,double xp2,double xp3,double xp4,double yp1,double yp2,double yp3,double yp4);
    double mapy(double xi,double eta,double xp1,double xp2,double xp3,double xp4,double yp1,double yp2,double yp3,double yp4);
#include "tensors.h"


// This module describes two of the most commonly executed blocks in the code:
//
//      (*) sampling the basis functions on a single element
//      (*) integrating the basis functions on a single element
//
// We declare these inner methods in-line to avoid function call overhead.
namespace L2ProjectInline2d
{

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
    //    (i,j)    : index of cell to be integrated.
    //
    //    mcomps_out : number of components to be written to.  This block
    //    writes to fout( i, j, 1:mcomps_out, 1:kmax_fout ).
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
            const int i, const int j,
            const int mcomps_out,
            const int kmax_fout, 
            const int mpoints,
            const dTensor2& wgt_phi_transpose,
            const dTensor2& fvals,
            dTensorBC4* fout)
    {
        for (int meq=1; meq<=mcomps_out; meq++)
        for (int k=1; k<=kmax_fout; k++)
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
                //
                //   nothing:           10.64012 9.82373 9.76175 9.79337 9.78146
                //   wgt_phi_transpose: 9.27372  9.37958 9.25974 9.25850 9.24290
                //   vget:              9.19254  9.26525 9.19923 9.15430 9.26599 9.28980 9.16268 9.16859
                //

                // "nothing":
                //tmp += wgt.get(mp)*phi.get(mp,k)*fvals.get(mp,meq);

                // wgt_phi_transpose:
                //tmp += wgt_phi_transpose.get(k,mp)*fvals.get(mp,meq);

                // vget with wgt_phi_transpose (fastest execution time):
                tmp += wgt_phi_transpose.get(k,mp)*fvals.vget(fvals_idx);

                fvals_idx += mcomps_out;
            }
            if(add)
            {  fout->fetch(i,j,meq,k) += tmp*0.25;  }
            else
            {  fout->set(i,j,meq,k, tmp*0.25);  }
        } 
    }


    // This routine samples the basis functions at a list of points provided by
    // spts.
    //
    // Input:
    // ------
    //
    //      i,j     - index into the grid describing what element to sample
    //      mpoints - number of points used for sampling (should be consistent 
    //                with spts)
    //      [others] 
    //
    // Returns:
    // --------
    //
    //     xpts( 1:mpts, 1:2 ) - in 2D, x = xpts( :, 1 ) and y = xpts( :, 2 ).
    //    qvals( 1:mpts, 1:meqn ) - the solution at each point.
    //  auxvals( 1:mpts, 1:maux ) - the auxiliary function at each point.
    //
    inline void set_vals_at_each_Gaussian_quadrature_point(
        const int i, const int j,
        const int mpoints,
        const int meqn,
        const int maux,
        const int kmax_qin,
        const int kmax_auxin,
        const double xc, const double yc,
        const double dx, const double dy,
        const dTensor2& spts, const dTensor2& phi,
        const dTensorBC4* qin, const dTensorBC4* auxin,
        dTensor2& xpts, dTensor2& qvals, dTensor2& auxvals)
    {
        for (int m=1; m<=(mpoints); m++)
        {
            // grid point (xi,eta) in the unit square
            const double& xi  = spts.get(m,1);
            const double& eta = spts.get(m,2);

            // location of each gaussian quadrature point in the square with
            // width dx and height dy. this will be passed to SourceTermFunc

            double xp1 = xc-0.5*dx;
            double yp1= yc-0.5*dy;
            double xp2 = xc+0.5*dx;
            double yp2= yc-0.5*dy;
            double xp3 = xc+0.5*dx;
            double yp3= yc+0.5*dy;
            double xp4 = xc-0.5*dx;
            double yp4= yc+0.5*dy;


            mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);

            double x1=mapx(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4);
            double y1=mapy(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4);


            xpts.set( m,1, x1  );
            xpts.set( m,2, y1 );

            // Solution values (q) at each grid point
            for (int me=1; me<=meqn; me++)
            {
                double val=qin->get(i,j,me,1);
                for (int k=2; k<=kmax_qin; k++)
                {
                    val += phi.get(m,k) * qin->get(i,j,me,k);
                }
                qvals.set(m,me, val );
            }

            // Auxiliary values (aux) at each grid point
            for (int ma=1; ma<=maux; ma++)
            {
                double val=auxin->get(i,j,ma,1);
                for (int k=2; k<=kmax_auxin; k++)
                {
                    val += phi.get(m,k) * auxin->get(i,j,ma,k);
                }
                auxvals.set(m,ma, val);
            }
        }
    }

}

#endif
