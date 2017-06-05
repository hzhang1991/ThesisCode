#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(
    const double t, 
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& source)
{

    const int numpts=xpts.getsize(1);
    const int meqn=source.getsize(2);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double r  = sqrt( pow(x,2) + pow(y,2) );

        double dr_dx = x / (r+1e-14);
        double dr_dy = y / (r+1e-14);

        // ------------------------------------------------------------- //
        // Exact solution is given by:
        // qex = (0.75 + 0.25*cos(2*pi*t) ) * ( 0.5 + 0.5*cos( pi*r ) )
        // ------------------------------------------------------------- //
        double qt = -0.5*pi*sin(2.0*pi*t)   * ( 0.5 + 0.5*cos( pi*r ) );
        double qx = (0.75+0.25*cos(2*pi*t)) * ( -0.5*pi*sin( pi*r )* dr_dx );
        double qy = (0.75+0.25*cos(2*pi*t)) * ( -0.5*pi*sin( pi*r )* dr_dy );

        for (int m=1; m<=meqn; m++)
        {
            double vx = vel_vec->get(m,1);
            double vy = vel_vec->get(m,2);

            double tmp = qt + vx*qx + vy*qy;
            source.set(i,m, tmp );

        }

    }

source.setall(0.);

}
