#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "Legendre2d.h"
#include "ApplyPosLimiter.h"
#include <iostream>
using namespace std;


// Routine that applies a modified Barth-Jespersen limiter to higher order moments of
// the conserved variables.  This limiter is locally applied, and conserves
// total mass by not adjusting cell averages.
//
// See: K. Michalak and C. Ollivier-Gooch, "Limiters for Unstructured Higher-Order
// Accurate Solutions of the Euler Equations"
inline double phi_func(double x)
{
   return (x<1.1)*Min(1.0/1.1*x,1.0)+(x>=1.1)*1.0;//(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}
void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);

void ProjectLeftEig(int ixy,
                    const dTensor1& Aux_ave,
                    const dTensor1& Q_ave,
                    const dTensor2& Qvals,
                    dTensor2& Wvals);

void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
    void (*ProjectRightEig)(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&),
    void (*ProjectLeftEig)(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&))
{

    ApplyPosLimiter(aux, q);

}


// Gaussian quadrature points.  However, an application will likely want to add in the edge
// points that are used for Riemann solves.
void SetPositivePointsGauss(const int& space_order, dTensor2& spts)
{

    dTensor1 s1d(space_order);

    // 1D Gaussian quadrature points
    switch ( space_order )
    {
        case 1:
            s1d.set(1, 0.0e0 );
            break;

        case 2:
            s1d.set(1, -1.0/sq3 );
            s1d.set(2,  1.0/sq3 );      
            break;

        case 3:
            s1d.set(1, -sq3/sq5 );
            s1d.set(2,  0.0e0 );
            s1d.set(3,  sq3/sq5 );
            break;

        case 4:
            s1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            s1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );      
            break;

        case 5:
            s1d.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            s1d.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(3,  0.0 );
            s1d.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );      
            break;
    }

    // This region has been commented out because we are no longer applying
    // the limiter at the corner points

    // 2D points -- corners
  spts.set(1,1, -1.0e0 );
  spts.set(1,2, -1.0e0 );

  spts.set(2,1,  1.0e0 );
  spts.set(2,2, -1.0e0 );

  spts.set(3,1, -1.0e0 );
  spts.set(3,2,  1.0e0 );

  spts.set(4,1,  1.0e0 );
  spts.set(4,2,  1.0e0 );

    // 2D points -- left, right, bottom and top edges
  for (int m=1; m<=space_order; m++)
  {
      double s = s1d.get(m);

      // left edge
      spts.set(4+m,1, -1.0e0 );
      spts.set(4+m,2,  s     );

      // right edge
      spts.set(4+space_order+m,1,  1.0e0 );
      spts.set(4+space_order+m,2,  s     );

      // bottom edge
      spts.set(4+2*space_order+m,1,  s     );
      spts.set(4+2*space_order+m,2, -1.0e0 );

      // top edge
      spts.set(4+3*space_order+m,1,  s     );
      spts.set(4+3*space_order+m,2,  1.0e0 );
  }

    // 2D points -- all interior points
    int z = 4+4*space_order;
    //z = 0;
    for (int m=1; m<=space_order; m++)
        for (int k=1; k<=space_order; k++)
        {
            double s1 = s1d.get(m);
            double s2 = s1d.get(k);

            z = z+1;

            spts.set(z,1, s1d.get(m) );
            spts.set(z,2, s1d.get(k) );
        }

}
