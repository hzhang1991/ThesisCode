#include "LegendrePolys1d.h"

// Evaluate the Legendre polynomials at the list of points given by spts.
// 
// Input:
// ------
//
//      spts( 1:numpts )   - quadrature points ( in the interval [-1,1] )
//
// Returns:
// --------
//
//       phi( 1:numpts, 1:MAX_ORDER ) - Legendre polynomials evaluated at each
//                                      point.
//
// For simplicity, we evaluate more polynomials than are really necessary.  (The
// 1D code runs fast enough to handle this.) 
//
// See also: $DOGPACK/lib/Quadrature.cpp
void evaluateLegendrePolys(const dTensor1& spts, dTensor2& phi )
{

    const int numpts = spts.getsize();

    for (int m=1; m<=numpts; m++)
    {
        // Legendre basis functions at each grid point
        phi.set( m,1, 1.0 );
        phi.set( m,2, sq3*spts.get(m) );
        phi.set( m,3, 0.5*sq5*( 3.0*pow(spts.get(m),2) - 1.0 ) );
        phi.set( m,4, 0.5*sq7*spts.get(m)
                *(5.0*pow(spts.get(m),2) - 3.0) );
        phi.set( m,5, (105.0/8.0)*pow(spts.get(m),4) 
                - (45.0/4.0)*pow(spts.get(m),2) + (9.0/8.0) );
        phi.set( m, 6, (63.0/8.0)*sq11 * ( pow(spts.get(m),5)  -
                    (10.0/9.0)*(pow(spts.get(m),3) ) + (5.0/21.0)*spts.get(m) ) );
    }

}

// function that evaluates legendre polynomials and all their derivatives
// at the list of points given by spts
void evaluateLegendrePolys( 
        const double dx, const dTensor1& spts, dTensor2& phi, dTensor2& phi_x)
{

    // Regular Basis functions
    evaluateLegendrePolys(spts, phi );

    // 1st Derivative of Basis functions
    const int numpts = spts.getsize();
    for (int m=1; m<=numpts; m++)
    {
        phi_x.set( m,1, 0.0 );
        phi_x.set( m,2, 2.0*sq3/dx );
        phi_x.set( m,3, 6.0*sq5*spts.get(m)/dx );
        phi_x.set( m,4, 3.0*sq7*(5.0*pow(spts.get(m),2)-1.0)/dx );
        phi_x.set( m,5, 15.0*spts.get(m)*
                (7.0*pow(spts.get(m),2)-3.0)/dx );
        // TODO - include the sixth-order derivative here!
    }

}
