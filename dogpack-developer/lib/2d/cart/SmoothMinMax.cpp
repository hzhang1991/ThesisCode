#include <cmath>
#include "dog_math.h"
#include <iostream>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
using namespace std;

// Function to smooth out rough edges in the initial condition
void SmoothMinMax(dTensorBC4& q)
{

    // Parameters that need to be read-in
    int mcomp = 1;
    int mcomp1 = mcomp+1;
    int* complist = new int[mcomp1];

    double* xlow  = new double[mcomp1];
    double* ylow  = new double[mcomp1];
    double* xhigh = new double[mcomp1];
    double* yhigh = new double[mcomp1];

    int* istart = new int[mcomp1];
    int* iend   = new int[mcomp1];
    int* jstart = new int[mcomp1];
    int* jend   = new int[mcomp1];

    double* qmax_allowed = new double[mcomp1];
    double* qmin_allowed = new double[mcomp1];

    int msx = 11;
    int msy = 11;

    complist[1] = 1;

    xlow[1]  = 0.0;  
    ylow[1]  = 0.25;
    xhigh[1] = 0.5;
    yhigh[1] = 0.75;

    qmax_allowed[1] = 1.0;
    qmin_allowed[1] = 0.0;

    // Find indeces of affected elements
    double dx = dogParamsCart2.get_dx();
    double dy = dogParamsCart2.get_dy();
    double xlm = dogParamsCart2.get_xlow();
    double ylm = dogParamsCart2.get_ylow();

    for (int k=1; k<=mcomp; k++)
    {
        istart[k] = 1 + floor( (xlow[k] - xlm)/dx + 1.0e-8 );
        iend[k] =        ceil( (xhigh[k] - xlm)/dx - 1.0e-8 );

        jstart[k] = 1 + floor( (ylow[k] - ylm)/dy + 1.0e-8 );
        jend[k] =        ceil( (yhigh[k] - ylm)/dy - 1.0e-8 );
    }

    // sub-grid spacing
    double dx1 = 2.0/double(msx-1);
    double dy1 = 2.0/double(msy-1);
    dTensor3 phi(msx,msy,dogParams.get_kmax());

    // Sample Legendre polynomials on each element msx*msy times
    for (int i1=1; i1<=msx; i1++)
        for (int j1=1; j1<=msy; j1++)
        {	
            // grid point (xi,eta)
            const double xi   = -1.0 + double(i1-1)*dx1;
            const double eta  = -1.0 + double(j1-1)*dy1;
            const double xi2  = xi*xi;
            const double xi3  = xi*xi2;
            const double xi4  = xi*xi3;
            const double eta2 = eta*eta;
            const double eta3 = eta*eta2;
            const double eta4 = eta*eta3;

            // Legendre basis functions at each gaussian quadrature point in the
            // interval [-1,1]x[-1,1].
            switch( dogParams.get_space_order() )
            {
                case 5:  // fifth order                                 
                    phi.set(i1,j1,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                    phi.set(i1,j1,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
                    phi.set(i1,j1,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
                    phi.set(i1,j1,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                    phi.set(i1,j1,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

                case 4:  // fourth order
                    phi.set(i1,j1,10, sq7*(2.5*eta3 - 1.5*eta) );
                    phi.set(i1,j1,9,  sq7*(2.5*xi3 - 1.5*xi) );
                    phi.set(i1,j1,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                    phi.set(i1,j1,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

                case 3:  // third order
                    phi.set(i1,j1,6,  sq5*(1.5*eta2 - 0.5) );
                    phi.set(i1,j1,5,  sq5*(1.5*xi2 - 0.5) );
                    phi.set(i1,j1,4,  3.0*xi*eta );                  

                case 2:  // second order                
                    phi.set(i1,j1,3, sq3*eta );
                    phi.set(i1,j1,2, sq3*xi  );

                case 1:  // first order
                    phi.set(i1,j1,1, 1.0 );

                    break;                
            }
        }

    // Find min and max on each element
    //     if qmax>qmax_allowed or qmin<qmin_allowed
    //     then zero-out highest moments and repeat  
    for (int m=1; m<=mcomp; m++)    
    for (int i=istart[m]; i<=iend[m]; i++)
    for (int j=jstart[m]; j<=jend[m]; j++)	
    {

        int mstop = 0;
        int morder = dogParams.get_space_order();
        int kmax  = (morder*(morder+1))/2;
        int kless = (morder*(morder-1))/2;

        while (mstop==0 && morder>1)
        {
            double qmax = -1e8;
            double qmin =  1e8;
            for (int i1=1; i1<=msx; i1++)
            for (int j1=1; j1<=msy; j1++)
            {
                double qtmp = 0.0;
                for (int k=1; k<=kmax; k++)
                {
                    qtmp = qtmp + q.get(i,j,complist[m],k)*phi.get(i1,j1,k);
                }

                qmax = Max(qmax,qtmp);
                qmin = Min(qmin,qtmp);
            }

            if (qmax>qmax_allowed[m] || qmin<qmin_allowed[m])
            {
                for (int k=(kless+1); k<=kmax; k++)
                {
                    q.set(i,j,complist[m],k,  0.0e0 );
                }
                morder = morder - 1;
                kmax  = (morder*(morder+1))/2;
                kless = (morder*(morder-1))/2;
            }
            else
            {
                mstop = 1;
            }
        }

    }

    // Cleanup allocated arrays
    delete[] complist;
    delete[] xlow;
    delete[] ylow;
    delete[] xhigh;
    delete[] yhigh;

    delete[] istart;
    delete[] iend;
    delete[] jstart;
    delete[] jend;

    delete[] qmax_allowed;
    delete[] qmin_allowed;

}
