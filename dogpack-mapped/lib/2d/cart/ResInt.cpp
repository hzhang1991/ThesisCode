#include "dogdefs.h"

// -------------------------------------------------------------------------- //
// This module describes the integrtion process on the residual which is used
// for SDC time stepping.
//
// Interpolate-then-integrate L0, L1, .... over several 
//    sub-elements of a time element of width dt
//
// Other options for a selection of time points exist, but are not implemented
// in this code.  See, SetSDCtimePoints for which time points are used for the
// integration process.
//
// TODO: what's the difference between this function and the one that's later
// in this file?  It would be nice to move this into a separate
// time-integrator, along with all of the RK time stepping options in order to
// reuse this code for the different parts (e.g. 1d, 2d/cart, 2d/unst, ... ).
// ( -DS )
//
// -------------------------------------------------------------------------- //
void ResInt(double dt, dTensorBC4 const*const*const L, dTensorBC5& ILout)
{
    // convenience aliases
    //
    const dTensorBC4& L0 = *L[0];
    const dTensorBC4& L1 = *L[1];
    const dTensorBC4& L2 = *L[2];
    const dTensorBC4& L3 = *L[3];
    const dTensorBC4& L4 = *L[4];

    const int    mbc = ILout.getmbc();
    const int     mx = ILout.getsize(1);
    const int     my = ILout.getsize(2);
    const int   meqn = ILout.getsize(3);
    const int   kmax = ILout.getsize(4);
    const int  meth2 = 1+ILout.getsize(5);

    // Choose order of accuracy in integration
    switch( meth2 )
    {
        case 2:  // 2nd order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	      
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = L0.get(i,j,m,k);
                            const double f1 = L1.get(i,j,m,k);

                            double tmp = 0.5*dt*( f0 + f1 );
                            ILout.set(i,j,m,k,1, tmp );
                        }

            break;

        case 3:  // 3rd order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = L0.get(i,j,m,k);
                            const double f1 = L1.get(i,j,m,k);
                            const double f2 = L2.get(i,j,m,k);

                            double tmp;
                            tmp = dt/24.0*( 5.0*f0 + 8.0*f1 - f2 );
                            ILout.set(i,j,m,k,1, tmp );

                            tmp = dt/24.0*( 5.0*f2 + 8.0*f1 - f0 );
                            ILout.set(i,j,m,k,2, tmp );
                        }

            break;

        case 4:  // 4th order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = L0.get(i,j,m,k);
                            const double f1 = L1.get(i,j,m,k);
                            const double f2 = L2.get(i,j,m,k);
                            const double f3 = L3.get(i,j,m,k);

                            double tmp;
                            tmp =  (dt/576.0)*(59.0*f0 - 14.0*f2 + 5.0*f3 + 94.0*f1);
                            ILout.set(i,j,m,k,1, tmp );

                            tmp = -(dt/36.0)*(2.0*f0 - 11.0*f2 + 2.0*f3 - 11.0*f1);
                            ILout.set(i,j,m,k,2, tmp );

                            tmp = (dt/576.0)*(5.0*f0 + 94.0*f2 + 59.0*f3 - 14.0*f1);
                            ILout.set(i,j,m,k,3, tmp );
                        }

            break;

        case 5:  // 5th order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = L0.get(i,j,m,k);
                            const double f1 = L1.get(i,j,m,k);
                            const double f2 = L2.get(i,j,m,k);
                            const double f3 = L3.get(i,j,m,k);
                            const double f4 = L4.get(i,j,m,k);

                            double tmp;
                            tmp = (dt/480.0)*((23.0+4.0*sq2)*f0 + (-13.0*sq2+64.0)*f1 
                                    + (-72.0*sq2+96.0)*f2 + (-43.0*sq2+64.0)*f3
                                    + (-7.0+4.0*sq2)*f4);
                            ILout.set(i,j,m,k,1, tmp );

                            tmp = (sq2*dt/960.0)*((-8.0-15.0*sq2)*f0 + 146.0*f1
                                    + 144.0*f2 - 34.0*f3 
                                    + (-8.0+15.0*sq2)*f4);
                            ILout.set(i,j,m,k,2, tmp );

                            tmp = (sq2*dt/960.0)*((-8.0+15.0*sq2)*f0 - 34.0*f1 
                                    + 144.0*f2 + 146.0*f3 
                                    + (-8.0-15.0*sq2)*f4);
                            ILout.set(i,j,m,k,3, tmp );

                            tmp = (dt/480.0)*((-7.0+4.0*sq2)*f0 + (-43.0*sq2+64.0)*f1 
                                    + (-72.0*sq2+96.0)*f2 + (-13.0*sq2+64.0)*f3
                                    + (23.0+4.0*sq2)*f4);
                            ILout.set(i,j,m,k,4, tmp );
                        }

            break;

    }

}

// -------------------------------------------------------------------------- //
//
// Interpolate-then-integrate L0, L1, .... over several 
//    sub-elements of a time element of width dt
//
// -------------------------------------------------------------------------- //
void ResInt(double dt, dTensorBC4* L[], dTensorBC5& ILout)
{
    const int    mbc = ILout.getmbc();
    const int     mx = ILout.getsize(1);
    const int     my = ILout.getsize(2);
    const int   meqn = ILout.getsize(3);
    const int   kmax = ILout.getsize(4);
    const int  meth2 = 1+ILout.getsize(5);

    // Choose order of accuracy in integration
    switch( meth2 )
    {
        case 2:  // 2nd order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	      
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = (*L[0]).get(i,j,m,k);
                            const double f1 = (*L[1]).get(i,j,m,k);

                            double tmp = 0.5*dt*( f0 + f1 );
                            ILout.set(i,j,m,k,1, tmp );
                        }

            break;


        case 3:  // 3rd order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = (*L[0]).get(i,j,m,k);
                            const double f1 = (*L[1]).get(i,j,m,k);
                            const double f2 = (*L[2]).get(i,j,m,k);

                            double tmp;
                            tmp = dt/24.0*( 5.0*f0 + 8.0*f1 - f2 );
                            ILout.set(i,j,m,k,1, tmp );

                            tmp = dt/24.0*( 5.0*f2 + 8.0*f1 - f0 );
                            ILout.set(i,j,m,k,2, tmp );
                        }

            break;

        case 4:  // 4th order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = (*L[0]).get(i,j,m,k);
                            const double f1 = (*L[1]).get(i,j,m,k);
                            const double f2 = (*L[2]).get(i,j,m,k);
                            const double f3 = (*L[3]).get(i,j,m,k);

                            double tmp;
                            tmp =  (dt/576.0)*(59.0*f0 - 14.0*f2 + 5.0*f3 + 94.0*f1);
                            ILout.set(i,j,m,k,1, tmp );

                            tmp = -(dt/36.0)*(2.0*f0 - 11.0*f2 + 2.0*f3 - 11.0*f1);
                            ILout.set(i,j,m,k,2, tmp );

                            tmp = (dt/576.0)*(5.0*f0 + 94.0*f2 + 59.0*f3 - 14.0*f1);
                            ILout.set(i,j,m,k,3, tmp );
                        }

            break;

        case 5:  // 5th order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)	
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            const double f0 = (*L[0]).get(i,j,m,k);
                            const double f1 = (*L[1]).get(i,j,m,k);
                            const double f2 = (*L[2]).get(i,j,m,k);
                            const double f3 = (*L[3]).get(i,j,m,k);
                            const double f4 = (*L[4]).get(i,j,m,k);

                            double tmp;
                            tmp = (dt/480.0)*((23.0+4.0*sq2)*f0 + (-13.0*sq2+64.0)*f1 
                                    + (-72.0*sq2+96.0)*f2 + (-43.0*sq2+64.0)*f3
                                    + (-7.0+4.0*sq2)*f4);
                            ILout.set(i,j,m,k,1, tmp );

                            tmp = (sq2*dt/960.0)*((-8.0-15.0*sq2)*f0 + 146.0*f1
                                    + 144.0*f2 - 34.0*f3 
                                    + (-8.0+15.0*sq2)*f4);
                            ILout.set(i,j,m,k,2, tmp );

                            tmp = (sq2*dt/960.0)*((-8.0+15.0*sq2)*f0 - 34.0*f1 
                                    + 144.0*f2 + 146.0*f3 
                                    + (-8.0-15.0*sq2)*f4);
                            ILout.set(i,j,m,k,3, tmp );

                            tmp = (dt/480.0)*((-7.0+4.0*sq2)*f0 + (-43.0*sq2+64.0)*f1 
                                    + (-72.0*sq2+96.0)*f2 + (-13.0*sq2+64.0)*f3
                                    + (23.0+4.0*sq2)*f4);
                            ILout.set(i,j,m,k,4, tmp );
                        }
            break;
    }
}
// -------------------------------------------------------------------------- //



// -------------------------------------------------------------------------- //
//
// These are the time points used for the integration of the Residual.
//
// (moved from DogSolveSDC, because these need to be compatible with the
// integration formulas presented in this module).
//
// t     = current time value
// dt    = time step size
// tvec  = vector containing intermediate time values for sdc time stepping
// dtvec = vector containing intermediate time step values for sdc 
//		time stepping
//
// -------------------------------------------------------------------------- //
void SetSDCtimePoints(int morder, double t, double dt, 
        dTensor1& dtvec, dTensor1& tvec)
{
    // --------------------------------------------------------------
    // Select Gauss-Lobatto points in time and
    // construct the initial right-hand side for SDC
    //
    // These time points need to be compatible with
    //
    //
    switch(morder)
    {
        case 2: 
            dtvec.set(1, dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );	  
            break;

        case 3:
            dtvec.set(1, 0.5*dt );
            dtvec.set(2, 0.5*dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );
            tvec.set(3, tvec.get(2) + dtvec.get(2) );
            break;

        case 4:
            dtvec.set(1, 0.25*dt );
            dtvec.set(2, 0.50*dt );
            dtvec.set(3, 0.25*dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );
            tvec.set(3, tvec.get(2) + dtvec.get(2) );
            tvec.set(4, tvec.get(3) + dtvec.get(3) );
            break;

        case 5:
            dtvec.set(1, (0.5 - 0.25*sq2) * dt );
            dtvec.set(2,       (0.25*sq2) * dt );
            dtvec.set(3,       (0.25*sq2) * dt );
            dtvec.set(4, (0.5 - 0.25*sq2) * dt );

            tvec.set(1, t );
            tvec.set(2, tvec.get(1) + dtvec.get(1) );
            tvec.set(3, tvec.get(2) + dtvec.get(2) );
            tvec.set(4, tvec.get(3) + dtvec.get(3) );
            tvec.set(5, tvec.get(4) + dtvec.get(4) );
            break;
    }
}
// -------------------------------------------------------------------------- //
