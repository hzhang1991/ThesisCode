#include "constants.h"
#include "tensors.h"

// Interpolate-then-integrate L0, L1, .... over several 
//    sub-elements of a time element of width dt
//
void ResInt(double dt, 
        const dTensorBC3& L0, 
        const dTensorBC3& L1, 
        const dTensorBC3& L2, 
        const dTensorBC3& L3, 
        const dTensorBC3& L4, 
        const dTensorBC3& L5,
        dTensorBC4& ILout)
{
    const int    mbc = ILout.getmbc();
    const int melems = ILout.getsize(1);
    const int   meqn = ILout.getsize(2);
    const int   kmax = ILout.getsize(3);
    const int  meth2 = 1+ILout.getsize(4);

    // Choose order of accuracy in integration
    switch( meth2 )
    {
        case 2:  // 2nd order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(melems+mbc); i++)
                for (int m=1; m<=meqn; m++)
                    for (int k=1; k<=kmax; k++)
                    {
                        double f0 = L0.get(i,m,k);
                        double f1 = L1.get(i,m,k);

                        double tmp = 0.5*dt*( f0 + f1 );
                        ILout.set(i,m,k,1, tmp );
                    }

            break;

        case 3:  // 3rd order in time

#pragma omp parallel for      
            for (int i=(1-mbc); i<=(melems+mbc); i++)
                for (int m=1; m<=meqn; m++)
                    for (int k=1; k<=kmax; k++)
                    {
                        double f0 = L0.get(i,m,k);
                        double f1 = L1.get(i,m,k);
                        double f2 = L2.get(i,m,k);

                        double tmp = dt/24.0*( 5.0*f0 + 8.0*f1 - f2 );
                        ILout.set(i,m,k,1, tmp );

                        tmp = dt/24.0*( 5.0*f2 + 8.0*f1 - f0 );
                        ILout.set(i,m,k,2, tmp );
                    }

            break;

        case 4:  // 4th order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(melems+mbc); i++)
                for (int m=1; m<=meqn; m++)
                    for (int k=1; k<=kmax; k++)
                    {
                        double f0 = L0.get(i,m,k);
                        double f1 = L1.get(i,m,k);
                        double f2 = L2.get(i,m,k);
                        double f3 = L3.get(i,m,k);

                        double tmp =  (dt/576.0)*(59.0*f0 - 14.0*f2 + 5.0*f3 + 94.0*f1);
                        ILout.set(i,m,k,1, tmp );

                        tmp = -(dt/36.0)*(2.0*f0 - 11.0*f2 + 2.0*f3 - 11.0*f1);
                        ILout.set(i,m,k,2, tmp );

                        tmp = (dt/576.0)*(5.0*f0 + 94.0*f2 + 59.0*f3 - 14.0*f1);
                        ILout.set(i,m,k,3, tmp );
                    }

            break;

        case 5:  // 5th order in time

#pragma omp parallel for
            for (int i=(1-mbc); i<=(melems+mbc); i++)
                for (int m=1; m<=meqn; m++)
                    for (int k=1; k<=kmax; k++)
                    {
                        double f0 = L0.get(i,m,k);
                        double f1 = L1.get(i,m,k);
                        double f2 = L2.get(i,m,k);
                        double f3 = L3.get(i,m,k);
                        double f4 = L4.get(i,m,k);

                        double tmp = (dt/480.0)*((23.0+4.0*sq2)*f0 + (-13.0*sq2+64.0)*f1 
                                + (-72.0*sq2+96.0)*f2 + (-43.0*sq2+64.0)*f3
                                + (-7.0+4.0*sq2)*f4);
                        ILout.set(i,m,k,1, tmp );

                        tmp = (sq2*dt/960.0)*((-8.0-15.0*sq2)*f0 + 146.0*f1
                                + 144.0*f2 - 34.0*f3 
                                + (-8.0+15.0*sq2)*f4);
                        ILout.set(i,m,k,2, tmp );

                        tmp = (sq2*dt/960.0)*((-8.0+15.0*sq2)*f0 - 34.0*f1 
                                + 144.0*f2 + 146.0*f3 
                                + (-8.0-15.0*sq2)*f4);
                        ILout.set(i,m,k,3, tmp );

                        tmp = (dt/480.0)*((-7.0+4.0*sq2)*f0 + (-43.0*sq2+64.0)*f1 
                                + (-72.0*sq2+96.0)*f2 + (-13.0*sq2+64.0)*f3
                                + (23.0+4.0*sq2)*f4);
                        ILout.set(i,m,k,4, tmp );
                    }

            break;

        case 6:

            // integrate the right hand side for each subinterval:
            dTensor2 M(6,5);

            M.set(1,1, 19.0/288.0);
            M.set(1,2, -3.0/800.0);
            M.set(1,3, 11.0/7200.0);
            M.set(1,4, -11.0/7200.0);
            M.set(1,5,    3.0/800.0); 
            M.set(2,1, 1427.0/7200.0);
            M.set(2,2, 637.0/7200.0);
            M.set(2,3, -31.0/2400.0);
            M.set(2,4, 77.0/7200.0);
            M.set(2,5, -173.0/7200.0); 
            M.set(3,1, -133.0/1200.0);
            M.set(3,2, 511.0/3600.0);
            M.set(3,3, 401.0/3600.0);
            M.set(3,4, -43.0/1200.0);
            M.set(3,5, 241.0/3600.0);
            M.set(4,1, 241.0/3600.0);
            M.set(4,2, -43.0/1200.0);
            M.set(4,3, 401.0/3600.0);
            M.set(4,4, 511.0/3600.0);
            M.set(4,5, -133.0/1200.0); 
            M.set(5,1, -173.0/7200.0);
            M.set(5,2, 77.0/7200.0);
            M.set(5,3, -31.0/2400.0);
            M.set(5,4, 637.0/7200.0);
            M.set(5,5, 1427.0/7200.0); 
            M.set(6,1, 3.0/800.0);
            M.set(6,2, -11.0/7200.0);
            M.set(6,3, 11.0/7200.0);
            M.set(6,4, -3.0/800.0);
            M.set(6,5, 19.0/288.0);


            // integration based on uniformly chosen quadrature points
#pragma omp parallel for
            for (int i=(1-mbc); i<=(melems+mbc); i++)
                for (int m=1; m<=meqn; m++)
                    for (int k=1; k<=kmax; k++)
                    {
                        double f0 = L0.get(i,m,k);
                        double f1 = L1.get(i,m,k);
                        double f2 = L2.get(i,m,k);
                        double f3 = L3.get(i,m,k);
                        double f4 = L4.get(i,m,k);
                        double f5 = L5.get(i,m,k);
                        double tmp = 0.;

                        // in this case, dt = large time step dt ...
                        tmp = dt*( M.get(1,1) * f0 + M.get(2,1) * f1 + M.get(3,1) * f2 + M.get(4,1) * f3 + M.get(5,1) * f4  + M.get(6,1) * f5 );
                        ILout.set(i,m,k,1, tmp );

                        tmp = dt*( M.get(1,2) * f0 + M.get(2,2) * f1 + M.get(3,2) * f2 + M.get(4,2) * f3 + M.get(5,2) * f4  + M.get(6,2) * f5 );
                        ILout.set(i,m,k,2, tmp );

                        tmp = dt*( M.get(1,3) * f0 + M.get(2,3) * f1 + M.get(3,3) * f2 + M.get(4,3) * f3 + M.get(5,3) * f4  + M.get(6,3) * f5 );
                        ILout.set(i,m,k,3, tmp );

                        tmp = dt*( M.get(1,4) * f0 + M.get(2,4) * f1 + M.get(3,4) * f2 + M.get(4,4) * f3 + M.get(5,4) * f4  + M.get(6,4) * f5 );
                        ILout.set(i,m,k,4, tmp );

                        tmp = dt*( M.get(1,5) * f0 + M.get(2,5) * f1 + M.get(3,5) * f2 + M.get(4,5) * f3 + M.get(5,5) * f4  + M.get(6,5) * f5 );
                        ILout.set(i,m,k,5, tmp );
                    }

            break;

    }

}
