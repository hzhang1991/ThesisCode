#include <cmath>
#include "dogdefs.h"
#include "DogParamsCart2.h"
#include <iostream>
using namespace std;

// Euler equations - shock-diffraction problem.  This routine is used to set
// the boundary conditions for a problem with "Geometry".  See
//
// D.C. Seal, Q. Tang, Z. Xu and A.J. Christlieb, "An explicit high-order
// single-stage single-step positivity-preserving finite difference WENO method
// for the compressible Euler equations", (http://arxiv.org/abs/1411.0328) and 
// references therein for a description of the problem.
//
// This requires setting a number of different boundary conditions:
//
//     * Hard surface BC's for every point along the wedge
//     * Either zero-extrapolation, or apply the IC's for every other
//       boundary.
//
// Moreover, this routine is used for defining positivity-preserving schemes,
// so of these routines, we also set BCs for the low-order problem as well.
//
// See also: $DOGPACK/apps/2d/euler/mach_reflection for a simpler case of
// reflective boundary conditions.

// Initial conditions
void qinitfuncs(double x,double y,double &rho,double &press,double &u1,double &u2,double &u3,double &energy)
{

    const double gamma = 1.4;
    const double gm1   = 1.4 - 1.0;
    const double gp1   = 1.4 + 1.0;

    const double Mach  = 5.09;
    const double M2    = Mach*Mach;

    const double vs = Mach*sqrt( gamma*press/rho );

    // Correct with the Mach number for incoming shock
    {
        press = press*(2.0*gamma*M2- gm1)/gp1;
        rho   = rho*gp1*M2 / (gm1*M2 + 2.0);
        u1    = vs*(1.0 - 1.4 / rho );
        energy = press/(gamma-1.0e0)
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);
    }
    if(press<0.0){press=1.0e-12;}
}

// Vector defining sign of (even/odd) polynomials.
// That is, for the polynomials,
//
//    phi = {1, xi, eta, xi*eta, xi^2, eta^2 },
//
// Which of these define even/odd functions of xi/eta.  This can be used to
// extract hard-surface boundary conditions for the unkowns.
double xi_sign[]  = {1.0, -1.0, 1.0, -1.0, 1.0, 1.0};
double eta_sign[] = {1.0, 1.0, -1.0, -1.0, 1.0, 1.0};

// This is a user-supplied routine that sets the the boundary conditions
//
// Euler Equations - shock-diffraction problem.
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  
    const int kmax = q.getsize(4);

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    double gamma=1.4;

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;

    double energy = press/(gamma-1.0e0)
        + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }


    // Compute index where the step is located.  
    // q(istep,:), and q(:,jstep) are inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // ********************************************************************* //
    // Reset every value inside the wedge to a junk value.
    // In AfterFullTimeStep, this will be reset to zero for plotting purposes.
    // ********************************************************************* //
    for( int i = 1-mbc; i < istep; i++)
    for( int j = 1-mbc; j < jstep; j++ )
    {

        q.set(i,j,1,1, rho );
        q.set(i,j,2,1, u1 );
        q.set(i,j,3,1, u2 );
        q.set(i,j,4,1, u3 );
        q.set(i,j,5,1, energy );

        for(int m=1; m<=meqn; m++)
        for(int k=2; k<=kmax; k++)
        { q.set(i,j,m,k, 0.0 ); }

    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // RIGHT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    {
        q.set(i,j,1,1, rho );
        q.set(i,j,2,1, u1 );
        q.set(i,j,3,1, u2 );
        q.set(i,j,4,1, u3 );
        q.set(i,j,5,1, energy );

        // Initial conditions are constants, so higher modes can be zeroed out
        for (int m=1; m<=meqn; m++)
        for (int k=2; k<=kmax; k++)
        { q.set(i,j,m,k, 0.0 ); }

    }


    // ********************************************************************* //

    // ********************************************************************* //
    // TOP BOUNDARY "outflow" - zeroth-order extrapolation
    // ********************************************************************* //
    for (int i=1-mbc; i<= mx+mbc; i++)
    for (int j=my+1; j<=my+mbc; j++)
    for (int m=1; m<=meqn; m++)
    {
        for(int k=1;k<=kmax;k++)
        {
            double tmp = q.get(i,my,m,k);                    
            q.set(i,j,m,k, tmp );
        }
    }
    // ********************************************************************* //

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // The first two of these regions are hard-surface boundaries.  The third
    // will be zeroth-order extrapolation.
    //
    // ********************************************************************* //


    // Segment 1
    for( int i = 1;       i <= istep-1; i++)
    for( int j = jstep-1; j >= jstep-mbc; j-- )
    {
        for(int k=1; k<=kmax; k++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = eta_sign[k-1]*q.get(i, 2*jstep-1-j, m,k);                    
                q.set(i,j,m,k, tmp );
            }
            // Flip the momemtum u2:
            q.set(i,j,3,k, -q.get(i,j,3,k) );
        }
    }

    // Segment 2
    for (int i = istep-1; i >= istep-mbc; i--)
    for (int j = 1; j < jstep; j++ )
    {
        for (int m=1; m<=meqn; m++)
        {
            for(int k=1;k<=kmax;k++)
            {

                double tmp = xi_sign[k-1]*q.get(2*istep+1-i, j, m,k);                    
                q.set(i,j,m,k, tmp );
            }
        }
        for(int k=1;k<=kmax;k++)
        {
            // Flip the momemtum u1:
            q.set(i,j,2,k, -q.get(i,j,2,k) );
        }
    }

    // Segment 3 (Changed on 12/23/2014 -DS)
    for (int i=istep; i<= mx+mbc; i++)
    for (int j=0; j>=(1-mbc); j--)
    for (int k=1; k<=kmax; k++)
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(i,1, m,k);                    
            q.set(i,j,m,k, tmp );
        }
    }

    // ********************************************************************* //
    // LEFT BOUNDARY (the part above the wall)
    //
    //     We consider the section x=0, 6 < y < 11 here.
    //     the other part of the "left" boundary has already been dealt with.
    // ********************************************************************* //
    qinitfuncs(0.0, 7.0, rho, press, u1, u2, u3, energy);
    for (int i=1-mbc; i<= 0; i++)
    for (int j=jstep; j<=my+mbc; j++)
    {

        q.set(i,j,1,1, rho    );
        q.set(i,j,2,1, rho*u1 );
        q.set(i,j,3,1, rho*u2 );
        q.set(i,j,4,1, rho*u3 );
        q.set(i,j,5,1, energy );

        // Initial conditions are constants, so higher modes can be zeroed out
        for(int m=1; m<=meqn; m++)
        for(int k=2; k<=kmax; k++)
        { q.set(i,j,m,k, 0.0 ); }

    }


}      

void SetBndValuesX(dTensorBC4& q, dTensorBC4& aux)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  
    const int kmax = q.getsize(4);

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    double gamma=1.4;

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;

    double energy = press/(gamma-1.0e0)
        + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // RIGHT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    {
        q.set(i,j,1,1, rho );
        q.set(i,j,2,1, u1);
        q.set(i,j,3,1, u2 );
        q.set(i,j,4,1, u3 );
        q.set(i,j,5,1, energy );

        // Initial conditions are constants, so higher modes can be zeroed out
        for (int m=1; m<=meqn; m++)
        for (int k=2; k<=kmax; k++)
        { q.set(i,j,m,k, 0.0 ); }

    }


    // ********************************************************************* //

    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:), and q(:,jstep) are inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 2
    for (int i = istep-1; i >= istep-mbc; i--)
    for (int j = 1; j < jstep; j++ )
    {
        for (int m=1; m<=meqn; m++)
        {
            for(int k=1;k<=kmax;k++)
            {

                double tmp = xi_sign[k-1]*q.get(2*istep+1-i, j, m,k);                    
                q.set(i,j,m,k, tmp );
            }
        }
        for(int k=1;k<=kmax;k++)
        {
            // Flip the momemtum u1:
            q.set(i,j,2,k, -q.get(i,j,2,k) );
        }
    }


    // ********************************************************************* //
    // LEFT BOUNDARY (the part above the wall)
    //
    //     We consider the section x=0, 6 < y < 11 here.
    //     the other part of the "left" boundary has already been dealt with.
    // ********************************************************************* //
    qinitfuncs(0.0, 7.0, rho, press, u1, u2, u3, energy);
    for (int i=1-mbc; i<= 0; i++)
    for (int j=jstep; j<=my+mbc; j++)
    {
        q.set(i,j,1,1, rho );
        q.set(i,j,2,1, rho*u1 );
        q.set(i,j,3,1, rho*u2 );
        q.set(i,j,4,1, rho*u3 );
        q.set(i,j,5,1, energy );

        // Initial conditions are constants, so higher modes can be zeroed out
        for(int m=1; m<=meqn; m++)
        for(int k=2; k<=kmax; k++)
        { q.set(i,j,m,k, 0.0 ); }

    }

}      

void SetBndValuesY(dTensorBC4& q, dTensorBC4& aux)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  
    const int kmax = q.getsize(4);

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    double gamma=1.4;

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;

    double energy = press/(gamma-1.0e0)
        + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //


    // ********************************************************************* //
    // TOP BOUNDARY "outflow" - zeroth-order extrapolation
    // ********************************************************************* //
    for (int i=1-mbc; i<= mx+mbc; i++)
    for (int j=my+1; j<=my+mbc; j++)
    for (int m=1; m<=meqn; m++)
    {
        for(int k=1;k<=kmax;k++)
        {
            double tmp = q.get(i,my,m,k);                    
            q.set(i,j,m,k, tmp );
        }
    }
    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:), and q(:,jstep) are inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 1
    for( int i = 1;       i <= istep-1; i++)
    for( int j = jstep-1; j >= jstep-mbc; j-- )
    {
        for (int m=1; m<=meqn; m++)
        {
            for(int k=1;k<=kmax;k++)
            {
                double tmp = eta_sign[k-1]*q.get(i, 2*jstep-1-j, m,k);                    
                q.set(i,j,m,k, tmp );
            }
        }
        for(int k=1;k<=kmax;k++)
        {
            // Flip the momemtum u2:
            q.set(i,j,3,k, -q.get(i,j,3,k) );}
    }

    // Segment 3 (Changed on 12/23/2014 -DS)
    for (int i=istep; i<= mx+mbc; i++)
    for (int j=0; j>=(1-mbc); j--)
    for (int k=1; k<=kmax; k++)
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(i,1, m,k);                    
            q.set(i,j,m,k, tmp );
        }
    }

}      


// This is a user-supplied routine that sets the the boundary conditions
//
// Euler Equations - shock-diffraction problem.
//
// This "duplicate" of the above routine is intended to be used with
// ConstructLFL, which only observes the cell averages of the polynomial.
void SetBndValues(dTensorBC3& q, dTensorBC3& aux)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    double gamma=1.4;

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;

    double energy = press/(gamma-1.0e0)
        + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }


    // Compute index where the step is located.  
    // q(istep,:), and q(:,jstep) are inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    for( int i = 1-mbc; i < istep; i++)
    for( int j = 1-mbc; j < jstep; j++ )
    {
        q.set(i,j,1, rho );
        q.set(i,j,2, u1 );
        q.set(i,j,3, u2 );
        q.set(i,j,4, u3 );
        q.set(i,j,5, energy );
    }


    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // RIGHT BOUNDARY (copy initial conditions)
    // ********************************************************************* //

    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    {
        q.set(i,j,1, rho );
        q.set(i,j,2, u1 );
        q.set(i,j,3, u2 );
        q.set(i,j,4, u3 );
        q.set(i,j,5, energy );
    }


    // ********************************************************************* //

    // ********************************************************************* //
    // TOP BOUNDARY "outflow" - zeroth-order extrapolation
    // ********************************************************************* //
    for (int i=1-mbc; i<= mx+mbc; i++)
    for (int j=my+1; j<=my+mbc; j++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp = q.get(i,my,m);                    
        q.set(i,j,m, tmp );
    }
    // ********************************************************************* //

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 1
    for( int i = 1;       i <= istep-1; i++)
    for( int j = jstep-1; j >= jstep-mbc; j-- )
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(i, 2*jstep-1-j, m);                    
            q.set(i,j,m, tmp );
        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );
    }

    // Segment 2
    for (int i = istep-1; i >= istep-mbc; i--)
    for (int j = 1; j < jstep; j++ )
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(2*istep+1-i, j, m);                    
            q.set(i,j,m, tmp );
        }
        // Flip the momemtum u1:
        q.set(i,j,2, -q.get(i,j,2) );
    }

    // Segment 3
    for (int i=istep; i<= mx+mbc; i++)
    for (int j=0; j>=(1-mbc); j--)
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i,1-j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );
    }

    // ********************************************************************* //
    // LEFT BOUNDARY (the part above the wall)
    //
    //     We consider the section x=0, 6 < y < 11 here.
    //     the other part of the "left" boundary has already been dealt with.
    // ********************************************************************* //
    qinitfuncs(0.0, 7.0, rho, press, u1, u2, u3, energy);
    for (int i=1-mbc; i<= 0; i++)
    for (int j=jstep; j<=my+mbc; j++)
    {

        q.set(i,j,1, rho    );
        q.set(i,j,2, rho*u1 );
        q.set(i,j,3, rho*u2 );
        q.set(i,j,4, rho*u3 );
        q.set(i,j,5, energy );

    }


}      

void SetBndValuesX(dTensorBC3& q, dTensorBC3& aux)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    double gamma=1.4;

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;

    double energy = press/(gamma-1.0e0)
        + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // RIGHT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    {
        q.set(i,j,1, rho );
        q.set(i,j,2, u1);
        q.set(i,j,3, u2 );
        q.set(i,j,4, u3 );
        q.set(i,j,5, energy );
    }
    // ********************************************************************* //

    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:), and q(:,jstep) are inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 2
    for (int i = istep-1; i >= istep-mbc; i--)
    for (int j = 1; j < jstep; j++ )
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(2*istep+1-i, j, m);                    
            q.set(i,j,m, tmp );
        }
            // Flip the momemtum u1:
        q.set(i,j,2, -q.get(i,j,2) );
    }


    // ********************************************************************* //
    // LEFT BOUNDARY (the part above the wall)
    //
    //     We consider the section x=0, 6 < y < 11 here.
    //     the other part of the "left" boundary has already been dealt with.
    // ********************************************************************* //
    qinitfuncs(0.0, 7.0, rho, press, u1, u2, u3, energy);
    for (int i=1-mbc; i<= 0; i++)
    for (int j=jstep; j<=my+mbc; j++)
    {
        q.set(i,j,1, rho );
        q.set(i,j,2, rho*u1 );
        q.set(i,j,3, rho*u2 );
        q.set(i,j,4, rho*u3 );
        q.set(i,j,5, energy );

    }

}      

void SetBndValuesY(dTensorBC3& q, dTensorBC3& aux)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    double gamma=1.4;

    double rho   =  1.4;
    double u1    =  0.0;
    double u2    =  0.0;
    double u3    =  0.0;
    double press =  1.0;

    double energy = press/(gamma-1.0e0)
        + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //


    // ********************************************************************* //
    // TOP BOUNDARY "outflow" - zeroth-order extrapolation
    // ********************************************************************* //
    for (int i=1-mbc; i<= mx+mbc; i++)
    for (int j=my+1; j<=my+mbc; j++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp = q.get(i,my,m);                    
        q.set(i,j,m, tmp );
    }
    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:), and q(:,jstep) are inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 1
    for( int i = 1;       i <= istep-1; i++)
    for( int j = jstep-1; j >= jstep-mbc; j-- )
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(i, 2*jstep-1-j, m);                    
            q.set(i,j,m, tmp );
        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );
    }

    // Segment 3
    for (int i=istep; i<= mx+mbc; i++)
    for (int j=0; j>=(1-mbc); j--)
    {

        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i,1-j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );

    }


}
