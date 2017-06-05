#include <cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"
#include "EulerParams.h"  

// This is a user-supplied routine that sets the the boundary conditions
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{

    void L2Project(int istart, int iend, int jstart, int jend,
        const dTensorBC4& q,
        const dTensorBC4& aux, dTensorBC4& Fout,
        void (*Func)(const dTensor2& xpts,
            const dTensor2& qvals,
            const dTensor2& auxvals,
            dTensor2& source));

// This call no longer exists, and broke this application 4-1-2013 (-DS)
//  void L2Project(const int istart, 
//          const int iend, 
//          const int jstart, 
//          const int jend,
//          const int QuadOrder, 
//          const int BasisOrder_fout,
//          dTensorBC4* fout,
//          void (*Func)(const dTensor2&,dTensor2&));
    const int space_order = dogParams.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  
    const double dx=dogParamsCart2.get_dx();
    const double xlow=dogParamsCart2.get_xlow();
    const double x0 = eulerParams.x0;

// TODO: how am I supposed to read in the time!? (-DS)
//  const double t = dogParams.get_time();
// printf("In SetBndValues, t = %2.15e\n", dogParams.get_time() );


    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // LEFT BOUNDARY
    // ********************************************************************* //
//  void LeftFunc(const dTensor2& xpts,
//          dTensor2& qvals);
    void LeftFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals);
    L2Project(1-mbc, 0, 1-mbc, my+mbc,  q, aux, q, &LeftFunc);     
    // ********************************************************************* //


    // ********************************************************************* //
    // RIGHT BOUNDARY
    // ********************************************************************* //
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    for (int m=1; m<=meqn; m++)
    for (int ell=1; ell<=kmax; ell++)
    {
        double tmp = q.get(mx,j,m,ell);                    
        q.set(i,j,m,ell, tmp );
    }
    // ********************************************************************* //


    // ********************************************************************* //
    // BOTTOM BOUNDARY (Reflective boundary conditions for right half)
    //
    // This is defined by extrapolating an even function for all components,
    // except the momentum u2 in the y-direction.  This function is extended
    // as an odd function.
    //
    // In order to accomodate this, we consider the form of the polynomials,
    // and extend them as appropriate by looking at whether or not they are
    // even or odd functions of nu, the variable that corresponds to y.
    //
    // This method is worked out up to third-order in space only.
    // ********************************************************************* //
    for (int i=1; i<=mx; i++)
    for (int j=0; j>=(1-mbc); j--)
    {
        double x = xlow + (double(i)-1.0)*dx;      

        if (x>=(x0-1.0e-12))
        {

            // Extrapolate an even function of the currenct cell into the
            // ghost cells.
            for (int m=1; m<=meqn; m++)
            for (int ell=1; ell<=kmax; ell++)
            {

                double tmp = q.get(i,1-j,m,1);                    
                q.set(i,j,m,1, tmp );

                if (dogParams.get_space_order()>1)
                {
                    tmp = q.get(i,1-j,m,2);                    
                    q.set(i,j,m,2,  tmp );

                    tmp = q.get(i,1-j,m,3);                    
                    q.set(i,j,m,3, -tmp );
                }

                if (dogParams.get_space_order()>2)
                {
                    tmp = q.get(i,1-j,m,4);                    
                    q.set(i,j,m,4, -tmp );

                    tmp = q.get(i,1-j,m,5);                    
                    q.set(i,j,m,5,  tmp );

                    tmp = q.get(i,1-j,m,6);                    
                    q.set(i,j,m,6,  tmp );
                }
            }

            // Extrapolate an "odd" function for the momentum, rho*u2.
            for (int ell=1; ell<=kmax; ell++)
            {
                double tmp;
                tmp = -q.get(i,j,3,ell);
                q.set(i,j,3,ell, tmp );
            }
        }

    }
    // Use the initial conditions in front of the wedge
    void BotFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
            dTensor2& qvals);
    L2Project(1-mbc, mx+mbc, 1-mbc, 0,  q, aux, q, &BotFunc);     
    // ********************************************************************* //


    // ***********************************************
    // TOP BOUNDARY
    //************************************************  
    void TopFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
            dTensor2& qvals);
    L2Project(1-mbc, mx+mbc, my+1, my+mbc,  q, aux, q, &TopFunc);     
    //************************************************

}      

// This function is idential to the "left hand" initial conditions of the
// problem.
void LeftFunc(const dTensor2& xpts,
    const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);
    const double gamma = eulerParams.gamma;

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double rho,u1,u2,u3,press;

        rho   =  8.0;
        u1    =  8.25*cos(pi/6.0);
        u2    = -8.25*sin(pi/6.0);
        u3    =  0.0;
        press =  116.5;

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );      
    }
}

// We use the same "initial conditions"
void BotFunc(const dTensor2& xpts,
    const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals)
{
    const int numpts   = xpts.getsize(1);
    const double gamma = eulerParams.gamma;
    const double x0    = eulerParams.x0;

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double rho,u1,u2,u3,press;

        if ( x < x0 )
        {
            rho   =  8.0;
            u1    =  8.25*cos(pi/6.0);
            u2    = -8.25*sin(pi/6.0);
            u3    =  0.0;
            press =  116.5;
        }
        else
        {
            rho   =  1.4;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  1.0;
        }	

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );      
    }
}

void TopFunc(const dTensor2& xpts,
    const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);
    const double gamma = eulerParams.gamma;
    const double x0    = eulerParams.x0;
    const double t     = dogParams.get_time();

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double rho,u1,u2,u3,press;
        if ( x < x0+(20.0*t+1.0)*osq3)
        {
            rho   =  8.0;
            u1    =  8.25*cos(pi/6.0);
            u2    = -8.25*sin(pi/6.0);
            u3    =  0.0;
            press =  116.5;
        }
        else
        {
            rho   =  1.4;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  1.0;
        }	

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );      
    }
}
