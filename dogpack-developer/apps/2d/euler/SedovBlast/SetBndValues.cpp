#include "dogdefs.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
//              TODO - I do not believe this comment - I think this looks like
//              reflective on left and bottom, not zeroth-order extrapolation
//              everywhere.  (12/23/2014) -DS.
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for (int i=0; i>=(1-mbc); i--)
    for (int j=1; j<=my; j++)
    for (int ell=1; ell<=kmax; ell++)
    {
        double tmp = q.get(1-i,j,1,ell);
        q.set(i,j,1,ell, tmp );

        tmp = q.get(1-i,j,2,ell);
        q.set(i,j,2,ell, -tmp );

        tmp = q.get(1-i,j,3,ell);
        q.set(i,j,3,ell, tmp );

        tmp = q.get(1-i,j,4,ell);
        q.set(i,j,4,ell, tmp );

        tmp = q.get(1-i,j,5,ell);
        q.set(i,j,5,ell, tmp );
    }
    // ***********************************************


    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    for (int m=1; m<=meqn; m++)
    for (int ell=1; ell<=kmax; ell++)
    {
        double tmp = q.get(mx,j,m,ell);                    
        q.set(i,j,m,ell, tmp );
    }
    // ***********************************************


    // ***********************************************
    // BOTTOM BOUNDARY
    // ***********************************************
    for (int j=0; j>=(1-mbc); j--)
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    for (int ell=1; ell<=kmax; ell++)
    {
        double tmp = q.get(i,1-j,1,ell);
        q.set(i,j,1,ell, tmp );

        tmp = q.get(i,1-j,2,ell);
        q.set(i,j,2,ell, tmp );

        tmp = q.get(i,1-j,3,ell);
        q.set(i,j,3,ell, -tmp );

        tmp = q.get(i,1-j,4,ell);
        q.set(i,j,4,ell, tmp );

        tmp = q.get(i,1-j,5,ell);
        q.set(i,j,5,ell, tmp );

    }
    // ***********************************************


    // ***********************************************
    // TOP BOUNDARY
    // ***********************************************
    for (int j=(my+1); j<=(my+mbc); j++)
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)      
    for (int m=1; m<=meqn; m++)
    for (int ell=1; ell<=kmax; ell++)
    {
        double tmp = q.get(i,my,m,ell);                    
        q.set(i,j,m,ell, tmp );
    }
    // ***********************************************

}
