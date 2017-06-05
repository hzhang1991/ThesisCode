#include "tensors.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{
    int i,j,ell,m;
    double tmp;
    int mx   = q.getsize(1);
    int my   = q.getsize(2);
    int meqn = q.getsize(3);
    int kmax = q.getsize(4);
    int mbc  = q.getmbc();
    int maux = aux.getsize(3);

double xi_sign[]  = {1.0, -1.0, 1.0, -1.0, 1.0, 1.0,1.0,-1.0,-1.0,1.0};
double eta_sign[] = {1.0, 1.0, -1.0, -1.0, 1.0, 1.0,-1.0,1.0,1.0,-1.0};

        
    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------
   
    //for (ell=1; ell<=kmax; ell++)
    for (ell=1; ell<=kmax; ell++)
    {
    
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (i=0; i>=(1-mbc); i--)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(1,j,m,ell);
                    q.set(i,j,m,ell, 0.0+(ell<=1)*tmp );
                    //q.set(i,j,m,ell, (xi_sign[ell-1])*tmp );
                    //q.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (i=(mx+1); i<=(mx+mbc); i++)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(mx,j,m,ell);                    
                    q.set(i,j,m,ell, 0.0+(ell<=1)*tmp );
                    //q.set(i,j,m,ell, (xi_sign[ell-1])*tmp );
                    //q.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // BOTTOM BOUNDARY
        // ***********************************************
        for (j=0; j>=(1-mbc); j--)
        {
            for (i=(1-mbc); i<=(mx+mbc); i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,1,m,ell);                    
                    q.set(i,j,m,ell, 0.0+(ell<=1)*tmp );
                    //q.set(i,j,m,ell, (eta_sign[ell-1])*tmp );
                    //q.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // TOP BOUNDARY
        // ***********************************************
        for (j=(my+1); j<=(my+mbc); j++)
        {
            for (i=(2-mbc); i<=(mx+mbc); i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,my,m,ell);                    
                    q.set(i,j,m,ell, 0.0+(ell<=1)*tmp );
                    //q.set(i,j,m,ell, (eta_sign[ell-1])*tmp );
                    //q.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
    }
    
}
