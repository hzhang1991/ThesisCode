#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"
#include "DogParamsCart2.h"

///////////////////////////////////////////////////////////////////////////////
//
//  Function to Convert 1d legendre weights into 2d legendre weights.  I.e.
//  This routine takes a representation q = q(x) and fills it out as a
//  representation q = q(x,y).  Every polynomial weight involving the second
//  variable is zero.
//
//  Parameters:
//
//     mopt == 1:  Convert q(x) -> q(x,y):  melems1d = mx.
//     mopt == 2:  Convert q(y) -> q(x,y):  melems1d = my.
//
//     q2d(mx,my,meqn,kmax, mbc)  - 2d Legendre weights  
//
//     q1d(melems1d, meqn, kmax1d, mbc ) - 1d Legendre weights 
//
//     q2d(mx, my, meqn, kmax, mbc ) - 2d Legendre weights after extension
//
///////////////////////////////////////////////////////////////////////////////
void ConvertQ1dToQ2d(const int &mopt, int istart, int iend, 
                     int jstart, int jend, 
                     const dTensorBC3& q1d, dTensorBC4& q2d)
{

    const int mx   = q2d.getsize(1);
    const int my   = q2d.getsize(2);
    const int meqn = q2d.getsize(3);  // assume meqn == 1
    const int kmax = q2d.getsize(4);
    const int mbc  = q2d.getmbc();
    
    const int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d = int(sqrt(mpoints));
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double xlength = dogParamsCart2.get_xlength();
    const double ylength = dogParamsCart2.get_ylength();

    for(int i=istart; i<= iend; i++)
    for(int j=jstart; j<= jend; j++)
    for(int me = 1; me <= meqn; me++)
    {
        // determine polynomial that determines current weight
        switch( kmax1d + (mopt-1)*100 )
        {
            // q = q(x) - extend to q = q(x,y)
            case 5:
            q2d.set(i, j, me, 15, 0.0 );
            q2d.set(i, j, me, 14, q1d.get(i, me, 5) );
            q2d.set(i, j, me, 13, 0.0 );
            q2d.set(i, j, me, 12, 0.0 );
            q2d.set(i, j, me, 11, 0.0 );

            case 4:
            q2d.set(i, j, me, 10, 0.0 );
            q2d.set(i, j, me, 9, q1d.get(i, me, 4)  );
            q2d.set(i, j, me, 8, 0.0 );
            q2d.set(i, j, me, 7, 0.0 );

            case 3:
            q2d.set(i, j, me, 6, 0.0 );
            q2d.set(i, j, me, 5, q1d.get(i, me, 3)  );
            q2d.set(i, j, me, 4, 0.0 );

            case 2:
            q2d.set(i, j, me, 3, 0.0 );
            q2d.set(i, j, me, 2, q1d.get(i, me, 2)  );

            case 1:
            q2d.set(i, j, me, 1, q1d.get(i, me, 1) );
            break;

            default:
            unsupported_value_error( kmax1d + (mopt-1)*100 );
            break;
        }
            //integration in x-direction:  
//            case 101: k2d = 1; break;
//            case 102: k2d = 3; break;
//            case 103: k2d = 6; break;
//            case 104: k2d = 10; break;
//            case 105: k2d = 15; break;
            
    }

}//end of function ConvertQ1dToQ2d 

void ConvertQ2dToQ1d(const int &mopt, int istart, int iend, 
                     const dTensorBC4& qin, dTensorBC3& qout)
{

    const int mx   = qin.getsize(1);
    const int my   = qin.getsize(2);
    const int meqn = qin.getsize(3);  // assume meqn == 1
    const int kmax = qin.getsize(4);
    const int mbc  = qin.getmbc();
    
    const int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d = int(sqrt(mpoints));
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double xlength = dogParamsCart2.get_xlength();
    const double ylength = dogParamsCart2.get_ylength();

    // qin is a function of x only, so only need to look at first row //
    const int j = 1;
    for(int i=istart; i<= iend; i++)
    for(int me = 1; me <= meqn; me++)
    {
        // determine polynomial that determines current weight
        switch( kmax1d + (mopt-1)*100 )
        {
            // q = q(x) - lower to q = q(x)
            case 5:
            qout.set(i, me, 5, qin.get(i, j, me, 14) );

            case 4:
            qout.set(i, me, 4, qin.get(i, j, me, 9) );

            case 3:
            qout.set(i, me, 3, qin.get(i, j, me, 5) );

            case 2:
            qout.set(i, me, 2, qin.get(i, j, me, 2) );

            case 1:
            qout.set(i, me, 1, qin.get(i, j, me, 1) );
            break;

            default:
            unsupported_value_error( kmax1d + (mopt-1)*100 );
            break;
        }
            //integration in x-direction:  
//            case 101: k2d = 1; break;
//            case 102: k2d = 3; break;
//            case 103: k2d = 6; break;
//            case 104: k2d = 10; break;
//            case 105: k2d = 15; break;
            
    }

}//end of function ConvertQ1dToQ2d 
