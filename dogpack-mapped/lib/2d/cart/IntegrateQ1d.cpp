#include "dogdefs.h"
#include <cmath>
#include "DogParamsCart2.h"

///////////////////////////////////////////////////////////////////////////////
//
//  Function to integrate Q in one variable only.
//
//  I.e. evaluate the integral \int q(x,y)\ dx OR \int q(x,y)\ dy.
//  depending on parameter mopt
//
//  Parameters:
//
//     mopt == 1:  integrate in y-direction:  melems1d = mx.
//     mopt == 2:  integrate in x-direction:  melems1d = my.
//
//     q2d(mx,my,meqn,kmax, mbc)  - 2d Legendre weights  
//
//     q1d(melems1d, meqn, kmax1d, mbc ) - 1d Legendre weights (after
//     integration)
//
///////////////////////////////////////////////////////////////////////////////
void IntegrateQ1d(const int mopt, 
		  const dTensorBC4& q2d, 
		  dTensorBC3& q1d)
{
  void IntegrateQ1d(const int mopt, 
		    const dTensorBC4* q2d, 
		    dTensorBC3* q1d);
  IntegrateQ1d(mopt,&q2d,&q1d);
}

void IntegrateQ1d(const int mopt, 
		  const dTensorBC4* q2d, 
		  dTensorBC3* q1d)
{
    const int mx   = q2d->getsize(1);
    const int my   = q2d->getsize(2);
    const int meqn = q2d->getsize(3);
    const int kmax = q2d->getsize(4);
    const int mbc  = q2d->getmbc();
    
    const int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d = int(sqrt(mpoints));
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    
    // quick error check
    assert_eq(mopt,1);
    assert_eq(kmax1d,q1d->getsize(3));
    assert_eq(mbc,q1d->getmbc());

    if( mopt == 1)
    {  // integrate in y-direction


        iTensor1 k2d(5);
        k2d.set(1,  1);
        k2d.set(2,  2);
        k2d.set(3,  5);
        k2d.set(4,  9);
        k2d.set(5, 14);

        for( int i = (1-mbc); i<= mx+mbc; i++ )
        for( int me = 1; me <= meqn; me++ )
        for( int k = 1; k <= kmax1d; k++ )
        {
            double tmp = 0.0;
            // \sum_j dy * Q^(k2d)_{i,j}
            for( int j = 1; j<= my; j++) 
            {
                tmp += dy * q2d->get(i,j,me, k2d.get(k) );    
            }

            // save integrated value into q
            q1d->set(i,me,k, tmp);

        } // end of integrating cell i


    }

    else if( mopt == 2 )
    { //integrate in x-direction

        iTensor1 k2d(5);
        k2d.set(1,  1);
        k2d.set(2,  3);
        k2d.set(3,  6);
        k2d.set(4, 10);
        k2d.set(5, 15);

        for( int j = (1-mbc); j<= my+mbc; j++ )
        for( int me = 1; me <= meqn; me++ )
        for( int k = 1; k <= kmax1d; k++ )
        {
            double tmp = 0.0;
            // \sum_j dy * Q^(k2d)_{i,j}
            for( int i = 1; i<= mx; i++) 
            {
                tmp += dx * q2d->get(i,j,me, k2d.get(k) );    
            }

            // save integrated value into q
            q1d->set(j,me,k, tmp);

        }

    }
    else
    {
        assert_eq(1,0);
    }

}//end of function IntegrateQ1d

///////////////////////////////////////////////////////////////////////////////
//
// Function to compute first moment of q.  
// Currently will only support computing
//
//     \int_y  y * q(x,y) dy
//
///////////////////////////////////////////////////////////////////////////////
void IntegrateQ1dMoment1(const dTensorBC4& q2d, 
			 dTensorBC3& q1d)
{
  void IntegrateQ1dMoment1(const dTensorBC4* q2d, 
			   dTensorBC3* q1d);
  IntegrateQ1dMoment1(&q2d,&q1d);
}

void IntegrateQ1dMoment1(const dTensorBC4* q2d, 
			 dTensorBC3* q1d)
{
    int mx   = q2d->getsize(1);
    int my   = q2d->getsize(2);
    int meqn = q2d->getsize(3);
    int kmax = q2d->getsize(4);
    int mbc  = q2d->getmbc();
    
    int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    int kmax1d = int(sqrt(mpoints));
    double dx = dogParamsCart2.get_dx();
    double dy = dogParamsCart2.get_dy();
    const double ylow = dogParamsCart2.get_ylow();
    
#pragma omp parallel for
    for( int i = (1-mbc); i<= (mx+mbc); i++ )
    for( int me=1; me <= meqn; me++ )
    {

        dTensor1 qi(kmax1d);
        qi.setall(0.);

        // \sum_j dy * Q^(k2d)_{i,j}
        for( int j = 1; j<= my; j++) 
        {

            // polynomial weights ...
            dTensor1 tmp(kmax1d);  

            const double yj = ylow + (double(j)-0.5)*dy;

            switch( kmax1d )
            {
                case 5:

                    tmp.set(5, yj * q2d->get(i,j,me,14) );

                    tmp.set(4, yj * q2d->get(i,j,me,9) +
                        dy*sq3/6*q2d->get(i,j,me,11) );

                    tmp.set(3, yj * q2d->get(i,j,me,5) + 
                        dy * sq3 / 6.0 * q2d->get(i,j,me,7) );

                    tmp.set(2, yj * q2d->get(i,j,me,2) + 
                        dy *sq3 / 6.0 * q2d->get(i,j,me,4) );

                    tmp.set(1, yj * q2d->get(i,j,me,1) + 
                        dy *sq3 / 6.0 * q2d->get(i,j,me,3) );
                    break;

                case 4:

                    tmp.set(4, yj * q2d->get(i,j,me,9) );

                    tmp.set(3, yj * q2d->get(i,j,me,5) + 
                        dy * sq3 / 6.0 * q2d->get(i,j,me,7) );

                    tmp.set(2, yj * q2d->get(i,j,me,2) + 
                        dy *sq3 / 6.0 * q2d->get(i,j,me,4) );

                    tmp.set(1, yj * q2d->get(i,j,me,1) + 
                        dy *sq3 / 6.0 * q2d->get(i,j,me,3) );
                    break;

                case 3:

                    tmp.set(3, yj * q2d->get(i,j,me,5) );

                    tmp.set(2, yj * q2d->get(i,j,me,2) + 
                        dy *sq3 / 6.0 * q2d->get(i,j,me,4) );

                    tmp.set(1, yj * q2d->get(i,j,me,1) + 
                        dy *sq3 / 6.0 * q2d->get(i,j,me,3) );
                    break;

                case 2:
                    tmp.set(2, yj * q2d->get(i,j,me,2) );

                    tmp.set(1, yj * q2d->get(i,j,me,1) + 
                        dy *sq3 / 6.0 * q2d->get(i,j,me,3) ) ;
                    break;

                case 1:
                    tmp.set(1, yj * q2d->get(i,j,me,1) );

                    break;

                default:
                    perror("  bad spatial order = chosen\n");
                    exit(1);

            }

            for( int know=1; know <= kmax1d; know++ )
            { qi.set(know, qi.get(know) + dy * tmp.get(know) ); }

        } // end of integrating cell i

        // save integrated value into q
        for( int k = 1; k <= kmax1d; k++ )
        { q1d->set(i,meqn,k, qi.get(k) ); }

    }

}

///////////////////////////////////////////////////////////////////////////////
//
// Function to compute second moment of q.  
// Currently will only support computing
//
//     1/2 \int y^2 * q(x,y) dy
//
///////////////////////////////////////////////////////////////////////////////
void IntegrateQ1dMoment2(const dTensorBC4& q2d, 
			 dTensorBC3& q1d)
{
  void IntegrateQ1dMoment2(const dTensorBC4* q2d, 
			   dTensorBC3* q1d);
  IntegrateQ1dMoment2(&q2d,&q1d);
}

void IntegrateQ1dMoment2(const dTensorBC4* q2d, 
			 dTensorBC3* q1d)
{

    int mx   = q2d->getsize(1);
    int my   = q2d->getsize(2);
    int meqn = q2d->getsize(3);
    int kmax = q2d->getsize(4);
    int mbc  = q2d->getmbc();
   // int maux = aux.getsize(3);
    
    int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    int kmax1d = int(sqrt(mpoints));
    double dx = dogParamsCart2.get_dx();
    double dy = dogParamsCart2.get_dy();
    const double ylow = dogParamsCart2.get_ylow();
    
#pragma omp parallel for
    for( int i = (1-mbc); i<= (mx+mbc); i++ )
    for( int me=1; me <= meqn; me++ )
    {

        dTensor1 qi(kmax1d);
        qi.setall(0.);

        // \sum_j dy * Q^(k2d)_{i,j}
        for( int j = 1; j<= my; j++) 
        {
            const double yj = ylow + (double(j)-0.5)*dy;

            // polynomial weights ...
            dTensor1 tmp(kmax1d);  

            switch( kmax1d )
            {


                case 5:

                    tmp.set(5, q2d->get(i,j,me,14) * (dy*dy/12.0+yj*yj) );

                    tmp.set(4, 
                        q2d->get(i,j,me,11) * (dy*yj*sq3/3.0) + 
                        q2d->get(i,j,me,9 ) * (dy*dy/12.0+pow(yj,2)) );

                    tmp.set(3, 
                        q2d->get(i,j,me,13) * sq5/30.0*pow(dy,2) + 
                        q2d->get(i,j,me,5) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,7) );

                    tmp.set(2, (
                        q2d->get(i,j,me,2) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,4) +
                        pow(dy,2)*sq5/30.0*q2d->get(i,j,me,8) ) );

                    tmp.set(1, (
                        q2d->get(i,j,me,1) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,3) +
                        pow(dy,2)*sq5/30.0*q2d->get(i,j,me,6) ) );
                    break;
 
                case 4:

                    tmp.set(4, q2d->get(i,j,me,9) * (dy*dy/12.0+pow(yj,2)) );

                    tmp.set(3, (q2d->get(i,j,me,5) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,7) ) );

                    tmp.set(2, (
                        q2d->get(i,j,me,2) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,4) +
                        pow(dy,2)*sq5/30.0*q2d->get(i,j,me,8) ) );

                    tmp.set(1, (
                        q2d->get(i,j,me,1) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,3) +
                        pow(dy,2)*sq5/30.0*q2d->get(i,j,me,6) ) );
                    
                    break;

                case 3:
                    tmp.set(3, q2d->get(i,j,me,5) * (dy*dy/12.0+pow(yj,2)) );

                    tmp.set(2, (
                        q2d->get(i,j,me,2) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,4) ) );

                    tmp.set(1, (
                        q2d->get(i,j,me,1) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,3) +
                        pow(dy,2)*sq5/30.0*q2d->get(i,j,me,6) )
                        );
                    break;

                case 2:

                    tmp.set(2, q2d->get(i,j,me,2) * (dy*dy/12.0+pow(yj,2)) );

                    tmp.set(1, (q2d->get(i,j,me,1) * (dy*dy/12.0+pow(yj,2)) +
                        yj*dy/sq3*q2d->get(i,j,me,3) ) );
                    break;
                case 1:

                    tmp.set(1, q2d->get(i,j,me,1) * (dy*dy/12.0+pow(yj,2)) );
                    break;
               
                default:
                    perror(" not implemented yet\n");

            }

            // note: factor of 1/2 shows up here ...
            for( int know=1; know <= kmax1d; know++ )
            { qi.set(know, qi.get(know) + 0.5 * dy * tmp.get(know) ); }

        } // end of integrating cell i

        // save integrated value into q
        for( int k = 1; k <= kmax1d; k++ )
        { q1d->set(i,meqn,k, qi.get(k) ); }

    }

}

///////////////////////////////////////////////////////////////////////////////
//
// Function to compute third moment of q.  
// Currently will only support computing
//
//     \int y^3 * q(x,y) dy
//
///////////////////////////////////////////////////////////////////////////////
void IntegrateQ1dMoment3(const dTensorBC4& q2d, 
			 dTensorBC3& q1d)
{
  void IntegrateQ1dMoment3(const dTensorBC4* q2d, 
			   dTensorBC3* q1d);
  IntegrateQ1dMoment3(&q2d,&q1d);
}

void IntegrateQ1dMoment3(const dTensorBC4* q2d, 
			 dTensorBC3* q1d)
{

    const int mx   = q2d->getsize(1);
    const int my   = q2d->getsize(2);
    const int meqn = q2d->getsize(3);
    const int kmax = q2d->getsize(4);
    const int mbc  = q2d->getmbc();
    
    const int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d = int(sqrt(mpoints));
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double ylow = dogParamsCart2.get_ylow();
   
    const double dy2 = dy*dy;
    const double dy3 = dy*dy2;

#pragma omp parallel for
    for( int i = (1-mbc); i<= (mx+mbc); i++ )
    for( int me=1; me <= meqn; me++ )
    {

        dTensor1 qi(kmax1d);
        qi.setall(0.);

        // \sum_j dy * Q^(k2d)_{i,j}
        for( int j = 1; j<= my; j++) 
        {
            const double yj  = ylow + (double(j)-0.5)*dy;
            const double yj2 = yj*yj;
            const double yj3 = yj*yj2;

            // polynomial weights ...
            dTensor1 tmp(kmax1d);  

            switch( kmax1d )
            {

                case 5:

                    tmp.set(5,
                        q2d->get(i,j,me,14) * (
                            yj*dy2/4.0 + yj3 )
                    );

                    tmp.set(4, 
                        q2d->get(i,j,me,9) * (
                            yj*(dy2+4.0*yj2)/4.0)
                        + q2d->get(i,j,me,11) * (
                            dy*sq3/40.0*(dy2+20.0*yj2) )
                    );

                    tmp.set(3, 
                        q2d->get(i,j,me,5) * (
                            yj*(dy2+4.0*yj2)/4.0)
                        + q2d->get(i,j,me,7) * (
                            dy*sq3/40.0*(dy2+20.0*yj2) )
                        + q2d->get(i,j,me,13) * (
                            yj*dy2*sq5/10.0 )
                    );

                    tmp.set(2, 
                        q2d->get(i,j,me,2) * (
                            yj*(dy2+4.0*yj2)/4.0)
                        + q2d->get(i,j,me,4) * (
                            dy*sq3/40.0*(dy2+20.0*yj2) )
                        + q2d->get(i,j,me,8) * (
                            yj*dy2*sq5/10.0 )
                        + q2d->get(i,j,me,12) * (
                            sq7*dy3/140.0 )
                    );

                    tmp.set(1, q2d->get(i,j,me,1) * (
                                0.25*dy*dy*yj + pow(yj,3)) 
                            + q2d->get(i,j,me,3) * (
                                dy3/40.0*sq3 + 0.5*yj2*dy*sq3 )
                            + q2d->get(i,j,me,6) * (
                                yj*dy2*sq5/10.0 )
                            + q2d->get(i,j,me,10) * (
                                sq7/140.0*dy3 )
                           );


                    break;

                case 4:

                    tmp.set(4, 
                        q2d->get(i,j,me,9) * (
                            yj*(dy2+4.0*yj2)/4.0)
                    );

                    tmp.set(3, 
                        q2d->get(i,j,me,5) * (
                            yj*(dy2+4.0*yj2)/4.0)
                        + q2d->get(i,j,me,7) * (
                            dy*sq3/40.0*(dy2+20.0*yj2) )
                    );

                    tmp.set(2, 
                        q2d->get(i,j,me,2) * (
                            yj*(dy2+4.0*yj2)/4.0)
                        + q2d->get(i,j,me,4) * (
                            dy*sq3/40.0*(dy2+20.0*yj2) )
                        + q2d->get(i,j,me,8) * (
                            yj*dy2*sq5/10.0 )
                    );

                    tmp.set(1, q2d->get(i,j,me,1) * (
                                0.25*dy*dy*yj + pow(yj,3)) 
                            + q2d->get(i,j,me,3) * (
                                dy3/40.0*sq3 + 0.5*yj2*dy*sq3 )
                            + q2d->get(i,j,me,6) * (
                                yj*dy2*sq5/10.0 )
                            + q2d->get(i,j,me,10) * (
                                sq7/140.0*dy3 )
                           );

                    break;

                case 3:

                    tmp.set(3, 
                        q2d->get(i,j,me,5) * (
                            yj*(dy2+4.0*yj2)/4.0)
                    );

                    tmp.set(2, 
                        q2d->get(i,j,me,2) * (
                            yj*(dy2+4.0*yj2)/4.0)
                        + q2d->get(i,j,me,4) * (
                            dy*sq3/40.0*(dy2+20.0*yj2) )
                    );

                    tmp.set(1, q2d->get(i,j,me,1) * (
                        0.25*dy*dy*yj + pow(yj,3)) 
                        + q2d->get(i,j,me,3) * (
                            dy3/40.0*sq3 + 0.5*yj2*dy*sq3 )
                        + q2d->get(i,j,me,6) * (
                            yj*dy2*sq5/10.0 )
                    );

                    break;

                case 2:

                    tmp.set(2, 
                        q2d->get(i,j,me,2) * (
                            yj*(dy2+4.0*yj2)/4.0)
                    );

                    tmp.set(1, 
                        q2d->get(i,j,me,1) * (
                            0.25*dy*dy*yj + pow(yj,3) )
                        + q2d->get(i,j,me,3) * (
                            dy3/40.0*sq3 + 0.5*yj2*dy*sq3 )
                    );
                    break;

                case 1:

                    tmp.set(1, q2d->get(i,j,me,1) * (
                        0.25*dy*dy*yj + pow(yj,3)) 
                    );
                    break;
               
                default:
                    perror(" not implemented yet\n");

            }

            for( int know=1; know <= kmax1d; know++ )
            { qi.set(know, qi.get(know) + dy * tmp.get(know) ); }

        } // end of integrating cell i

        // save integrated value into q
        for( int k = 1; k <= kmax1d; k++ )
        { q1d->set(i,meqn,k, qi.get(k) ); }

    }

}
