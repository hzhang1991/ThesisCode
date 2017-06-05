#include <cmath>
#include "constants.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "tensors.h"
#include "assert.h"
using namespace std;

//<<<<<<< .working
//void ArtificialViscosity(const dTensorBC4& aux, 
//        const dTensorBC4& q, dTensorBC4& Lstar)
//=======
void ArtificialViscosity(const dTensorBC4& aux,
			 const dTensorBC4& q, 
			 dTensorBC4& Lstar)
//>>>>>>> .merge-right.r1600
{
  void ArtificialViscosity(const dTensorBC4* aux,
			   const dTensorBC4* q, 
			   dTensorBC4* Lstar);
  ArtificialViscosity(&aux,&q,&Lstar);
}


void ArtificialViscosity(const dTensorBC4* aux,
			 const dTensorBC4* q, 
			 dTensorBC4* Lstar)
{
    int     mx = q->getsize(1);
    int     my = q->getsize(2);
    int   meqn = q->getsize(3);
    int   kmax = q->getsize(4);
    int    mbc = q->getmbc();
    double dx  = dogParamsCart2.get_dx();
    double dy  = dogParamsCart2.get_dx();
    int morder = dogParams.get_space_order();
    int kless = int((morder*(morder-1))/2);
    double dxymax = Max(dx,dy);
//<<<<<<< .working
//    void CreateLdiff(const dTensorBC2& eps, const dTensorBC4& q, 
//            dTensorBC4& Lstar);
//=======
    void CreateLdiff(const dTensorBC2& eps, const dTensorBC4* q, 
		     dTensorBC4* Lstar);
//>>>>>>> .merge-right.r1600

    assert_printf(meqn<=1, "Currently only supported for scalar equations: "
            "meqn = %d", meqn);

    // Indicator for artificial viscosity limiter
    int m = 1;
    double s0   = 0.00005e0/pow(double(morder),4);
    double eps0 = 0.001;
    dTensorBC2 eps(mx,my,mbc);
    int mz = 0;

    for (int i=(1-mbc); i<=(mx+mbc); i++)
        for (int j=(1-mbc); j<=(my+mbc); j++)
        {	
            double Q2sum = 0.0;

            for (int k=1; k<=kless; k++)
            {  Q2sum = Q2sum + pow(q->get(i,j,m,k),2);  }

            double Q2kmax = 0.0;

            for (int k=(kless+1); k<=kmax; k++)
            {  Q2kmax = Q2kmax + pow(q->get(i,j,m,k),2);  }

            Q2sum = Q2sum + Q2kmax;

            double se;
            if (Q2kmax>1.0e-10)
            {  se =  (Q2kmax/Q2sum);  }
            else
            {  se = 0.0e0;  }

            if (se>=s0)
            {  
                mz = mz + 1;

                // Set the artificial diffusion parameter
                // in this element
                eps.set(i,j, eps0 * pow(dxymax, morder-1) );
            }
            else
            {  eps.set(i,j, 0.0 );  }
        }

    if (mz>0)
      {  CreateLdiff(eps,q,Lstar);  }

}


// Compute explicit time-stepping diffusion operator
//<<<<<<< .working
//void CreateLdiff(const dTensorBC2& eps, const dTensorBC4& q, 
//        dTensorBC4& Lstar)
//=======
void CreateLdiff(const dTensorBC2& eps, const dTensorBC4* q, 
		 dTensorBC4* Lstar)
//>>>>>>> .merge-right.r1600
{
    int mx   = q->getsize(1);
    int my   = q->getsize(2);
    int meqn = q->getsize(3);
    int kmax = q->getsize(4);
    int mbc  = q->getmbc();
    dTensorBC4 V(mx,my,meqn,kmax,mbc);
    dTensorBC4 W(mx,my,meqn,kmax,mbc);

    // Compute gradient
    void ComputeGradient_AV(const dTensorBC4* q, dTensorBC4& V, 
			    dTensorBC4& W);
    ComputeGradient_AV(q,V,W);

    // Add divergence of gradient to Lstar
    void ComputeRHS_AV(const dTensorBC2& eps, const dTensorBC4& V, 
		       const dTensorBC4& W, dTensorBC4* Lstar);
    ComputeRHS_AV(eps,V,W,Lstar);
}


// -------------------------------------------------------------
// Function to compute gradient
// -------------------------------------------------------------
void ComputeGradient_AV(const dTensorBC4* q, dTensorBC4& V, 
        dTensorBC4& W)
{
    int mx   = q->getsize(1);
    int my   = q->getsize(2);
    int meqn = q->getsize(3);
    int kmax = q->getsize(4);
    int mbc  = q->getmbc();
    double dx = dogParamsCart2.get_dx();
    double dy = dogParamsCart2.get_dy();
    int m = 1;

    for (int i=(1-mbc); i<=(mx+mbc-1); i++)
        for (int j=(1-mbc); j<=(my+mbc-1); j++)
        {
            dTensor1 SU(kmax);
            dTensor1 TU(kmax);
            dTensor1 FUV_p(kmax);
            dTensor1 FUV_m(kmax);
            dTensor1 FUW_p(kmax);
            dTensor1 FUW_m(kmax);

            // Various contributions from integration-by-parts
            switch( dogParams.get_space_order() )
            {
                case 5: // fifth order

                    // Integral over element for x-derivative
                    SU.set(15, 0.0 );
                    SU.set(14, 6.0*(sq3*q->get(i,j,m,2)+sq7*q->get(i,j,m,9)) );
                    SU.set(13, 2.0*sq3*sq5*q->get(i,j,m,8) );
                    SU.set(12, 2.0*sq3*q->get(i,j,m,10) );
                    SU.set(11, 2.0*sq7*(q->get(i,j,m,3)+sq5*q->get(i,j,m,7)) );
                    SU.set(10, 0.0 );
                    SU.set(9,  2.0*sq7*(q->get(i,j,m,1)+sq5*q->get(i,j,m,5)) );
                    SU.set(8,  2.0*sq3*q->get(i,j,m,6) );
                    SU.set(7,  2.0*sq3*sq5*q->get(i,j,m,4) );
                    SU.set(6,  0.0 );
                    SU.set(5,  2.0*sq3*sq5*q->get(i,j,m,2) );
                    SU.set(4,  2.0*sq3*q->get(i,j,m,3) );
                    SU.set(3,  0.0 );
                    SU.set(2,  2.0*sq3*q->get(i,j,m,1) );
                    SU.set(1,  0.0 );

                    // Boundary contribution for x-derivative (right side)
                    FUV_p.set(15, q->get(i+1,j,m,15) );
                    FUV_p.set(14, 3.0*q->get(i+1,j,m,1) - 3.0*sq3*q->get(i+1,j,m,2) + 3.0*sq5*q->get(i+1,j,m,5) 
                            - 3.0*sq7*q->get(i+1,j,m,9) + 9.0*q->get(i+1,j,m,14) );
                    FUV_p.set(13, sq5*q->get(i+1,j,m,6) - sq3*sq5*q->get(i+1,j,m,8) + 5.0*q->get(i+1,j,m,13) );
                    FUV_p.set(12, sq3*q->get(i+1,j,m,10) - 3.0*q->get(i+1,j,m,12) );
                    FUV_p.set(11, sq7*q->get(i+1,j,m,3) - sq3*sq7*q->get(i+1,j,m,4) + sq5*sq7*q->get(i+1,j,m,7) 
                            - 7.0*q->get(i+1,j,m,11) );
                    FUV_p.set(10, q->get(i+1,j,m,10) - sq3*q->get(i+1,j,m,12) );
                    FUV_p.set(9,  sq7*q->get(i+1,j,m,1) - sq3*sq7*q->get(i+1,j,m,2) + sq5*sq7*q->get(i+1,j,m,5) 
                            - 7.0*q->get(i+1,j,m,9) + 3.0*sq7*q->get(i+1,j,m,14) );
                    FUV_p.set(8,  sq3*q->get(i+1,j,m,6) - 3.0*q->get(i+1,j,m,8) + sq3*sq5*q->get(i+1,j,m,13) );
                    FUV_p.set(7,  sq5*q->get(i+1,j,m,3) - sq3*sq5*q->get(i+1,j,m,4) + 5.0*q->get(i+1,j,m,7) 
                            - sq5*sq7*q->get(i+1,j,m,11) );
                    FUV_p.set(6,  q->get(i+1,j,m,6) - sq3*q->get(i+1,j,m,8) + sq5*q->get(i+1,j,m,13) );
                    FUV_p.set(5,  sq5*q->get(i+1,j,m,1) - sq3*sq5*q->get(i+1,j,m,2) + 5.0*q->get(i+1,j,m,5) 
                            - sq5*sq7*q->get(i+1,j,m,9) + 3.0*sq5*q->get(i+1,j,m,14) );
                    FUV_p.set(4,  sq3*q->get(i+1,j,m,3) - 3.0*q->get(i+1,j,m,4) + sq3*sq5*q->get(i+1,j,m,7) 
                            - sq3*sq7*q->get(i+1,j,m,11) );
                    FUV_p.set(3,  q->get(i+1,j,m,3) - sq3*q->get(i+1,j,m,4) + sq5*q->get(i+1,j,m,7)
                            - sq7*q->get(i+1,j,m,11) );
                    FUV_p.set(2,  sq3*q->get(i+1,j,m,1) - 3.0*q->get(i+1,j,m,2) + sq3*sq5*q->get(i+1,j,m,5)
                            - sq3*sq7*q->get(i+1,j,m,9) + 3.0*sq3*q->get(i+1,j,m,14) );
                    FUV_p.set(1,  q->get(i+1,j,m,1) - sq3*q->get(i+1,j,m,2) + sq5*q->get(i+1,j,m,5)
                            - sq7*q->get(i+1,j,m,9) + 3.0*q->get(i+1,j,m,14) );

                    // Boundary contribution for x-derivative (left side)
                    FUV_m.set(15,  q->get(i,j,m,15) );
                    FUV_m.set(14,  3.0*q->get(i,j,m,1) - 3.0*sq3*q->get(i,j,m,2) + 3.0*sq5*q->get(i,j,m,5)
                            - 3.0*sq7*q->get(i,j,m,9) + 9.0*q->get(i,j,m,14) );
                    FUV_m.set(13,  sq5*q->get(i,j,m,6) - sq3*sq5*q->get(i,j,m,8) + 5.0*q->get(i,j,m,13) );
                    FUV_m.set(12, -sq3*q->get(i,j,m,10) + 3.0*q->get(i,j,m,12) );
                    FUV_m.set(11, -sq7*q->get(i,j,m,3) + sq3*sq7*q->get(i,j,m,4) - sq5*sq7*q->get(i,j,m,7) 
                            + 7.0*q->get(i,j,m,11) );
                    FUV_m.set(10,  q->get(i,j,m,10) - sq3*q->get(i,j,m,12) );
                    FUV_m.set(9,  -sq7*q->get(i,j,m,1) + sq3*sq7*q->get(i,j,m,2) - sq5*sq7*q->get(i,j,m,5) 
                            + 7.0*q->get(i,j,m,9) - 3.0*sq7*q->get(i,j,m,14) );
                    FUV_m.set(8,  -sq3*q->get(i,j,m,6) + 3.0*q->get(i,j,m,8) - sq3*sq5*q->get(i,j,m,13) );
                    FUV_m.set(7,   sq5*q->get(i,j,m,3) - sq3*sq5*q->get(i,j,m,4) + 5.0*q->get(i,j,m,7) 
                            - sq5*sq7*q->get(i,j,m,11) );
                    FUV_m.set(6,   q->get(i,j,m,6) - sq3*q->get(i,j,m,8) + sq5*q->get(i,j,m,13) );
                    FUV_m.set(5,   sq5*q->get(i,j,m,1) - sq3*sq5*q->get(i,j,m,2) + 5.0*q->get(i,j,m,5) 
                            - sq5*sq7*q->get(i,j,m,9) + 3.0*sq5*q->get(i,j,m,14) );
                    FUV_m.set(4,  -sq3*q->get(i,j,m,3) + 3.0*q->get(i,j,m,4) - sq3*sq5*q->get(i,j,m,7) 
                            + sq3*sq7*q->get(i,j,m,11) );
                    FUV_m.set(3,   q->get(i,j,m,3) - sq3*q->get(i,j,m,4) + sq5*q->get(i,j,m,7)
                            - sq7*q->get(i,j,m,11) );
                    FUV_m.set(2,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,2) - sq3*sq5*q->get(i,j,m,5)
                            + sq3*sq7*q->get(i,j,m,9) - 3.0*sq3*q->get(i,j,m,14) );
                    FUV_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,2) + sq5*q->get(i,j,m,5)
                            - sq7*q->get(i,j,m,9) + 3.0*q->get(i,j,m,14) );

                    // Integral over element for y-derivative
                    TU.set(15, 6.0*(sq3*q->get(i,j,m,3)+sq7*q->get(i,j,m,10)) );
                    TU.set(14, 0.0 );
                    TU.set(13, 2.0*sq3*sq5*q->get(i,j,m,7) );
                    TU.set(12, 2.0*sq7*(q->get(i,j,m,2)+sq5*q->get(i,j,m,8)) );
                    TU.set(11, 2.0*sq3*q->get(i,j,m,9) );
                    TU.set(10, 2.0*sq7*(q->get(i,j,m,1)+sq5*q->get(i,j,m,6)) );
                    TU.set(9,  0.0 );
                    TU.set(8,  2.0*sq3*sq5*q->get(i,j,m,4) );
                    TU.set(7,  2.0*sq3*q->get(i,j,m,5) );
                    TU.set(6,  2.0*sq3*sq5*q->get(i,j,m,3) );
                    TU.set(5,  0.0 );
                    TU.set(4,  2.0*sq3*q->get(i,j,m,2) );	    
                    TU.set(3,  2.0*sq3*q->get(i,j,m,1) );
                    TU.set(2,  0.0 );	    
                    TU.set(1,  0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FUW_p.set(15,  3.0*q->get(i,j+1,m,1) - 3.0*sq3*q->get(i,j+1,m,3) + 3.0*sq5*q->get(i,j+1,m,6)
                            - 3.0*sq7*q->get(i,j+1,m,10) + 9.0*q->get(i,j+1,m,15) );
                    FUW_p.set(14,  q->get(i,j+1,m,14) );
                    FUW_p.set(13,  sq5*q->get(i,j+1,m,5) - sq3*sq5*q->get(i,j+1,m,7) + 5.0*q->get(i,j+1,m,13) );
                    FUW_p.set(12,  sq7*q->get(i,j+1,m,2) - sq3*sq7*q->get(i,j+1,m,4) + sq5*sq7*q->get(i,j+1,m,8) 
                            - 7.0*q->get(i,j+1,m,12) );
                    FUW_p.set(11,  sq3*q->get(i,j+1,m,9) - 3.0*q->get(i,j+1,m,11) );
                    FUW_p.set(10,  sq7*q->get(i,j+1,m,1) - sq3*sq7*q->get(i,j+1,m,3) + sq5*sq7*q->get(i,j+1,m,6)
                            - 7.0*q->get(i,j+1,m,10) + 3.0*sq7*q->get(i,j+1,m,15) );
                    FUW_p.set(9,   q->get(i,j+1,m,9) - sq3*q->get(i,j+1,m,11) );
                    FUW_p.set(8,   sq5*q->get(i,j+1,m,2) - sq3*sq5*q->get(i,j+1,m,4) + 5.0*q->get(i,j+1,m,8)
                            - sq5*sq7*q->get(i,j+1,m,12) );
                    FUW_p.set(7,   sq3*q->get(i,j+1,m,5) - 3.0*q->get(i,j+1,m,7) + sq3*sq5*q->get(i,j+1,m,13) );
                    FUW_p.set(6,   sq5*q->get(i,j+1,m,1) - sq3*sq5*q->get(i,j+1,m,3) + 5.0*q->get(i,j+1,m,6)
                            - sq5*sq7*q->get(i,j+1,m,10) + 3.0*sq5*q->get(i,j+1,m,15) );
                    FUW_p.set(5,   q->get(i,j+1,m,5) - sq3*q->get(i,j+1,m,7) + sq5*q->get(i,j+1,m,13) );
                    FUW_p.set(4,   sq3*q->get(i,j+1,m,2) - 3.0*q->get(i,j+1,m,4) + sq3*sq5*q->get(i,j+1,m,8)
                            - sq3*sq7*q->get(i,j+1,m,12) );
                    FUW_p.set(3,   sq3*q->get(i,j+1,m,1) - 3.0*q->get(i,j+1,m,3) + sq3*sq5*q->get(i,j+1,m,6)
                            - sq3*sq7*q->get(i,j+1,m,10) + 3.0*sq3*q->get(i,j+1,m,15) );
                    FUW_p.set(2,   q->get(i,j+1,m,2) - sq3*q->get(i,j+1,m,4) + sq5*q->get(i,j+1,m,8)
                            - sq7*q->get(i,j+1,m,12) );
                    FUW_p.set(1,   q->get(i,j+1,m,1) - sq3*q->get(i,j+1,m,3) + sq5*q->get(i,j+1,m,6)
                            - sq7*q->get(i,j+1,m,10) + 3.0*q->get(i,j+1,m,15) );

                    // Boundary contribution for y-derivative (bottom side)
                    FUW_m.set(15,  3.0*q->get(i,j,m,1) - 3.0*sq3*q->get(i,j,m,3) + 3.0*sq5*q->get(i,j,m,6)
                            - 3.0*sq7*q->get(i,j,m,10) + 9.0*q->get(i,j,m,15) );
                    FUW_m.set(14,  q->get(i,j,m,14) );
                    FUW_m.set(13,  sq5*q->get(i,j,m,5) - sq3*sq5*q->get(i,j,m,7) + 5.0*q->get(i,j,m,13) );
                    FUW_m.set(12, -sq7*q->get(i,j,m,2) + sq3*sq7*q->get(i,j,m,4) - sq5*sq7*q->get(i,j,m,8)
                            + 7.0*q->get(i,j,m,12) );
                    FUW_m.set(11, -sq3*q->get(i,j,m,9) + 3.0*q->get(i,j,m,11) );
                    FUW_m.set(10, -sq7*q->get(i,j,m,1) + sq3*sq7*q->get(i,j,m,3) - sq5*sq7*q->get(i,j,m,6)
                            + 7.0*q->get(i,j,m,10) - 3.0*sq7*q->get(i,j,m,15) );
                    FUW_m.set(9,   q->get(i,j,m,9) - sq3*q->get(i,j,m,11) );
                    FUW_m.set(8,   sq5*q->get(i,j,m,2) - sq3*sq5*q->get(i,j,m,4) + 5.0*q->get(i,j,m,8)
                            - sq5*sq7*q->get(i,j,m,12) );
                    FUW_m.set(7,  -sq3*q->get(i,j,m,5) + 3.0*q->get(i,j,m,7) - sq3*sq5*q->get(i,j,m,13) );
                    FUW_m.set(6,   sq5*q->get(i,j,m,1) - sq3*sq5*q->get(i,j,m,3) + 5.0*q->get(i,j,m,6)
                            - sq5*sq7*q->get(i,j,m,10) + 3.0*sq5*q->get(i,j,m,15) );
                    FUW_m.set(5,   q->get(i,j,m,5) - sq3*q->get(i,j,m,7) + sq5*q->get(i,j,m,13) );
                    FUW_m.set(4,  -sq3*q->get(i,j,m,2) + 3.0*q->get(i,j,m,4) - sq3*sq5*q->get(i,j,m,8)
                            + sq3*sq7*q->get(i,j,m,12) );
                    FUW_m.set(3,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,3) - sq3*sq5*q->get(i,j,m,6)
                            + sq3*sq7*q->get(i,j,m,10) - 3.0*sq3*q->get(i,j,m,15) );
                    FUW_m.set(2,   q->get(i,j,m,2) - sq3*q->get(i,j,m,4) + sq5*q->get(i,j,m,8)
                            - sq7*q->get(i,j,m,12) );
                    FUW_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,3) + sq5*q->get(i,j,m,6)
                            - sq7*q->get(i,j,m,10) + 3.0*q->get(i,j,m,15) );

                    break;

                case 4: // fourth order

                    // Integral over element for x-derivative
                    SU.set(10, 0.0 );
                    SU.set(9,  2.0*sq7*(q->get(i,j,m,1)+sq5*q->get(i,j,m,5)) );
                    SU.set(8,  2.0*sq3*q->get(i,j,m,6) );
                    SU.set(7,  2.0*sq3*sq5*q->get(i,j,m,4) );
                    SU.set(6,  0.0 );
                    SU.set(5,  2.0*sq3*sq5*q->get(i,j,m,2) );
                    SU.set(4,  2.0*sq3*q->get(i,j,m,3) );
                    SU.set(3,  0.0 );
                    SU.set(2,  2.0*sq3*q->get(i,j,m,1) );
                    SU.set(1,  0.0 );

                    // Boundary contribution for x-derivative (right side)
                    FUV_p.set(10, q->get(i+1,j,m,10) );
                    FUV_p.set(9,  sq7*q->get(i+1,j,m,1) - sq3*sq7*q->get(i+1,j,m,2) + sq5*sq7*q->get(i+1,j,m,5) 
                            - 7.0*q->get(i+1,j,m,9) );
                    FUV_p.set(8,  sq3*q->get(i+1,j,m,6) - 3.0*q->get(i+1,j,m,8) );
                    FUV_p.set(7,  sq5*q->get(i+1,j,m,3) - sq3*sq5*q->get(i+1,j,m,4) + 5.0*q->get(i+1,j,m,7) );
                    FUV_p.set(6,  q->get(i+1,j,m,6) - sq3*q->get(i+1,j,m,8) );
                    FUV_p.set(5,  sq5*q->get(i+1,j,m,1) - sq3*sq5*q->get(i+1,j,m,2) + 5.0*q->get(i+1,j,m,5) 
                            - sq5*sq7*q->get(i+1,j,m,9) );
                    FUV_p.set(4,  sq3*q->get(i+1,j,m,3) - 3.0*q->get(i+1,j,m,4) + sq3*sq5*q->get(i+1,j,m,7) );
                    FUV_p.set(3,  q->get(i+1,j,m,3) - sq3*q->get(i+1,j,m,4) + sq5*q->get(i+1,j,m,7) );
                    FUV_p.set(2,  sq3*q->get(i+1,j,m,1) - 3.0*q->get(i+1,j,m,2) + sq3*sq5*q->get(i+1,j,m,5)
                            - sq3*sq7*q->get(i+1,j,m,9) );
                    FUV_p.set(1,  q->get(i+1,j,m,1) - sq3*q->get(i+1,j,m,2) + sq5*q->get(i+1,j,m,5)
                            - sq7*q->get(i+1,j,m,9) );

                    // Boundary contribution for x-derivative (left side)
                    FUV_m.set(10,  q->get(i,j,m,10) );
                    FUV_m.set(9,  -sq7*q->get(i,j,m,1) + sq3*sq7*q->get(i,j,m,2) - sq5*sq7*q->get(i,j,m,5) 
                            + 7.0*q->get(i,j,m,9) );
                    FUV_m.set(8,  -sq3*q->get(i,j,m,6) + 3.0*q->get(i,j,m,8) );
                    FUV_m.set(7,   sq5*q->get(i,j,m,3) - sq3*sq5*q->get(i,j,m,4) + 5.0*q->get(i,j,m,7) );
                    FUV_m.set(6,   q->get(i,j,m,6) - sq3*q->get(i,j,m,8) );
                    FUV_m.set(5,   sq5*q->get(i,j,m,1) - sq3*sq5*q->get(i,j,m,2) + 5.0*q->get(i,j,m,5) 
                            - sq5*sq7*q->get(i,j,m,9) );
                    FUV_m.set(4,  -sq3*q->get(i,j,m,3) + 3.0*q->get(i,j,m,4) - sq3*sq5*q->get(i,j,m,7) );
                    FUV_m.set(3,   q->get(i,j,m,3) - sq3*q->get(i,j,m,4) + sq5*q->get(i,j,m,7) );
                    FUV_m.set(2,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,2) - sq3*sq5*q->get(i,j,m,5)
                            + sq3*sq7*q->get(i,j,m,9) );
                    FUV_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,2) + sq5*q->get(i,j,m,5)
                            - sq7*q->get(i,j,m,9) );

                    // Integral over element for y-derivative
                    TU.set(10, 2.0*sq7*(q->get(i,j,m,1)+sq5*q->get(i,j,m,6)) );
                    TU.set(9,  0.0 );
                    TU.set(8,  2.0*sq3*sq5*q->get(i,j,m,4) );
                    TU.set(7,  2.0*sq3*q->get(i,j,m,5) );
                    TU.set(6,  2.0*sq3*sq5*q->get(i,j,m,3) );
                    TU.set(5,  0.0 );
                    TU.set(4,  2.0*sq3*q->get(i,j,m,2) );	    
                    TU.set(3,  2.0*sq3*q->get(i,j,m,1) );
                    TU.set(2,  0.0 );	    
                    TU.set(1,  0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FUW_p.set(10,  sq7*q->get(i,j+1,m,1) - sq3*sq7*q->get(i,j+1,m,3) + sq5*sq7*q->get(i,j+1,m,6)
                            - 7.0*q->get(i,j+1,m,10) );
                    FUW_p.set(9,   q->get(i,j+1,m,9) );
                    FUW_p.set(8,   sq5*q->get(i,j+1,m,2) - sq3*sq5*q->get(i,j+1,m,4) + 5.0*q->get(i,j+1,m,8) );
                    FUW_p.set(7,   sq3*q->get(i,j+1,m,5) - 3.0*q->get(i,j+1,m,7) );
                    FUW_p.set(6,   sq5*q->get(i,j+1,m,1) - sq3*sq5*q->get(i,j+1,m,3) + 5.0*q->get(i,j+1,m,6)
                            - sq5*sq7*q->get(i,j+1,m,10) );
                    FUW_p.set(5,   q->get(i,j+1,m,5) - sq3*q->get(i,j+1,m,7) );
                    FUW_p.set(4,   sq3*q->get(i,j+1,m,2) - 3.0*q->get(i,j+1,m,4) + sq3*sq5*q->get(i,j+1,m,8) );
                    FUW_p.set(3,   sq3*q->get(i,j+1,m,1) - 3.0*q->get(i,j+1,m,3) + sq3*sq5*q->get(i,j+1,m,6)
                            - sq3*sq7*q->get(i,j+1,m,10) );
                    FUW_p.set(2,   q->get(i,j+1,m,2) - sq3*q->get(i,j+1,m,4) + sq5*q->get(i,j+1,m,8) );
                    FUW_p.set(1,   q->get(i,j+1,m,1) - sq3*q->get(i,j+1,m,3) + sq5*q->get(i,j+1,m,6)
                            - sq7*q->get(i,j+1,m,10) );

                    // Boundary contribution for y-derivative (bottom side)
                    FUW_m.set(10, -sq7*q->get(i,j,m,1) + sq3*sq7*q->get(i,j,m,3) - sq5*sq7*q->get(i,j,m,6)
                            + 7.0*q->get(i,j,m,10) );
                    FUW_m.set(9,   q->get(i,j,m,9) );
                    FUW_m.set(8,   sq5*q->get(i,j,m,2) - sq3*sq5*q->get(i,j,m,4) + 5.0*q->get(i,j,m,8) );
                    FUW_m.set(7,  -sq3*q->get(i,j,m,5) + 3.0*q->get(i,j,m,7) );
                    FUW_m.set(6,   sq5*q->get(i,j,m,1) - sq3*sq5*q->get(i,j,m,3) + 5.0*q->get(i,j,m,6)
                            - sq5*sq7*q->get(i,j,m,10) );
                    FUW_m.set(5,   q->get(i,j,m,5) - sq3*q->get(i,j,m,7) );
                    FUW_m.set(4,  -sq3*q->get(i,j,m,2) + 3.0*q->get(i,j,m,4) - sq3*sq5*q->get(i,j,m,8) );
                    FUW_m.set(3,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,3) - sq3*sq5*q->get(i,j,m,6)
                            + sq3*sq7*q->get(i,j,m,10) );
                    FUW_m.set(2,   q->get(i,j,m,2) - sq3*q->get(i,j,m,4) + sq5*q->get(i,j,m,8) );
                    FUW_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,3) + sq5*q->get(i,j,m,6)
                            - sq7*q->get(i,j,m,10) );

                    break;

                case 3: // third order

                    // Integral over element for x-derivative
                    SU.set(6,  0.0 );
                    SU.set(5,  2.0*sq3*sq5*q->get(i,j,m,2) );
                    SU.set(4,  2.0*sq3*q->get(i,j,m,3) );
                    SU.set(3,  0.0 );
                    SU.set(2,  2.0*sq3*q->get(i,j,m,1) );
                    SU.set(1,  0.0 );

                    // Boundary contribution for x-derivative (right side)
                    FUV_p.set(6,  q->get(i+1,j,m,6) );
                    FUV_p.set(5,  sq5*q->get(i+1,j,m,1) - sq3*sq5*q->get(i+1,j,m,2) + 5.0*q->get(i+1,j,m,5) );
                    FUV_p.set(4,  sq3*q->get(i+1,j,m,3) - 3.0*q->get(i+1,j,m,4) );
                    FUV_p.set(3,  q->get(i+1,j,m,3) - sq3*q->get(i+1,j,m,4) );
                    FUV_p.set(2,  sq3*q->get(i+1,j,m,1) - 3.0*q->get(i+1,j,m,2) + sq3*sq5*q->get(i+1,j,m,5) );
                    FUV_p.set(1,  q->get(i+1,j,m,1) - sq3*q->get(i+1,j,m,2) + sq5*q->get(i+1,j,m,5) );

                    // Boundary contribution for x-derivative (left side)
                    FUV_m.set(6,   q->get(i,j,m,6) );
                    FUV_m.set(5,   sq5*q->get(i,j,m,1) - sq3*sq5*q->get(i,j,m,2) + 5.0*q->get(i,j,m,5) );
                    FUV_m.set(4,  -sq3*q->get(i,j,m,3) + 3.0*q->get(i,j,m,4) );
                    FUV_m.set(3,   q->get(i,j,m,3) - sq3*q->get(i,j,m,4) );
                    FUV_m.set(2,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,2) - sq3*sq5*q->get(i,j,m,5) );
                    FUV_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,2) + sq5*q->get(i,j,m,5) );

                    // Integral over element for y-derivative
                    TU.set(6,  2.0*sq3*sq5*q->get(i,j,m,3) );
                    TU.set(5,  0.0 );
                    TU.set(4,  2.0*sq3*q->get(i,j,m,2) );	    
                    TU.set(3,  2.0*sq3*q->get(i,j,m,1) );
                    TU.set(2,  0.0 );	    
                    TU.set(1,  0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FUW_p.set(6,   sq5*q->get(i,j+1,m,1) - sq3*sq5*q->get(i,j+1,m,3) + 5.0*q->get(i,j+1,m,6) );
                    FUW_p.set(5,   q->get(i,j+1,m,5) );
                    FUW_p.set(4,   sq3*q->get(i,j+1,m,2) - 3.0*q->get(i,j+1,m,4) );
                    FUW_p.set(3,   sq3*q->get(i,j+1,m,1) - 3.0*q->get(i,j+1,m,3) + sq3*sq5*q->get(i,j+1,m,6) );
                    FUW_p.set(2,   q->get(i,j+1,m,2) - sq3*q->get(i,j+1,m,4) );
                    FUW_p.set(1,   q->get(i,j+1,m,1) - sq3*q->get(i,j+1,m,3) + sq5*q->get(i,j+1,m,6) );

                    // Boundary contribution for y-derivative (bottom side)
                    FUW_m.set(6,   sq5*q->get(i,j,m,1) - sq3*sq5*q->get(i,j,m,3) + 5.0*q->get(i,j,m,6) );
                    FUW_m.set(5,   q->get(i,j,m,5) );
                    FUW_m.set(4,  -sq3*q->get(i,j,m,2) + 3.0*q->get(i,j,m,4) );
                    FUW_m.set(3,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,3) - sq3*sq5*q->get(i,j,m,6) );
                    FUW_m.set(2,   q->get(i,j,m,2) - sq3*q->get(i,j,m,4) );
                    FUW_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,3) + sq5*q->get(i,j,m,6) );

                    break;

                case 2: // second order

                    // Integral over element for x-derivative
                    SU.set(3,  0.0 );
                    SU.set(2,  2.0*sq3*q->get(i,j,m,1) );
                    SU.set(1,  0.0 );

                    // Boundary contribution for x-derivative (right side)
                    FUV_p.set(3,  q->get(i+1,j,m,3) );
                    FUV_p.set(2,  sq3*q->get(i+1,j,m,1) - 3.0*q->get(i+1,j,m,2) );
                    FUV_p.set(1,  q->get(i+1,j,m,1) - sq3*q->get(i+1,j,m,2) );

                    // Boundary contribution for x-derivative (left side)
                    FUV_m.set(3,   q->get(i,j,m,3) );
                    FUV_m.set(2,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,2) );
                    FUV_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,2) );

                    // Integral over element for y-derivative
                    TU.set(3,  2.0*sq3*q->get(i,j,m,1) );
                    TU.set(2,  0.0 );	    
                    TU.set(1,  0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FUW_p.set(3,   sq3*q->get(i,j+1,m,1) - 3.0*q->get(i,j+1,m,3) );
                    FUW_p.set(2,   q->get(i,j+1,m,2) );
                    FUW_p.set(1,   q->get(i,j+1,m,1) - sq3*q->get(i,j+1,m,3) );

                    // Boundary contribution for y-derivative (bottom side)
                    FUW_m.set(3,  -sq3*q->get(i,j,m,1) + 3.0*q->get(i,j,m,3) );
                    FUW_m.set(2,   q->get(i,j,m,2) );
                    FUW_m.set(1,   q->get(i,j,m,1) - sq3*q->get(i,j,m,3) );

                    break;

                case 1: //first order

                    // Integral over element for x-derivative
                    SU.set(1, 0.0 );

                    // Boundary contribution for x-derivative (right side)
                    FUV_p.set(1, q->get(i+1,j,m,1) );

                    // Boundary contribution for x-derivative (left side)
                    FUV_m.set(1, q->get(i,j,m,1) );

                    // Integral over element for y-derivative
                    TU.set(1, 0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FUW_p.set(1, q->get(i,j+1,m,1) );

                    // Boundary contribution for y-derivative (bottom side)
                    FUW_m.set(1, q->get(i,j,m,1) );	    

                    break;
            }

            // ------------------------------
            // | Store gradient components  |
            // | ---------------------------|
            // |    V = q_x                 |
            // |    W = q_y                 |
            // ------------------------------
            for (int k=1; k<=kmax; k++)
            {
                V.set(i,j,m,k, ( FUV_p.get(k) - FUV_m.get(k) - SU.get(k) )/dx );
                W.set(i,j,m,k, ( FUW_p.get(k) - FUW_m.get(k) - TU.get(k) )/dy );
            }

        }
}


// -------------------------------------------------------------
// Add divergence of gradient to Lstar
// -------------------------------------------------------------
void ComputeRHS_AV(const dTensorBC2& eps, const dTensorBC4& V, 
		   const dTensorBC4& W, dTensorBC4* Lstar)
{
    int mx   = V.getsize(1);
    int my   = V.getsize(2);
    int meqn = V.getsize(3);
    int kmax = V.getsize(4);
    int mbc  = V.getmbc();
    double dx = dogParamsCart2.get_dx();
    double dy = dogParamsCart2.get_dy();
    int m = 1;

    for (int i=(2-mbc); i<=(mx+mbc-2); i++)
        for (int j=(2-mbc); j<=(my+mbc-2); j++)
        {	
            dTensor1 SV(kmax);
            dTensor1 TW(kmax);
            dTensor1 FV_p(kmax);
            dTensor1 FV_m(kmax);
            dTensor1 FW_p(kmax);
            dTensor1 FW_m(kmax);

            // Various contributions from integration-by-parts
            switch( dogParams.get_space_order() )
            {
                case 5: // fifth order

                    // Integral over element for x-derivative
                    SV.set(15, 0.0 );
                    SV.set(14, 6.0*(sq3*V.get(i,j,m,2)+sq7*V.get(i,j,m,9)) );
                    SV.set(13, 2.0*sq3*sq5*V.get(i,j,m,8) );
                    SV.set(12, 2.0*sq3*V.get(i,j,m,10) );
                    SV.set(11, 2.0*sq7*(V.get(i,j,m,3)+sq5*V.get(i,j,m,7)) );
                    SV.set(10, 0.0 );
                    SV.set(9, 2.0*sq7*(V.get(i,j,m,1)+sq5*V.get(i,j,m,5)) );
                    SV.set(8, 2.0*sq3*V.get(i,j,m,6) );
                    SV.set(7, 2.0*sq3*sq5*V.get(i,j,m,4) );
                    SV.set(6, 0.0 );
                    SV.set(5, 2.0*sq3*sq5*V.get(i,j,m,2) );
                    SV.set(4, 2.0*sq3*V.get(i,j,m,3) );
                    SV.set(3, 0.0 );
                    SV.set(2, 2.0*sq3*V.get(i,j,m,1) );
                    SV.set(1, 0.0 );	    

                    // Boundary contribution for x-derivative (right side)
                    FV_p.set(15, V.get(i,j,m,15) );
                    FV_p.set(14, 3.0*V.get(i,j,m,1) + 3.0*sq3*V.get(i,j,m,2) + 3.0*sq5*V.get(i,j,m,5) 
                            + 3.0*sq7*V.get(i,j,m,9) + 9.0*V.get(i,j,m,14) );
                    FV_p.set(13, sq5*V.get(i,j,m,6) + sq3*sq5*V.get(i,j,m,8) + 5.0*V.get(i,j,m,13) );
                    FV_p.set(12, sq3*V.get(i,j,m,10) + 3.0*V.get(i,j,m,12) );
                    FV_p.set(11, sq7*V.get(i,j,m,3) + sq3*sq7*V.get(i,j,m,4) 
                            + sq5*sq7*V.get(i,j,m,7) + 7.0*V.get(i,j,m,11) );
                    FV_p.set(10, V.get(i,j,m,10) + sq3*V.get(i,j,m,12) );
                    FV_p.set(9,  sq7*V.get(i,j,m,1) + sq3*sq7*V.get(i,j,m,2) 
                            + sq5*sq7*V.get(i,j,m,5) + 7.0*V.get(i,j,m,9) + 3.0*sq7*V.get(i,j,m,14) );
                    FV_p.set(8,  sq3*V.get(i,j,m,6) + 3.0*V.get(i,j,m,8) + sq3*sq5*V.get(i,j,m,13) );
                    FV_p.set(7,  sq5*V.get(i,j,m,3) + sq3*sq5*V.get(i,j,m,4) + 5.0*V.get(i,j,m,7) 
                            + sq5*sq7*V.get(i,j,m,11) );
                    FV_p.set(6,  V.get(i,j,m,6) + sq3*V.get(i,j,m,8) + sq5*V.get(i,j,m,13) );
                    FV_p.set(5,  sq5*V.get(i,j,m,1) + sq3*sq5*V.get(i,j,m,2) + 5.0*V.get(i,j,m,5) 
                            + sq5*sq7*V.get(i,j,m,9) + 3.0*sq5*V.get(i,j,m,14) );
                    FV_p.set(4,  sq3*V.get(i,j,m,3) + 3.0*V.get(i,j,m,4) + sq3*sq5*V.get(i,j,m,7) 
                            + sq3*sq7*V.get(i,j,m,11) );
                    FV_p.set(3,  V.get(i,j,m,3) + sq3*V.get(i,j,m,4) + sq5*V.get(i,j,m,7) 
                            + sq7*V.get(i,j,m,11) );
                    FV_p.set(2,  sq3*V.get(i,j,m,1) + 3.0*V.get(i,j,m,2) + sq3*sq5*V.get(i,j,m,5) 
                            + sq3*sq7*V.get(i,j,m,9) + 3.0*sq3*V.get(i,j,m,14) );	    
                    FV_p.set(1,  V.get(i,j,m,1) + sq3*V.get(i,j,m,2) + sq5*V.get(i,j,m,5) 
                            + sq7*V.get(i,j,m,9) + 3.0*V.get(i,j,m,14) );

                    // Boundary contribution for x-derivative (left side)
                    FV_m.set(15, V.get(i-1,j,m,15) );
                    FV_m.set(14, 3.0*V.get(i-1,j,m,1) + 3.0*sq3*V.get(i-1,j,m,2) 
                            + 3.0*sq5*V.get(i-1,j,m,5) + 3.0*sq7*V.get(i-1,j,m,9) + 9.0*V.get(i-1,j,m,14) );
                    FV_m.set(13, sq5*V.get(i-1,j,m,6) + sq3*sq5*V.get(i-1,j,m,8) + 5.0*V.get(i-1,j,m,13) );
                    FV_m.set(12, -sq3*V.get(i-1,j,m,10) - 3.0*V.get(i-1,j,m,12) );
                    FV_m.set(11, -sq7*V.get(i-1,j,m,3) - sq3*sq7*V.get(i-1,j,m,4) 
                            - sq5*sq7*V.get(i-1,j,m,7) - 7.0*V.get(i-1,j,m,11) );
                    FV_m.set(10, V.get(i-1,j,m,10) + sq3*V.get(i-1,j,m,12) );
                    FV_m.set(9, -sq7*V.get(i-1,j,m,1) - sq3*sq7*V.get(i-1,j,m,2) - sq5*sq7*V.get(i-1,j,m,5) 
                            - 7.0*V.get(i-1,j,m,9) - 3.0*sq7*V.get(i-1,j,m,14) );
                    FV_m.set(8, -sq3*V.get(i-1,j,m,6) - 3.0*V.get(i-1,j,m,8) - sq3*sq5*V.get(i-1,j,m,13) );
                    FV_m.set(7, sq5*V.get(i-1,j,m,3) + sq3*sq5*V.get(i-1,j,m,4) + 5.0*V.get(i-1,j,m,7) 
                            + sq5*sq7*V.get(i-1,j,m,11) );	    
                    FV_m.set(6, V.get(i-1,j,m,6) + sq3*V.get(i-1,j,m,8) + sq5*V.get(i-1,j,m,13) );
                    FV_m.set(5, sq5*V.get(i-1,j,m,1) + sq3*sq5*V.get(i-1,j,m,2) + 5.0*V.get(i-1,j,m,5)
                            + sq5*sq7*V.get(i-1,j,m,9) + 3.0*sq5*V.get(i-1,j,m,14) );
                    FV_m.set(4, -sq3*V.get(i-1,j,m,3) - 3.0*V.get(i-1,j,m,4) - sq3*sq5*V.get(i-1,j,m,7) 
                            - sq3*sq7*V.get(i-1,j,m,11) );	    
                    FV_m.set(3, V.get(i-1,j,m,3) + sq3*V.get(i-1,j,m,4) + sq5*V.get(i-1,j,m,7) 
                            + sq7*V.get(i-1,j,m,11) );
                    FV_m.set(2, -sq3*V.get(i-1,j,m,1) - 3.0*V.get(i-1,j,m,2) - sq3*sq5*V.get(i-1,j,m,5)
                            - sq3*sq7*V.get(i-1,j,m,9) - 3.0*sq3*V.get(i-1,j,m,14) );	    
                    FV_m.set(1, V.get(i-1,j,m,1) + sq3*V.get(i-1,j,m,2) + sq5*V.get(i-1,j,m,5)
                            + sq7*V.get(i-1,j,m,9) + 3.0*V.get(i-1,j,m,14) );

                    // Integral over element for y-derivative
                    TW.set(15, 6.0*(sq3*W.get(i,j,m,3)+sq7*W.get(i,j,m,10)) );
                    TW.set(14, 0.0 );
                    TW.set(13, 2.0*sq3*sq5*W.get(i,j,m,7) );
                    TW.set(12, 2.0*sq7*(W.get(i,j,m,2)+sq5*W.get(i,j,m,8)) );
                    TW.set(11, 2.0*sq3*W.get(i,j,m,9) );
                    TW.set(10, 2.0*sq7*(W.get(i,j,m,1)+sq5*W.get(i,j,m,6)) );
                    TW.set(9, 0 );
                    TW.set(8, 2.0*sq3*sq5*W.get(i,j,m,4) );
                    TW.set(7, 2.0*sq3*W.get(i,j,m,5) );	    
                    TW.set(6, 2.0*sq3*sq5*W.get(i,j,m,3) );
                    TW.set(5, 0.0 );
                    TW.set(4, 2.0*sq3*W.get(i,j,m,2) );	    
                    TW.set(3, 2.0*sq3*W.get(i,j,m,1) );
                    TW.set(2, 0.0 );	    
                    TW.set(1, 0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FW_p.set(15, 3.0*W.get(i,j,m,1) + 3.0*sq3*W.get(i,j,m,3) + 3.0*sq5*W.get(i,j,m,6) 
                            + 3.0*sq7*W.get(i,j,m,10) + 9.0*W.get(i,j,m,15) );
                    FW_p.set(14, W.get(i,j,m,14) );
                    FW_p.set(13, sq5*W.get(i,j,m,5) + sq3*sq5*W.get(i,j,m,7) + 5.0*W.get(i,j,m,13) );
                    FW_p.set(12, sq7*W.get(i,j,m,2) + sq3*sq7*W.get(i,j,m,4) + sq5*sq7*W.get(i,j,m,8)
                            + 7.0*W.get(i,j,m,12) );
                    FW_p.set(11, sq3*W.get(i,j,m,9) + 3.0*W.get(i,j,m,11) );
                    FW_p.set(10, sq7*W.get(i,j,m,1) + sq3*sq7*W.get(i,j,m,3) + sq5*sq7*W.get(i,j,m,6)
                            + 7.0*W.get(i,j,m,10) + 3.0*sq7*W.get(i,j,m,15) );
                    FW_p.set(9,  W.get(i,j,m,9) + sq3*W.get(i,j,m,11) );
                    FW_p.set(8,  sq5*W.get(i,j,m,2) + sq3*sq5*W.get(i,j,m,4) + 5.0*W.get(i,j,m,8)
                            + sq5*sq7*W.get(i,j,m,12) );
                    FW_p.set(7,  sq3*W.get(i,j,m,5) + 3.0*W.get(i,j,m,7) + sq3*sq5*W.get(i,j,m,13) );
                    FW_p.set(6,  sq5*W.get(i,j,m,1) + sq3*sq5*W.get(i,j,m,3) + 5.0*W.get(i,j,m,6) 
                            + sq5*sq7*W.get(i,j,m,10) + 3.0*sq5*W.get(i,j,m,15) );
                    FW_p.set(5,  W.get(i,j,m,5) + sq3*W.get(i,j,m,7) + sq5*W.get(i,j,m,13) );
                    FW_p.set(4,  sq3*W.get(i,j,m,2) + 3.0*W.get(i,j,m,4) + sq3*sq5*W.get(i,j,m,8)
                            + sq3*sq7*W.get(i,j,m,12) );
                    FW_p.set(3,  sq3*W.get(i,j,m,1) + 3.0*W.get(i,j,m,3) + sq3*sq5*W.get(i,j,m,6)
                            + sq3*sq7*W.get(i,j,m,10) + 3.0*sq3*W.get(i,j,m,15) );
                    FW_p.set(2,  W.get(i,j,m,2) + sq3*W.get(i,j,m,4) + sq5*W.get(i,j,m,8)
                            + sq7*W.get(i,j,m,12) );
                    FW_p.set(1,  W.get(i,j,m,1) + sq3*W.get(i,j,m,3) + sq5*W.get(i,j,m,6)
                            + sq7*W.get(i,j,m,10) + 3.0*W.get(i,j,m,15) );

                    // Boundary contribution for y-derivative (bottom side)    
                    FW_m.set(15, 3.0*W.get(i,j-1,m,1) + 3.0*sq3*W.get(i,j-1,m,3) 
                            + 3.0*sq5*W.get(i,j-1,m,6) + 3.0*sq7*W.get(i,j-1,m,10) + 9.0*W.get(i,j-1,m,15) );
                    FW_m.set(14, W.get(i,j-1,m,14) );
                    FW_m.set(13, sq5*W.get(i,j-1,m,5) + sq3*sq5*W.get(i,j-1,m,7) + 5.0*W.get(i,j-1,m,13) );
                    FW_m.set(12, -sq7*W.get(i,j-1,m,2) - sq3*sq7*W.get(i,j-1,m,4) - sq5*sq7*W.get(i,j-1,m,8) 
                            - 7.0*W.get(i,j-1,m,12) );
                    FW_m.set(11, -sq3*W.get(i,j-1,m,9) - 3.0*W.get(i,j-1,m,11) );
                    FW_m.set(10, -sq7*W.get(i,j-1,m,1) - sq3*sq7*W.get(i,j-1,m,3) - sq5*sq7*W.get(i,j-1,m,6) 
                            - 7.0*W.get(i,j-1,m,10) - 3.0*sq7*W.get(i,j-1,m,15) );
                    FW_m.set(9,  W.get(i,j-1,m,9) + sq3*W.get(i,j-1,m,11) );
                    FW_m.set(8,  sq5*W.get(i,j-1,m,2) + sq3*sq5*W.get(i,j-1,m,4) + 5.0*W.get(i,j-1,m,8)
                            + sq5*sq7*W.get(i,j-1,m,12) );
                    FW_m.set(7,  -sq3*W.get(i,j-1,m,5) - 3.0*W.get(i,j-1,m,7) - sq3*sq5*W.get(i,j-1,m,13) );
                    FW_m.set(6,  sq5*W.get(i,j-1,m,1) + sq3*sq5*W.get(i,j-1,m,3) + 5.0*W.get(i,j-1,m,6) 
                            + sq5*sq7*W.get(i,j-1,m,10) + 3.0*sq5*W.get(i,j-1,m,15) );
                    FW_m.set(5,  W.get(i,j-1,m,5) + sq3*W.get(i,j-1,m,7) + sq5*W.get(i,j-1,m,13) );
                    FW_m.set(4,  -sq3*W.get(i,j-1,m,2) - 3*W.get(i,j-1,m,4) - sq3*sq5*W.get(i,j-1,m,8)
                            - sq3*sq7*W.get(i,j-1,m,12) );
                    FW_m.set(3,  -sq3*W.get(i,j-1,m,1) - 3*W.get(i,j-1,m,3) - sq3*sq5*W.get(i,j-1,m,6) 
                            - sq3*sq7*W.get(i,j-1,m,10) - 3.0*sq3*W.get(i,j-1,m,15) );
                    FW_m.set(2,  W.get(i,j-1,m,2) + sq3*W.get(i,j-1,m,4) + sq5*W.get(i,j-1,m,8) 
                            + sq7*W.get(i,j-1,m,12) );
                    FW_m.set(1,  W.get(i,j-1,m,1) + sq3*W.get(i,j-1,m,3) + sq5*W.get(i,j-1,m,6) 
                            + sq7*W.get(i,j-1,m,10) + 3.0*W.get(i,j-1,m,15) );

                    break;

                case 4: // fourth order

                    // Integral over element for x-derivative
                    SV.set(10, 0.0 );
                    SV.set(9, 2.0*sq7*(V.get(i,j,m,1)+sq5*V.get(i,j,m,5)) );
                    SV.set(8, 2.0*sq3*V.get(i,j,m,6) );
                    SV.set(7, 2.0*sq3*sq5*V.get(i,j,m,4) );
                    SV.set(6, 0.0 );
                    SV.set(5, 2.0*sq3*sq5*V.get(i,j,m,2) );
                    SV.set(4, 2.0*sq3*V.get(i,j,m,3) );
                    SV.set(3, 0.0 );
                    SV.set(2, 2.0*sq3*V.get(i,j,m,1) );
                    SV.set(1, 0.0 );	    

                    // Boundary contribution for x-derivative (right side)
                    FV_p.set(10, V.get(i,j,m,10) );
                    FV_p.set(9,  sq7*V.get(i,j,m,1) + sq3*sq7*V.get(i,j,m,2) 
                            + sq5*sq7*V.get(i,j,m,5) + 7.0*V.get(i,j,m,9) );
                    FV_p.set(8,  sq3*V.get(i,j,m,6) + 3.0*V.get(i,j,m,8) );
                    FV_p.set(7,  sq5*V.get(i,j,m,3) + sq3*sq5*V.get(i,j,m,4) + 5.0*V.get(i,j,m,7) );
                    FV_p.set(6,  V.get(i,j,m,6) + sq3*V.get(i,j,m,8) );
                    FV_p.set(5,  sq5*V.get(i,j,m,1) + sq3*sq5*V.get(i,j,m,2) + 5.0*V.get(i,j,m,5) 
                            + sq5*sq7*V.get(i,j,m,9) );
                    FV_p.set(4,  sq3*V.get(i,j,m,3) + 3.0*V.get(i,j,m,4) + sq3*sq5*V.get(i,j,m,7) );
                    FV_p.set(3,  V.get(i,j,m,3) + sq3*V.get(i,j,m,4) + sq5*V.get(i,j,m,7) );
                    FV_p.set(2,  sq3*V.get(i,j,m,1) + 3.0*V.get(i,j,m,2) + sq3*sq5*V.get(i,j,m,5) 
                            + sq3*sq7*V.get(i,j,m,9) );	    
                    FV_p.set(1,  V.get(i,j,m,1) + sq3*V.get(i,j,m,2) + sq5*V.get(i,j,m,5) 
                            + sq7*V.get(i,j,m,9) );

                    // Boundary contribution for x-derivative (left side)
                    FV_m.set(10, V.get(i-1,j,m,10) );
                    FV_m.set(9, -sq7*V.get(i-1,j,m,1) - sq3*sq7*V.get(i-1,j,m,2) - sq5*sq7*V.get(i-1,j,m,5) 
                            - 7.0*V.get(i-1,j,m,9) );
                    FV_m.set(8, -sq3*V.get(i-1,j,m,6) - 3.0*V.get(i-1,j,m,8) );
                    FV_m.set(7, sq5*V.get(i-1,j,m,3) + sq3*sq5*V.get(i-1,j,m,4) + 5.0*V.get(i-1,j,m,7) );
                    FV_m.set(6, V.get(i-1,j,m,6) + sq3*V.get(i-1,j,m,8) );
                    FV_m.set(5, sq5*V.get(i-1,j,m,1) + sq3*sq5*V.get(i-1,j,m,2) + 5.0*V.get(i-1,j,m,5)
                            + sq5*sq7*V.get(i-1,j,m,9) );
                    FV_m.set(4, -sq3*V.get(i-1,j,m,3) - 3.0*V.get(i-1,j,m,4) - sq3*sq5*V.get(i-1,j,m,7) );	    
                    FV_m.set(3, V.get(i-1,j,m,3) + sq3*V.get(i-1,j,m,4) + sq5*V.get(i-1,j,m,7) );
                    FV_m.set(2, -sq3*V.get(i-1,j,m,1) - 3.0*V.get(i-1,j,m,2) - sq3*sq5*V.get(i-1,j,m,5)
                            - sq3*sq7*V.get(i-1,j,m,9) );	    
                    FV_m.set(1, V.get(i-1,j,m,1) + sq3*V.get(i-1,j,m,2) + sq5*V.get(i-1,j,m,5)
                            + sq7*V.get(i-1,j,m,9) );

                    // Integral over element for y-derivative
                    TW.set(10, 2.0*sq7*(W.get(i,j,m,1)+sq5*W.get(i,j,m,6)) );
                    TW.set(9, 0 );
                    TW.set(8, 2.0*sq3*sq5*W.get(i,j,m,4) );
                    TW.set(7, 2.0*sq3*W.get(i,j,m,5) );	    
                    TW.set(6, 2.0*sq3*sq5*W.get(i,j,m,3) );
                    TW.set(5, 0.0 );
                    TW.set(4, 2.0*sq3*W.get(i,j,m,2) );	    
                    TW.set(3, 2.0*sq3*W.get(i,j,m,1) );
                    TW.set(2, 0.0 );	    
                    TW.set(1, 0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FW_p.set(10, sq7*W.get(i,j,m,1) + sq3*sq7*W.get(i,j,m,3) + sq5*sq7*W.get(i,j,m,6)
                            + 7.0*W.get(i,j,m,10) );
                    FW_p.set(9,  W.get(i,j,m,9) );
                    FW_p.set(8,  sq5*W.get(i,j,m,2) + sq3*sq5*W.get(i,j,m,4) + 5.0*W.get(i,j,m,8) );
                    FW_p.set(7,  sq3*W.get(i,j,m,5) + 3.0*W.get(i,j,m,7) );
                    FW_p.set(6,  sq5*W.get(i,j,m,1) + sq3*sq5*W.get(i,j,m,3) + 5.0*W.get(i,j,m,6) 
                            + sq5*sq7*W.get(i,j,m,10) );
                    FW_p.set(5,  W.get(i,j,m,5) + sq3*W.get(i,j,m,7) );
                    FW_p.set(4,  sq3*W.get(i,j,m,2) + 3.0*W.get(i,j,m,4) + sq3*sq5*W.get(i,j,m,8) );
                    FW_p.set(3,  sq3*W.get(i,j,m,1) + 3.0*W.get(i,j,m,3) + sq3*sq5*W.get(i,j,m,6)
                            + sq3*sq7*W.get(i,j,m,10) );
                    FW_p.set(2,  W.get(i,j,m,2) + sq3*W.get(i,j,m,4) + sq5*W.get(i,j,m,8) );
                    FW_p.set(1,  W.get(i,j,m,1) + sq3*W.get(i,j,m,3) + sq5*W.get(i,j,m,6)
                            + sq7*W.get(i,j,m,10) );

                    // Boundary contribution for y-derivative (bottom side)    
                    FW_m.set(10, -sq7*W.get(i,j-1,m,1) - sq3*sq7*W.get(i,j-1,m,3) - sq5*sq7*W.get(i,j-1,m,6) 
                            - 7.0*W.get(i,j-1,m,10) );
                    FW_m.set(9,  W.get(i,j-1,m,9) );
                    FW_m.set(8,  sq5*W.get(i,j-1,m,2) + sq3*sq5*W.get(i,j-1,m,4) + 5.0*W.get(i,j-1,m,8) );
                    FW_m.set(7,  -sq3*W.get(i,j-1,m,5) - 3.0*W.get(i,j-1,m,7) );
                    FW_m.set(6,  sq5*W.get(i,j-1,m,1) + sq3*sq5*W.get(i,j-1,m,3) + 5.0*W.get(i,j-1,m,6) 
                            + sq5*sq7*W.get(i,j-1,m,10) );
                    FW_m.set(5,  W.get(i,j-1,m,5) + sq3*W.get(i,j-1,m,7) );
                    FW_m.set(4,  -sq3*W.get(i,j-1,m,2) - 3*W.get(i,j-1,m,4) - sq3*sq5*W.get(i,j-1,m,8) );
                    FW_m.set(3,  -sq3*W.get(i,j-1,m,1) - 3*W.get(i,j-1,m,3) - sq3*sq5*W.get(i,j-1,m,6) 
                            - sq3*sq7*W.get(i,j-1,m,10) );
                    FW_m.set(2,  W.get(i,j-1,m,2) + sq3*W.get(i,j-1,m,4) + sq5*W.get(i,j-1,m,8) );
                    FW_m.set(1,  W.get(i,j-1,m,1) + sq3*W.get(i,j-1,m,3) + sq5*W.get(i,j-1,m,6) 
                            + sq7*W.get(i,j-1,m,10) );

                    break;

                case 3: // third order

                    // Integral over element for x-derivative
                    SV.set(6, 0.0 );
                    SV.set(5, 2.0*sq3*sq5*V.get(i,j,m,2) );
                    SV.set(4, 2.0*sq3*V.get(i,j,m,3) );
                    SV.set(3, 0.0 );
                    SV.set(2, 2.0*sq3*V.get(i,j,m,1) );
                    SV.set(1, 0.0 );	    

                    // Boundary contribution for x-derivative (right side)
                    FV_p.set(6,  V.get(i,j,m,6) );
                    FV_p.set(5,  sq5*V.get(i,j,m,1) + sq3*sq5*V.get(i,j,m,2) + 5.0*V.get(i,j,m,5) );
                    FV_p.set(4,  sq3*V.get(i,j,m,3) + 3.0*V.get(i,j,m,4) );
                    FV_p.set(3,  V.get(i,j,m,3) + sq3*V.get(i,j,m,4) );
                    FV_p.set(2,  sq3*V.get(i,j,m,1) + 3.0*V.get(i,j,m,2) + sq3*sq5*V.get(i,j,m,5) );	    
                    FV_p.set(1,  V.get(i,j,m,1) + sq3*V.get(i,j,m,2) + sq5*V.get(i,j,m,5) );

                    // Boundary contribution for x-derivative (left side)
                    FV_m.set(6, V.get(i-1,j,m,6) );
                    FV_m.set(5, sq5*V.get(i-1,j,m,1) + sq3*sq5*V.get(i-1,j,m,2) + 5.0*V.get(i-1,j,m,5) );
                    FV_m.set(4, -sq3*V.get(i-1,j,m,3) - 3.0*V.get(i-1,j,m,4) );	    
                    FV_m.set(3, V.get(i-1,j,m,3) + sq3*V.get(i-1,j,m,4) );
                    FV_m.set(2, -sq3*V.get(i-1,j,m,1) - 3.0*V.get(i-1,j,m,2) - sq3*sq5*V.get(i-1,j,m,5) );	    
                    FV_m.set(1, V.get(i-1,j,m,1) + sq3*V.get(i-1,j,m,2) + sq5*V.get(i-1,j,m,5) );

                    // Integral over element for y-derivative
                    TW.set(6, 2.0*sq3*sq5*W.get(i,j,m,3) );
                    TW.set(5, 0.0 );
                    TW.set(4, 2.0*sq3*W.get(i,j,m,2) );	    
                    TW.set(3, 2.0*sq3*W.get(i,j,m,1) );
                    TW.set(2, 0.0 );	    
                    TW.set(1, 0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FW_p.set(6,  sq5*W.get(i,j,m,1) + sq3*sq5*W.get(i,j,m,3) + 5.0*W.get(i,j,m,6) );
                    FW_p.set(5,  W.get(i,j,m,5) );
                    FW_p.set(4,  sq3*W.get(i,j,m,2) + 3.0*W.get(i,j,m,4) );
                    FW_p.set(3,  sq3*W.get(i,j,m,1) + 3.0*W.get(i,j,m,3) + sq3*sq5*W.get(i,j,m,6) );
                    FW_p.set(2,  W.get(i,j,m,2) + sq3*W.get(i,j,m,4) );
                    FW_p.set(1,  W.get(i,j,m,1) + sq3*W.get(i,j,m,3) + sq5*W.get(i,j,m,6) );

                    // Boundary contribution for y-derivative (bottom side)    
                    FW_m.set(6,  sq5*W.get(i,j-1,m,1) + sq3*sq5*W.get(i,j-1,m,3) + 5.0*W.get(i,j-1,m,6) );
                    FW_m.set(5,  W.get(i,j-1,m,5) );
                    FW_m.set(4,  -sq3*W.get(i,j-1,m,2) - 3*W.get(i,j-1,m,4) );
                    FW_m.set(3,  -sq3*W.get(i,j-1,m,1) - 3*W.get(i,j-1,m,3) - sq3*sq5*W.get(i,j-1,m,6) );
                    FW_m.set(2,  W.get(i,j-1,m,2) + sq3*W.get(i,j-1,m,4) );
                    FW_m.set(1,  W.get(i,j-1,m,1) + sq3*W.get(i,j-1,m,3) + sq5*W.get(i,j-1,m,6) );

                    break;

                case 2: // second order

                    // Integral over element for x-derivative
                    SV.set(3, 0.0 );
                    SV.set(2, 2.0*sq3*V.get(i,j,m,1) );
                    SV.set(1, 0.0 );	    

                    // Boundary contribution for x-derivative (right side)
                    FV_p.set(3,  V.get(i,j,m,3) );
                    FV_p.set(2,  sq3*V.get(i,j,m,1) + 3.0*V.get(i,j,m,2) );	    
                    FV_p.set(1,  V.get(i,j,m,1) + sq3*V.get(i,j,m,2) );

                    // Boundary contribution for x-derivative (left side)
                    FV_m.set(3, V.get(i-1,j,m,3) );
                    FV_m.set(2, -sq3*V.get(i-1,j,m,1) - 3.0*V.get(i-1,j,m,2) );	    
                    FV_m.set(1, V.get(i-1,j,m,1) + sq3*V.get(i-1,j,m,2) );

                    // Integral over element for y-derivative
                    TW.set(3, 2.0*sq3*W.get(i,j,m,1) );
                    TW.set(2, 0.0 );	    
                    TW.set(1, 0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FW_p.set(3,  sq3*W.get(i,j,m,1) + 3.0*W.get(i,j,m,3) );
                    FW_p.set(2,  W.get(i,j,m,2) );
                    FW_p.set(1,  W.get(i,j,m,1) + sq3*W.get(i,j,m,3) );

                    // Boundary contribution for y-derivative (bottom side)    
                    FW_m.set(3, -sq3*W.get(i,j-1,m,1) - 3*W.get(i,j-1,m,3) );
                    FW_m.set(2,  W.get(i,j-1,m,2) );
                    FW_m.set(1,  W.get(i,j-1,m,1) + sq3*W.get(i,j-1,m,3) );

                    break;

                case 1: //first order

                    // Integral over element for x-derivative
                    SV.set(1, 0.0 );	    

                    // Boundary contribution for x-derivative (right side)
                    FV_p.set(1, V.get(i,j,m,1) );

                    // Boundary contribution for x-derivative (left side)
                    FV_m.set(1, V.get(i-1,j,m,1) );

                    // Integral over element for y-derivative
                    TW.set(1, 0.0 );

                    // Boundary contribution for y-derivative (top side)
                    FW_p.set(1, W.get(i,j,m,1) );

                    // Boundary contribution for y-derivative (bottom side)    
                    FW_m.set(1, W.get(i,j-1,m,1) );

                    break;
            }

            // Scale each piece by appropriate value of eps
            for (int k=1; k<=kmax; k++)
            {
                SV.set(k,   SV.get(k) * eps.get(i,j) );
                FV_p.set(k, FV_p.get(k) * 0.5*(eps.get(i+1,j)+eps.get(i,j)) );
                FV_m.set(k, FV_m.get(k) * 0.5*(eps.get(i-1,j)+eps.get(i,j)) );

                TW.set(k,   TW.get(k) * eps.get(i,j) );
                FW_p.set(k, FW_p.get(k) * 0.5*(eps.get(i,j+1)+eps.get(i,j)) );
                FW_m.set(k, FW_m.get(k) * 0.5*(eps.get(i,j-1)+eps.get(i,j)) );
            }

            // Add divergence of gradient to right-hand side "Lstar"
            for (int k=1; k<=kmax; k++)
            {
	      Lstar->set(i,j,m,k, Lstar->get(i,j,m,k) 
			 + ( FV_p.get(k) - FV_m.get(k) - SV.get(k) )/dx
			 + ( FW_p.get(k) - FW_m.get(k) - TW.get(k) )/dy );
            }

        }
}
