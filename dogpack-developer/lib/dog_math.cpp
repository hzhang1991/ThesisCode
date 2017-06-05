#include <cmath>
#include "dog_math.h"
#define NDIMS 0
#include "dogdefs.h"

//////////////////////////////////////////////////////////////////////////////
//  The purpose of this module is to describe common math functions
//  what we'd like to be able to perform in DogPack.
//
// See also: 
//
//   $DOGPACK/constants.h  for a list of commonly used variables such as pi.
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Returns the maximum value of a single vector
//////////////////////////////////////////////////////////////////////////////
double Max(const dTensor1& vec)
{
    double max = 0.0;
    for( int i=1; i <= vec.getsize(); i++)
    {
        max = Max(max, vec.get(i));
    }
    return max;
}

//////////////////////////////////////////////////////////////////////////////
// Returns the maximum value of a single vector
//////////////////////////////////////////////////////////////////////////////
/*double Max(const dTensorBC1& vec)
{
    double Max(double,double);
    double max = 0.0;
    double mbc = vec.getmbc();
    for( int i=1-mbc; i <= vec.getsize()+mbc; i++)
    {
        max = Max(max, vec.get(i));
    }
    return max;
}
*/

//////////////////////////////////////////////////////////////////////////////
// Returns a vector containing the maximum value of each component in vec1 and
// vec2
//////////////////////////////////////////////////////////////////////////////
void Max(const dTensor1& vec1, const dTensor1& vec2, dTensor1& max)
{
    assert_eq(vec1.getsize(),vec2.getsize());
    assert_eq(max.getsize(),vec1.getsize());

    for( int i=1; i <= vec1.getsize(); i++)
    {
        max.set(i,  Max(vec1.get(i), vec2.get(i)) );
    }
}

//////////////////////////////////////////////////////////////////////////////
// Returns the maximum value found in Matrix A
//////////////////////////////////////////////////////////////////////////////
double Max(const dTensorBC2& A)
{
    int mbc = A.getmbc();
    double tmp = A.get(1,1);
    for( int i=1-mbc; i <= A.getsize(1)+mbc; i++)
    for( int j=1-mbc; j <= A.getsize(2)+mbc; j++)
    {
        tmp = Max(tmp, A.get(i,j) );
    }
    return tmp;
}

// Compute the product: A*B. Save output in C.
//
// Tensors need to have the same size.  That is, if 
// A is mxn and B is nxl, then C should be mxl.
//
void TensorMultiply( const dTensor2& A, const dTensor2& B, dTensor2& C )
{
    // Quick error check:
    assert_eq( A.getsize( 2 ), B.getsize(1) );
    assert_eq( A.getsize( 1 ), C.getsize(1) );
    assert_eq( B.getsize( 2 ), C.getsize(2) );

    for( int k=1; k <= A.getsize(1); k++ )
    for( int l=1; l <= B.getsize(2); l++ )
    {

        // TODO: FOR LARGE MATRICES, IT MAY BE FASTER TO FIRST TAKE THE
        // TRANSPOSE OF B, AND THEN DO THE MULTIPLCIATION ... (-DS)
        double tmp = 0.;
        for( int j=1; j <= A.getsize(2); j++ )
        {
            tmp += A.get(k,j)*B.get(j, l );
        }
        
    }
}

// Compute the product: A*v. Save output in b.
//
void MatVecMultiply( const dTensor2& A, const dTensor1& v, dTensor1& b )
{
    // Quick error check:
    assert_eq( A.getsize( 2 ), v.getsize() );
    assert_eq( A.getsize( 1 ), b.getsize() );

    for( int k=1; k <= A.getsize(1); k++ )
    {

        double tmp = 0.;
        for( int j=1; j <= A.getsize(2); j++ )
        {
            tmp += A.get(k,j)*v.get(j);
        }
        b.set( k, tmp );
        
    }
}

//////////////////////////////////////////////////////////////////////////////
// Multiply two functions Q1 * Q2 and write to output qnew
//
// Q1(x) = \sum_{i,k} Q1_i^{(k)} \phi_i^{(k)}(x)
// Q2(x) = \sum_{i,k} Q2_i^{(k)} \phi_i^{(k)}(x)
//
// The resultant function qnew is after doing the projections.
//////////////////////////////////////////////////////////////////////////////
void MultiplyFunctions(const dTensorBC3& Q1, 
    const dTensorBC3& Q2, dTensorBC3& qout)
{
    const int melems = qout.getsize(1);
    const int meqn   = qout.getsize(2);
    const int kmax   = qout.getsize(3);
    const int mbc    = qout.getmbc();

    // TODO - I THINK THIS COULD BE MUCH MORE EFFICIENT IF WE WROTE THESE OUT BY
    // HAND RATHER THAN USING A MATRIX TO DO THE MULTIPLICATION.  THERE ARE A
    // LOT OF ZEROS IN THIS MATRIX THAT DO NOT NEED TO BE INCLUDED
    dTensor3 M(5,5,5);
    M.setall(0.);

    // k == 5
    M.set(5,1, 5, 1.0 );
    M.set(1,5, 5, 1.0 );

    M.set(2,4, 5, 4.0*sq3*sq7/21.0 );
    M.set(4,2, 5, M.get(2,4,5) );

    M.set(3,5, 5, 20.0/77.0*sq5 );
    M.set(5,3, 5, M.get(3,5,5) );

    M.set(3,3, 5, 6.0/7.0 );
    M.set(4,4, 5, 6.0/11.0 );
    M.set(5,5, 5, 486.0/1001.0 );

    // k == 4
    M.set(1,4, 4, 1.0 );
    M.set(4,1, 4, 1.0 );
    M.set(2,3, 4, sq3*sq5*sq7*3.0/35.0 );
    M.set(3,2, 4, M.get(2,3,4) );
    M.set(2,5, 4, 4.0*sq3*sq7/21.0 );
    M.set(5,2, 4, 4.0*sq3*sq7/21.0 );
    M.set(3,4, 4, 4.0*sq5/15.0 );
    M.set(4,3, 4, 4.0*sq5/15.0 );
    M.set(4,5, 4, 6.0/11.0 );
    M.set(5,4, 4, 6.0/11.0 );

    // k == 3
    M.set(1,3, 3, 1.0 );
    M.set(3,1, 3, 1.0 );
    M.set(2,2, 3, 2.0*sq5/5.0 );
    M.set(2,4, 3, 3.0*sq3*sq5*sq7/35.0 );
    M.set(4,2, 3, 3.0*sq3*sq5*sq7/35.0 );
    M.set(3,3, 3, 2.0*sq5/7.0 );
    M.set(3,5, 3, 6.0/7.0 );
    M.set(5,3, 3, 6.0/7.0 );
    M.set(4,4, 3, 4.0*sq5/15.0 );
    M.set(5,5, 3, 20.0*sq5/77.0 );

    // k == 2
    M.set(1,2, 2, 1.0 );
    M.set(2,1, 2, 1.0 );
    M.set(2,3, 2, 2.0*sq5/5.0 );
    M.set(3,2, 2, 2.0*sq5/5.0 );
    M.set(3,4, 2, sq3*sq5*sq7*3.0/35.0 );
    M.set(4,3, 2, sq3*sq5*sq7*3.0/35.0 );
    M.set(4,5, 2, sq3*sq7*4.0/21.0 );
    M.set(5,4, 2, sq3*sq7*4.0/21.0 );

    // k == 1
    M.set(1,1, 1, 1.0 );
    M.set(2,2, 1, 1.0 );
    M.set(3,3, 1, 1.0 );
    M.set(4,4, 1, 1.0 );
    M.set(5,5, 1, 1.0 );

//#pragma omp parallel for
    for(int i=1-mbc; i<=melems+mbc; i++ )
    for(int me=1; me <= meqn; me++)
    {
        double tmp = 0.0;
        for( int  k=1; k <= kmax; k++)
        {
            double tmp = 0.0;
            for( int k1=1; k1 <= kmax; k1++)
            for( int k2=1; k2 <= kmax; k2++)
            { tmp += Q1.get(i, me, k1)*Q2.get(i, me, k2)*M.get(k1,k2,k); }
            qout.set(i,me,k,tmp);
        }

    }

}

/*
        switch(kmax)
        {

            case 5:

                ////////////// k == 5 /////////////
                tmp = 0.0;

                tmp += 486.0/1001.0*Q1.get(i,me,5)*Q2.get(i,me,5);
                tmp += 6.0/11.0*Q1.get(i,me,4)*Q2.get(i,me,4);
                tmp += 6.0/7.0*Q1.get(i,me,3)*Q2.get(i,me,3);

                tmp += Q1.get(i,me,5) * Q2.get(i,me,3) * 
                    20.0 * sq5 / 77.0;
                tmp += Q2.get(i,me,5) * Q1.get(i,me,3) * 
                    20.0 * sq5 / 77.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,4) * 
                    4.0*sq3*sq7/21.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,4) * 
                    4.0*sq3*sq7/21.0;

                tmp += Q1.get(i,me,5) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,5) * Q1.get(i,me,1);
                qnew.set(i,me,5,tmp);


                ////////////// k == 4 /////////////
                tmp = 0.0;
                tmp += Q1.get(i,me,2) * Q2.get(i,me,5) * 
                    4.0 * sq3 * sq7 / 21.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,5) * 
                    4.0 * sq3 * sq7 / 21.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,5) * 
                    6.0 / 11.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,5) * 
                    6.0 / 11.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,3) * 
                    4.0 * sq5 / 15.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,3) * 
                    4.0 * sq5 / 15.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,4) * Q1.get(i,me,1);
                qnew.set(i,me,4,tmp);


                ////////////// k == 3 /////////////
                tmp = 0.0;

                tmp += 20.0*sq5/77.0*Q1.get(i,me,5)*Q2.get(i,me,5);
                tmp += Q1.get(i,me,3) * Q2.get(i,me,5) * 
                    6.0 / 7.0;
                tmp += Q2.get(i,me,3) * Q1.get(i,me,5) * 
                    6.0 / 7.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,4) * 
                        4.0 / (3.0*sq5);

                tmp += 2.0*sq5/7.0*Q1.get(i,me,3)*Q2.get(i,me,3);

                tmp += Q1.get(i,me,2) * Q2.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;

                tmp += 2.0 / sq5 * Q1.get(i,me,2)*Q2.get(i,me,2);

                tmp += Q1.get(i,me,3) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,3) * Q1.get(i,me,1);
                qnew.set(i,me,3,tmp);

                ////////////// k == 2 /////////////
                tmp = 0.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,5) * 4.0*sq3*sq7/21.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,5) * 4.0*sq3*sq7/21.0;

                tmp += Q1.get(i,me,3) * Q2.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;
                tmp += Q2.get(i,me,3) * Q1.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 2.0 / sq5;
                tmp += Q1.get(i,me,3) * Q2.get(i,me,2) * 2.0 / sq5;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,1);
                tmp += Q1.get(i,me,1) * Q2.get(i,me,2);
                qnew.set(i,me,2,tmp);

                ////////////// k == 1 /////////////
                tmp = 0.0;
                for( int k=1; k <= kmax; k++ )
                { tmp += Q1.get(i,me,k) * Q2.get(i,me,k); }
                qnew.set(i,me,1, tmp);
                break;


            case 4:

                ////////////// k == 4 /////////////
                tmp = 0.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,3) * 
                    4.0 * sq5 / 15.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,3) * 
                    4.0 * sq5 / 15.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,4) * Q1.get(i,me,1);
                qnew.set(i,me,4,tmp);


                ////////////// k == 3 /////////////
                tmp = 0.0;
                tmp += Q1.get(i,me,4) * Q2.get(i,me,4) * 
                        4.0 / (3.0*sq5);

                tmp += 2.0*sq5/7.0*Q1.get(i,me,3)*Q2.get(i,me,3);

                tmp += Q1.get(i,me,2) * Q2.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;

                tmp += 2.0 / sq5 * Q1.get(i,me,2)*Q2.get(i,me,2);

                tmp += Q1.get(i,me,3) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,3) * Q1.get(i,me,1);
                qnew.set(i,me,3,tmp);

                ////////////// k == 2 /////////////
                tmp = 0.0;

                tmp += Q1.get(i,me,3) * Q2.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;
                tmp += Q2.get(i,me,3) * Q1.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 2.0 / sq5;
                tmp += Q1.get(i,me,3) * Q2.get(i,me,2) * 2.0 / sq5;
                tmp += Q1.get(i,me,2) * Q2.get(i,me,1);
                tmp += Q1.get(i,me,1) * Q2.get(i,me,2);
                qnew.set(i,me,2,tmp);

                ////////////// k == 1 /////////////
                tmp = 0.0;
                for( int k=1; k <= kmax; k++ )
                {
                    tmp += Q1.get(i,me,k) * Q2.get(i,me,k);
                }
                qnew.set(i,me,1, tmp);
                break;

            default:
                perror(" method not implemented yet\n");
                exit(1);
        }
    }
    */
