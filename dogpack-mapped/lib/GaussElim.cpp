#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "tensors.h"

// Gaussian Elimination with Partial Pivoting Routine.
// Produces the solution to inMAT*outSOLN = inRHS,
// where inMAT is an NXN matrix, outSOLN and inRHS are 
// vectors of length N.
void GaussElim(dTensor2 inMAT, dTensor1 inRHS, dTensor1& outSOLN)
{
    int n1,m1,n2,n3;
    int n,i,j,k,p;
    double dm,maxA,sum,tmp;

    n1 = inMAT.getsize(1);
    m1 = inMAT.getsize(2);
    n2 = inRHS.getsize();
    n3 = outSOLN.getsize();

    if (n1!=m1) 
    { 
        printf(" ERROR in GaussElim.cpp:  coeff. matrix is not square \n");
        printf("        Size of coeff. matrix = %i X %i\n",n1,m1);
        printf("         Length of RHS vector = %i\n",n2);
        printf("     Length of outSOLN vector = %i\n",n3);
        exit(1);
    }

    if (n2!=n1) 
    { 
        printf(" ERROR in GaussElim.cpp:  incorrect size of RHS vector \n");
        printf("        Size of coeff. matrix = %i X %i\n",n1,m1);
        printf("         Length of RHS vector = %i\n",n2);
        printf("     Length of outSOLN vector = %i\n",n3);
        exit(1);
    }

    if (n3!=n1) 
    { 
        printf(" ERROR in GaussElim.cpp:  incorrect size of outSOLN vector \n");
        printf("        Size of coeff. matrix = %i X %i\n",n1,m1);
        printf("         Length of RHS vector = %i\n",n2);
        printf("     Length of outSOLN vector = %i\n",n3);
        exit(1);
    }  

    n = n1;

    for (i=1; i<=(n-1); i++)
    {
        // Select pivot
        maxA = -100.0e0;
        for (j=i; j<=n; j++)
        {
            if ( fabs(inMAT.get(j,i)) > maxA )
            {
                p = j;
                maxA = fabs(inMAT.get(j,i));
            }
        }

        // See if matrix is singular
        if (maxA <= 1.0e-14)
        {
            printf("  Error in GaussElim.cpp: cannot invert system\n");
        }

        // Pivot
        if (p!=i)
        {
            for (j=1; j<=n; j++)
            {
                tmp = inMAT.get(i,j);
                inMAT.set(i,j, inMAT.get(p,j));
                inMAT.set(p,j, tmp);
            }

            tmp = inRHS.get(i);
            inRHS.set(i, inRHS.get(p) );
            inRHS.set(p, tmp );
        }

        // Eliminate below diagonal
        for (j=(i+1); j<=n; j++)
        {
            dm = inMAT.get(j,i)/inMAT.get(i,i);
            for (k=1; k<=n; k++)
            {
                inMAT.set(j,k, (inMAT.get(j,k) - dm*inMAT.get(i,k)) );
            }
            inRHS.set(j, inRHS.get(j) - dm*inRHS.get(i) );
        }
    }

    // Backward substitution
    outSOLN.set(n, inRHS.get(n)/inMAT.get(n,n) );
    for (j=1; j<=(n-1); j++)
    {
        sum = 0.0e0;

        for (k=(n-j+1); k<=n; k++)
        { sum = sum + inMAT.get(n-j,k)*outSOLN.get(k); }

        outSOLN.set(n-j, (inRHS.get(n-j) - sum)/inMAT.get(n-j,n-j) );
    }

}

// Gaussian Elimination with Partial Pivoting Routine.
// Produces the inverse matrix outMat from the matrix inMat.
// Both of these matrices must be square and the same size.
//
// Note: this routine copies inMat when it's called.
//
void GaussElimMatrixInv(dTensor2 inMat, dTensor2& outMat)
{
    int n1,m1,n2,m2;
    int n,i,j,k,p;
    double dm,maxA,tmp;

    n1 = inMat.getsize(1);
    m1 = inMat.getsize(2);
    n2 = outMat.getsize(1);
    m2 = outMat.getsize(2);

    if (n1!=m1 || n2!=m2) 
    { 
        printf(" ERROR in GaussElim.cpp:  matrices are not square \n \n");
        exit(1);
    }

    if (n1!=n2 || m1!=m2) 
    { 
        printf(" ERROR in GaussElim.cpp:  input and output matrices are of different size  \n \n");
        exit(1);
    }

    n = n1;

    // Initialize outMat to the identity matrix
    for (i=1; i<=n; i++)
    {
        for (j=1; j<=n; j++)
        {
            if (i==j)
            { outMat.set(i,j,1.0e0); }
            else
            { outMat.set(i,j,0.0e0); }
        }
    }

    for (i=1; i<=n; i++)
    {
        // Select pivot
        maxA = -100.0e0;
        for (j=i; j<=n; j++)
        {
            if ( fabs(inMat.get(j,i)) > maxA )
            {
                p = j;
                maxA = fabs(inMat.get(j,i));
            }
        }

        // See if matrix is singular
        if (maxA <= 1.0e-14)
        {
            printf("  Error in GaussElim.cpp: cannot invert system \n \n");
            exit(1);
        }

        // Pivot matrices
        if (p!=i)
        {
            for (j=1; j<=n; j++)
            {
                tmp = inMat.get(i,j);
                inMat.set(i,j, inMat.get(p,j));
                inMat.set(p,j, tmp);

                tmp = outMat.get(i,j);
                outMat.set(i,j, outMat.get(p,j));
                outMat.set(p,j, tmp);
            }
        }

        // Eliminate below diagonal
        for (j=(i+1); j<=n; j++)
        {
            dm = inMat.get(j,i)/inMat.get(i,i);
            for (k=1; k<=n; k++)
            {
                inMat.set(j,k, (inMat.get(j,k) - dm*inMat.get(i,k)) );
                outMat.set(j,k, (outMat.get(j,k) - dm*outMat.get(i,k)) );
            }
        }

        // Eliminate above diagonal
        for (j=1; j<=(i-1); j++)
        {
            dm = inMat.get(j,i)/inMat.get(i,i);
            for (k=1; k<=n; k++)
            {
                inMat.set(j,k, (inMat.get(j,k) - dm*inMat.get(i,k)) );
                outMat.set(j,k, (outMat.get(j,k) - dm*outMat.get(i,k)) );
            }
        }

        // Rescale current row
        tmp = inMat.get(i,i);
        for (k=1; k<=n; k++)
        {
            inMat.set(i,k, (inMat.get(i,k)/tmp) );
            outMat.set(i,k, (outMat.get(i,k)/tmp) );
        }

    }

}
