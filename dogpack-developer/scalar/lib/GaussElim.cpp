#include "tensors.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

void GaussElim(dTensor2 inMAT, 
        dTensor1 inRHS, 
        dTensor1& outSOLN)
//
//    # Gaussian Elimination with Partial Pivoting Routine.
//    # Produces the solution to inMAT*outSOLN = inRHS,
//    # where inMAT is an NXN matrix, outSOLN and inRHS are 
//    # vectors of length N.
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
