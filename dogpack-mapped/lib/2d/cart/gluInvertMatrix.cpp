# include <stdio.h>
#include <stdlib.h>"
# include <math.h>
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;
typedef double ** Matrix;
typedef double * Row;
typedef double * Col;
typedef double Elem;
 
Matrix allocate_matrix(int n);
Col allocate_col(int n);
Row allocate_row(int n);
void free_matrix(Matrix M, int n);
 
void pivot_partial(Matrix A, Col S,Col B, int n);
void forward_elimination(Matrix A,Col B,int n);
Col back_substitution(Matrix A, Col B, int n);
Col scale_factor(Matrix A,int n);
void gauss(Matrix A, Col B,vector<double>& x, int n);
 
void swap_rows(Row *r1, Row*r2);
void print_matrix(Matrix M, int n, char * name);
void print_col(Col C, int n, char *name);
void print_row(Row R, int n, char *name);
 
 
 
bool gluInvertMatrix(vector<double>& m, vector<double>& b,vector<double>& x)
{
    int i,j;
    Matrix A;
    Col B;

    int n=b.size();
 
    A = allocate_matrix(n);
    for( i = 1; i <= n; ++i)
        for(j = 1; j <= n; ++j)
                A[i][j]=m[(i-1)*n+j-1];
    //print_matrix(A,n,"A");
    B = allocate_col(n);   
     
    for(j = 1; j <= n; ++j)
        B[j]=b[j-1];
        
    gauss(A,B,x,n);  
     


   double tmp=0.0;
for (i=0;i<=n-1;i++){double tmp3=0.0;
        for(j=0;j<n;j++){
           tmp3=tmp3+m[i*n+j]*x[j];     
        }
        double bmin=1.0;
        //if(abs(b[i])>1e-8){bmin=abs(b[i]);}
        //else{bmin=1.0;}
        tmp=max(tmp,abs(tmp3-b[i])/bmin);
}

if(tmp<1.0e-14){
    free_matrix(A,n);
    free(B + 1);
return 1;}
else{
cout<<"problem "<<tmp<<endl;
tmp=0.0;
for (i=0;i<=n-1;i++){double tmp3=0.0;
        for(j=0;j<n;j++){
           tmp3=tmp3+A[i+1][j+1]*x[j];     
        }
        tmp=max(tmp,abs(tmp3-B[i+1]));
        cout<<b[i]<<endl;
}
cout<<"problem2 "<<tmp<<endl;
print_matrix(A,n,"A");
    free_matrix(A,n);
    free(B + 1);
return 0;}
}
 
 
 
 
void print_matrix(Matrix M, int n, char * name)
{
    int i,j;
    printf("\n[%s] = ",name);
    printf("\n\n");
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; ++j)
            printf("%6lG ",M[i][j]);
        printf("\n");
    }
}
 
void print_col(Col C, int n, char * name)
{
    int j;
    printf("\n[%s] = ",name);
    printf("\n\n");
    for(j = 1; j <= n; ++j)
        printf("%6lg\n",C[j]);
             
}
void print_row(Row R, int n, char * name)
{
    int i;
    printf("\n[%s] = ",name);
    for(i = 1; i <= n; ++i)
        printf("%6lg ",R[i]);
    printf("\n");
}
 
Matrix allocate_matrix(int n)
{
    Matrix A;
    int i,j;
    A = (double **) malloc(n * sizeof(Row));
    if(!A)
    {
        printf("\nError : Could not allocate memory for matrix\n");
        exit(1);
    }      
    --A;
         
    for(i = 1; i <= n; ++i)
    {
        A[i] = (double *)malloc(n * sizeof(Elem));
        if(!A[i])
        {
            printf("\nError : Could not allocate memory for matrix\n");
            exit(1);
        }
        --A[i];
    }
    return A;
}
 
void free_matrix(Matrix M, int n)
{
    int i;
    for(i = 1; i <= n; ++i)
            free(M[i] + 1);
    free(M + 1);
}
 
Col allocate_col(int n)
{
    Col B;
     
    B = (double *)malloc(n * sizeof(Elem));
     
    if(!B)
    {
        printf("\nError : could not allocate memory\n");
        exit(1);
    }
    --B;
    return B;
}
 
Row allocate_row(int n)
{
    Row B;
    B = (double *)malloc(n * sizeof(Elem));
    if(!B)
    {
        printf("\nError : could not allocate memory\n");
        exit(1);
    }
    --B;
    return B;
}
 
 
Col scale_factor(Matrix A, int n)
{
    int i,j;
    Col S ;
    S = allocate_col(n);
     
    for(i = 1; i <= n; ++i)
    {
        S[i] = A[i][1];
        for(j = 2; j <= n; ++j)
        {
            if(S[i] < fabs(A[i][j]))
                S[i] = fabs(A[i][j]);
        }
         
    }
    return S;
         
}
 
void pivot_partial(Matrix A, Col S,Col B, int n)
{
    int i,j;
    Elem temp;
    for(j = 1; j <= n; ++j)
    {
        for(i = j + 1; i <= n; ++i)
        {
            if(S[i] == 0)
            {
                if(B[i] == 0)
                    printf("\nSystem doesnt have a unique solution");
                else
                    printf("\nSystem is inconsistent");
                exit(1);
            }
            if(fabs(A[i][j]/S[i]) > fabs(A[j][j]/S[j]))
            {
                swap_rows(&A[i],&A[j]);
                temp = B[i];
                B[i] = B[j];
                B[j] = temp;
            }
        }
         
        if(A[j][j] == 0)
        {
            printf("\nSingular System Detected\n");
            exit(1);
        }
    }
     
}
 
void swap_rows(Row *r1, Row*r2)
{
    Row temp;
    temp = *r1;
    *r1 = *r2;
    *r2 = temp;
}
 
void forward_elimination(Matrix A,Col B,int n)
{
    int i,j,k;
    double m;
     
    for(k = 1; k <= n-1; ++k)
    {
        for(i = k + 1; i <= n; ++i)
        {
            m =  A[i][k] / A[k][k];
            for(j = k + 1; j <= n; ++j)
            {
                A[i][j] -= m * A[k][j];
                if(i == j && A[i][j] == 0)
                {
                    printf("\nSingular system detected");
                    exit(1);
                }                  
            }
            B[i] -= m * B[k];
             
        }
    }
             
}
 
Col back_substitution(Matrix A, Col B, int n)
{
    int i,j;
    Elem sum;
    Col X = allocate_col(n);
    X[n] = B[n]/A[n][n];
    for(i = n - 1; i >= 1; --i)
    {
        sum = 0;
        for(j = i + 1; j <= n; ++j)
            sum += A[i][j] * X[j];
        X[i] = (B[i] - sum) / A[i][i];
    }
    return X;
}
 
void gauss(Matrix A, Col B,vector<double>& x, int n)
{
    int i,j;
    Col S, X;
    S = scale_factor(A,n);
    pivot_partial(A,S,B,n);
    forward_elimination(A,B,n);
    X = back_substitution(A,B,n);
    for(j = 1; j <= n; ++j)
        x[j-1]=X[j];
     
    free(S + 1);
    free(X + 1);
}
