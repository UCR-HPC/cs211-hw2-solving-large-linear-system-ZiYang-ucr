#ifndef __MY_C__
#define __MY_C__

#include "include.h"

int mydgetrf(double *A,int *ipiv,int n)
{
    //TODO
    //The return value (an integer) can be 0 or 1
    //If 0, the matrix is irreducible and the result will be ignored
    //If 1, the result is valid

     int i, j, k, maxind;
    double max, temp;
    
    for (i = 0; i < n - 1; i++) {
        // Pivoting step
        maxind = i;
        max = fabs(A[i * n + i]);
        
        for (j = i + 1; j < n; j++) {
            if (fabs(A[j * n + i]) > max) {
                max = fabs(A[j * n + i]);
                maxind = j;
            }
        }

        // Check if matrix is singular
        if (max == 0) {
            printf("LU factorization failed: coefficient matrix is singular\n");
            return;
        }

        // Swap rows if necessary
        if (maxind != i) {
            int temp_pvt = pvt[i];
            pvt[i] = pvt[maxind];
            pvt[maxind] = temp_pvt;

            for (k = 0; k < n; k++) {
                temp = A[i * n + k];
                A[i * n + k] = A[maxind * n + k];
                A[maxind * n + k] = temp;
            }
        }

        // Factorization step
        for (j = i + 1; j < n; j++) {
            A[j * n + i] /= A[i * n + i];  // Divide by pivot element

            for (k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];  // Update elements
            }
        }
    }
    return 1;
}

void mydtrsv(char UPLO,double *A,double *B,int n,int *ipiv)
{
    //TODO

    int i, j;

    // Forward substitution (Solving Ly = b)
    for (i = 0; i < n; i++) {
        x[i] = b[pvt[i]];

        for (j = 0; j < i; j++) {
            x[i] -= A[i * n + j] * x[j];
        }
    }

    // Backward substitution (Solving Ux = y)
    for (i = n - 1; i >= 0; i--) {
        for (j = i + 1; j < n; j++) {
            x[i] -= A[i * n + j] * x[j];
        }
        x[i] /= A[i * n + i];
    }
}

void my_f(double *A,double *B,int n)
{
    int *ipiv=(int*)malloc(n*sizeof(int));
    for (int i=0;i<n;i++)
        ipiv[i]=i;
    if (mydgetrf(A,ipiv,n)==0) 
    {
        printf("LU factoration failed: coefficient matrix is singular.\n");
        return;
    }
    mydtrsv('L',A,B,n,ipiv);
    mydtrsv('U',A,B,n,ipiv);
}

#endif