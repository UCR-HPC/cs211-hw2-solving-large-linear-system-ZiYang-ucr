#ifndef __MY_C__
#define __MY_C__

// #include "include.h"

int mydgetrf(double *A,int *ipiv,int n)
{
    //TODO
    //The return value (an integer) can be 0 or 1
    //If 0, the matrix is irreducible and the result will be ignored
    //If 1, the result is valid

 int i, j, k, maxind;
    double max, temp;

    for (i = 0; i < n; i++) {
        ipiv[i] = i;
    }

    for (i = 0; i < n - 1; i++) {
        maxind = i;
        max = fabs(A[i * n + i]);
        for (j = i + 1; j < n; j++) {
            if (fabs(A[j * n + i]) > max) {
                max = fabs(A[j * n + i]);
                maxind = j;
            }
        }

        if (max == 0) {
            printf("LU factorization failed: matrix is singular\n");
            return 0;
        }

       
        if (maxind != i) {
       
            int temp_piv = ipiv[i];
            ipiv[i] = ipiv[maxind];
            ipiv[maxind] = temp_piv;

            
            for (k = 0; k < n; k++) {
                temp = A[i * n + k];
                A[i * n + k] = A[maxind * n + k];
                A[maxind * n + k] = temp;
            }
        }

        for (j = i + 1; j < n; j++) {
            A[j * n + i] /= A[i * n + i];  

            for (k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];  
            }
        }
    }

    return 1;
}

void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv) 
{
    double *z = (double *)malloc(n * sizeof(double));

    if (UPLO == 'L') {
        // Forward substitution (solve Ly = B)
        for (int i = 0; i < n; i++) {
            z[i] = B[ipiv[i]];  // Apply the pivoting to B
            for (int j = 0; j < i; j++) {
                z[i] -= A[i * n + j] * z[j];
            }
        }

        // Copy the result back to B
        for (int i = 0; i < n; i++) {
            B[i] = z[i];
        }
    } 
    else if (UPLO == 'U') {
        // Backward substitution (solve Ux = y, where y is the result from forward substitution)
        for (int i = n - 1; i >= 0; i--) {
            for (int j = i + 1; j < n; j++) {
                B[i] -= A[i * n + j] * B[j];
            }
            B[i] /= A[i * n + i];
        }
    } 
    else {
        printf("Invalid UPLO value. Use 'L' for lower triangular (forward substitution) or 'U' for upper triangular (backward substitution).\n");
    }

    free(z);
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