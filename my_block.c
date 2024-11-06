#ifndef __MY_BLOCK_C__
#define __MY_BLOCK_C__

#include "include.h"

void mydgemm(double *A,double *B,int n,int bid,int b)
{
    register int i, j, k;
    int A_row_start = (bid + 1) * b;
    int A_col_start = bid * b;
    int B_row_start = bid * b;
    int B_col_start = (bid + 1) * b;

    int A_row_end = n;
    int A_col_end = (bid + 1) * b;
    int B_row_end = (bid + 1) * b;
    int B_col_end = n;

    for (i = A_row_start; i < A_row_end; i++) {
        for (j = B_col_start; j < B_col_end; j++) {
           register double sum = 0.0;
            for (k = A_col_start; k < A_col_end; k++) {
                sum += A[i * n + k] * B[k * n + j];
            } 
            A[i * n + j] = sum;  // Store result in place within A or a designated output matrix
        }
    }  
}
// void mydgemm(double *A, double *B, int n, int bid, int b) {
//     register int i, j, k;

//     int A_row_start = (bid + 1) * b;
//     int A_col_start = bid * b;
//     int B_row_start = bid * b;
//     int B_col_start = (bid + 1) * b;

//     int A_row_end = (A_row_start + b > n) ? n : A_row_start + b;
//     int A_col_end = (A_col_start + b > n) ? n : A_col_start + b;
//     int B_row_end = (B_row_start + b > n) ? n : B_row_start + b;
//     int B_col_end = (B_col_start + b > n) ? n : B_col_start + b;

//     double C_block[b][b];
//     for (i = 0; i < b; i++) {
//         for (j = 0; j < b; j++) {
//             C_block[i][j] = 0.0;  
//         }
//     }

//     for (i = A_row_start; i < A_row_end; i++) {
//         for (j = B_col_start; j < B_col_end; j++) {
//             double sum = 0.0;
//             double c = 0.0; 
//             for (k = A_col_start; k < A_col_end; k++) {
//                 double product = A[i * n + k] * B[k * n + j];
//                 double y = product - c;
//                 double t = sum + y;
//                 c = (t - sum) - y;
//                 sum = t;
//             }
//             C_block[i - A_row_start][j - B_col_start] = sum;
//         }
//     }

//     for (i = 0; i < (A_row_end - A_row_start); i++) {
//         for (j = 0; j < (B_col_end - B_col_start); j++) {
//             A[(A_row_start + i) * n + (B_col_start + j)] += C_block[i][j];
//         }
//     }
// }
// void mydgemm(double *A, double *B, int n, int bid, int b) {
//   for (int A_row_start = 0; A_row_start < n; A_row_start += b) {
//         int A_row_end = (A_row_start + b > n) ? n : A_row_start + b;
        
//         for (int B_col_start = 0; B_col_start < n; B_col_start += b) {
//             int B_col_end = (B_col_start + b > n) ? n : B_col_start + b;
            
//             for (int A_col_start = 0; A_col_start < n; A_col_start += b) {
//                 int A_col_end = (A_col_start + b > n) ? n : A_col_start + b;

//                 for (int i = A_row_start; i < A_row_end; i++) {
//                     for (int j = B_col_start; j < B_col_end; j++) {
//                         double sum = 0.0;
                        
//                         for (int k = A_col_start; k < A_col_end; k++) {
//                             if (A[i * n + k] < TMP_MAX / B[k * n + j]) {
//                                 sum += A[i * n + k] * B[k * n + j];
//                             }
//                         }
                        
//                         A[i * n + j] += sum;
//                     }
//                 }
//             }
//         }
//     }
// }
// void mydgemm(double *A, double *B, int n, int bid, int b) {
//     register int i, j, k;

//     // Define starting and ending indices for A and B blocks
//     int A_row_start = (bid + 1) * b;
//     int A_col_start = bid * b;
//     int B_row_start = bid * b;
//     int B_col_start = (bid + 1) * b;

//     // Adjust row/column end indices based on matrix size to avoid overflow
//     int A_row_end = (A_row_start + b > n) ? n : A_row_start + b;
//     int A_col_end = (A_col_start + b > n) ? n : A_col_start + b;
//     int B_row_end = (B_row_start + b > n) ? n : B_row_start + b;
//     int B_col_end = (B_col_start + b > n) ? n : B_col_start + b;

//     // Allocate C_block dynamically based on block size b
//     double **C_block = (double **)malloc(b * sizeof(double *));
//     for (i = 0; i < b; i++) {
//         C_block[i] = (double *)malloc(b * sizeof(double));
//         for (j = 0; j < b; j++) {
//             C_block[i][j] = 0.0;  // Initialize to zero
//         }
//     }

//     // Perform the matrix multiplication on the block level
//     for (i = A_row_start; i < A_row_end; i++) {
//         for (j = B_col_start; j < B_col_end; j++) {
//             register double sum = 0.0;
//             for (k = A_col_start; k < A_col_end; k++) {
//                 sum += A[i * n + k] * B[k * n + j];
//             }
//             C_block[i - A_row_start][j - B_col_start] = sum;
//         }
//     }

//     // Copy C_block values back to A (or designated output matrix)
//     for (i = 0; i < (A_row_end - A_row_start); i++) {
//         for (j = 0; j < (B_col_end - B_col_start); j++) {
//             A[(A_row_start + i) * n + (B_col_start + j)] = C_block[i][j];
//         }
//     }

//     // Free dynamically allocated memory
//     for (i = 0; i < b; i++) {
//         free(C_block[i]);
//     }
//     free(C_block);
// }


// void mydgemm(double *A, double *B, int n, int bid, int b) {
//     register int i, j, k;

//     // Define starting and ending indices for A and B blocks
//     int A_row_start = (bid + 1) * b;
//     int A_col_start = bid * b;
//     int B_row_start = bid * b;
//     int B_col_start = (bid + 1) * b;

//     // Adjust row/column end indices based on matrix size to avoid overflow
//     int A_row_end = (A_row_start + b > n) ? n : A_row_start + b;
//     int A_col_end = (A_col_start + b > n) ? n : A_col_start + b;
//     int B_row_end = (B_row_start + b > n) ? n : B_row_start + b;
//     int B_col_end = (B_col_start + b > n) ? n : B_col_start + b;

//     // Allocate C_block dynamically to avoid stack overflow
//       double C_block[32][32] = {0.0};


//     for (i = A_row_start; i < A_row_end; i++) {
//         for (j = B_col_start; j < B_col_end; j++) {
//             register double sum = 0.0;
//             for (k = A_col_start; k < A_col_end; k++) {
//                 sum += A[i * n + k] * B[k * n + j];
//             }
//             C_block[i - A_row_start][j - B_col_start] = sum;
//         }
//     }

//     for (i = 0; i < (A_row_end - A_row_start); i++) {
//         for (j = 0; j < (B_col_end - B_col_start); j++) {
//             A[(A_row_start + i) * n + (B_col_start + j)] = C_block[i][j];
//         }
//     }
//     // // Free dynamically allocated memory
//     // for (i = 0; i < b; i++) {
//     //     free(C_block[i]);
//     // }
   
// }
// int mydgetrf_block(double *A, int *ipiv, int n) {
//     int b = 128;  // Adjust block size based on performance needs
//     for (int i = 0; i < n; i += b) {
//         int end = (i + b > n) ? n : (i + b);

//         // LU Decomposition on the current block
//         for (int j = i; j < end; j++) {
//             int pivot = j;

//             // Find the pivot
//             for (int k = j + 1; k < n; k++) {
//                 if (fabs(A[k * n + j]) > fabs(A[pivot * n + j])) {
//                     pivot = k;
//                 }
//             }

//             // Row swapping for partial pivoting
//             if (pivot != j) {
//                 for (int k = 0; k < n; k++) {
//                     double temp = A[j * n + k];
//                     A[j * n + k] = A[pivot * n + k];
//                     A[pivot * n + k] = temp;
//                 }
//                 int temp = ipiv[j];
//                 ipiv[j] = ipiv[pivot];
//                 ipiv[pivot] = temp;
//             }

//             // LU Decomposition for this block
//             for (int k = j + 1; k < end; k++) {
//                 A[k * n + j] /= A[j * n + j];
//                 for (int l = j + 1; l < end; l++) {
//                     A[k * n + l] -= A[k * n + j] * A[j * n + l];
//                 }
//             }
//         }

//         // Update the trailing submatrix using the block
//         for (int j = i + b; j < n; j += b) {
//             mydgemm(A, A, n, i, b);  // Update submatrix based on mydgemm results
//         }
//     }
//     return 1;
// }
int mydgetrf_block(double *A, int *ipiv, int n) {
    int b = 128;  // Block size, set based on cache size or through testing for best performance

    for (int i = 0; i < n; i += b) {
        int end = (i + b > n) ? n : (i + b);

        // Step 1: LU Decomposition on the current block
        for (int j = i; j < end; j++) {
            int pivot = j;

            // Partial pivoting (optimized)
            double max_val = fabs(A[j * n + j]);
            for (int k = j + 1; k < n; k++) {
              register double temp_val = fabs(A[k * n + j]);
              if (fabs(A[j * n + j]) < 1e-12) {
    printf("Warning: Possible numerical instability at element %d\n", j);
}
                if (temp_val > max_val) {
                    max_val = temp_val;
                    pivot = k;
                }
            }

            // Row swapping for partial pivoting
            if (pivot != j) {
                for (int k = 0; k < n; k++) {
                  register double temp = A[j * n + k];
                    A[j * n + k] = A[pivot * n + k];
                    A[pivot * n + k] = temp;
                }
                register int temp = ipiv[j];
                ipiv[j] = ipiv[pivot];
                ipiv[pivot] = temp;
            }

            // LU Decomposition for this block
            register double reciprocal = 1.0 / A[j * n + j]; 
            for (int k = j + 1; k < end; k++) {
                A[k * n + j] *= reciprocal;  // Optimized division
                register double lu_factor = A[k * n + j];
                for (int l = j + 1; l < end; l++) {
                    A[k * n + l] -= lu_factor * A[j * n + l];
                }
            }
        }
        

        // Step 2: Update the trailing submatrix using the block
        for (int j = i + b; j < n; j += b) {
            mydgemm(A, A, n, i, b);  // Update submatrix based on mydgemm results
        }
    }
    return 1;
}


void my_block_f(double *A,double *B,int n)
{
    int *ipiv=(int*)malloc(n*sizeof(int));
    for (int i=0;i<n;i++)
        ipiv[i]=i;
    if (mydgetrf_block(A,ipiv,n)==0) 
    {
        printf("LU factoration failed: coefficient matrix is singular.\n");
        return;
    }
    mydtrsv('L',A,B,n,ipiv);
    mydtrsv('U',A,B,n,ipiv);
}

#endif