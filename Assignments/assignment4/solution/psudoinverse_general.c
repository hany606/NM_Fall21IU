// Source: https://fractalytics.io/moore-penrose-matrix-optimization-cuda-c
// Confirmation on results: https://www.wolframalpha.com/input/?i=pseudoinverse+%5B+%5B5%2C+-2%2C+2%2C+7%2C9+%5D%2C++%5B1%2C+0%2C+0%2C+3%2C1+%5D%2C++%5B-3%2C+1%2C+5%2C+0%2C2+%5D%2C++%5B3%2C+-1%2C+-9%2C+4%2C6+%5D%2C++%5B1%2C+0%2C+4%2C+4%2C1+%5D%5D
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// C++ program to find Moore-Penrose inverse  matrix
 
 
#define N 4
 
void Trans_2D_1D(double matrix_2D[N][N], double *matrix) {
 
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            matrix[i * N + j] = matrix_2D[i][j];
        }
        printf("\n");
    }
 
    return;
 
}
 
 
 
 
void Transpose(double *matrix, double *t_matrix) {
    for (int i = 0; i<N; i++) {
        for (int j = 0; j<N; j++) {
             t_matrix[j * N + i]= matrix[i * N + j];
        }
        printf("\n");
    }
 
    return;
 
}
 
 
 
void MatrixMult(double *matrix_1, double *matrix_2, double *matrix_product) {
    int k;
    for (int i = 0; i<N; i++) {
        for (int j = 0; j<N; j++) {             // not j<M
            matrix_product[i * N + j] = 0;
            for (k = 0; k<N; k++) {
                matrix_product[i * N + j] += matrix_1[i * N + k] * matrix_2[k * N + j];
            }
        }
    }
     return;
 
}
 
 
 
// Function to get cofactor 
void getCofactor(double *A, double *temp, int p, int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            // Copying into temporary matrix only those element
            // which are not in given row and column
            if (row != p && col != q)
            {
                temp[i * N + j++] = A[row * N + col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
 
// Recursive function for finding determinant of matrix.
int determinant(double *A, int n)
{
    int D = 0; // Initialize result
 
    // Base case : if matrix contains single element
    if (n == 1)
        return A[0];
 
    double temp[N*N]; // To store cofactors
 
    int sign = 1; // To store sign multiplier
 
    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0 * N + f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}
 
// Function to get adjoint
void adjoint(double *A,double *adj)
{
    if (N == 1)
    {
        adj[0] = 1;
        return;
    }
 
    // temp is used to store cofactors 
    int sign = 1;
    double temp[N*N];
 
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Get cofactor
            getCofactor(A, temp, i, j, N);
 
            // sign of adj positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j * N + i] = (sign)*(determinant(temp, N-1));
        }
    }
}
 
// Function to calculate and store inverse, returns false if
// matrix is singular
int inverse(double *A, double *inverse)
{
    // Find determinant of A[][]
    int det = determinant(A, N); 
    if (det == 0)
    {
        printf("Singular matrix, can't find its inverse");

        return 0;
    }
 
    // Find adjoint
    double adj[N*N];
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            inverse[i * N + j]= adj[i * N + j]/((double) det);
 
    return 1;
}
 
// Generic function to display the matrix. We use it to display
// both adjoin and inverse. adjoin is integer matrix and inverse
// is a double.
void display(double *A)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
            printf("%.12f ", A[i * N + j]);
        printf("\n");
    }
}
 
// Driver program
int main()
{
    double matrix[N][N] = { 
                            // {5, -2, 2, 7,9},
                            // {1, 0, 0, 3,1},
                            // {-3, 1, 5, 0,2},
                            // {3, -1, -9, 4,6},
                            // {1, 0, 4, 4,1},
                            {45.000000000000, -16.000000000000, -28.000000000000, 54.000000000000},
                            {-16.000000000000, 6.000000000000, 10.000000000000, -18.000000000000 },
                            {-28.000000000000, 10.000000000000, 126.000000000000, -6.000000000000 },
                            {54.000000000000, -18.000000000000, -6.000000000000, 90.000000000000 }
                            // {8, 0, 3, 8,6,5,2},
                            // {5, 6, 4, 1,3,2,0}
            };
 
    // double *matrix = (double *) malloc(N*N*sizeof(double));
    double *t_matrix= (double *) malloc(N*N*sizeof(double));
    double *matrix_mult= (double *) malloc(N*N*sizeof(double));
    double *pseudoinverse= (double *) malloc(N*N*sizeof(double));
    double *adj= (double *) malloc(N*N*sizeof(double)); // To store adjoint 
    double *inv= (double *) malloc(N*N*sizeof(double)); // To store inverse 
 
    display(*matrix);

    Transpose(*matrix, t_matrix);
    printf("\nThe Transpose is :\n");

    display(t_matrix);
 
    printf("The product of the matrix is: \n");

    MatrixMult(t_matrix, *matrix, matrix_mult);
    display(matrix_mult);

    printf("\nThe Inverse is :\n");

    if (inverse(matrix_mult, inv))
        display(inv);
 
    MatrixMult(inv,t_matrix,pseudoinverse);
    
    printf("\nThe Monroe-penrose inverse is :\n");
    
    display(pseudoinverse);
 
    return 0;
}