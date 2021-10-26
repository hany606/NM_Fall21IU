// Source: https://fractalytics.io/moore-penrose-matrix-optimization-cuda-c
// Confirmation on results: https://www.wolframalpha.com/input/?i=pseudoinverse+%5B+%5B5%2C+-2%2C+2%2C+7%2C9+%5D%2C++%5B1%2C+0%2C+0%2C+3%2C1+%5D%2C++%5B-3%2C+1%2C+5%2C+0%2C2+%5D%2C++%5B3%2C+-1%2C+-9%2C+4%2C6+%5D%2C++%5B1%2C+0%2C+4%2C+4%2C1+%5D%5D
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// C program to find Moore-Penrose inverse  matrix
 
 
#define N_ 5
#define M_ 4 
 
void Transpose(double *matrix, double *t_matrix, int n) {
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<M_; j++) {
             t_matrix[j * n + i]= matrix[i * M_ + j];
        }
        printf("\n");
    }

    return;
 
}
 
 
 
void MatrixMult(double *matrix_1, double *matrix_2, double *matrix_product, int n1, int n2, int n3) {
    int k;
    for (int i = 0; i<n1; i++) {
        for (int j = 0; j<n3; j++) {             // not j<M
            matrix_product[i * n1 + j] = 0;
            for (k = 0; k<n2; k++) {
                matrix_product[i * n1 + j] += matrix_1[i * n2 + k] * matrix_2[k * n3 + j];
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
                temp[i * M_ + j++] = A[row * M_ + col];
 
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
 
    double temp[M_*M_]; // To store cofactors
 
    int sign = 1; // To store sign multiplier
 
    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0 * M_ + f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}
 
// Function to get adjoint
void adjoint(double *A,double *adj)
{
    if (M_ == 1)
    {
        adj[0] = 1;
        return;
    }
 
    // temp is used to store cofactors 
    int sign = 1;
    double temp[M_*M_];
 
    for (int i=0; i<M_; i++)
    {
        for (int j=0; j<M_; j++)
        {
            // Get cofactor
            getCofactor(A, temp, i, j, M_);
 
            // sign of adj positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j * M_ + i] = (sign)*(determinant(temp, M_-1));
        }
    }
}
 
// Function to calculate and store inverse, returns false if
// matrix is singular
int inverse(double *A, double *inverse)
{
    // Find determinant of A[][]
    int det = determinant(A, M_); 
    if (det == 0)
    {
        printf("Singular matrix, can't find its inverse");

        return 0;
    }
 
    // Find adjoint
    double adj[M_*M_];
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<M_; i++)
        for (int j=0; j<M_; j++)
            inverse[i * M_ + j]= adj[i * M_ + j]/((double) det);
 
    return 1;
}

// Generic function to display the matrix. We use it to display
// both adjoin and inverse. adjoin is integer matrix and inverse
// is a double.

void displayM_(double *A){
    for (int i=0; i<M_; i++)
    {
        for (int j=0; j<M_; j++)
            printf("%.12f ", A[i * M_ + j]);
        printf("\n");
    }
}
void displayN_M_(double *A, int n){
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<M_; j++)
            printf("%.12f ", A[i * M_ + j]);
        printf("\n");
    }
}

void displayM_N_(double *A, int n){
    for (int i=0; i<M_; i++)
    {
        for (int j=0; j<n; j++)
            printf("%.12f ", A[i * M_ + j]);
        printf("\n");
    }
}
 
// Driver program
int main()
{
    const int n = N_;
    double matrix[n][M_] = { {5, -2, 2, 7},
                              {1,  0, 0, 3},
                              {-3, 1, 5, 0},
                              {3, -1,-9, 4},
                              {1,  0, 4, 4},
                            // {8, 0, 3, 8,6,5,2},
                            // {5, 6, 4, 1,3,2,0}
            };
 
    // double *matrix = (double *) malloc(N*N*sizeof(double));
    double *t_matrix= (double *) malloc(M_*n*sizeof(double));  // To store the transpose
    double *matrix_mult= (double *) malloc(M_*M_*sizeof(double)); // To store 4xn @ nx4
    double *pseudoinverse= (double *) malloc(M_*n*sizeof(double)); // To store the pseudoinverse
    double *adj= (double *) malloc(M_*M_*sizeof(double)); // To store adjoint 
    double *inv= (double *) malloc(M_*M_*sizeof(double)); // To store inverse 
 
    for(int i = 0; i < 1; i++){

        displayN_M_(*matrix, n);

        Transpose(*matrix, t_matrix, n);
        printf("\nThe Transpose is :\n");

        displayM_N_(t_matrix, n);
    
        printf("The product of the matrix is: \n");

        MatrixMult(t_matrix, *matrix, matrix_mult, M_, n, M_);
        displayM_(matrix_mult);

        printf("\nThe Inverse is :\n");

        if (inverse(matrix_mult, inv))
            displayM_(inv);
    
        MatrixMult(inv,t_matrix,pseudoinverse, M_, M_, n);
        
        printf("\nThe Monroe-penrose inverse is :\n");
        
        displayM_N_(pseudoinverse, n);
    }
 
    return 0;
}