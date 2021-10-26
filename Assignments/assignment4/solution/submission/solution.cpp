// Newton's method
// Source: 
// https://personal.math.ubc.ca/~pwalls/math-python/roots-optimization/newton/
// https://www.lakeheadu.ca/sites/default/files/uploads/77/docs/RemaniFinal.pdf

// Gradient matrix (Jacobian)  J: 4xn -> s.t. n is the number of the satalities that we recieved their data
// for n = 4 -> J: 4x4, we will just calculate the inverse
// otherwise; n > 4 -> J: 4xn, we will calculate the pseudo inverse (J^T@J)^{-1}@J^T
//      Thus, we will need inverse of (J^T@J)

// Method:
// F = [f_1; f_2; f_3; f_4; f_5, ...., f_n].T
// \phi = \sum_{i=0}^{n} r_i(x,y,z,t)       (Lecture notations)
// r_i(x,y,z,t) = f_i
// f_i = (x-x_i)^2 + (y-y_i)^2 + (z-z_i)^2 - (t-t_i)^2 = 0 
// Convert the problem to optimization min_{x,y,z,t} (\phi) -> this means we will satisfy the equations and get the values of x,y,z,t that will make the summation at least as possible (satisfy the equality)
// Solve it using Gradient based optimization method for unconstrained optimization, we might add a constraint on time to be positive but we can change it back to unconstrained by adding this constraint as extra cost (Lagrange multiplier)
// Moreover, we can solve using simple steepest gradient descent or Newton's method
// 1. Get the jacobian matrix (1st order derivative) (nx4) such that the row: [\partial f_i / \partial x, \partial f_i / \partial y, \partial f_i / \partial z, \partial f_i / \partial t] \forall i \in {0..n}
// 2. Get the inverse of the jacobian using pseudo inverse as it is not squared matrix
// 3. Multiple the inverse of the jacobian with F   (Update step)
// 4. Subtract the update step vector from the previous iteration result
// And repeat these steps for each iteration till convergence with error margin for each new incoming data


#include <stdio.h>
#include <math.h>
#include "blackbox.h"


#define N 30
#define alpha 1//0.5

static const double eps = 1e-7;

double calc_error(double *xi, double *yi, double *zi, double *ti,
                  double x, double y, double z, double t,
                  int n) {
    double sum = 0.0;
    for(int i = 0; i < n; i++){
        double error = (x-xi[i])*(x-xi[i]) + (y-yi[i])*(y-yi[i]) + (z-zi[i])*(z-zi[i]) - ((t-ti[i])*(t-ti[i]));
        sum += error * error;
    }
    return sum;
}

// https://stackoverflow.com/questions/33058848/generate-a-random-double-between-1-and-1/33058967
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}


// ------------------------------------------------------------------------------
// C code to find Moore-Penrose inverse  matrix
// Modified from https://fractalytics.io/moore-penrose-matrix-optimization-cuda-c

#define N_ 5
#define M_ 4 
 
void Transpose(double *matrix, double *t_matrix, int n) {
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<M_; j++) {
             t_matrix[j * n + i]= matrix[i * M_ + j];
        }
        // printf("\n");
    }

    return;
 
}
 
 
 
void MatrixMult(double *matrix_1, double *matrix_2, double *matrix_product, int n1, int n2, int n3) {
    int k;
    for (int i = 0; i<n1; i++) {
        for (int j = 0; j<n3; j++) {
            matrix_product[i * n1 + j] = 0;
            for (k = 0; k<n2; k++) {
                matrix_product[i * n1 + j] += matrix_1[i * n2 + k] * matrix_2[k * n3 + j];
            }
        }
    }
     return;
 
}


// Just a quick solution for the bug
void tmpMatrixMult(double *matrix_1, double *matrix_2, double *matrix_product, int n1, int n2) {

    for(int k = 0; k < n1; k++){
        matrix_product[k] = 0;
            for (int j = 0; j<n2; j++) {
                matrix_product[k] += matrix_1[k*n2+j]* matrix_2[j];
            }

    }

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
        // printf("Singular matrix, can't find its inverse");

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

void displayN_1_(double *A, int n){
    for (int i=0; i<n; i++)
    {
        printf("%.12f ", A[i]);
    }
    printf("\n");

}
// ------------------------------------------------------------------------------

void compute_grad(double *xi, double *yi, double *zi, double *ti,
                  double x, double y, double z, double t,
                  int n, double *grad_matrix)
{
    for(int i = 0; i < n; i++){
        grad_matrix[i*n+0] = 2*(x-xi[i]);
        grad_matrix[i*n+1] = 2*(y-yi[i]);
        grad_matrix[i*n+2] = 2*(z-zi[i]);
        grad_matrix[i*n+3] = -2*(t-ti[i]);
    }
}

void compute_function(double *xi, double *yi, double *zi, double *ti,
                  double x, double y, double z, double t,
                  int n, double *function_matrix)
{
    for(int i = 0; i < n; i++){
        double error = (x-xi[i])*(x-xi[i]) + (y-yi[i])*(y-yi[i]) + (z-zi[i])*(z-zi[i]) - ((t-ti[i])*(t-ti[i]));
        function_matrix[i] = error;
    }
}

int main() {
    gps_init();
    double x = randfrom(0,1), y = randfrom(0,1), z = randfrom(0,1), t = randfrom(0,1);
    double xi[N], yi[N], zi[N], ti[N];

    double *matrix_mult= (double *) malloc(M_*M_*sizeof(double)); // To store 4xn @ nx4
    double *adj= (double *) malloc(M_*M_*sizeof(double)); // To store adjoint 
    double *inv= (double *) malloc(M_*M_*sizeof(double)); // To store inverse 
    double *update_matrix= (double *) malloc(M_*sizeof(double)); // To store 4xn @ nx4

    while (true) {
        const int n = gps_read(xi, yi, zi, ti);
        if (n < 4) break; // If gps_read returns any value less than four, it means that we have lost the signal, and the program should quit (return 0;).
        // if(n==0) break;
        double *t_matrix= (double *) malloc(M_*n*sizeof(double));  // To store the transpose
        double *pseudoinverse= (double *) malloc(M_*n*sizeof(double)); // To store the pseudoinverse
        double *grad_matrix = (double *) malloc(n*M_*sizeof(double));
        double *function_matrix = (double *) malloc(n*sizeof(double));


        int steps = 0;
        while (calc_error(xi, yi, zi, ti, x,y,z,t, n) > eps){
            compute_grad(xi, yi, zi, ti, x,y,z,t, n, grad_matrix);

            Transpose(grad_matrix, t_matrix, n);
            MatrixMult(t_matrix, grad_matrix, matrix_mult, M_, n, M_);
            if (inverse(matrix_mult, inv) > 0){    
                MatrixMult(inv,t_matrix,pseudoinverse, M_, M_, n);
                compute_function(xi, yi, zi, ti, x,y,z,t, n, function_matrix);
                tmpMatrixMult(pseudoinverse, function_matrix, update_matrix, M_, n);

                x = x - alpha*update_matrix[0];
                y = y - alpha*update_matrix[1];
                z = z - alpha*update_matrix[2];
                t = t - alpha*update_matrix[3];
            }

            steps +=1;
            if(steps > 10)
               break;
        }
        // printf("%d %.12f\n", steps, calc_error(xi, yi, zi, ti, x,y,z,t, n));
        free(t_matrix);
        free(pseudoinverse);
        free(grad_matrix);
        gps_submit(x, y, z, t);
    }
    free(update_matrix);
    free(inv);
    free(adj);
    free(matrix_mult);
    return 0;
}

