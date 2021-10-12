
// Solving System of Linear Algebric equations using simple iteration method:
// Here, we don't have direct access to A matrix, neither x vector, we have only a multiplaction of A with our vector
// x_{k+1} = x_{k} - tau*A@x_{k} + tau*b
// tau = 2/(L+l) -> but we don't have the eigen values

// Minimal residual method:
// x_{k+1} = x_{k} - tau_{k}(A@x_{k} - b)
// r_{k+1} = r_{k} - tau_{k}(A@r_{k})
// tau_{k} = (A@r_{k}, r_{k}) / (A@r_{k}, A@r_{k})


// Note: 
//      '@' means matrix multiplication, 'b' is 'f' in the slides
//      '(a,b)' means that it is equal to a dot product of 'a' and 'b' vectors

#include <math.h>
#include <stdlib.h>
#include "blackbox.h"

static const double eps = 1e-12;

double calculateTau(double *rk, int n, double *Ark) {

    blackbox_mult(rk, Ark); // Get A@r_{k}

    double numerator = 0;
    double denominator = 0;
    // Calculate the dot proudcts -> dot prouct = sum of element-wise multiplication
    for (int i = 0; i < n; i++) {
        numerator += Ark[i] * rk[i];    
        denominator += Ark[i] * Ark[i];
    }
    double tau = numerator / denominator;
    return tau;
}

double residual(const double *b, const double *x, const int n, double *res) {
    double ret = 0;
    blackbox_mult(x, res);
    for (int r = 0; r < n; ++r) {
        res[r] -= b[r];
        ret += res[r] * res[r];
    }
    return sqrt(ret);
}

int main() {
    blackbox_init();
    
    // A@x_{original_sol} = b
    // x_{current_sol} = x
    // x_{original_sol} is not accessable
    // A is not accessable

    // Initialize the variables
    const int n = blackbox_size();                        // number of equations -> number of unknowns
    double *b  = (double *) malloc(n*sizeof(double));     // Store the result of multiplication of the original x with A
    double *xk  = (double *) malloc(n*sizeof(double));    // Store the final solutions (x_{k})
    double *rk = (double *) malloc(n*sizeof(double));     // Store the residual
    double *d  = (double *) malloc(n*sizeof(double));     // Store the result of multiplication of the solution with A and calculate error
    double *Ark = (double *) malloc(n*sizeof(double));    // Store A@r_{k}
    double *Axk  = (double *) malloc(n*sizeof(double));   // Store A@x_{k}

    double alpha = 0.8; // just for intialization of r_{k}
    double tau = 1; // initial value for tau

    // Get the RHS of the equation Ax=b -> Getting b
    blackbox_rhs(b);

    // Initialize the vectors of rk and xk
    for (int i = 0; i < n; i++){
        xk[i] = 0;
        // r_{k} = A@x_{k} - b , as A@x_{k} = 0, then:
        rk[i] = -b[i]; // -alpha * b[i] both of them works fine  // Also, intializing it with 0 is bad as it makes the denominator = 0
    }


    int counter = 0;
    // Loop until the residual error is less than eps
    while (residual(b, xk, n, d) > eps) {
        double convergence_error = 0;
        blackbox_mult(xk, Axk); // Get A@x_{k} = d
        for (int i = 0; i < n; i++) {
            // x_{k+1} = x_{k} - tau_{k}(A@x_{k} - b)
            xk[i] = xk[i] + tau * (b[i] - Axk[i]);
            // r_{k+1} = r_{k} - tau_{k}(A@r_{k})
            rk[i] = rk[i] - tau*Ark[i];

        }
        // tau_{k} = (A@r_{k}, r_{k}) / (A@r_{k}, A@r_{k})
        tau = calculateTau(rk, n, Ark);
        counter += 1;
    }
    // printf("%d\n", counter);

    blackbox_submit(xk);
    return 0;
}

