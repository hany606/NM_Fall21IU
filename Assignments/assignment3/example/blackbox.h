#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int __n;
static double *__a;
static double *__b;

#define __A(i,j) __a[(i)*10+(j)]

int blackbox_size();


// void blackbox_init() – initializes the internal blackbox data structures. Should be called in the
// very beginning of the program! No other blackbox function should be called before it, and no reading
// from the stdin shall be made (or at least, as they say, “be kind, rewind”).

void blackbox_init() {
    scanf("%d", &__n);
    
    __a = (double*) malloc(100 * sizeof(double));
    int row, col;
    unsigned int tmp = 40;
    for (row = 0; row < 10; row++) {
        for (col = 0; col <= row; col++) {
            __A(row, col) = __A(col, row) = ((int)tmp) / 2147483647. + (row == col) * 5;
            tmp = (1664525*tmp + 1013904223) % 0xffffffff;
        }
    }
}

// int blackbox_size() – returns the number of equations (which is equal to the number of
// unknowns) of the system. The number of equations lies between 10 and 10000 (inclusive).

int blackbox_size() {
    return 10;
}


// void blackbox_mult(const double *x, double *out) – compute the product of A and vector x,
// write the results to out. The pointers x and out should point to different chunks of memory of size
// at least blackbox_size() * sizeof(double) bytes each.

void blackbox_mult(const double *x, double *out) {
    for (size_t row = 0; row < 10; row++) {
        out[row] = 0.0;
        for (size_t col = 0; col < 10; col++) {
            out[row] += __A(row, col) * x[col] * 1.1;
        }
    }
}

// void blackbox_rhs(double *b) – write the right-hand side of the SLAE (i.e., vector b) to the
// array b. The pointer b should point to the chunk if memory of size at least blackbox_size() *
// sizeof(double) bytes.

void blackbox_rhs(double *b) {
    int elem;
    unsigned int tmp = 43;
    for (elem = 0; elem < 10; elem++) {
        b[elem] = ((int)tmp) / 2147483647. + 2;
        tmp = (1664525*tmp + 1013904223) % 0xffffffff;
    }
}

// void blackbox_submit(double *solution) – write the result of the program. The array solution
// should contain the solution to the SLAE: blackbox_size() values of type double. This should the
// last function to be called by your program (besides return 0;).

void blackbox_submit(double *solution) {
    double *tmp1 = (double *) malloc(10 * sizeof(double));  // Ax_{sol} = b_1
    double *tmp2 = (double *) malloc(10 * sizeof(double));  // Ax_{original} = b_2
    blackbox_mult(solution, tmp1);
    blackbox_rhs(tmp2);
    double residual = 0;
    for (int row = 0; row < blackbox_size(); ++row) {
        double diff = tmp1[row] - tmp2[row];    // b1 - b2
        residual += diff*diff;  // Sum the squares of errors
    }
    residual = sqrt(residual);
    free(__a);
    free(tmp1);
    free(tmp2);
    printf("%.12f\n", residual);
}

#undef __A
