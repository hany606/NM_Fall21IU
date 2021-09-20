// Simpson's 3/8 Rule
#include <stdio.h>
#include "blackbox.h"
#include <cmath>        // std::abs

static double eps = 1e-9;

int main() {
    int n;
    scanf("%d", &n);
    blackbox_init(n);

    const double N = 1e2;
    const int b = 1;
    const int a = -1;
    const double dx = b-a;

    double sum = 0;
    const double xi_1 = a;

    for (size_t i = 1; i < N; i++) {
        double xi = a + (b-a)/N * i;
        double t1 = blackbox(xi_1)/8.0;
        double t2 = 3*blackbox((2*xi_1+xi)/3.0)/8.0;
        double t3 = 3*blackbox((xi_1+2*xi)/3.0)/8.0;
        double t4 = blackbox(xi)/8.0;
    
        double s = (t1+t2+t3+t4)*(b-a)/N;

        sum += std::abs(s); 
        printf("%.12lf\t%.12lf\n", xi, sum);
    }

    printf("%.12lf\n", sum);
}

