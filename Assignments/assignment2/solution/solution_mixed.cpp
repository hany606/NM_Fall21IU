// Multiple Segments for Simpson 3/8 Rule (Mixed Simpson (Composite) 1/3 and 3/8)
// Source: 
// http://mathforcollege.com/nm/mws/gen/07int/mws_gen_int_txt_simpson3by8.pdf
// https://www.freecodecamp.org/news/simpsons-rule/

#include <stdio.h>
#include "blackbox.h"
#include <cmath>        // std::abs


static double eps = 1e-9;

int main() {
    int c;
    scanf("%d", &c);
    blackbox_init(c);

    const size_t n1 = 2e5;
    const size_t n2 = 3e5;
    
    const double b = 1;
    const double a = -1;

    const size_t n = n1+n2;

    double x[n+1];

    const double h = (b-a) / (n);

    double sum = 0;
    const double xi_1 = a;
    for (size_t i = 0; i <= n; i++)
        x[i] = a+i*h;
    
    double I1 = blackbox(x[0]) + blackbox(x[n1]);
    for (size_t i = 1; i <= n1-1; i+=2)
        I1 += 4*(blackbox(x[i]));

    for (size_t i = 2; i <= n1-2; i+=2)
        I1 += 2*(blackbox(x[i]));

    I1 *= (h/3.0);

    double I2 = blackbox(x[0]) + blackbox(x[n2]);
    for (size_t i = 1; i <= n2-2; i+=3)
        I2 += 3*(blackbox(x[i]));

    for (size_t i = 2; i <= n2-1; i+=3)
        I2 += 3*(blackbox(x[i]));

    for (size_t i = 3; i <= n2-3; i+=3)
        I2 += 2*(blackbox(x[i]));


    I2 *= (3.0*h/8.0);

    sum = I1 + I2;
    // printf("%.12lf\n", I1);
    // printf("%.12lf\n", I2);

    printf("%.12lf\n", sum);
}

