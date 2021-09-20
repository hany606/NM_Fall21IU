// Simpson 3/8 Rule
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

    const size_t n1 = 3e5;
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

    for (size_t i = 1; i <= n/3; i++)
	{
		sum += blackbox(x[3*i - 3]) + 3*blackbox(x[3*i - 2]) + 3*blackbox(x[3*i - 1]) + blackbox(x[3*i]);
	}
	
    sum *= 3.0*h/8.0;
    printf("%.12lf\n", sum);
}

