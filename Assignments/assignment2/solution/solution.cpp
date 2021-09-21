// Composite Simpson 3/8 Rule
// Sources: 
// 3rd Lecture
// http://mathforcollege.com/nm/mws/gen/07int/mws_gen_int_txt_simpson3by8.pdf
// https://www.freecodecamp.org/news/simpsons-rule/
// https://www.wikiwand.com/en/Simpson%27s_rule
// https://opentextbc.ca/calculusv2openstax/chapter/numerical-integration/
// https://sites.und.edu/timothy.prescott/apex/web/apex.Ch8.S7.html#SSx3

#include <stdio.h>
#include "blackbox.h"


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

    // Error calculation
    // double mx_val = blackbox_df(4);
    // double error = -h*h*h*h/80*(b-a)*mx_val;    // Source: https://www.wikiwand.com/en/Simpson%27s_rule#/Composite_Simpson's_3/8_rule
    // printf("Error: %.12lf\n", error);

 

}


// Used for testing:
// double blackbox(double x) {
//     switch (__n_method) {
//         case 1: //2.000000000000
//             return 2*x + 1;
//         case 2: //1.682941969616
//             return std::cos(x);
//         case 3: //0.000000000000
//             return std::sin(x);
//         case 4: //0.000000000000
//             return std::tan(x);
//         case 5: //0.666666666667
//             return x*x;
//         case 6: //2.350402387288
//             return std::exp(x);
//         case 7: //0.562774800352
//             return std::log10(x+2);
//         case 8: //1.885618082563
//             return std::sqrt(x+1);
//         case 9: //0.000000000000
//             return x*std::tan(x)*std::tan(x);
//         case 10: //29.681284231115
//             return std::exp(-x*5);
//         case 11: //1.493648265625
//             return std::exp(-x*x);
//         case 12: //0.0
//             return std::sin(x*x*x);
//         case 13: //0.620536603447
//             return std::sin(x*x);
//         case 14: //0.181818181818
//             return std::pow(x,10);

//     }
// }