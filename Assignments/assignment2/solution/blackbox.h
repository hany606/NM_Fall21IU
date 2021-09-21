#include <math.h>
#include <cmath>


static int __n_method;

void blackbox_init(int n) {
    __n_method = n;
}

double blackbox(double x) {
    switch (__n_method) {
        case 1: //2.000000000000
            return 2*x + 1;
        case 2: //1.682941969616
            return std::cos(x);
        case 3: //0.000000000000
            return std::sin(x);
        case 4: //0.000000000000
            return std::tan(x);
        case 5: //0.666666666667
            return x*x;
        case 6: //2.350402387288
            return std::exp(x);
        case 7: //0.562774800352
            return std::log10(x+2);
        case 8: //1.885618082563
            return std::sqrt(x+1);
        case 9: //0.000000000000
            return x*std::tan(x)*std::tan(x);
        case 10: //29.681284231115
            return std::exp(-x*5);
        case 11: //1.493648265625
            return std::exp(-x*x);
        case 12: //0.0
            return std::sin(x*x*x);
        case 13: //0.620536603447
            return std::sin(x*x);
        case 14:
            return std::pow(x,10);

    }
}

double blackbox_df(int k) {
    switch(__n_method) {
        case 1: // 2*x + 1
            if (k==1)
                return 2;
            else
                return 0;
        case 14: // 2*x + 1
            if (k==4)
                return 5040;
            else
                return 0;
    }

}

double blackbox_period() {
    switch(__n_method) {
        case 1: // 2*x + 1
            return 0;
    }
}



