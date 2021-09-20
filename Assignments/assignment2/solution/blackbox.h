#include <math.h>
#include <cmath>


static int __n_method;

void blackbox_init(int n) {
    __n_method = n;
}

double blackbox(double x) {
    switch (__n_method) {
        case 1:
            return 2*x + 1;
        case 2:
            return std::cos(x);
    }
}

double blackbox_df(int k) {
    switch(__n_method) {
        case 1: // 2*x + 1
            if (k==1)
                return 2;
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



