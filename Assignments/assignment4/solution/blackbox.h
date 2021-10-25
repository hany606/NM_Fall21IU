#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int __time, __method;

// initializes the internal blackbox data structures. Should be called at the very
// beginning of the program! No other blackbox function should be called before it, and no reading
// from the stdin shall occur.
void gps_init() {
    scanf("%d", &__method);
    __time = 0;
}

double _ti(double x, double y, double z, double t, double xi, double yi, double zi) {
    const double dist2 = (x-xi)*(x-xi) + (y-yi)*(y-yi) + (z-zi)*(z-zi);
    return t - sqrt(dist2);
}


// reads data from the receiver.
// The data is written to the arrays x, y, z, t, and the function returns the number of satellites for
// which we have the data. The arrays x, y, z, t must be allocated by you, and contain at least 30
// elements each. For example, if the function returns 5, that means that first five elements of x now
// contain x coordinates of the five satellites, etc. If gps_read returns any value less than four, it means
// that we have lost the signal, and the program should quit (return 0;).
int gps_read(double *x, double *y, double *z, double *t) {
    if (__method == 1) {
        if (__time > 20)
            return 0;
        __time += 1;
        x[0]=1;
        y[0]=0;
        z[0]=0;
        x[1]=0;
        y[1]=1;
        z[1]=0;
        x[2]=0;
        y[2]=0;
        z[2]=1;
        x[3]=1;
        y[3]=1;
        z[3]=1;
        for (int i = 0; i < 4; i++) {
            t[i] = _ti(0.5, 0.5, 0.5, __time/10., x[i], y[i], z[i]);
        }
        return 4;
    }
    return 0;
}


// accepts calculated coordinates
// of the receiver. Must be called after every successful (return value ≥ 4) call to gps_read.
void gps_submit(double x, double y, double z, double t) {
    printf("%.12f %.12f %.12f %.12f\n", x, y, z, t);
}

