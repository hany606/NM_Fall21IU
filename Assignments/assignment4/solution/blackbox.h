#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int __time, __method;

void gps_init() {
    scanf("%d", &__method);
    __time = 0;
}

double _ti(double x, double y, double z, double t, double xi, double yi, double zi) {
    const double dist2 = (x-xi)*(x-xi) + (y-yi)*(y-yi) + (z-zi)*(z-zi);
    return t - sqrt(dist2);
}

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

void gps_submit(double x, double y, double z, double t) {
    printf("%.12f %.12f %.12f %.12f\n", x, y, z, t);
}

