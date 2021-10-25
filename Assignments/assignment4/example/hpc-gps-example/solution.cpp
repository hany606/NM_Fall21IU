#include <stdio.h>
#include <math.h>
#include "blackbox.h"

#define N 30

int main() {
    gps_init();
    double x = 0, y = 0, z = 0, t = 0;
    double xi[N], yi[N], zi[N], ti[N];
    while (true) {
        const int n = gps_read(xi, yi, zi, ti);
        if (n < 4) break;
        // TODO: find x, y, z, t
        gps_submit(x, y, z, t);
    }
    return 0;
}

