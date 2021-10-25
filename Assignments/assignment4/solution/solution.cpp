// Source: Newton's method for single variable
// https://personal.math.ubc.ca/~pwalls/math-python/roots-optimization/newton/

#include <stdio.h>
#include <math.h>
#include "blackbox.h"


#define N 30
static const double eps = 1e-6;

double calc_error(double *xi, double *yi, double *zi, double *ti,
                  double x, double y, double z, double t,
                  int n) {
    double sum = 0.0;
    for(int i = 0; i < n; i++){
        double error = (x-xi[i])*(x-xi[i]) + (y-yi[i])*(y-yi[i]) + (z-zi[i])*(z-zi[i]) - ((t-ti[i])*(t-ti[i]));
        sum += error * error;
    }
    return sum;
}




int main() {
    gps_init();
    double x = 0, y = 0, z = 0, t = 0;
    double xi[N], yi[N], zi[N], ti[N];
    while (true) {
        const int n = gps_read(xi, yi, zi, ti);
        if (n < 4) break; // If gps_read returns any value less than four, it means that we have lost the signal, and the program should quit (return 0;).
        // if(n==0) break;
        // TODO: find x, y, z, t
        double dFdx = 0.0;
        double dFdy = 0.0;
        double dFdz = 0.0;
        double dFdt = 0.0;
        
        for(int i = 0; i < n; i++){
            dFdx += (x-xi[i]);
            dFdy += (y-yi[i]);
            dFdz += (z-zi[i]);
            dFdt += (t-ti[i]);
            
        }

        int steps = 0;
        while (calc_error(xi, yi, zi, ti, x,y,z,t, n) > eps){
            // x = ?
            // y = ?
            // z = ?
            // t = ?
            // break;
            if(steps > 10)
               break;
        }
        gps_submit(x, y, z, t);
    }
    return 0;
}

