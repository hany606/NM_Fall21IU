// Source: Gradient descent
// http://fa.bianp.net/teaching/2018/eecs227at/gradient_descent.html

#include <stdio.h>
#include <math.h>
#include "blackbox.h"


#define N 30
#define alpha 0.001
#define lambda 2 // extra cost for the time if it was negative

static const double eps = 1e-6;

double calc_error(double *xi, double *yi, double *zi, double *ti,
                  double x, double y, double z, double t,
                  int n) 
{
    double sum = 0.0;
    for(int i = 0; i < n; i++){
        double error = (x-xi[i])*(x-xi[i]) + (y-yi[i])*(y-yi[i]) + (z-zi[i])*(z-zi[i]) - ((t-ti[i])*(t-ti[i]));
        sum += error * error;
    }
    return sum;
}



// https://stackoverflow.com/questions/33058848/generate-a-random-double-between-1-and-1/33058967
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}


int main() {
    gps_init();
    double x = randfrom(0,1), y = randfrom(0,1), z = randfrom(0,1), t = randfrom(0,1);
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
        
        int sign_t = -1;    // if t is positive
        if(t < 0) sign_t = 1;

        for(int i = 0; i < n; i++){
            dFdx += 2.0*(x-xi[i]);
            dFdy += 2.0*(y-yi[i]);
            dFdz += 2.0*(z-zi[i]);
            dFdt += 2.0*(t-ti[i]) + sign_t*lambda;
        }

        int steps = 0;
        while (true){
            double error = calc_error(xi, yi, zi, ti, x,y,z,t, n);
            if(error <= eps) break;
            // printf("%.12lf\n", error);
            x = x - alpha*dFdx;
            y = y - alpha*dFdy;
            z = z - alpha*dFdz;
            t = t - alpha*dFdt;
            steps += 1;
            // break;
            if(steps > 100)
               break;
        }
        gps_submit(x, y, z, t);
    }
    return 0;
}

