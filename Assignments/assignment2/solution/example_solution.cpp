// Left side integration
#include <stdio.h>
#include "blackbox.h"

static double eps = 1e-9;

int main() {
    int n;
    scanf("%d", &n);
    blackbox_init(n);

    const size_t num_steps = 1e8;
    const double h = 2.0 / num_steps;

    double sum = 0;
    const double x0 = -1 + h/2;
    for (size_t i = 0; i < num_steps; i++) {
        sum += blackbox(-1 + 0.5*h + i*h) * h;  // -1 is shifting as the integration from -1 to 1 and multiplied by h -> the area
    }

    printf("%.12lf\n", sum);
}

