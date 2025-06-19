#include <stdio.h>
#include "PSolver.h"

int main() {

    double a[NX][NY] = {0};

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {

            a[i][j] = 1.0;

        }
    }

    double rms = rms_function(NX, NY, a);

    printf("%.6f\n", rms);
}