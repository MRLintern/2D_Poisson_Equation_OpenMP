#include <stdio.h>
#include "PSolver.h"

int main() {

    double u[NX][NY] = {0};

    // Fill with nonzero values
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {

            u[i][j] = 1.0;
        }
    }

    apply_dirichlet_bc(NX, NY, u);

    // Check edges
    for (int i = 0; i < NX; i++) {

        if (u[i][0] != 0 || u[i][NY-1] != 0) {

            printf("FAIL\n");

            return 1;
        }
    }

    for (int j = 0; j < NY; j++) {

        if (u[0][j] != 0 || u[NX-1][j] != 0) {

            printf("FAIL\n");

            return 1;
        }
    }

    printf("PASS\n");
}
