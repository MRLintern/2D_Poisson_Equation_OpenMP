// PSolver.c

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "PSolver.h"

// Compute root-mean-square (RMS) of a 2D array
double rms_function(int nx, int ny, double a[NX][NY]) {

    double v = 0.0;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {

            v += a[i][j] * a[i][j];
        }
    }
    
    return sqrt(v / (double)(nx * ny));
}

// Compute source term f(x,y) = pi*pi*(x*x + y*y) * sin(pi*x*y)
void source_function(int nx, int ny, double f[NX][NY]) {

    double x, y;

    for (int j = 0; j < ny; j++) {

        y = (double)j / (double)(ny - 1);

        for (int i = 0; i < nx; i++) {

            x = (double)i / (double)(nx - 1);

            f[i][j] = pow(PI, 2) * (x * x + y * y) * sin(PI * x * y);
        }
    }
}

// Jacobi iterative solver for the 2D Poisson equation with Dirichlet BCs
void jacobi(int nx, int ny, double dx, double dy, double f[NX][NY], int input_iteration, int output_iteration, double u[NX][NY], double u_new[NX][NY]) {

    // Precompute constants for stencil computation
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    // Iterative loop from input_iteration to output_iteration
    for (int step = input_iteration; step < output_iteration; step++) {

        // Perform Jacobi update over interior points (excluding boundaries)
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {

                u[i][j] = (

                    dy2 * (u_new[i + 1][j] + u_new[i - 1][j]) +
                    dx2 * (u_new[i][j + 1] + u_new[i][j - 1]) -
                    dx2 * dy2 * f[i][j]
                ) / denom;
            }
        }

        // Copy updated values to u_new for the next iteration
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {

                u_new[i][j] = u[i][j];
            }
        }
            
        // Explicitly enforce zero Dirichlet boundary conditions on u_new
        apply_dirichlet_bc(nx, ny, u_new);
    }
}

// Exact analytical solution for comparison: u(x,y) = sin(pi*x*y)
double potential_exact_function(double x, double y) {

    return sin(PI * x * y);
}

// Laplacian of the exact potential: del (del u) = pi*pi*(x*x + y*y) * sin(pi*x*y)
double laplacian_potential(double x, double y) {

    return pow(PI, 2) * (x * x + y * y) * sin(PI * x * y);
}

// Utility function to print a timestamp
void timestamp() {

    time_t now = time(NULL);
    printf("Timestamp: %s", ctime(&now));
}

// Apply homogeneous Dirichlet boundary conditions: u = 0 on the domain boundary
void apply_dirichlet_bc(int nx, int ny, double u[NX][NY]) {

    // Set left and right boundary values to 0
    for (int j = 0; j < ny; j++) {

        u[0][j] = 0.0;           // left boundary
        u[nx - 1][j] = 0.0;      // right boundary
    }

    // Set bottom and top boundary values to 0
    for (int i = 0; i < nx; i++) {

        u[i][0] = 0.0;           // bottom boundary
        u[i][ny - 1] = 0.0;      // top boundary
    }
}