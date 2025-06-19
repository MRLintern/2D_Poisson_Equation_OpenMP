// Main driver for solving the 2D Poisson equation on a unit square
// using the Jacobi method. The source term and exact solution are
// analytically known
// Results are saved to CSV.

#include <stdio.h>
#include <math.h>
#include <omp.h> // for OpenMP
#include <sys/stat.h> // for mkdir()
#include <sys/types.h> // for mode_t
#include "PSolver.h"


int main() {

    // create results directory if it doesn't exist
    mkdir("results", 0777);

    // Grid size (internal + boundaries)
    int nx = NX;
    int ny = NY;

    // Grid spacing (assuming unit square [0,1] x [0,1])
    double dx = 1.0 / (nx - 1);
    double dy = 1.0 / (ny - 1);

    // Number of Jacobi iterations
    int num_iterations = 5000;

    // Declare arrays for potential, new potential, source, error
    double u[NX][NY] = {0};             // Numerical solution
    double u_new[NX][NY] = {0};         // Working buffer for Jacobi update
    double f[NX][NY] = {0};             // Source term
    double u_exact[NX][NY] = {0};       // Exact solution
    double error[NX][NY] = {0};         // Error = u - u_exact

    // Print start timestamp
    timestamp();

    printf("Initializing source term...\n");

    source_function(nx, ny, f);

    printf("Running Jacobi solver with %d iterations...\n", num_iterations);
    jacobi(nx, ny, dx, dy, f, 0, num_iterations, u, u_new);

    printf("Computing exact solution and error...\n");

    // Calculate the exact solution and error
    for (int j = 0; j < ny; j++) {

        double y = (double)j / (ny - 1);

        for (int i = 0; i < nx; i++) {

            double x = (double)i / (nx - 1);

            u_exact[i][j] = potential_exact_function(x, y);

            error[i][j] = u[i][j] - u_exact[i][j];
        }
    }

    // Compute and print RMS error
    double rms_error = rms_function(nx, ny, error);

    printf("RMS error = %e\n", rms_error);

    // Output numerical solution and error to CSV files
    FILE *file_sol = fopen("results/solution.csv", "w");
    FILE *file_err = fopen("results/error.csv", "w");

    if (!file_sol || !file_err) {

        fprintf(stderr, "Error opening output files.\n");

        return 1;
    }

    printf("Writing output to solution.csv and error.csv...\n");

    // Write headers
    fprintf(file_sol, "x,y,u_numeric\n");
    fprintf(file_err, "x,y,error\n");

    for (int j = 0; j < ny; j++) {

        double y = (double)j / (ny - 1);

        for (int i = 0; i < nx; i++) {

            double x = (double)i / (nx - 1);

            fprintf(file_sol, "%f,%f,%f\n", x, y, u[i][j]);
            fprintf(file_err, "%f,%f,%f\n", x, y, error[i][j]);
        }
    }

    // close files
    fclose(file_sol);
    fclose(file_err);

    // Print finish timestamp
    timestamp();

    // end of main()
}