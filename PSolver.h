// Header file for the 2D Poisson equation solver using the Jacobi method.

// Provides declarations for all major functions used in solving
// the Poisson equation numerically, including initialization,
// exact solution, source term, iterative solver, and boundary conditions.


#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Problem Size --- //

#define NX 100        // Number of grid points in x-direction
#define NY 100        // Number of grid points in y-direction
#define PI 3.14159265358979323846  // Definition of pi

// -----------------------------------------------------------------------------
// Function: rms_function
// Purpose:  Compute the root-mean-square (RMS) of a 2D array
// Inputs:   nx, ny - grid dimensions
//           a      - 2D array to compute RMS on
// Returns:  Root-mean-square of all elements in array
// -----------------------------------------------------------------------------

double rms_function(int nx, int ny, double a[NX][NY]);

// -----------------------------------------------------------------------------
// Function: source_function
// Purpose:  Compute the source term f(x, y) = pi*pi*(x*x + y*y) * sin(pi*x*y)
// Inputs:   nx, ny - grid dimensions
// Outputs:  f      - 2D array containing source term values
// -----------------------------------------------------------------------------

void source_function(int nx, int ny, double f[NX][NY]);

// -----------------------------------------------------------------------------
// Function: jacobi
// Purpose:  Perform Jacobi iterations to solve the Poisson equation
// Inputs:   nx, ny           - grid dimensions
//           dx, dy           - spatial grid spacings
//           f                - source term array
//           input_iteration  - starting iteration index
//           output_iteration - ending iteration index
// In/Out:   u                - current solution estimate (updated in-place)
//           u_new            - copy of previous iteration's solution
// -----------------------------------------------------------------------------

void jacobi(int nx, int ny, double dx, double dy, double f[NX][NY],
            int input_iteration, int output_iteration,
            double u[NX][NY], double u_new[NX][NY]);

// -----------------------------------------------------------------------------
// Function: potential_exact_function
// Purpose:  Evaluate exact solution u(x,y) = sin(pi*x*y) at a given point
// Inputs:   x, y - coordinates
// Returns:  Exact solution u(x, y)
// -----------------------------------------------------------------------------

double potential_exact_function(double x, double y);

// -----------------------------------------------------------------------------
// Function: laplacian_potential
// Purpose:  Evaluate Laplacian of the exact solution (i.e., ∇²u)
// Inputs:   x, y - coordinates
// Returns:  Laplacian of u(x, y)
// -----------------------------------------------------------------------------

double laplacian_potential(double x, double y);

// -----------------------------------------------------------------------------
// Function: timestamp
// Purpose:  Print a timestamp to standard output
// -----------------------------------------------------------------------------

void timestamp();

// -----------------------------------------------------------------------------
// Function: apply_dirichlet_bc
// Purpose:  Explicitly apply zero Dirichlet boundary conditions (u = 0)
// Inputs:   nx, ny - grid dimensions
// In/Out:   u      - 2D solution array to enforce boundary values
// Notes:    Enforces u = 0 on the domain boundary for robustness
// -----------------------------------------------------------------------------

void apply_dirichlet_bc(int nx, int ny, double u[NX][NY]);

#endif  // POISSON_SOLVER_H
