/*
This is an example of how we can employ the OpenMP API to solve the Poisson equation.

Note to the reader/user: 

It will be helpful if you have an understanding of numerical partial differential equations (NPDE)
and numerical linear algebra (NLA).
There are many online resources for OpenMP, NPDE & NLA; YouTube is a good place to go.


The Poisson equation: 

	-div(grad(u)) = f, 

is discretized using the finite difference method.
The system of algebraic equations produced by this method is then iterated through 
using Jacobi's method to seek solutions at points in a [0,1]x[0,1] unit square.
OpenMP will be used to parallelise the Jacobi iterative solver.
This method will be used until convergence has been detected.
This is where OpenMP is useful; it will help speed up the computation of the 
system of algebraic equations.

This system of algebraic equations can be represented as a matrix equation:

	Au = f

div(grad(u)) reprents the Laplace operator; this is A in the matrix equation.
u is the potential field, something we try find; in the case of an inverse problem.
f represents the source term; the data we already have.

Example scenario: geophysics):

Note: this isn't represented in this program.

f = gravitational field. 
u = mass density distribution which produces the gravitational field, the thing we want to find.

In the case that follows, we actually know both u and f.
We're not trying to solve an inverse problem, but look at how OpenMP speeds up computation.


	u = sin(pi*x*y)
	
	f = pow(pi,2)*(pow(x,2)+pow(y,2))*sin(pi*x*y)
	
	Note: pow(x,2) means x^2 (x squared); pow() comes from the math library.

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

//number of interior points in the grid
#define NX 161
#define NY 161

#define PI 3.141592653589793

//function prototypes

//root-mean-square function; used for error analysis. a represents elemets of the Matrix, A
double rms_function(int nx, int ny, double a[NX][NY]);
//source function represents the rhs of the equation; the source term. f represents elements of the source matrix
void source_function(int nx, int ny, double f[NX][NY]);
//jacobi iteration function
void jacobi(int nx, int ny, double dx, double dy, double f[NX][NY], int input_iteration, int output_iteration, double potential_array[NX][NY], double potential_new_array[NX][NY]);
//used for time aspects of the program
void timestamp();
//exact potential solution function
double potential_exact_function(double x, double y);
//laplacian operator; applied to the potential
double laplacian_potential(double x, double y);


int main ()
{
  //for convergence testing
  int converged;
  //used for differences in solution
  double diff;
  //for holding spacing values between grid points
  double dx, dy;
  double error;
  //elements of the source matrix
  double f[NX][NY];
  //counters
  int i, j;
  //thread I.D.
  int id;
  //output iteration on a point in the grid
  int output_iteration;
  //input iteration on a point in the grid
  int input_iteration;
  //interior points in the grid; 2D Cartesian plane; non-complex geometry
  int nx = NX;
  int ny = NY;
  double tolerance = 0.000001;
  //holds potential solution as an array
  double potential_array[NX][NY];
  //normed potential
  double potential_norm;
  //difference in potential solutions as an array
  double potential_diff_array[NX][NY];
  //exact potential solution array
  double potential_exact_array[NX][NY];
  //new potential solution array
  double potential_new_array[NX][NY];
  //new potential normed
  double potential_new_norm;
  //wall time
  double wtime;
  //points in the square
  double x, y;
  
  //calculate grid spacing 
  dx = 1.0 / ( double ) ( nx - 1 );
  dy = 1.0 / ( double ) ( ny - 1 );
  
  timestamp();
  printf ("\n");
  printf("Program Uses the OpenMP API to parallelize the Jacobi Iterative Solver to Increase Computation Time to Solve the 2D Poisson's Equation.\n");
  printf ("\n");
  //tell users how many processors are being used.
  printf ("The number of processors is %d\n", omp_get_num_procs ( ) );
# pragma omp parallel
{
  //how many threads are being used
  id = omp_get_thread_num ( );
  if (id == 0)
  {
    printf ("The maximum number of threads is %d\n", omp_get_num_threads()); 
  }
} 
  //provide user with some basic info r.e. geometry
  printf ("\n");
  printf ("Dimensions of square: 0 <= x <= 1, 0 <= y <= 1.\n");
  printf ("\n");
  printf ("The number of interior x grid points is %d\n", nx);
  printf ("The number of interior y grid points is %d\n", ny);
  printf ("The x grid spacing is %f\n", dx);
  printf ("The y grid spacing is %f\n", dy);

  //source function for Poissons Equation
  //f represents elements of the source matrix
  source_function(nx, ny, f);

  for (j = 0; j < ny; j++)
  {
    for (i = 0; i < nx; i++)
    {
      if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
      {
        potential_new_array[i][j] = f[i][j];
      }
      else
      {
        potential_new_array[i][j] = 0.0;
      }
    }
  }
  //calculating the norm
  potential_new_norm = rms_function(nx, ny, potential_new_array);
/*
  //Set up the exact solution potential_exact_array
*/
  for ( j = 0; j < ny; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny - 1 );
    for ( i = 0; i < nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx - 1 );
      potential_exact_array[i][j] = potential_exact_function( x, y );
    }
  }
  potential_norm = rms_function(nx, ny, potential_exact_array);
  printf ( "RMS of exact solution = %g\n", potential_norm );

  //Do the iteration.

  converged = 0;

  printf ( "\n" );
  printf ( "  Step    ||potential_new_array||     ||potential_new_array - potential_array||     ||potential_new_array - Exact||\n" );
  printf ( "\n" );

  for (j = 0; j < ny; j++)
  {
    for (i = 0; i < nx; i++)
    {
      potential_diff_array[i][j] = potential_new_array[i][j] - potential_exact_array[i][j];
    }
  }
  error = rms_function(nx, ny, potential_diff_array);
  
  printf ( "  %4d  %14g                  %14g\n", 0, potential_newnorm, error);
  
  //get wall time
  wtime = omp_get_wtime ( );

  output_iteration = 0;

  for (;;)
  {
    input_iteration = output_iteration;
    output_iteration = input_iteration + 500;

    //jacobi function carries out 500 steps in parallel before coming  back to check for convergence
  
    jacobi(nx, ny, dx, dy, f, input_iteration, output_iteration, potential_array, potential_new_array);

    //Check for convergence.

    poential_norm = potential_new_norm;
    potential_new_norm = rms_function(nx, ny, potential_new_array);

    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        potential_diff_array[i][j] = potential_new_array[i][j] - potential_array[i][j];
      }
    }
    //differences in rms for the potential difference in potential solutions
    diff = rms_function(nx, ny, potential_diff_array);

    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        potential_diff_array[i][j] = potential_new_array[i][j] - potential_exact_array[i][j];
      }
    }
    error = rms_function(nx, ny, potential_diff_array);

    printf ( "  %4d  %14g  %14g  %14g\n", output_iteration, potential_new_norm, diff, error);

    //if values greater than tolerance, solutions not converged and blow-up
    if ( diff <= tolerance )
    {
      converged = 1;
      break;
    }

  }

  if ( converged )
  {
    printf ( "  The iteration has converged.\n" );
  }
  else
  {
    printf ( "  The iteration has NOT converged.\n" );
  }

  //time taken
  wtime = omp_get_wtime ( ) - wtime;
  printf ( "\n" );
  printf ( "  Elapsed seconds = %g\n", wtime );

  printf ( "\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp();

  return 0;
}
//function definitions

//rms returns the root-mean-square norm of a vector stored as a matrix.
//NX & NY are the number of rows and columns in the matrix A; note: a[NX][NY] are the matrix A elements
//rms will return the root-mean-square value of matrix A

double rms_function(int nx, int ny, double a[NX][NY])
{
  //counters
  int i, j;
  //approximate vector for potential (u)
  double v;

  v = 0.0;

  for (j = 0; j < ny; j++)
  {
    for (i = 0; i < nx; i++)
    {
      v = v + a[i][j] * a[i][j];
    }
  }
  v = sqrt(v / (double)(nx*ny));

  return v;
}

//source function; RHS of Poissons Equation.
//initialized the RHS of the matrix equation.
//f[NX][NY] are the elements of the source matrix.
void source_function(int nx, int ny, double f[NX][NY] )
/*
   Recall: Au = f
   
   When u(i,j) is a boundary value, then f(i,j) holds the boundary data.
   
   i.e. 
   
   	u(i,j) = f(i,j)

  When not at boundary points, Poissons Equation has the form (after Finte Differences have been applied):

        (1/pow(h,2))*(u(i+1,j) + u(i-1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j)) = f(i,j)

  Where h is the spacing between the points (or nodes) in the grid.

*/
{
  double fnorm; //normed source elements
  int i, j; //counters
  double x, y; //points in the square
  
  for (j = 0; j < ny; j++)
  {
    y = ( double ) ( j ) / ( double ) ( ny - 1 );
    for ( i = 0; i < nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx - 1 );
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        f[i][j] = potential_exact_function(x, y);
      }
      else
      {
        f[i][j] = - laplacian_potential(x, y);
      }
    }
  }

  fnorm = rms_function(nx, ny, f);

  printf ( "  RMS of f = %g\n", fnorm );

  return;
}
//jacobi function
//function to iterate through the system of equations. Parallized using OpenMP
void jacobi(int nx, int ny, double dx, double dy, double f[NX][NY], int input_iteration, int output_iteration, double potential_array[NX][NY], double potential_new_array[NX][NY] )

/*

  The geometry of the unit square is simple.
  Partial derivatives in x and y assumed the same. With that in mind, 

      - ( d/dx d/dx + d/dy d/dy ) u(x,y) 

  can be simplified to:

      
      (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j))/dx/dy

    
   double potential_array[NX][NY] is the solution estimate on iteration iteration_output - 1
   double potential_new_array[NX][NY] is the solution estimate on iteration iteration_input
   for the output, potential_new_array[NX][NY] is the solution estimate on iteration iteration_output
   
*/
{
  int i, j, iter; //counters
  
# pragma omp parallel \
  shared(dx, dy, f, output_iteration, input_iteration, nx, ny, potential_array, potential_new_array) \
  private(i, iter, j)

  for (iter = input_iteration + 1; it <= output_iteration; iter++)
  {

  //Save the current estimate.

# pragma omp for
    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        potential_array[i][j] = potential_new_array[i][j];
      }
    }

  //Compute a new estimate.

# pragma omp for
    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1)
        {
          potential_new_array[i][j] = f[i][j];
        }
        else
        { 
          potential_new_array[i][j] = 0.25*( potential_array[i-1][j] + potential_array[i][j+1] + potential_array[i][j-1] + potential_array[i+1][j] + f[i][j]*dx*dy);
        }
      }
    }

  }
  return;
}

//function prints the current YMDHMS date as a time stamp
void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf ("%s\n", time_buffer);

  return;
  
# undef TIME_SIZE
}

//evalutes the exact solution at points x & y
double potential_exact_function(double x, double y)
{
  
  double value;

  value = sin(PI*x*y);

  return value;
}

//laplacian_potential evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.
double laplacian_potential(double x, double y)
{
  
  double value;

  value = -pow(PI,2)*(pow(x,2) + pow(y,2))*sin(PI*x*y);

  return value;
}

# undef NX
# undef NY
