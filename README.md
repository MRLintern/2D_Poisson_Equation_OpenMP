## 2D Poisson's Equation (OpenMP Version)

## Introduction
2D Poisson's Equation. Discretized using the Finite Difference Method &amp; Solved by Parallelising the Jacobi Iterative Method via the OpenMP API.

This is an example of how we can employ the OpenMP API to solve the Poisson equation.

It will help if you have an understanding of numerical partial differential equations (NPDE)
and numerical linear algebra (NLA).
There are many online resources for OpenMP, NPDE & NLA; YouTube is a good place to go.

The Poisson equation: 

	-div(grad(u)) = f, 

is discretized using the finite difference method.
The system of algebraic equations produced by this method is then iterated through 
using Jacobi's method to seek solutions at points in a [0,1]x[0,1] unit square.
OpenMP will be used to parallelize the Jacobi iterative solver.
This method will be used until convergence has been detected.
This is where OpenMP is useful; it will help speed up the computation of the 
system of algebraic equations.

This system of algebraic equations can be represented as a matrix equation:

	Au = f

and,

	div(grad),
	
represents the Laplace operator; this is A in the matrix equation.
u is the potential field, something we try find; in the case of an inverse problem.
f represents the source term; the data we already have.

Example scenario: geophysics):

Note: this isn't represented in this program.

f = gravitational field. 
u = mass density distribution which produces the gravitational field, the thing we want to find.

In the case that follows, we actually know both u and f.
We're not trying to solve an inverse problem, but look at how OpenMP speeds up computation.

The potential and source for the example application:

	u = sin(pi*x*y)
	
	f = pow(pi,2)*(pow(x,2)+pow(y,2))*sin(pi*x*y)
	
	Note: pow(x,2) means x^2 (x squared); pow() comes from the math library.

    Note: in the code, u = potential and f = source 
    

## Do I have OpenMP?
You should do if you have an up to date compiler. e.g. `gcc`.
You don't install OpenMP. OpenMP is a feature of the compiler. Check that the compiler you are using implements OpenMP, which is an API for parallel programming in `C/C++/Fortran`. 
The `gcc` compiler version used for this program is: `9.4.0`.

## Instructions
1. Clone the repo: `https://github.com/MRLintern/2D_Poisson_Equation_OpenMP.git`
2. `$make`
3. `$./main

