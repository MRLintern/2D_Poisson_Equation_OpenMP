## 2D Poisson's Equation

### Introduction: Poisson's Equation
* ___Poisson's Equation___ is a __Partial Differential Equation (PDE)__ that describes how a __Scalar Field__ (like __Electric Potential__, __Temperature__, or __Gravitational Potential__) varies in __Space__ due to a __Source Term__.
* The mathematical form of this PDE is:
  	* `del * (del u) = f(x, y)`.
  	* `del * del` is the __Laplace Operator__.
  	* `u` is the unknown __Scalar Field__.
  	* `f(x, y)` is the __Source Term__ in `2D Space`; note: `3D` case: `f(x, y, z)`.
* ___What does this equation tell us?___
* __Poisson’s Equation__ tells you how sources (e.g., __charge__, __mass__, __heat__) generate or influence a __Potential Field__.
* In terms of __Fluid Dynamics__, `u` represents the ___Pressure Field___ inside an __Incompressible Flow__ and `f` represents __Negative Divergence__ of __Velocity__.
* What this means is How __Velocity Divergence__ relates to __Pressure__.
* Note: this project/software doesn't look any physical scenario.
### Poisson's Equation and OpenMP
* `OpenMP` helps solve Poisson’s equation more efficiently by `Parallelizing` the `computational workload`, which is often the most time-consuming part in numerical methods.
* In methods like __Finite Difference__ or __Finite Element__, you Discretize the Domain into a Grid and iteratively solve for `u(x,y)` using an algorithm like `Jacobi`, `Gauss-Seidel`, or `SOR`.
* For a `2D grid` of size `N×N`, each iteration involves updating every interior grid point using neighboring values.
* This is computationally expensive, especially for:
  	* Fine grids (large `N`).
  	* Many iterations (to converge).
#### How OpenMP Helps
##### 1. Parallelizes Grid Updates
* Each point in the grid update is `independent` (especially in `Jacobi iteration`), so `OpenMP` can safely `split the workload` among `multiple threads`.
* This means `multiple points are computed simultaneously`, reducing total computation time.
##### 2. Speeds Up Iterations
* Since each iteration involves updating the entire grid, `OpenMP` allows `multiple CPU cores` to work in parallel, `reducing the time per iteration`. This is especially effective for high-resolution grids.
##### 3. Improves Scalability
* OpenMP is scalable across many cores, so performance improves with hardware — making it a simple yet powerful way to `scale up simulations`.
### Requirements
* __Compiler__: `gcc 13.1.0`; needs to support `OpenMP`.
* __OS__: `Ubuntu 20.04`.
* A `text editor for any changes.
* [Make](https://www.gnu.org/software/make/).
### Getting and Using the Software
* `$ git clone git@github.com:MRLintern/2D_Poisson_Equation_OpenMP.git`
* `$ cd 2D_Poisson_Equation_OpenMP`
* If you want to generate the results and produce the plot yourself:
* `$ make clean`
* __Unit Testing__:
* `$ make test`
* All test cases (should) pass. Now create the executable:
* `$ make`
* `$ ./main`
* __Plot the Results__:
* You don't need to be in the `results` directory create the plot; just stay inside `2D_Poisson_Equation_OpenMP`. Note: the plot, `solution_and_error_surface.png`, will be saved inside `results`.
* `$ python3 plot_results.py`
* An image will be created and saved in the `results` directory.
