# acse-5-assignment-nan

# Group Member
* Sokratis Anagnostopoulos
* Lingaona Zhu
* Hao Lu

# BLAS Pre-requisites
* Have OpenBLAS installed to recogonise <cblas.h>
* To make sure it compiles, the path of installed OpenBLAS need to be exported first:
* Below is the example used in macOS
```
	brew install OpenBLAS
	export LDFLAGS="-L/usr/local/opt/openblas/lib"
  	export CPPFLAGS="-I/usr/local/opt/openblas/include"
```
* Then in the command line compile using: 
```
	gcc-9 -lstdc++ -g -I/usr/local/opt/openblas/include -L/usr/local/opt/openblas/lib -lopenblas main.cpp 
```
* If using Windows please have OpenBLAS installed first, or other open libraries that include cblas.h

# Linear Solvers
* Dense Matrix Solver
	1. Gaussian Elimination
	2. LU Decomposition
	3. Gauss-Seidel Iteration
	4. Jacobi Iteration
	5. Cholesky Factorisation
* Sparse Matrix Solver
	1. LU Decomposition
	2. Gauss-Seidel Iteration
	3. Jacobi Iteration
	4. Cholesky Factorisation
* Tridiagonal (banded) Matrix
	* Thomas Algorithm

## Parameters
$$$Ax = b$
* A the left hand side matrix
* b the right hand side constants
* For iterative solvers:
	1. maxit: maximum iteration time
	2. tolerance: tolerance/criterion to stop the iteration
	3. relaxation factor: for more stable performance
