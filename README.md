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
* **Dense Matrix Solver**
	* Gaussian Elimination
	* LU Decomposition
	* Gauss-Seidel Iteration
	* Jacobi Iteration
	* Cholesky Factorisation
* **Sparse Matrix Solver**
	* LU Decomposition
	* Gauss-Seidel Iteration
	* Jacobi Iteration
	* Cholesky Factorisation
* **Tridiagonal (banded) Matrix**
	* Thomas Algorithm

## Parameters
* **Regular Parameters:**
	* A the left hand side matrix
	* b the right hand side constants
* **Parameters for iterative solvers:**
	* maxit: maximum iteration time
	* tolerance: tolerance/criterion to stop the iteration
	* relaxation factor: for more stable performance
