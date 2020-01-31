#include "solver.h"

template<class T>
void swap_rows(Matrix<T>& A, int& j, int& k)
{
	const int N = A.cols;
	auto* jline = new double[N];
	for (int i = 0; i < N; i++)
	{
		jline[i] = A.values[j * A.rows + i];
	}

	auto* kline = new double[N];
	for (int i = 0; i < N; i++)
	{
		kline[i] = A.values[k * A.rows + i];
	}

	for (int i = 0; i < N; i++)
	{
		A.values[k * A.rows + i] = jline[i];
		A.values[j * A.rows + i] = kline[i];
	}

	delete[] jline;
	delete[] kline;
}

// find the index of the largest value
// and return it as the index of row
template<class T>
int argmax(int& k, Matrix<T>& A)
{
	int index = k * A.cols + k;
	for (int i = k; i < A.cols; i++)
	{
		if (abs(A.values[i * A.cols + k]) > abs(A.values[index]))
		{
			index = i * A.cols + k;
		}
	}
	return index/A.cols;
}

template<class T>
void LU_decomposition_pp(Matrix<T>& A, Matrix<T>& L, Matrix<T>& P_)
{
	const int m = A.cols;

	// create an empty matrix L
	for (int i = 0; i < m * m; i++)
	{
		L.values[i] = 0;
	}

	// create eye matrix P_

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) { P_.values[i * m + j] = 1; }
			else P_.values[i * m + j] = 0;
		}
	}

	for (int k = 0; k < m - 1; k++)
	{
		int j = argmax(k, A);
		swap_rows(A, j, k);
		swap_rows(P_, j, k);
		swap_rows(L, j, k);
		for (int i = k + 1; i < m; i++)
		{
			const double s = A.values[i * m + k] / A.values[k * m + k];
			L.values[i * m + k] = s;
			for (j = k; j < m; j++)
			{
				A.values[i * m + j] -= A.values[k * m + j] * s;
			}
		}
	}

	// add diagnal eye into L
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) { L.values[i * m + j] = 1; }
		}
	}

	// get the transpose of P_
	P_.transpose(P_);
}

template<class T>
void backward_substitution(Matrix<T>& U, T* b, T* output)
{
	const int m = U.rows;

	for (int k = m - 1; k >= 0; --k)
	{
		T s(0);
		for (int j = k + 1; j < m; j++)
		{
			s += U.values[k * m + j] * output[j];
		}
		output[k] = (b[k] - s) / U.values[k * m + k];
	}
}

template<class T>
void forward_substitution(Matrix<T>& L, T* b, T* output)
{
	const int m = L.rows;

	for (int k = 0; k < m; k++)
	{
		T s(0);
		for (int j = 0; j < k; j++)
		{
			s += L.values[k * m + j] * output[j];
		}
		output[k] = (b[k] - s) / L.values[k * m + k];
	}
	/*cout << "Printing forward sub vector output: " << endl;
	for (int i = 0; i < m; i++)
	{
		cout << " " << output[i];
	}*/
}

template<class T>
void LU_solver(Matrix<T>& A, T* x, T* b)
{
	if (A.rows != A.cols)
	{
		cerr << "Cannot decomposite non-square matrix into LU. \n";
		return;
	}

	const int m = A.cols;
	auto* L = new Matrix<T>(m, m, true);
	auto* P = new Matrix<T>(m, m, true);
	LU_decomposition_pp(A, *L, *P);
	// print out P, L, U
	//L->printMatrix();
	
	// A is now the upper triangle
	//A.printMatrix();
	// get the inverse i.e. the transpose of P
	P->transpose(*P);
	//P->printMatrix();

	auto* pinvb = new double[m * 1];
	P->matVecMult(b, pinvb);
	auto* y = new double[m * 1];
	forward_substitution(*L, pinvb, y);
	auto* solution = new double[m * 1];
	backward_substitution(A, y, solution);

	cout << "LU Solution: \n";
	for (int i = 0; i < m; i++)
	{
		cout << "x" << i << ": " << solution[i] << endl;
	}

	delete L;
	delete P;
	delete[] pinvb;
	delete[] y;
	delete[] solution;
}

template<class T>
void gauss_elimination(Matrix<T>& A, T* x, T* b)
{
	if (A.rows != A.cols)
	{
		cerr << "Cannot apply Gaussian Elimination on non-square matrix. \n";
		return;
	}

	const int m = A.cols;
	int i, j, k;
    for (i = 0; i < m; i++) {                   //Pivotisation
        for (k = i + 1; k < m; k++) {
            if (abs(A.values[i * m + i]) < abs(A.values[k * m + i])) {
                for (j = 0; j < m; j++) {
                    auto* temp = new double;
                    *temp = A.values[i * m + j];
                    A.values[i * m + j] = A.values[k * m + j];
                    A.values[k * m + j] = *temp;
                    delete temp;
                }
                auto* temp2 = new double;
                *temp2 = b[i];
                *(b + i) = b[k];
                *(b + k) = *temp2;
                delete temp2;
            }
        }
    }

	// Perform Gauss-Elimination
    for (i = 0; i < m - 1; i++) {
        for (k = i + 1; k < m; k++) {
            auto* t = new double;
            *t = A.values[k * m + i] / A.values[i * m + i];
            for (j = 0; j < m; j++) {

				//make the elements below the pivot elements equal to zero or eliminate the variables
                A.values[k * m + j] = A.values[k * m + j] - *t * A.values[i * m + j];
            }
            *(b + k) = b[k] - *t * b[i];
            delete t;
        }
    }

	backward_substitution(A, b, x);
}

template<class T>
void gauss_seidel(Matrix<T>& A, T* x, T* b, T er, T urf) {
	// gauss_idel(Matrix<T> & A, T * x, T * b, er, urf, nmax, rmx)

	// solves systems of linear eqs., with the Gauss-Seidel method

	/*
	// input

	n = number of eqs.
	a  = matrix of coefficients
	b  = matrix of constants
	x  = initial values for the unknown (nq values)
	er = termination criterion
	urf = relaxation factor
	nmax = maximum number of iterations

	// output

	x = solution
	rmx = maximum residual (verification)
	*/

	if (A.rows != A.cols)
	{
		cerr << "Cannot apply Gauss-Seidel on non-square matrix. \n";
		return;
	}

	const int n = A.cols;
	const int m = A.cols;
	int i, j, k, niter;
	double sum, xold;
	int nmax = 10000;
	double rmx = 0.0;

	// initialisation of x
	for (i = 0; i < n; i++) {
		*(x + i) = 0.0;
	}

	//Pivotisation
	for (i = 0; i < m; i++) {
		for (k = i + 1; k < m; k++) {
			if (abs(A.values[i * m + i]) < abs(A.values[k * m + i])) {
				for (j = 0; j < m; j++) {
					auto* temp = new double;
					*temp = A.values[i * m + j];
					A.values[i * m + j] = A.values[k * m + j];
					A.values[k * m + j] = *temp;
					delete temp;
				}
				auto* temp2 = new double;
				*temp2 = b[i];
				*(b + i) = b[k];
				*(b + k) = *temp2;
				delete temp2;
			}
		}
		cout << endl;
	}

	// start of iterations
	for (niter = 0; niter < nmax; niter++) {

		double ea = 0.0;

		// cout << "iter " << niter << endl;  // print iterations until convergence

		// new x's calculation 
		for (i = 0; i < n; i++) {
			xold = x[i];
			sum = b[i];

			// cout << "i " << i << endl;

			for (j = 0; j <= i - 1; j++) {
				sum += - A.values[i * n + j] * x[j];
				// cout << "iter j1: " << j << ", A: " << A.values[i * n + j] << ", x: " << x[j] << ", sum: " << sum << endl;
			}

			for (j = i + 1; j < n; j++) {
				sum += - A.values[i * n + j] * x[j];
				// cout << "iter j2: " << j << ", A: " << A.values[i * n + j] << ", x: " << x[j] << ", sum: " << sum << endl;
			}

			*(x + i) = sum / A.values[i * n + i];

			// error and underelaxation
			double ern = abs(*(x + i) - xold);

			// cout << i << " " <<  ern << endl;
			// cout << endl;
			ea = __max(ea, ern);
			*(x + i) = *(x + i) * urf + xold * (1. - urf);
		}

		// Checks for exit
		if (ea < er) {
			break;
		}

	}

	if (niter >= nmax) {
		cout << "Warning! Iterations' limit";
	}

	//// Verification
	//for (i = 0; i < n; i++) {
	//	sum = - b[i];
	//	for (j = 0; j < n; j++) {
	//		sum += A.values[i * n + j] * x[j];
	//	}
	//	rmx = __max(abs(sum), rmx);
	//}
}

template<class T>
void thomas(Matrix<T>& A, T* x, T* r) {
	int i;
	double term;
	const int n = A.cols;
	A.values[2 * n] = A.values[2 * n] / A.values[1 * n];
	r[0] = r[0] / A.values[1 * n];

	//-- forward elimination
	for (i = 1; i < n; i++) {
		term = A.values[1 * n + i] - A.values[0 * n + i] * A.values[2 * n + i - 1];
		A.values[2 * n + i] = A.values[2 * n + i] / term;
		r[i] = (r[i] - A.values[0 * n + i] * r[i - 1]) / term;
	}
	//-- back substitution
	x[n - 1] = r[n - 1];
	for (i = n - 2; i >= 0; i--) {
			x[i] = r[i] - A.values[2 * n + i] * x[i + 1];
		}
}

template<class T>
void cholesky_fact(Matrix<T>& A, T* x, T* b) {

	int i, j, k;
	const int n = A.rows;
	//-- elements of [L]^T  ------------------------------
	for (i = 0; i < n - 1; i++) {   // !!!!!!!!!!!!
		for (j = i + 1; j < n; j++) {
			A.values[i * n + j] = A.values[j * n + i];
		}
	}

	auto* au = new Matrix<T>(n, n, true);
	auto* al = new Matrix<T>(n, n, true);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			al->values[i * n + j] = 0.;
			au->values[i * n + j] = 0.;
			if (i >= j) {
				al->values[i * n + j] = A.values[i * n + j];
			}
			if (i <= j) {
				au->values[i * n + j] = A.values[i * n + j];
			}
		}
	}

	auto* y = new double[n];
	*y = b[0] / al->values[0];

	forward_substitution(*al, b, y);
	backward_substitution(*au, y, x);

	delete al;
	delete au;
	delete[] y;
}

template<class T>
void cholesky_dense(Matrix<T>& A, T* x, T* b) {

	int i, j, k;
	const int n = A.rows;

	// Dense matrix
	for (k = 0; k < n; k++) {
		A.values[k * n + k] = sqrt(A.values[k * n + k]);
		for (i = k + 1; i < n; i++) {
			A.values[i * n + k] = A.values[i * n + k] / A.values[k * n + k];
		}

		for (j = k + 1; j < n; j++) {
			for (i = j; i < n; i++) {
				A.values[i * n + j] -= A.values[i * n + k] * A.values[j * n + k];
			}
		}
	}
	cholesky_fact(A, x, b);
}


template<class T>
void cholesky_sparse(Matrix<T>& A, T* x, T* b) {

	int i, j, k;
	int n = A.rows;

	int* sk = new int[n + 1];
	// int sk[n];

	// Sparse matrix
	for (k = 0; k < n; k++) {
		A.values[k * n + k] = sqrt(A.values[k * n + k]);

		int nsk = 0;
		for (i = k + 1; i < n; i++) {
			A.values[i * n + k] = A.values[i * n + k] / A.values[k * n + k];

			if (A.values[i * n + k] != 0.) {
				
				nsk += 1;
				sk[nsk] = i;
				
			}
		}

		if (nsk > 0); {
			for (j = sk[1]; j <= sk[nsk]; j++) {
				for (i = j; i <= sk[nsk]; i++) {
					A.values[i * n + j] = A.values[i * n + j] - A.values[i * n + k] * A.values[j * n + k];
				}
			}
		}
	}

	delete[] sk;

	cholesky_fact(A, x, b);
}
