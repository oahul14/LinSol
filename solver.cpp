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

	//back-substitution
	for (i = m - 1; i >= 0; i--) {

		// make the variable to be calculated equal to the rhs of the last equation
		*(x + i) = b[i];
		for (j = i + 1; j < m; j++) {

			// subtract all the lhs values except the coefficient of the variable whose value is being calculated
			if (j != i) {
				*(x + i) = x[i] - A.values[i * m + j] * x[j];
			}
		}
		//now finally divide the rhs by the coefficient of the variable to be calculated
		*(x + i) = x[i] / A.values[i * m + i];
	}
}
