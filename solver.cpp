#include "solver.h"


template<class T>
void swap_rows(Matrix<T>& A, int& j, int& k)
{
	const int N = A.cols;
	for (int i = 0; i < A.cols; i++)
	{
		T* num = new T;
		*num = A.values[k * A.rows + i];
		A.values[k * A.rows + i] = A.values[j * A.rows + i];
		A.values[j * A.rows + i] = *num;
		delete num;
	}
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
			for (int j = k; j < m; j++)
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
	L->printMatrix();
	
	// A is now the upper triangle
	A.printMatrix();
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