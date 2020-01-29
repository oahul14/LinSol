#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "solver.h"
#include "solver.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <vector>
#include <typeinfo>
#include <type_traits>

using namespace std;

int main()
{
	const int no = 5;
	const int rows(no), cols(no);
	const int m(no);
	// const int m(5);
	//// testing our Matrix class
	// double array[25] = { 1, 0, 3, 7, 2, 1, 0, 4, 5, 4, 1, -2, 4, 1, 6, 2, 6, 9, 1, 5, 2, 3, 6, 8, 0 };

	// diagonal dominance
	double array[25] = { 22, 0, 3, 7, 2, 1, 20, 4, 5, 4, 1, -2, 24, 1, 6, 2, 6, 9, 21, 5, 2, 3, 6, 8, 20 };

	// tridiagonal
	// double array[12] = { 0, 6, 2, 4, 2, 3, 5, 3, 3, 9, 2, 0 };

	// Intialisation of A, b and x
	Matrix<double>* A = new Matrix<double>(rows, cols, array);
	//std::cout << typeid(A).name() << std::endl;
	A->printMatrix();

	vector<double> barray = { 1, 2, -3, 8, 3 };
	// vector<double> barray = { 21, 69, 34, 22};
	auto* b = new double[m];
	double* x = new double[m];


	// Printing
	cout << "\nRHS: " << endl;
	for (int i = 0; i < m; i++)
	{
		*(b+i) = barray[i];
		cout << *(b + i) << endl; 
	}
	cout << endl;

	

	///// SOLVERS /////

	// Just uncomment the method to be used

	// LU_solver(*A, x, b);
	// gauss_elimination(*A, x, b);
	gauss_seidel(*A, x, b, 1e-3, 0.9);
	// cholesky(*A, x, b, 1);
	// thomas(*A, x, b);

	cout << "Solution: \n";
	for (int i = 0; i < m; i++)
	{
		cout << "x" << i << ": " << x[i] << endl;
	}
	
	delete A;
	delete[] b;
	delete[] x;

	return 0;
}

