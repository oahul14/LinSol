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

using namespace std;

//int main()
//{
//	const int rows(5), cols(5);
//	const int m(5);
//	//// testing our Matrix class
//	double array[25] = { 1,0,3,7,2,1,0,4,5,4,1,-2,4,1,6,2, 6, 9,1,5,2,3,6,8,0 };
//
//	auto* A = new Matrix<double>(rows, cols, array);
//
//	A->printMatrix();
//
//	// double barray[] = { 1,2,-3,8,3 };
//	vector<double> barray = { 1, 2, -3, 8, 3};
//	auto* b = new double[m*1];
//	cout << "\nRHS: " << endl;
//	for (int i = 0; i < m; i++)
//	{
//		*(b+i) = barray[i];
//		cout << *(b + i) << endl; 
//	}
//	cout << endl;
//	auto* x = new double[m * 1];
//
//	LU_solver(*A, x, b);
//	
//	delete A;
//	delete[] b;
//	delete[] x;
//
//	return 0;
//}

using namespace std;

int& assignInt()
{
	int *i = new int;
	return *i;
}

int main()
{
	int i = assignInt();
	cout << "values of i: " << i << endl;
	delete i;
}