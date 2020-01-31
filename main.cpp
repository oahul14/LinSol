#include <iostream>
#include <math.h>
#include <ctime>
#include "solver.h"
#include "solver.cpp"
#include <vector>

using namespace std;



int main()
{
	const int rows(10), cols(10);
	const int m(10);
	//// testing our Matrix class
	double array[100] = {	6,0,3,0,0,0,0,3,0,0,
							0,8,0,0,0,0,4,0,0,1,
							3,0,9,0,0,2,0,0,0,0,
							0,0,0,7,0,0,3,0,0,0,
							0,0,0,0,9,0,0,5,0,4,
							0,0,2,0,0,8,0,0,6,0,
							0,4,0,3,0,0,6,0,0,0,
							3,0,0,0,5,0,0,7,0,2,
							0,0,0,0,0,6,0,0,9,1,
                            0,1,0,0,4,0,0,2,1,8     };

	auto* A = new Matrix<double>(rows, cols, array);

	A->printMatrix();

	// double barray[] = { 1,2,-3,8,3 };
	vector<double> barray = { 1, 2, -3, 8, 3, 6, 9, 4, 2, 3};
	auto* b = new double[m*1];
	cout << "\nRHS: " << endl;
	for (int i = 0; i < m; i++)
	{
		*(b+i) = barray[i];
		cout << *(b + i) << endl; 
	}
	cout << endl;
	auto* x = new double[m * 1];

	LU_solver(*A, x, b);
	
	delete A;
	delete[] b;
	delete[] x;

	return 0;
}

/*// using namespace std;

// int& assignInt()
// {
// 	int *i = new int;
// 	return *i;
// }

// int main()
// {
// 	int i = assignInt();
// 	cout << "values of i: " << i << endl;
// 	delete i;
// }*/
