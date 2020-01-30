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
#include <memory>

using namespace std;

int main()
{
    int rows(9), cols(9);

    /*int rp[6] = { 0,4,7,15,20,25 };
    int ci[25] = {0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4};
    int v[25] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
    int nnzs = 25;*/

    /*int rp[6] = {0,2,4,6,10,13};
    int ci[13] = {0,4,1,3,2,3,1,2,3,4,0,3,4};
    int v[13] = {3,5,3,7,4,5,7,5,1,2,5,2,9};
    int nnzs = 13;*/

    int rp[10] = { 0,3,6,10,15,20,23,27,30,31 };
    int ci[31] = { 0,2,4,1,3,6,0,2,3,7,1,2,3,4,6,0,3,4,5,7,4,5,6,1,3,5,6,2,4,7,8 };
    int v[31] = { 1,6,3,2,4,5,6,3,9,1,4,9,4,7,1,3,7,5,3,4,3,3,2,5,1,2,9,1,4,6,3 };
    int nnzs = 31;

    int *row_position = new int[rows+1];
    int *col_index = new int[nnzs];
    double *values = new double[nnzs];

    for (int i = 0; i < nnzs; i++)
    {
        col_index[i] = ci[i];
        values[i] = v[i];
    }

    //cout << "-------col_index, values assignment done-------" << endl;

    for (int i = 0; i < rows+1; i++)
    {
        row_position[i] = rp[i];
    }

    //cout << "-------row position assignment done-------" << endl;
    
    //initialise x and b
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, values, row_position, col_index);
    cout << "Linear system: Matrix" << endl;
    sparse_mat->printMatrix();
    sparse_mat->printDense();
    
    vector<double> barray = { 1, 2, -3, 8, 3, 6, 9, 4, 2 };
    auto* b = new double[rows * 1];
    cout << "\nLinear system: RHS [";
    for (int i = 0; i < rows; i++)
    {
        *(b + i) = barray[i];
        cout << *(b + i) << " ";
    }
    cout << "]" << endl;

    auto* x = new double[rows * 1];

    sparse_mat->LU_solver(b, x);

    delete[] row_position;
    delete[] col_index;
    delete[] values;
    delete sparse_mat;

    return 0;
}