#include <iostream>
#include <algorithm>
#include <memory>
#include "Matrix.h"
#include "CSRMatrix.h"
#include "solver.h"
#include "solver.cpp"
#include <vector>


using namespace std;

int main()
{
    int rows(15), cols(15);
    auto* dense_mat = new Matrix<double>(rows, cols, true);
    dense_mat->genRanSparse(0.7, true);
    dense_mat->printMatrix();
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, 1, true);
    sparse_mat->fromDense(*dense_mat);
    sparse_mat->printMatrix();

    auto* b = new double[rows * 1];
    cout << "\nLinear system: RHS [";
    for (int i = 0; i < rows; i++)
    {
        b[i] = rand() % 10 + 1;
        cout << b[i] << " ";
    }
    cout << "]" << endl;

    auto* x = new double[rows * 1];

    //********************** LU decomposition ***********************
    // LU decomposition
    sparse_mat->LU_solver(b, x);
    cout << endl;
    LU_dense(*dense_mat, b, x);
    cout << endl;

//    int maxit = 10000;
//    double tol = 1e-4;
////    //*********************** Gauss-Seidel **************************
////    gauss_seidel_dense(*dense_mat, x, b, maxit, tol);
//    gauss_seidel_sparse(*sparse_mat, x, b, maxit, tol);
    
    delete dense_mat;
    delete sparse_mat;
    delete[] x;
    delete[] b;
    
//    int size = 10;
//    vector<int> temp_col{1,4,2,7,3,8,6,11,13,5};
//    vector<double> temp_value{60,9,24,35,150,2,67,32,44,55};
//
//    int tempc, tempv, k;
//    for (int i = 1; i < size; i++)
//    {
//        for (int j = 1; j < size-i; j++)
//        {
//            if (temp_col[j] > temp_col[j+1])
//            {
//                tempc = temp_col[j];
//                tempv = temp_value[j];
//
//                temp_col[j] = temp_col[j+1];
//                temp_value[j] = temp_value[j+1];
//
//                temp_col[j+1] = tempc;
//                temp_value[j+1] = tempv;
//            }
//        }
//    }
//    for (int i = 0; i < size; i++)
//    {
//        cout << temp_col[i] << "  " << temp_value[i] << endl;
//    }

    
    
    return 0;
}
