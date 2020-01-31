#include <iostream>
#include <algorithm>
#include <memory>
#include "Matrix.h"
#include "CSRMatrix.h"
#include "solver.h"
#include "solver.cpp"


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
//    //*********************** Gauss-Seidel **************************
//    gauss_seidel_dense(*dense_mat, x, b, maxit, tol);
//    gauss_seidel_sparse(*sparse_mat, x, b, maxit, tol);
    
    
    delete dense_mat;
    delete sparse_mat;
    
    
    return 0;
}
