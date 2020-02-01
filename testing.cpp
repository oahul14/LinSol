#include "testing.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <Accelerate/Accelerate.h>
#include <memory>
#include <ctime>
#include <string>

using namespace std;

template<class T>
void xinit(T* x, int n)
{
    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
    }
}

template<class T>
void printVec(T* x, int size, bool print)
{
    if (print)
    {
        int act_size = 10;
        if (size <= 10) act_size = size;
        else cout << "(1 - 10)" << endl;
        cout << "Solution: [";
        for (int i = 0; i < act_size; i++)
        {
            cout << setw(8) << x[i] << " ";
        }
        cout << "]" << endl;
    }
}

template <class T>
void test_gauss_elimination(Matrix<T>& A, T* x, T* b, bool print)
{
    gauss_elimination(A, x, b);
    cout << "\nGauss Elimination Solution Dense: ";
    printVec(x, A.cols, print);
}

template<class T>
void test_LU_dense(Matrix<T>& A, T* x, T* b, bool print)
{
    LU_dense(A, x, b);
    cout << "\nLU Solution Dense: ";
    printVec(x, A.cols, print);
}

template<class T>
void test_LU_sparse(CSRMatrix<T>& A, T* x, T* b, bool print)
{
    LU_sparse(A, x, b);
    cout << "\nLU Solution Sparse: ";
    printVec(x, A.cols, print);
}

template<class T>
void test_jacobi_dense(Matrix<T>& A, T* x, T* b, int maxit, double tol, bool blas, bool print)
{
    if (!blas)
    {
        jacobi_dense(A, x, b, maxit, tol);
        cout << "\nJacobi Solution Dense: ";
    }
    else
    {
        jacobi_dense_blas(A, x, b, maxit, tol);
        cout << "\nJacobi Solution Dense BLAS: ";
    }
    printVec(x, A.cols, print);
}

template <class T>
void test_jacobi_sparse(CSRMatrix<T>& A, T* x, T* b, int maxit, double tol, bool print)
{
    jacobi_sparse(A, x, b, maxit, tol);
    cout << "\nJacobi Solution Sparse: ";
    printVec(x, A.cols, print);
}

template<class T>
void test_gauss_seidel_dense(Matrix<T>& A,T* x, T* b, int maxit, T er, T urf, bool blas, bool print)
{
    if (!blas)
    {
        gauss_seidel_dense(A, x, b, maxit, er, urf);
        cout << "\nGauss-Seidel Solution Dense: ";
    }
    else
    {
        gauss_seidel_dense_blas(A, x, b, maxit, er, urf);
        cout << "\nGauss-Seidel Solution Dense BLAS: ";
    }
    printVec(x, A.cols, print);
}

template <class T>
void test_gauss_seidel_sparse(CSRMatrix<T>& A, T* x, T* b, int maxit, double tol, bool print)
{
    gauss_seidel_sparse(A, x, b, maxit, tol);
    cout << "\nGauss Seidel Solution Sparse: ";
    printVec(x, A.cols, print);
}

template <class T>
void test_thomas_tri(Matrix<T>& A, T* x, T* b, bool print)
{
    thomas(A, x, b);
    cout << "\nThomas Solution Tridiagonal: ";
    printVec(x, A.cols, print);
}

template <class T>
void test_cholesky_dense(Matrix<T>& A, T* x, T* b, bool print)
{
    cholesky_dense(A, x, b);
    cout << "\nCholesky Solution Dense: ";
    printVec(x, A.cols, print);
}

template <class T>
void test_cholesky_sparse(Matrix<T>& A, T* x, T* b, bool print)
{
    cholesky_sparse(A, x, b);
    cout << "\nCholesky Solution Sparse: ";
    printVec(x, A.cols, print);
}

void time_all_rand(int n, string mat_type, bool print)
{
    int rows(n), cols(n);
    Matrix<double>* dense_mat = new Matrix<double>(rows, cols, true);
    if (mat_type == "dense") dense_mat->genRanDense(true);
    else if (mat_type == "sparse") dense_mat->genRanSparse(0.7, true);
    else if (mat_type == "tridiagonal") dense_mat->genRanTri(true);
    else
    {
        cout << "Matrix type not implemented. ";
        return;
    }
    if (rows <= 10) dense_mat->printMatrix();
    
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, 1, true);
    sparse_mat->fromDense(*dense_mat);
    if (rows <= 10) sparse_mat->printMatrix();

    auto* b = new double[rows * 1];
    //cout << "\nLinear system: RHS [";
    for (int i = 0; i < rows; i++)
    {
        b[i] = rand() % 10 + 1;
        //cout << b[i] << " ";
    }
    //cout << "]" << endl;
    auto* x = new double[rows * 1];
    
    // ***************  Tests  ****************
    // gauss elimination
    xinit(x, rows);
    clock_t time_ge;
    time_ge = clock();
    test_gauss_elimination(*dense_mat, x, b, print);
    time_ge = clock() - time_ge;

    // LU dense
    xinit(x, rows);
    clock_t time_LU_dense;
    time_LU_dense = clock();
    test_LU_dense(*dense_mat, x, b, print);
    time_LU_dense = clock() - time_LU_dense;

    // LU sparse
    xinit(x, rows);
    clock_t time_LU_sparse;
    time_LU_sparse = clock();
    test_LU_sparse(*sparse_mat, x, b, print);
    time_LU_sparse = clock() - time_LU_sparse;

    int maxit = 10000;
    double tol = 1e-7;

    // jacobi dense
    xinit(x, rows);
    clock_t time_j_dense;
    time_j_dense = clock();
    test_jacobi_dense(*dense_mat, x, b, maxit, tol, false, print);
    time_j_dense = clock() - time_j_dense;

    // jacobi dense BLAS
    xinit(x, rows);
    clock_t time_j_dense_blas;
    time_j_dense_blas = clock();
    test_jacobi_dense(*dense_mat, x, b, maxit, tol, true, print);
    time_j_dense_blas = clock() - time_j_dense_blas;

    // jacobi sparse
    xinit(x, rows);
    clock_t time_j_sparse;
    time_j_sparse = clock();
    test_jacobi_sparse(*sparse_mat, x, b, maxit, tol, print);
    time_j_sparse = clock() - time_j_sparse;

    // gauss_seidel dense
    xinit(x, rows);
    clock_t time_gs_dense;
    time_gs_dense = clock();
    test_gauss_seidel_dense(*dense_mat, x, b, maxit, tol, 1., false, print);
    time_gs_dense = clock() - time_gs_dense;

    // gauss_seidel dense blas
    xinit(x, rows);
    clock_t time_gs_dense_blas;
    time_gs_dense_blas = clock();
    test_gauss_seidel_dense(*dense_mat, x, b, maxit, tol, 1., true, print);
    time_gs_dense_blas = clock() - time_gs_dense_blas;

    // gauss_seidel sparse
    xinit(x, rows);
    clock_t time_gs_sparse;
    time_gs_sparse = clock();
    test_gauss_seidel_sparse(*sparse_mat, x, b, maxit, tol, print);
    time_gs_sparse = clock() - time_gs_sparse;
    
    // thomas tridiagonal
    xinit(x, rows);
    clock_t time_th_tri;
    time_th_tri = clock();
    test_thomas_tri(*dense_mat, x, b, print);
    time_th_tri = clock() - time_th_tri;
    
    //cholesky dense
    xinit(x, rows);
    clock_t time_cho_dense;
    time_cho_dense = clock();
    test_cholesky_dense(*dense_mat, x, b, print);
    time_cho_dense = clock() - time_cho_dense;
    
    //cholesky sparse
    xinit(x, rows);
    clock_t time_cho_sparse;
    time_cho_sparse = clock();
    test_cholesky_sparse(*dense_mat, x, b, print);
    time_cho_sparse = clock() - time_cho_sparse;
    
    // show timings
    cout << "\nTimings: " << endl;
    cout << setw(40) << "Gauss Elimination dense: " << setw(10) << (float)time_ge / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU dense: " << setw(10) << (float)time_LU_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU sparse: " << setw(10) << (float)time_LU_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Jacobi dense: " << setw(10) << (float)time_j_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Jacobi dense blas: " << setw(10) << (float)time_j_dense_blas / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Jacobi sparse: " << setw(10) << (float)time_j_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Gauss-Seidel dense: " << setw(10) << (float)time_gs_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Gauss-Seidel dense blas: " << setw(10) << (float)time_gs_dense_blas / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Gauss-Seidel sparse: " << setw(10) << (float)time_gs_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Thomas tridiagonal: " << setw(10) << (float)time_th_tri / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Cholesky dense: " << setw(10) << (float)time_gs_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Cholesky sparse: " << setw(10) << (float)time_gs_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    
    delete dense_mat;
    delete sparse_mat;
    delete[] b;
    delete[] x;
}

template<class T>
void time_all_given(int n, T* mat_array, T* b, bool print)
{
    int rows(n), cols(n);
    Matrix<T>* dense_mat = new Matrix<T>(rows, cols, mat_array);
    
    if (rows <= 20) dense_mat->printMatrix();
    
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, 1, true);
    sparse_mat->fromDense(*dense_mat);
    if (rows <= 20) sparse_mat->printMatrix();
    
    auto* x = new double[rows * 1];
    
    // ***************  Tests  ****************
    // gauss elimination
    xinit(x, rows);
    clock_t time_ge;
    time_ge = clock();
    test_gauss_elimination(*dense_mat, x, b, print);
    time_ge = clock() - time_ge;

    // LU dense
    xinit(x, rows);
    clock_t time_LU_dense;
    time_LU_dense = clock();
    test_LU_dense(*dense_mat, x, b, print);
    time_LU_dense = clock() - time_LU_dense;

    // LU sparse
    xinit(x, rows);
    clock_t time_LU_sparse;
    time_LU_sparse = clock();
    test_LU_sparse(*sparse_mat, x, b, print);
    time_LU_sparse = clock() - time_LU_sparse;

    int maxit = 10000;
    double tol = 1e-7;

    // jacobi dense
    xinit(x, rows);
    clock_t time_j_dense;
    time_j_dense = clock();
    test_jacobi_dense(*dense_mat, x, b, maxit, tol, false, print);
    time_j_dense = clock() - time_j_dense;

    // jacobi dense BLAS
    xinit(x, rows);
    clock_t time_j_dense_blas;
    time_j_dense_blas = clock();
    test_jacobi_dense(*dense_mat, x, b, maxit, tol, true, print);
    time_j_dense_blas = clock() - time_j_dense_blas;

    // jacobi sparse
    xinit(x, rows);
    clock_t time_j_sparse;
    time_j_sparse = clock();
    test_jacobi_sparse(*sparse_mat, x, b, maxit, tol, print);
    time_j_sparse = clock() - time_j_sparse;

    // gauss_seidel dense
    xinit(x, rows);
    clock_t time_gs_dense;
    time_gs_dense = clock();
    test_gauss_seidel_dense(*dense_mat, x, b, maxit, tol, 1., false, print);
    time_gs_dense = clock() - time_gs_dense;

    // gauss_seidel dense blas
    xinit(x, rows);
    clock_t time_gs_dense_blas;
    time_gs_dense_blas = clock();
    test_gauss_seidel_dense(*dense_mat, x, b, maxit, tol, 1., true, print);
    time_gs_dense_blas = clock() - time_gs_dense_blas;

    // gauss_seidel sparse
    xinit(x, rows);
    clock_t time_gs_sparse;
    time_gs_sparse = clock();
    test_gauss_seidel_sparse(*sparse_mat, x, b, maxit, tol, print);
    time_gs_sparse = clock() - time_gs_sparse;
    
    // thomas tridiagonal
    xinit(x, rows);
    clock_t time_th_tri;
    time_th_tri = clock();
    test_thomas_tri(*dense_mat, x, b, print);
    time_th_tri = clock() - time_th_tri;
    
    //cholesky dense
    xinit(x, rows);
    clock_t time_cho_dense;
    time_cho_dense = clock();
    test_cholesky_dense(*dense_mat, x, b, print);
    time_cho_dense = clock() - time_cho_dense;
    
    //cholesky sparse
    xinit(x, rows);
    clock_t time_cho_sparse;
    time_cho_sparse = clock();
    test_cholesky_sparse(*dense_mat, x, b, print);
    time_cho_sparse = clock() - time_cho_sparse;
    
    // show timings
    cout << "\nTimings: " << endl;
    cout << setw(40) << "Gauss Elimination dense: " << setw(10) << (float)time_ge / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU dense: " << setw(10) << (float)time_LU_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU sparse: " << setw(10) << (float)time_LU_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Jacobi dense: " << setw(10) << (float)time_j_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Jacobi dense blas: " << setw(10) << (float)time_j_dense_blas / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Jacobi sparse: " << setw(10) << (float)time_j_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Gauss-Seidel dense: " << setw(10) << (float)time_gs_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Gauss-Seidel dense blas: " << setw(10) << (float)time_gs_dense_blas / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Gauss-Seidel sparse: " << setw(10) << (float)time_gs_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Thomas tridiagonal: " << setw(10) << (float)time_th_tri / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Cholesky dense: " << setw(10) << (float)time_gs_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "Cholesky sparse: " << setw(10) << (float)time_gs_sparse / CLOCKS_PER_SEC << " seconds" << endl;
    
    delete dense_mat;
    delete sparse_mat;
    delete[] x;
}

