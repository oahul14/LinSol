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
void init_condition(T* x, T* b, T* bog, int n)
{
    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
        b[i] = bog[i];
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
void test_LU_dense(Matrix<T>& A, T* x, T* b, bool blas, bool print)
{
    if (!blas)
    {
        LU_dense(A, x, b);
        cout << "\nLU Solution Dense: ";
    }
    else
    {
        LU_dense_blas(A, x, b);
        cout << "\nLU Solution Dense BLAS: ";
    }
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
    int maxit = 10000;
    double tol = 1e-7;
    
    Matrix<double>* dense_mat = new Matrix<double>(rows, cols, true);
    if (mat_type == "dense") dense_mat->genRanDense(true);
    else if (mat_type == "sparse") dense_mat->genRanSparse(0.7, true);
    else if (mat_type == "tridiagonal") dense_mat->genRanTri(true);
    else
    {
        cout << "Matrix type not implemented. ";
        return;
    }
    Matrix<double>* dense_mat_hard = new Matrix<double>(*dense_mat);
    if (rows <= 20) dense_mat->printMatrix();
    
    auto* b = new double[rows * 1];
    //cout << "\nLinear system: RHS [";
    for (int i = 0; i < rows; i++)
    {
        b[i] = rand() % 10 + 1;
        //cout << b[i] << " ";
    }
    //cout << "]" << endl;
    
    auto* bog = new double[rows * 1];
    for (int i = 0; i < rows; i++)
    {
        bog[i] = b[i];
    }
    
    auto* x = new double[rows * 1];
    
    // ***************  Tests  ****************
    // gauss elimination
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c0 = new Matrix<double>(*dense_mat_hard);
    clock_t time_ge;
    time_ge = clock();
    test_gauss_elimination(*dense_mat_c0, x, b, print);
    time_ge = clock() - time_ge;
    delete dense_mat_c0;

    // LU dense
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c1 = new Matrix<double>(*dense_mat_hard);
    clock_t time_LU_dense;
    time_LU_dense = clock();
    test_LU_dense(*dense_mat_c1, x, b, false, print);
    time_LU_dense = clock() - time_LU_dense;
    delete dense_mat_c1;
    
    // LU dense blas
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c1b = new Matrix<double>(*dense_mat_hard);
    clock_t time_LU_dense_blas;
    time_LU_dense_blas = clock();
    test_LU_dense(*dense_mat_c1b, x, b, true, print);
    time_LU_dense_blas = clock() - time_LU_dense_blas;
    delete dense_mat_c1b;
    
    // LU sparse
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c2 = new Matrix<double>(*dense_mat_hard);
    auto* sparse_mat_c2 = new CSRMatrix<double>(rows, cols, 1, true);
    sparse_mat_c2->fromDense(*dense_mat_c2);
    clock_t time_LU_sparse;
    time_LU_sparse = clock();
    test_LU_sparse(*sparse_mat_c2, x, b, print);
    time_LU_sparse = clock() - time_LU_sparse;
    delete dense_mat_c2;
    delete sparse_mat_c2;

    // jacobi dense
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c3 = new Matrix<double>(*dense_mat_hard);
    clock_t time_j_dense;
    time_j_dense = clock();
    test_jacobi_dense(*dense_mat_c3, x, b, maxit, tol, false, print);
    time_j_dense = clock() - time_j_dense;
    delete dense_mat_c3;

    // jacobi dense BLAS
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c4 = new Matrix<double>(*dense_mat_hard);
    clock_t time_j_dense_blas;
    time_j_dense_blas = clock();
    test_jacobi_dense(*dense_mat_c4, x, b, maxit, tol, true, print);
    time_j_dense_blas = clock() - time_j_dense_blas;
    delete dense_mat_c4;
    
    // jacobi sparse
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c5 = new Matrix<double>(*dense_mat_hard);
    auto* sparse_mat_c5 = new CSRMatrix<double>(rows, cols, 1, true);
    sparse_mat_c5->fromDense(*dense_mat_c5);
    clock_t time_j_sparse;
    time_j_sparse = clock();
    test_jacobi_sparse(*sparse_mat_c5, x, b, maxit, tol, print);
    time_j_sparse = clock() - time_j_sparse;
    delete dense_mat_c5;
    delete sparse_mat_c5;

    // gauss_seidel dense
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c6 = new Matrix<double>(*dense_mat_hard);
    clock_t time_gs_dense;
    time_gs_dense = clock();
    test_gauss_seidel_dense(*dense_mat_c6, x, b, maxit, tol, 1., false, print);
    time_gs_dense = clock() - time_gs_dense;
    delete dense_mat_c6;

    // gauss_seidel dense blas
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c7 = new Matrix<double>(*dense_mat_hard);
    clock_t time_gs_dense_blas;
    time_gs_dense_blas = clock();
    test_gauss_seidel_dense(*dense_mat_c7, x, b, maxit, tol, 1., true, print);
    time_gs_dense_blas = clock() - time_gs_dense_blas;
    delete dense_mat_c7;
    
    // gauss_seidel sparse
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c8 = new Matrix<double>(*dense_mat_hard);
    auto* sparse_mat_c8 = new CSRMatrix<double>(rows, cols, 1, true);
    sparse_mat_c8->fromDense(*dense_mat_c8);
    clock_t time_gs_sparse;
    time_gs_sparse = clock();
    test_gauss_seidel_sparse(*sparse_mat_c8, x, b, maxit, tol, print);
    time_gs_sparse = clock() - time_gs_sparse;
    delete dense_mat_c8;
    delete sparse_mat_c8;
    
    // thomas tridiagonal
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c9 = new Matrix<double>(*dense_mat_hard);
    clock_t time_th_tri;
    time_th_tri = clock();
    test_thomas_tri(*dense_mat_c9, x, b, print);
    time_th_tri = clock() - time_th_tri;
    delete dense_mat_c9;
    
    //cholesky dense
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c10 = new Matrix<double>(*dense_mat_hard);
    clock_t time_cho_dense;
    time_cho_dense = clock();
    test_cholesky_dense(*dense_mat_c10, x, b, print);
    time_cho_dense = clock() - time_cho_dense;
    delete dense_mat_c10;
    
    //cholesky sparse
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c11 = new Matrix<double>(*dense_mat_hard);
    clock_t time_cho_sparse;
    time_cho_sparse = clock();
    test_cholesky_sparse(*dense_mat_c11, x, b, print);
    time_cho_sparse = clock() - time_cho_sparse;
    delete dense_mat_c11;
    
    // show timings
    cout << "\nTimings: " << endl;
    cout << setw(40) << "Gauss Elimination dense: " << setw(10) << (float)time_ge / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU dense: " << setw(10) << (float)time_LU_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU dense blas: " << setw(10) << (float)time_LU_dense_blas / CLOCKS_PER_SEC << " seconds" << endl;
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
    delete dense_mat_hard;
    delete[] bog;
    delete[] x;
}

template<class T>
void time_all_given(int n, T* mat_array, T* b, bool print)
{
    int rows(n), cols(n);
    int maxit = 10000;
    double tol = 1e-7;
    Matrix<T>* dense_mat = new Matrix<T>(rows, cols, mat_array);
    Matrix<T>* dense_mat_hard = new Matrix<T>(rows, cols, mat_array);
    if (rows <= 20) dense_mat->printMatrix();
    
    auto* x = new T[rows * 1];
    auto* bog = new T[rows * 1];
    for (int i = 0; i < n; i++)
    {
        bog[i] = b[i];
    }
    
    // ***************  Tests  ****************
    // gauss elimination
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c0 = new Matrix<T>(*dense_mat_hard);
    clock_t time_ge;
    time_ge = clock();
    test_gauss_elimination(*dense_mat_c0, x, b, print);
    time_ge = clock() - time_ge;
    delete dense_mat_c0;

    // LU dense
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c1 = new Matrix<T>(*dense_mat_hard);
    clock_t time_LU_dense;
    time_LU_dense = clock();
    test_LU_dense(*dense_mat_c1, x, b, print);
    time_LU_dense = clock() - time_LU_dense;
    delete dense_mat_c1;
    
    // LU dense blas
    init_condition(x, b, bog, n);
    Matrix<double>* dense_mat_c1b = new Matrix<double>(*dense_mat_hard);
    clock_t time_LU_dense_blas;
    time_LU_dense_blas = clock();
    test_LU_dense(*dense_mat_c1b, x, b, true, print);
    time_LU_dense_blas = clock() - time_LU_dense_blas;
    delete dense_mat_c1b;
    
    // LU sparse
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c2 = new Matrix<T>(*dense_mat_hard);
    auto* sparse_mat_c2 = new CSRMatrix<T>(rows, cols, 1, true);
    sparse_mat_c2->fromDense(*dense_mat_c2);
    clock_t time_LU_sparse;
    time_LU_sparse = clock();
    test_LU_sparse(*sparse_mat_c2, x, b, print);
    time_LU_sparse = clock() - time_LU_sparse;
    delete dense_mat_c2;
    delete sparse_mat_c2;

    // jacobi dense
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c3 = new Matrix<T>(*dense_mat_hard);
    clock_t time_j_dense;
    time_j_dense = clock();
    test_jacobi_dense(*dense_mat_c3, x, b, maxit, tol, false, print);
    time_j_dense = clock() - time_j_dense;
    delete dense_mat_c3;

    // jacobi dense BLAS
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c4 = new Matrix<T>(*dense_mat_hard);
    clock_t time_j_dense_blas;
    time_j_dense_blas = clock();
    test_jacobi_dense(*dense_mat_c4, x, b, maxit, tol, true, print);
    time_j_dense_blas = clock() - time_j_dense_blas;
    delete dense_mat_c4;
    
    // jacobi sparse
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c5 = new Matrix<T>(*dense_mat_hard);
    auto* sparse_mat_c5 = new CSRMatrix<T>(rows, cols, 1, true);
    sparse_mat_c5->fromDense(*dense_mat_c5);
    clock_t time_j_sparse;
    time_j_sparse = clock();
    test_jacobi_sparse(*sparse_mat_c5, x, b, maxit, tol, print);
    time_j_sparse = clock() - time_j_sparse;
    delete dense_mat_c5;
    delete sparse_mat_c5;

    // gauss_seidel dense
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c6 = new Matrix<T>(*dense_mat_hard);
    clock_t time_gs_dense;
    time_gs_dense = clock();
    test_gauss_seidel_dense(*dense_mat_c6, x, b, maxit, tol, 1., false, print);
    time_gs_dense = clock() - time_gs_dense;
    delete dense_mat_c6;

    // gauss_seidel dense blas
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c7 = new Matrix<T>(*dense_mat_hard);
    clock_t time_gs_dense_blas;
    time_gs_dense_blas = clock();
    test_gauss_seidel_dense(*dense_mat_c7, x, b, maxit, tol, 1., true, print);
    time_gs_dense_blas = clock() - time_gs_dense_blas;
    delete dense_mat_c7;
    
    // gauss_seidel sparse
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c8 = new Matrix<T>(*dense_mat_hard);
    auto* sparse_mat_c8 = new CSRMatrix<T>(rows, cols, 1, true);
    sparse_mat_c8->fromDense(*dense_mat_c8);
    clock_t time_gs_sparse;
    time_gs_sparse = clock();
    test_gauss_seidel_sparse(*sparse_mat_c8, x, b, maxit, tol, print);
    time_gs_sparse = clock() - time_gs_sparse;
    delete dense_mat_c8;
    delete sparse_mat_c8;
    
    // thomas tridiagonal
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c9 = new Matrix<T>(*dense_mat_hard);
    clock_t time_th_tri;
    time_th_tri = clock();
    test_thomas_tri(*dense_mat_c9, x, b, print);
    time_th_tri = clock() - time_th_tri;
    delete dense_mat_c9;
    
    //cholesky dense
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c10 = new Matrix<T>(*dense_mat_hard);
    clock_t time_cho_dense;
    time_cho_dense = clock();
    test_cholesky_dense(*dense_mat_c10, x, b, print);
    time_cho_dense = clock() - time_cho_dense;
    delete dense_mat_c10;
    
    //cholesky sparse
    init_condition(x, b, bog, n);
    Matrix<T>* dense_mat_c11 = new Matrix<T>(*dense_mat_hard);
    clock_t time_cho_sparse;
    time_cho_sparse = clock();
    test_cholesky_sparse(*dense_mat_c11, x, b, print);
    time_cho_sparse = clock() - time_cho_sparse;
    delete dense_mat_c11;
    
    // show timings
    cout << "\nTimings: " << endl;
    cout << setw(40) << "Gauss Elimination dense: " << setw(10) << (float)time_ge / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU dense: " << setw(10) << (float)time_LU_dense / CLOCKS_PER_SEC << " seconds" << endl;
    cout << setw(40) << "LU dense: " << setw(10) << (float)time_LU_dense_blas / CLOCKS_PER_SEC << " seconds" << endl;
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
    delete dense_mat_hard;
    delete[] x;
    delete[] bog;
}

