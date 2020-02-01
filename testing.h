#pragma once
#include <iostream>
#include <Accelerate/Accelerate.h>
#include <memory>
#include <ctime>
#include "Matrix.h"
#include "CSRMatrix.h"
#include "solver.h"
#include "solver.cpp"
#include <algorithm>
#include <iomanip>

using namespace std;

template<class T>
void printVec(T* x, int size, bool print);

template<class T>
void test_LU_dense(Matrix<T>& A, T* b, bool print);

template<class T>
void test_LU_sparse(CSRMatrix<T>& A, T* b, bool print);

template <class T>
void test_gauss_elimination(Matrix<T>& A, T* b, bool print);

template<class T>
void test_jacobi_dense(Matrix<T>& A, T* b, int maxit, double tol, bool blas, bool print);

template <class T>
void test_jacobi_sparse(CSRMatrix<T>& A, T* b, int maxit, double tol, bool print);

template<class T>
void test_gauss_seidel_dense(Matrix<T>& A, T* b, int maxit, T er, T urf, bool blas, bool print);

template <class T>
void test_gauss_seidel_sparse(CSRMatrix<T>& A, T* b, int maxit, double tol, bool print);

template<class T>
void test_thomas_tri(Matrix<T>& A, T* b, bool print);

template <class T>
void test_cholesky_dense(Matrix<T>& A, T* b, bool print);

template <class T>
void test_cholesky_sparse(Matrix<T>& A, T* b, bool print);

void time_all_rand(int n, string mat_type, bool print);

template<class T>
void time_all_given(T* mat_array, T* b, bool);

