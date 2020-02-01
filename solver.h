#pragma once
#include <iostream>
#include <tuple>
#include <math.h>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

using namespace std;

template<class T>
void LU_decomposition_pp(Matrix<T>& A, Matrix<T>& L, Matrix<T>& P_);

template<class T>
void backward_substitution(Matrix<T>& U, T* b, T* output);

template<class T>
void forward_substitution(Matrix<T>& L, T* b, T* output);

template<class T>
void LU_dense(Matrix<T>& A, T* x, T* b);

template<class T>
void gauss_elimination(Matrix<T>& A, T* x, T* b);

template<class T>
void gauss_seidel_dense(Matrix<T>& A, T* x, T* b, int maxit, double er, double urf);

template<class T>
void jacobi_dense(Matrix<T>& A, T* x, T* b, int maxit, double tolerance);

template<class T>
void jacobi_dense_blas(Matrix<T>& A, T* x, T* b, int maxit, double tolerance);

template<class T>
void gauss_seidel_sparse(CSRMatrix<T>& A, T* x, T* b, int maxit, double tolerance);

template<class T>
void jacobi_sparse(CSRMatrix<T>& A, T* x, T* b, int maxit, double tolerance);

template<class T>
void thomas(Matrix<T>& A, T* x, T* r);

template<class T>
void cholesky_fact(Matrix<T>& A, T* x, T* b);

template<class T>
void cholesky_dense(Matrix<T>& A, T* x, T* b);

template<class T>
void cholesky_sparse(Matrix<T>& A, T* x, T* b);
