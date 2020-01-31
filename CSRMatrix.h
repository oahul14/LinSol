#pragma once
#include "Matrix.h"
#include <memory>

template <class T>
class CSRMatrix: public Matrix<T>
{
public:

    // constructor where we want to preallocate ourselves
    CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
    // constructor where we already have allocated memory outside
    CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index);
    // destructor
    ~CSRMatrix();

    // Print out the values in our matrix
    virtual void printMatrix();
    void printDense();

    // Perform some operations with our matrix
    void matVecMult(double *input, double *output);
    // Perform some operations with our matrix
    void matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output);
    void fromDense(Matrix<T>& dense_in);
    
    void LU_decomposition(int* p_col);
    void forward_substitution(T* b, T* output);
    void backward_substitution(T* b, T* output);
    void LU_solver(T* b, T* output);

    // Explicitly using the C++11 nullptr here
    unique_ptr<int[]> row_position;
    unique_ptr<int[]> col_index;

    // How many non-zero entries we have in the matrix
    int nnzs=-1;

// Private variables - there is no need for other classes 
// to know about these variables
private:
   
};
