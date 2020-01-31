#include <iostream>
#include <algorithm>
#include <memory>
#include <iomanip>
#include "Matrix.h"

// constructor 1
template<class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate): 
	rows(rows), cols(cols), size_of_values(rows * cols), preallocated(preallocate)
{
	// if we want to handle memory ourselves
	if (this->preallocated)
	{
		this->values.reset(new T[this->size_of_values]);
	 }
};

//constructor 2
template<class T>
Matrix<T>::Matrix(int rows, int cols, T* values_ptr): 
	rows(rows), cols(cols), size_of_values(rows * cols) 
{
	this->values.reset(values_ptr);
};

// copy constructor
template<class T>
Matrix<T>::Matrix(Matrix& B)
{
	this->cols = B.cols;
	this->rows = B.rows;
	this->values.reset(B.values);
	this->size_of_values = B.size_of_values;
}

// destructor
template<class T>
Matrix<T>::~Matrix()
{
	// if (this->preallocated){
	// 	delete[] this->values;
	// }
};

template<class T>
void Matrix<T>::printValues() 
{
	std::cout << "\nPrinting values" << std::endl;
	for (int i = 0; i < this->size_of_values; i++)
	{
		std::cout << this->values[i] << " ";
	}
	std::cout << std::endl;
};

template<class T>
void Matrix<T>::printMatrix() 
{
	std::cout << "\nPrinting matrix" << std::endl;
	for (int j = 0; j < this->rows; j++)
	{
		std::cout << std::endl;
		for (int i = 0; i < this->cols; i++)
		{
			// row-major ordering here
			std::cout << setw(3) << this->values[i + j * this->cols] << " ";
		}
	}
	std::cout << std::endl;
};

template<class T>
void Matrix<T>::transpose(Matrix<T> &itself)
{
	const int m = this->rows;
	const int n = this->cols;  
	auto* data_cm = new T[m * n];
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
		{
			const T rm = this->values[i * n + j];
			data_cm[j * m + i] = rm;
		}
	}
	
	itself.cols = m;
	itself.rows = n;
	for (int i = 0; i < m * n; i++)
	{
		itself.values[i] = data_cm[i];
	}
	delete[] data_cm;
}

template<class T>
void Matrix<T>::genRanDense(bool dom)
{
    //generate a dominant dense matrix
    // if dom, means dominant on the diagonal
    int rows = this->rows;
    int cols = this->cols;
    shared_ptr<double[]> ran_mat(new double[rows * cols]); //cannot use unique_ptr since that the pointer will points to other directions
    srand((unsigned) time(NULL)); //generate different random values
    if (dom) ran_mat[0] = rand() % 100+100;
    else ran_mat[0] = rand() % 10+1;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (j > i)
            {
                ran_mat[i * cols + j] = rand() % 5 + 0;//set the upper triangle's value from 1 to 5
            }
            else if (i > j)
            {
                ran_mat[i * cols + j] = ran_mat[j * cols + i];
            }
            else if (i=j) //the diagonal values
            {
                if (dom) ran_mat[0] = ran_mat[i * cols + j] = rand() % 100 + 100;
                else ran_mat[0] = ran_mat[i * cols + j] = rand() % 10 + 1;
            }
        }
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            this->values[i * cols + j] = ran_mat[i * cols + j];
//            cout << this->values[i * cols + j] << "  ";
        }
//        cout << endl;
    }
}

template<class T>
void Matrix<T>::genRanSparse(double sparsity, bool dom)
{
    //generate sparse matrix
    // if dom, means dominant on the diagonal
    int rows = this->rows;
    int cols = this->cols;
    //sparsity < 1
    if (sparsity >= 1)
    {
        cerr << "Sparsity should be less than 1." << endl;
        return;
    }
    else if (sparsity < 0.5)
    {
        cerr << "Warning: Sparsity too small." << endl;
    }
    
    int non_zeros = (rows * cols - rows) / 2 * (1 - sparsity);
    
    shared_ptr<double[]> ran_mat(new double[rows * cols]); //cannot use unique_ptr since that the pointer will points to other directions
    srand((unsigned)time(NULL)); //generate different random values
    
    for (int i = 0; i < cols * rows; i++)
    {
        ran_mat[i] = 0;
    }
    
    int k = 0;
    while (k < non_zeros)    //put the random nonzeros into random positions
    {
        int randi = rand() % rows + 0;
        int randj = rand() % cols + 0;
        if (randi != randj)
        {
            ran_mat[randi * cols + randj] = rand() % 10 + 5;
            k++;
        }
    }
    unique_ptr<double[]> ran_mat_T(new double[rows * cols]);
    unique_ptr<double[]> ran_sparse(new double[rows * cols]);
    
    for (int j = 0; j < cols; j++)
    {
        for (int i = 0; i < rows; i++)
        {
            const double rm = ran_mat[i * cols + j];
            ran_mat_T[j * rows + i] = rm;
        }
    }
    for (int i = 0; i < rows; i ++)
    {
        for (int j = 0; j < rows; j++)
        {
            ran_sparse[i * cols + j] = ran_mat[i * cols + j] + ran_mat_T[i * cols + j];
//            cout << setw(4) << ran_sparse[i * cols + j] << " ";
        }
//        cout << endl;
    }
    
    if (dom) ran_sparse[0] = rand() % 100 + 100;
    else ran_sparse[0] = rand() % 10 + 1;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (i = j) //the diagonal values
            {
                if (dom) ran_sparse[i * cols + j] = rand() % 100 + 50; //set the dominant matrix
                else ran_sparse[i * cols + j] = rand() % 10 + 1;
            }
        }
    }
//    for (int i = 0; i < rows; i++)
//    {
//        for (int j = 0; j < rows; j++)
//        {
//            cout << ran_sparse[i * cols + j] << "  ";
//        }
//        cout << endl;
//    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            this->values[i * cols + j] = ran_sparse[i * cols + j];
        }
    }
    
}

template<class T>
void Matrix<T>::matMatMult(Matrix& mat_right, Matrix& output)
{
	if (this->cols != mat_right.rows)
	{
		std::cerr << "input dimensions dont match" << std::endl;
		return;
	}
	else
	{
		output.values = new T[this->rows * this->cols];
	}

	for (int i = 0; i < output.size_of_values; i++)
	{
		output.values[i] = 0;
	}

	for (int i = 0; i < this->rows; i++)
	{
		for (int k = 0; k < this->cols; k++)
		{
			for (int j = 0; j < mat_right.cols; j++)
			{
				output.values[i * output.cols + j] += 
					this->values[i * this->cols + k] * 
					mat_right.values[k * mat_right.cols + j];
			}
		}
	}
};

template<class T>
void Matrix<T>::matVecMult(T* vec, T* output)
{
	for (int i = 0; i < rows; i++)
	{
		output[i] = 0;
	}

	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < this->cols; j++)
		{
			output[i] += this->values[i * cols + j] * vec[j];
		}
	}
};

template<class T>
void Matrix<T>::matMatMult_colMajor(Matrix& mat_right, Matrix& output)
{
	if (this->cols != mat_right.rows)
	{
		std::cerr << "input dimensions dont match" << std::endl;
		return;
	}

	for (int i = 0; i < output.size_of_values; i++)
	{
		output.values[i] = 0;
	}

	for (int j = 0; j < mat_right.cols; j++)
	{
		for (int k = 0; k < this->cols; k++)
		{
			for (int i = 0; i < this->rows; i++)
			{
				output.values[i * output.cols + j] +=
					this->values[i * this->cols + k] *
					mat_right.values[k * mat_right.cols + j];
			}
		}
	}
};

