#include <iostream>
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
			std::cout << this->values[i + j * this->cols] << " ";
		}
	}
	std::cout << std::endl;
};

template<class T>
T Matrix<T>::det()
{
	T* det_result(0);

	return *det_result;
}

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

