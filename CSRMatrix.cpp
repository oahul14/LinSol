#include <iostream>
#include "CSRMatrix.h"

// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
      this->values = new T[this->nnzs];
      this->row_position = new int[this->rows+1];
      this->col_index = new int[this->nnzs];
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // Delete the values array
   if (this->preallocated){
      delete[] this->row_position;
      delete[] this->col_index;
   }
   // The super destructor is called after we finish here
   // This will delete this->values if preallocated is true
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
   // std::cout << "Printing matrix" << std::endl;
   // std::cout << "Values: ";
   // for (int j = 0; j< this->nnzs; j++)
   // {  
   //    std::cout << this->values[j] << " ";      
   // }
   // std::cout << std::endl;
   // std::cout << "row_position: ";
   // for (int j = 0; j< this->rows+1; j++)
   // {  
   //    std::cout << this->row_position[j] << " ";      
   // }
   // std::cout << std::endl;   
   // std::cout << "col_index: ";
   // for (int j = 0; j< this->nnzs; j++)
   // {  
   //    std::cout << this->col_index[j] << " ";      
   // }
   // std::cout << std::endl;  

   auto* values = new T[this->rows * this->cols];
   for (int i = 0; i < this->rows; i++)
   {
      for (int j = 0; j < this->cols; j++)
      {
         values[i * this->rows + j] = 0;
      }
   }

   auto* row_index = new int[nnzs];

   int counter(0);
   for (int j = 1; j < this->rows + 1; j++)
   {
      int num = this->row_position[j] - this->row_position[j-1];
      for (int i = counter; i < counter + num; i++)
      {
         row_index[i] = j - 1;
      }
      counter += num;
   }
   
   for (int i = 0; i < this->nnzs; i++)
   {
      values[row_index[i] * this->rows + this->col_index[i]] = this->values[i];
   }
   delete[] row_index;

   std::cout << "\nPrinting sparse matrix:" << std::endl;
	for (int i = 0; i < this->rows; i++)
	{
		std::cout << std::endl;
		for (int j = 0; j < this->cols; j++)
		{
			// row-major ordering here
			std::cout << values[j + i * this->cols] << " ";
		}
	}
	std::cout << std::endl;
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(double *input, double *output)
{
   if (input == nullptr || output == nullptr)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   int val_counter = 0;
   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * input[this->col_index[val_index]];

      }
   }
}


// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output)
{

   // Check our dimensions match
   if (this->cols != mat_right.rows)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (output.values != nullptr) 
   {
      // Check our dimensions match
      if (this->rows != output.rows || this->cols != output.cols)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }      
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      std::cerr << "OUTPUT HASN'T BEEN ALLOCATED" << std::endl;

   }

   // HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??
}