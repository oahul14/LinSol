#pragma once
#include <iostream>
#include "CSRMatrix.h"
#include <cmath>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <vector>

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
      this->values.reset(new T[this->nnzs]);
      this->row_position.reset(new int[this->rows+1]);
      this->col_index.reset(new int[this->nnzs]);
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs)
{
    this->row_position.reset(row_position);
    this->col_index.reset(col_index);
}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix: sparse" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->values[j] << " ";      
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {  
      std::cout << this->row_position[j] << " ";      
   }
   std::cout << std::endl;   
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->col_index[j] << " ";      
   }
   std::cout << std::endl << std::endl;  
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

template<class T>
void CSRMatrix<T>::fromDense(Matrix<T>& dense_in)
{
    int m = dense_in.rows;
    int n = dense_in.cols;
    this->rows = m;
    this->cols = n;
    this->row_position.reset(new int[m + 1]);

    int counter = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (dense_in.values[i * n + j] != 0) counter++;
        }
        this->row_position[i + 1] = counter;
    }
    this->row_position[0] = 0;
    this->nnzs = this->row_position[m];
    
    vector<T> temp_value;
    vector<int> temp_col;

    this->values.reset(new T[this->nnzs]);
    this->col_index.reset(new int[this->nnzs]);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (dense_in.values[i * n + j] != 0)
            {
                //cout << dense_in.values[i * n + j] << endl;
                temp_value.push_back(dense_in.values[i * n + j]);
                temp_col.push_back(j);
            }
        }
    }
    for (int i = 0; i < this->nnzs; i++)
    {
        //cout << temp_value[i] << endl;
        this->values[i] = temp_value[i];
        this->col_index[i] = temp_col[i];
    }
}

template<class T>
void CSRMatrix<T>::LU_decomposition(int* p_col)
{
    for (int i = 0; i < this->rows; i++)
    {
        cout << "Loop over rows+1: i = " << i << endl;
        cout << "---------------" << endl;
        // get the non-zeros number for i-1th row
        int nnzs_current = this->row_position[i + 1] - this->row_position[i];

        // for each diagonal entry:
        // if possible, find the largest row index below
        // so should loop over col_index, which is also the loop
        // for CSRMatrix values
        double max_value(0);
        int max_value_index(i);
        int nnzs_max_row;
        int max_value_start_index;
        int current_start_index = this->row_position[i];

        // find max value and its index in the value array
        //cout << "this->row_position[i]: " << this->row_position[i] << endl;
        for (int k = this->row_position[i]; k < this->nnzs; k++)
        {
            //cout << "loop over nnzs: k = " << k << endl;
            // find the MAX value index in values 
            // for i - 1th col, and check if the col index if the same as rows'
            if (this->col_index[k] == i)
            {
                //cout << "--------: " << max_value << " at col:" << this->col_index[k] << endl;
                //cout << "current row: " << i << endl;
                //cout << "Find max value: " << this->values[k] << " and its col: " << this->col_index[k] << endl;
                if (abs(this->values[k]) > abs(max_value))
                {
                    max_value_index = k;
                    max_value = this->values[k];
                }
            }
        }

        //cout << "max value index: " << max_value_index << ", ";
        //cout << "Max value: " << max_value << endl;

        for (int m = 0; m < this->rows; m++)
        {
            if ((max_value_index >= this->row_position[m]) && (max_value_index < this->row_position[m + 1]))
            {
                //cout << "i: " << i << "  m: " << m << endl;
                //cout << p_col[i] << "   " << p_col[m] << endl;
                unique_ptr<int> temp_p_col(new int);
                *temp_p_col = p_col[m];
                p_col[m] = p_col[i];
                p_col[i] = *temp_p_col;
                //cout << p_col[i] << "   " << p_col[m] << endl;
            }
        }


        //cout << "Getting nnzs for max value row: " << endl;

        // as the matrix is diagonal positive
        // so dont need to worry about 0s on the diagnal after partial pivoting
        for (int j = 0; j < this->rows; j++)
        {
            if ((max_value_index >= this->row_position[j]) & (max_value_index <= this->row_position[j + 1]))
            {
                max_value_start_index = this->row_position[j];
                //cout << this->row_position[j] << " - " << this->row_position[j+1] << endl;
                nnzs_max_row = this->row_position[j + 1] - this->row_position[j];
            }
        }
        //cout << "How many non-zeros in the current row: " << nnzs_current << endl;
        //cout << "How many non-zeros in the max value row: " << nnzs_max_row << endl;
        //cout << endl;

        // create pointers to store current row and max value row non-zeros and their col indices;
        /*vector<double> temp_current_value(nnzs_current);
        vector<int> temp_current_col(nnzs_current);
        vector<double> temp_max_value(nnzs_max_row);
        vector<int> temp_max_col(nnzs_max_row);*/
        unique_ptr<double[]> temp_current_value(new double[nnzs_current]);
        unique_ptr<int[]> temp_current_col(new int[nnzs_current]);
        unique_ptr<double[]> temp_max_value(new double[nnzs_max_row]);
        unique_ptr<int[]> temp_max_col(new int[nnzs_max_row]);

        int nnzs_in_between = max_value_start_index - current_start_index - nnzs_current;
        int in_between_start_index = current_start_index + nnzs_current;
        //cout << "Current row start index: " << current_start_index << endl;
        //cout << "Max value row start index: " << max_value_start_index << endl;
        //cout << "This is how many non-zeros in between two rows: " << nnzs_in_between << endl;
        //cout << "And the in between start index is: " << in_between_start_index << endl;
        
        //cout << "Below is the current row non-zeros: " << endl;
        // store current row non-zeros value/col
        for (int in = 0; in < nnzs_current; in++)
        {
            temp_current_value[in] = this->values[current_start_index + in];
            temp_current_col[in] = this->col_index[current_start_index + in];
            //cout << "Value: " << temp_current_value[in] << " Col: " << temp_current_col[in] << endl;
        }

        //cout << "Below is the max value row non-zeros: " << endl; 
        // store max value row non-zeros value/col
        for (int in = 0; in < nnzs_max_row; in++)
        {
            //cout << "Start index for max value row: " << max_value_start_index << endl;
            temp_max_value[in] = this->values[max_value_start_index + in];
            temp_max_col[in] = this->col_index[max_value_start_index + in];
            //cout << "Value: " << temp_max_value[in] << " Col: " << temp_max_col[in] << endl;
        }

        // check if the rows need to be changed
        // if needed the swapped row should be ahead of the swapping row
        if (max_value_start_index > current_start_index)
        {
            // check if there is in between value, if there is:
            if (in_between_start_index != max_value_start_index)
            {
                // create pointers to store temp value/col for in between elements
                unique_ptr<double[]> temp_between_value(new double[nnzs_in_between]);
                unique_ptr<int[]> temp_between_col(new int[nnzs_in_between]);
                //cout << "Below is the in between non-zeros: " << endl;
                // store in between value/col
                for (int in = 0; in < nnzs_in_between; in++)
                {
                    temp_between_value[in] = this->values[in_between_start_index + in];
                    temp_between_col[in] = this->col_index[in_between_start_index + in];
                    //cout << "Value: " << temp_between_value[in] << " Col: " << temp_between_col[in] << endl;
                }

                // refill value/col: max value row -> current row
                // cout << endl << "Now the matrix value is changed: " << endl;
                for (int ivc = 0; ivc < nnzs_max_row; ivc++)
                {
                    this->values[current_start_index + ivc] = temp_max_value[ivc];
                    this->col_index[current_start_index + ivc] = temp_max_col[ivc];
                    //cout << "New max value/col position: " << current_start_index + ivc << endl;
                }

                // refill value/col: change pos of in between ones
                for (int ivc = 0; ivc < nnzs_in_between; ivc++)
                {
                    this->values[current_start_index + nnzs_max_row + ivc] = temp_between_value[ivc];
                    this->col_index[current_start_index + nnzs_max_row + ivc] = temp_between_col[ivc];
                    //cout << "New in between value/col position: " << current_start_index + nnzs_max_row + ivc << endl;
                }

                // refill value/col: current now behind in between ones
                for (int ivc = 0; ivc < nnzs_current; ivc++)
                {
                    this->values[current_start_index + nnzs_max_row + nnzs_in_between + ivc] = temp_current_value[ivc];
                    this->col_index[current_start_index + nnzs_max_row + nnzs_in_between + ivc] = temp_current_col[ivc];
                    //cout << "New current value/col position: " << current_start_index + nnzs_max_row + nnzs_in_between + ivc << endl;
                }

                // if there is a difference between nnzs of swapping two rows
                // if not, no need to change row position
                
                if (nnzs_current != nnzs_max_row)
                {   
                    int end_index;
                    for (int o = 0; o < this->rows; o++)
                    {
                        if ((this->row_position[o + 1] > max_value_index) && (this->row_position[o] <= max_value_index))
                        {
                            end_index = this->row_position[o + 1];
                            //cout << end_index;
                        }
                    }
                    int nnzs_diff = nnzs_max_row - nnzs_current;
                    for (int irp = i + 1; irp < this->rows; irp++)
                    {
                        if (this->row_position[irp] < end_index)
                        {
                            this->row_position[irp] += nnzs_diff;
                        }
                        else break;
                    }
                }//cout << "find 5.5 and -1.5: " << this->values[8] << "  " << this->values[9] << endl;
            }
            // if no elements in between
            else
            {
                // refill values/col: max value -> current 
                for (int ivc = 0; ivc < nnzs_max_row; ivc++)
                {
                    this->values[current_start_index + ivc] = temp_max_value[ivc];
                    this->col_index[current_start_index + ivc] = temp_max_col[ivc];
                    //cout << "New max value/col: " << this->values[current_start_index + ivc] << " " << this->col_index[current_start_index + ivc] << " at position: " << current_start_index + ivc << endl;
                }

                // refill values/col: current -> max value
                for (int ivc = 0; ivc < nnzs_current; ivc++)
                {
                    this->values[current_start_index + nnzs_max_row + ivc] = temp_current_value[ivc];
                    this->col_index[current_start_index + nnzs_max_row + ivc] = temp_current_col[ivc];
                    //cout << "New current value/col: " << this->values[current_start_index + nnzs_max_row + ivc] << " " << this->col_index[current_start_index + nnzs_max_row + ivc] << " at position: " << current_start_index + nnzs_max_row + ivc << endl;
                }
                // check if row position needs to be changed
                // like before, only need to change if the swapping rows have difference nnzs
                if (nnzs_current != nnzs_max_row)
                {
                    int nnzs_diff = nnzs_max_row - nnzs_current;
                    this->row_position[i + 1] += nnzs_diff;
                }
            }
        }
        //cout << "\nSwapped matrix: " << endl;
        //this->printMatrix();
        //cout << endl;
        //this->printDense();
        //cout << endl;

        // if (max_value != 0)
        double s;
        int found_start_index;
        int found_next_start_index;
        int nnzs_found = 0;
        int nnzs_found_iright;
        int found_iright_start_index;
        int nnzs_current_iright = nnzs_max_row;
        int current_iright_start_index = current_start_index;
        int counter = 0;
        for (int k = this->row_position[i + 1]; k < this->nnzs; k++)
        {
            //cout << "This->row_position[i + 1]: " << this->row_position[i+1] << endl;
            if (this->col_index[k] == i)
            {
                cout << "\n\nK loop: " << k << " check col at k: " << this->col_index[k] << " and its value: " << this->values[k] << " with max nnzs: " << this->row_position[this->rows] << endl;
                this->values[k] /= max_value;
                // from now on the lower part will become L except for the diagonal 1s
                // store the constant divisor
                s = this->values[k];
                //cout << "s: " << s << " which is now at: " << k << "th position" << endl;

                for (int irp = 0; irp < this->rows; irp++)
                {
                    //cout << "***" << this->row_position[irp] << "***" << endl;
                    if ((this->row_position[irp] <= k) && (this->row_position[irp + 1] > k))
                    {
                        nnzs_found = this->row_position[irp + 1] - this->row_position[irp];
                        found_start_index = this->row_position[irp];
                        found_next_start_index = this->row_position[irp + 1];
                        //cout << "nnzs found: " << nnzs_found << " and found start index: " << found_start_index << endl;
                    }
                }

                found_iright_start_index = k + 1;
                nnzs_found_iright = found_next_start_index - found_iright_start_index;
                nnzs_current_iright = nnzs_max_row;
                cout << "\nFound iright start index: " << found_iright_start_index << " with value: " << this->values[found_iright_start_index] << " and nnzs found iright: " << nnzs_found_iright << endl;

                current_iright_start_index = current_start_index;

                for (int ivc = current_start_index; ivc < current_start_index + nnzs_max_row; ivc++)
                {

                    if (this->col_index[ivc] <= i)
                    {
                        current_iright_start_index++;
                        nnzs_current_iright--;
                    }
                }
                cout << "Current iright start index: " << current_iright_start_index << " with value: " << this->values[current_iright_start_index] << " and nnzs current iright: " << nnzs_current_iright << endl;
                
                //this->printMatrix();
//                cout << endl;
                // if the current row has 0 non-zeros, no change would be made
                if (nnzs_current_iright > 0)
                {
                    // find what will be the nnzs of found row
                    int nnzs_found_iright_new = nnzs_current_iright + nnzs_found_iright;
                    for (int na = current_iright_start_index; na < current_iright_start_index + nnzs_current_iright; na++)
                    {
                        for (int nb = found_iright_start_index; nb < found_iright_start_index + nnzs_found_iright; nb++)
                        {
                            if (this->col_index[na] == this->col_index[nb]) nnzs_found_iright_new--;
                        }
                    }
                    cout << "Temp store length: " << nnzs_found_iright_new << endl;
                    unique_ptr<double[]> temp_modified_value(new double[nnzs_found_iright_new]);
                    unique_ptr<int[]> temp_modified_col(new int[nnzs_found_iright_new]);

                    int num_add = nnzs_found_iright_new - nnzs_found_iright;
                    cout << "num_add: " << num_add << endl;

                    vector<double> temp_value(0);
                    vector<int> temp_col(0);
                    cout << "vector sizes 1: " << temp_col.size() << "  " << temp_value.size() << endl;
                    for (int c = current_iright_start_index; c < current_iright_start_index + nnzs_current_iright; c++)
                    {
                        int col_c = this->col_index[c];
                        bool flag = false;
                        int equal_index;
                        for (int f = found_iright_start_index; f < found_iright_start_index + nnzs_found_iright; f++)
                        {
                            if (this->col_index[f] == col_c)
                            {
                                equal_index = f;
                                flag = true;
                            }
                        }
                        if (flag)
                        {
                            temp_value.push_back(this->values[equal_index] - this->values[c] * s);
                            temp_col.push_back(this->col_index[equal_index]);
                        }
                        else
                        {
                            temp_value.push_back(0 - this->values[c] * s);
                            temp_col.push_back(this->col_index[c]);
                        }
                    }
                    cout << "vector sizes 2: " << temp_col.size() << "  " << temp_value.size() << endl;
                    
//                    for (int f = found_iright_start_index + nnzs_found_iright - 1; f >= found_iright_start_index; f--)
//                    {
//                        if (!count(temp_col.begin(), temp_col.end(), this->col_index[f]))
//                        {
//                            //cout << this->values[f] << endl;
//                            for (int c = current_iright_start_index; c < current_iright_start_index + nnzs_current_iright; c++)
//                            {
//                                //cout << this->col_index[c-1] << "-----" << this->col_index[c] << endl;
//                                if (this->col_index[f] > this->col_index[c])
//                                {
//                                    if (!count(temp_col.begin(), temp_col.end(), this->col_index[f]))
//                                    {
//                                        temp_value.push_back(this->values[f]);
//                                        temp_col.push_back(this->col_index[f]);
//                                    }
//                                }
//                                else if ((this->col_index[f] > this->col_index[c-1]) && (this->col_index[f] < this->col_index[c]))
//                                {
//                                    if (!count(temp_col.begin(), temp_col.end(), this->col_index[f]))
//                                    {
//                                        temp_value.insert(temp_value.begin() + (c - current_iright_start_index ), this->values[f]);
//                                        temp_col.insert(temp_col.begin() + (c - current_iright_start_index ), this->col_index[f]);
//                                        cout << "Middle Temp value / col: " << this->values[f] << " / " << this->col_index[f] << endl;
//                                    }
//                                }
//                            }
//                        }
//                    }
                    
                    // get full record of found non-zeros at right of the row NO. i
                    // first get the overlapped points and calculate it
                    // then for
                    for (int f = found_iright_start_index; f < found_iright_start_index + nnzs_found_iright; f++)
                    {
                        if (!count(temp_col.begin(), temp_col.end(), this->col_index[f]))
                        {
                            temp_col.push_back(this->col_index[f]);
                            temp_value.push_back(this->values[f]);
                        }
                    }
                    cout << "vector sizes 3: " << temp_col.size() << "  " << temp_value.size() << endl;
                    if (temp_col.size() != nnzs_found_iright_new)
                    {
                        cerr << "Invalid new length of found nnzs iright." << endl;
                        return;
                    }
                    
                    unique_ptr<int> tempc(new int);
                    unique_ptr<T> tempv(new T);
                    for (int ik = 1; ik < nnzs_found_iright_new; ik++)
                    {
                        for (int jk = 1; jk < nnzs_found_iright_new - ik; jk++)
                        {
                            if (temp_col[jk] > temp_col[jk + 1])
                            {
                                *tempc = temp_col[jk];
                                *tempv = temp_value[jk];
                                
                                temp_col[jk] = temp_col[jk + 1];
                                temp_value[jk] = temp_value[jk + 1];
                                
                                temp_col[jk + 1] = *tempc;
                                temp_value[jk + 1] = *tempv;
                            }
                        }
                    }
                    
                    if (num_add < 0)
                    {
                        cerr << "Invalid num_add." << endl;
                        return;
                    }
                    
                    if (found_start_index + nnzs_found_iright)
                    
                    for (int ii = 0; ii < temp_value.size(); ii++)
                    {
                        cout << "Temp value / col: " << temp_value[ii] << " / " << temp_col[ii] << endl;
                    }

                    // if num_add == 0: no nnzs change, only change in values for iright at found row
                    if (num_add == 0)
                    {
                        for (int neg = 0; neg < nnzs_found_iright; neg++)
                        {
                            this->values[found_iright_start_index + neg] = temp_value[neg];
                            this->col_index[found_iright_start_index + neg] = temp_col[neg];
                        }
                    }
                    else
                    {
                        vector<double> new_value(this->nnzs + num_add);
                        vector<int> new_col(this->nnzs + num_add);

                        // assign unchanged values/col numbers
                        for (int ni = 0; ni < found_iright_start_index; ni++)
                        {
                            new_value[ni] = this->values[ni];
                            new_col[ni] = this->col_index[ni];
                            cout << "Begin Part: value / col: " << new_value[ni] << " / " << new_col[ni] << endl;
                        }

                        // insert temp modified values and col from the start of the found iright non-zeros
                        for (int temp_ni = 0; temp_ni < nnzs_found_iright_new; temp_ni++)
                        {
                            new_value[found_iright_start_index + temp_ni] = temp_value[temp_ni];
                            new_col[found_iright_start_index + temp_ni] = temp_col[temp_ni];
                            cout << "Add Part: value / col: " << temp_value[temp_ni] << " / " << temp_col[temp_ni] << endl;
                        }

                        // push back the rest of non-zeros
                        for (int ni = found_iright_start_index + nnzs_found_iright; ni < this->nnzs; ni++)
                        {
                            new_value[ni + num_add] = this->values[ni];
                            new_col[ni + num_add] = this->col_index[ni];
                            cout << "Rest Part: value / col: " << new_value[ni+num_add] << " / " << new_col[ni+num_add] << endl;
                        }

                        for (int irp = 0; irp < this->rows; irp++)
                        {
                            if ((found_iright_start_index >= this->row_position[irp]) && (found_iright_start_index <= this->row_position[irp + 1]))
                            {
                                for (int irp_ = irp + 1; irp_ < this->rows + 1; irp_++)
                                {
                                    this->row_position[irp_] += num_add;
                                }
                            }
                        }

                        this->nnzs += num_add;
                        this->values.reset(new double[this->nnzs + num_add]);
                        this->col_index.reset(new int[this->nnzs + num_add]);

                        for (int pp = 0; pp < this->nnzs; pp++)
                        {
                            this->values[pp] = new_value[pp];
                            this->col_index[pp] = new_col[pp];
                            //cout << "New value: " << new_value[pp] << " and new col: " << new_col[pp] << endl;
                        }
                    }
                }
            
                this->printMatrix();
                cout << endl;
                counter++;
                if (counter == this->rows)
                {
                    cerr << "Unable to decompose, please try another matrix" << endl;
                    return;
                }
            }
        
        
        }
        //cout << endl;
//        this->printMatrix();
        //cout << endl;
        //cout << endl;
        //this->printDense();
        //cout << endl << endl;
        if (i == this->rows - 2) break;
    }
    //this->printDense();
//    cout << endl;
//    this->printMatrix();
    
}

template<class T>
void CSRMatrix<T>::forward_substitution(T* b, T* output)
{   
    //cout << "Calculating forward: " << endl;
    for (int r = 0; r < this->rows; r++)
    {
        //cout << r << ": ";
        int start_index = this->row_position[r];
        int end_index;
        for (int i = start_index; i < this->row_position[r + 1]; i++)
        {
            if (this->col_index[i] == r)
            {
                end_index = i;
            }
        }
        T val;
        int col;
        double s(0);
        for (int k = start_index; k < end_index; k++)
        {
            col = this->col_index[k];
            if (this->col_index[k] != r) val = this->values[k];
            else val = 1;
            for (int c = 0; c < this->rows; c++)
            {
                if (c == col) s += val * output[c];
            }
            //cout << s << " ";
            //cout << val << "(" << col << ")" << " ";
        }
        //cout << "b[r]: " << b[r] << " [k, k]: " << this->values[end_index] << " at index: " << end_index;
        output[r] = (b[r] - s) / 1;
        //cout << "output: " << output[r];
        //cout << endl;
    }

//    cout << "Printing forward sub vector output: " << endl;
//    for (int i = 0; i < this->rows; i++)
//    {
//        cout << " " << output[i];
//    }
//    cout << endl;
}

template<class T>
void CSRMatrix<T>::backward_substitution(T* b, T* output)
{
    for (int r = this->rows; r > 0; r--)
    {
        int start_index;
        int end_index = this->row_position[r] - 1;
        //cout << end_index << " with value: " << this->values[end_index] << endl;
        //cout << "row_position[r-1]: " << this->row_position[r-1];
        for (int i = end_index; i >= this->row_position[r - 1]; i--)
        {
            if (this->col_index[i] == r - 1)
            {
                start_index = i;
                //cout << this->col_index[i] << endl;
            }
        }
        T val;
        int col;
        double s(0);
        for (int k = start_index + 1; k <= end_index; k++)
        {
            col = this->col_index[k];
            val = this->values[k];
            for (int c = 0; c < this->rows; c++)
            {
                if (c == col) s += val * output[c];
            }
            //cout << s << " ";
            //cout << val << "(" << col << ")" << " ";
        }
        //cout << b[r - 1] << endl;
        output[r - 1] = (b[r - 1] - s) / this->values[start_index];
        //cout << endl;
    }

//    cout << "Printing backward sub vector output: " << endl;
//    for (int i = 0; i < this->rows; i++)
//    {
//        cout << " " << output[i];
//    }
//    cout << endl;
}

template<class T>
void CSRMatrix<T>::LU_solver(T* b, T* output)
{
    // initialising an identity matrix to record row swapping
    auto* p_col = new int[this->rows * 1];
    for (int c = 0; c < this->rows; c++)
    {
        p_col[c] = c;
    }
    this->LU_decomposition(p_col);

    // get P @ b
    auto* patb = new double[this->rows * 1];
    for (int i = 0; i < this->rows; i++)
    {
        patb[i] = b[p_col[i]];
    }
    auto* y = new double[this->rows * 1];
    this->forward_substitution(patb, y);


    this->backward_substitution(y, output);

    cout << "\nLU Solution Sparse: " << endl;
    for (int i = 0; i < this->rows; i++)
    {
        cout << "x" << i << ": " << output[i] << endl;
    }

    delete[] p_col;
    delete[] patb;
    delete[] y;
}
