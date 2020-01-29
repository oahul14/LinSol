#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "solver.h"
#include "solver.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <vector>
#include <memory>

using namespace std;

int main()
{
    int rows(5), cols(5);

    /*int rp[6] = { 0,4,7,15,20,25 };
    int ci[25] = {0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4};
    int v[25] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
    int nnzs = 25;*/

    /*int rp[6] = {0,2,4,6,10,13};
    int ci[13] = {0,4,1,3,2,3,1,2,3,4,0,3,4};
    int v[13] = {3,5,3,7,4,5,7,5,1,2,5,2,9};
    int nnzs = 13;*/

    int rp[6] = { 0,2,4,7,11,13 };
    int ci[13] = { 0,2,1,3,0,2,3,1,2,3,4,3,4 };
    int v[13] = { 1,6,2,4,6,3,9,4,9,4,7,7,5 };
    int nnzs = 13;

    int *row_position = new int[rows+1];
    int *col_index = new int[nnzs];
    double *values = new double[nnzs];

    for (int i = 0; i < nnzs; i++)
    {
        col_index[i] = ci[i];
        values[i] = v[i];
    }

    cout << "-------col_index, values assignment done-------" << endl;

    for (int i = 0; i < rows+1; i++)
    {
        row_position[i] = rp[i];
    }

    cout << "-------row position assignment done-------" << endl;


    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, values, row_position, col_index);

    //CSRMatrix<double> sparse_mat(rows, cols, nnzs, values, row_position, col_index);
    
    for (int i = 0; i < sparse_mat->rows; i++)
    {
        //cout << "Loop over rows+1: i = " << i << endl;
        //cout << "---------------" << endl;
        // get the non-zeros number for i-1th row
        int nnzs_current = sparse_mat->row_position[i+1] - sparse_mat->row_position[i];
        
        // for each diagonal entry:
        // if possible, find the largest row index below
        // so should loop over col_index, which is also the loop
        // for CSRMatrix values
        double max_value(sparse_mat->values[sparse_mat->row_position[i]]);
        int max_value_index(i);
        int nnzs_max_row;
        int max_value_start_index;
        int current_start_index = sparse_mat->row_position[i];

        // find max value and its index in the value array
        for (int k = sparse_mat->row_position[i]; k < sparse_mat->nnzs; k++)
        {
            //cout << "loop over nnzs: k = " << k << endl;
            // find the MAX value index in values 
            // for i - 1th col, and check if the col index if the same as rows'
            if (sparse_mat->col_index[k] == i)
            {
                //cout << "current row: " << i << endl;
                //cout << "Find max value: " << sparse_mat->values[k] << " and its col: " << sparse_mat->col_index[k] << endl;
                if (sparse_mat->values[k] > max_value)
                {
                    max_value_index = k;
                    max_value = sparse_mat->values[k];
                }
            }
        }
        //cout << "Max value index: " << max_value_index << ", ";
        //cout << "Max value: " << max_value << endl;
        //cout << "Getting nnzs for max value row: " << endl;

        // as the matrix is diagonal positive
        // so dont need to worry about 0s on the diagnal after partial pivoting
        for (int j = 0; j < sparse_mat->rows; j++)
        {
            if ((max_value_index >= sparse_mat->row_position[j]) & (max_value_index <= sparse_mat->row_position[j+1]))
            {
                max_value_start_index = sparse_mat->row_position[j];
                //cout << sparse_mat->row_position[j] << " - " << sparse_mat->row_position[j+1] << endl;
                nnzs_max_row = sparse_mat->row_position[j+1] - sparse_mat->row_position[j];
            }
        }
        //cout << "How many non-zeros in the current row: " << nnzs_current << endl;
        //cout << "How many non-zeros in the max value row: " << nnzs_max_row << endl;
        //cout << endl;

        // create two pointers to store current row and max value row non-zeros and their col indices;
        vector<double> temp_current_value(nnzs_current);
        vector<int> temp_current_col(nnzs_current);
        vector<double> temp_max_value(nnzs_max_row); 
        vector<int> temp_max_col(nnzs_max_row);

        int nnzs_in_between = max_value_start_index - current_start_index - nnzs_current;
        int in_between_start_index = current_start_index + nnzs_current;
        //cout << "Current row start index: " << current_start_index << endl;
        //cout << "This is how many non-zeros in between two rows: " << nnzs_in_between << endl;
        //cout << "And the in between start index is: " << in_between_start_index << endl;

        //cout << "Below is the current row non-zeros: " << endl;
        // store current row non-zeros value/col
        for (int in = 0; in < nnzs_current; in++)
        {
            temp_current_value[in] = sparse_mat->values[current_start_index + in];
            temp_current_col[in] = sparse_mat->col_index[current_start_index + in];
            //cout << "Value: " << temp_current_value[in] << " Col: " << temp_current_col[in] << endl;
        }

        //cout << "Below is the max value row non-zeros: " << endl; 
        // store max value row non-zeros value/col
        for (int in = 0; in < nnzs_max_row; in++)
        {
            //cout << "Start index for max value row: " << max_value_start_index << endl;
            temp_max_value[in] = sparse_mat->values[max_value_start_index + in];
            temp_max_col[in] = sparse_mat->col_index[max_value_start_index + in];
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
                    temp_between_value[in] = sparse_mat->values[in_between_start_index + in];
                    temp_between_col[in] = sparse_mat->col_index[in_between_start_index + in];
                    //cout << "Value: " << temp_between_value[in] << " Col: " << temp_between_col[in] << endl;
                }

                // refill value/col: max value row -> current row
                // cout << endl << "Now the matrix value is changed: " << endl;
                for (int ivc = 0; ivc < nnzs_max_row; ivc++)
                {
                    sparse_mat->values[current_start_index + ivc] = temp_max_value[ivc];
                    sparse_mat->col_index[current_start_index + ivc] = temp_max_col[ivc];
                    //cout << "New max value/col position: " << current_start_index + ivc << endl;
                }

                // refill value/col: change pos of in between ones
                for (int ivc = 0; ivc < nnzs_in_between; ivc++)
                {
                    sparse_mat->values[current_start_index + nnzs_max_row + ivc] = temp_between_value[ivc];
                    sparse_mat->col_index[current_start_index + nnzs_max_row + ivc] = temp_between_col[ivc];
                    //cout << "New in between value/col position: " << current_start_index + nnzs_max_row + ivc << endl;
                }

                // refill value/col: current now behind in between ones
                for (int ivc = 0; ivc < nnzs_current; ivc++)
                {
                    sparse_mat->values[current_start_index + nnzs_max_row + nnzs_in_between + ivc] = temp_current_value[ivc];
                    sparse_mat->col_index[current_start_index + nnzs_max_row + nnzs_in_between + ivc] = temp_current_col[ivc];
                    //cout << "New current value/col position: " << current_start_index + nnzs_max_row + nnzs_in_between + ivc << endl;
                }

                // if there is a difference between nnzs of swapping two rows
                // if not, no need to change row position
                if (nnzs_current != nnzs_max_row)
                {
                    int nnzs_diff = nnzs_max_row - nnzs_current;
                    for (int irp = i + 1; irp < sparse_mat->rows; irp++)
                    {
                        if (sparse_mat->row_position[irp] < current_start_index + nnzs_max_row + nnzs_in_between)
                        {
                            sparse_mat->row_position[irp] += nnzs_diff;
                        }
                        else break;
                    }
                }
            }
            // if no elements in between
            else
            {
                // refill values/col: max value -> current 
                for (int ivc = 0; ivc < nnzs_max_row; ivc++)
                {
                    sparse_mat->values[current_start_index + ivc] = temp_max_value[ivc];
                    sparse_mat->col_index[current_start_index + ivc] = temp_max_col[ivc];
                    //cout << "New max value/col position: " << current_start_index + ivc << endl;
                }

                // refill values/col: current -> max value
                for (int ivc = 0; ivc < nnzs_current; ivc++)
                {
                    sparse_mat->values[current_start_index + nnzs_max_row + ivc] = temp_current_value[ivc];
                    sparse_mat->col_index[current_start_index + nnzs_max_row + ivc] = temp_current_col[ivc];
                    //cout << "New current value/col position: " << current_start_index + nnzs_max_row + ivc << endl;
                }
                // check if row position needs to be changed
                // like before, only need to change if the swapping rows have difference nnzs
                if (nnzs_current != nnzs_max_row)
                {
                    int nnzs_diff = nnzs_max_row - nnzs_current;
                    sparse_mat->row_position[i] += nnzs_diff;
                }
            }
        }
        //cout << "Swapped matrix: " << endl;
        //sparse_mat->printMatrix();
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
        for (int k = sparse_mat->row_position[i + 1]; k < sparse_mat->nnzs; k++)
        {
            //cout << "loop over nnzs: k = " << k << endl;
            // find the MAX value index in values 
            // for i - 1th col, and check if the col index if the same as rows'

            if (sparse_mat->col_index[k] == i)
            {
                sparse_mat->values[k] /= max_value;
                // from now on the lower part will become L except for the diagonal 1s
                // store the constant divisor
                s = sparse_mat->values[k];
                //cout << "s: " << s << " which is now at: " << k << "th position" << endl;

                for (int irp = 0; irp < sparse_mat->rows; irp++)
                {
                    //cout << "***" << sparse_mat->row_position[irp] << "***" << endl;
                    if ((sparse_mat->row_position[irp] <= k) && (sparse_mat->row_position[irp + 1] > k)) 
                    {
                        nnzs_found = sparse_mat->row_position[irp + 1] - sparse_mat->row_position[irp];
                        found_start_index = sparse_mat->row_position[irp];
                        found_next_start_index = sparse_mat->row_position[irp + 1];
                        //cout << "nnzs found: " << nnzs_found << " and found start index: " << found_start_index << endl;
                    }
                }

                found_iright_start_index = k + 1;
                nnzs_found_iright = found_next_start_index - found_iright_start_index;
                //cout << "Found iright start index: " << found_iright_start_index << " and nnzs found iright: " << nnzs_found_iright << endl;

                for (int ivc = current_start_index; ivc < current_start_index + nnzs_max_row; ivc++)
                {

                    if (sparse_mat->col_index[ivc] <= i)
                    {
                        current_iright_start_index++;
                        nnzs_current_iright--;
                    }
                }
                //cout << "Current iright start index: " << current_iright_start_index << " and nnzs current iright: " << nnzs_current_iright << endl;
            }
        }
        // if the current row has 0 non-zeros, no change would be made
        if (nnzs_current_iright > 0)
        {
            vector<double> temp_modified_value(nnzs_current_iright);
            vector<int> temp_modified_col(nnzs_current_iright);

            int num_add = nnzs_current_iright - nnzs_found_iright;
            //cout << "num_add: " << num_add << endl;
            if (num_add >= 0)
            {
                for (int t = 0; t < nnzs_current_iright; t++)
                {
                    // check when there is valid value at current row, if there is also a valid value
                    // at the row where non-zero was found at the begin of the row

                    if (sparse_mat->col_index[current_iright_start_index + t] == sparse_mat->col_index[found_iright_start_index])
                    {
                        temp_modified_value[t] = sparse_mat->values[found_iright_start_index] - sparse_mat->values[current_iright_start_index + t] * s;
                        temp_modified_col[t] = sparse_mat->col_index[current_iright_start_index + t];
                        //cout << "Pos If: " << "temp modified value: " << temp_modified_value[t] << " and col: " << temp_modified_col[t] << endl;
                    }
                    else
                    {
                        //cout << sparse_mat->values[current_iright_start_index + t] << endl;
                        temp_modified_value[t] = 0 - sparse_mat->values[current_iright_start_index + t] * s;
                        temp_modified_col[t] = sparse_mat->col_index[current_iright_start_index + t];
                        //cout << "Pos Else: " << "temp modified value: " << temp_modified_value[t] << " and col: " << temp_modified_col[t] << endl;
                    }
                }
            }
            else
            {
                for (int t = 0; t < nnzs_found_iright; t++)
                {
                    if (sparse_mat->col_index[current_iright_start_index + t] == sparse_mat->col_index[found_iright_start_index])
                    {
                        temp_modified_value[t] = sparse_mat->values[found_iright_start_index] - sparse_mat->values[current_iright_start_index + t] * s;
                        temp_modified_col[t] = sparse_mat->col_index[found_iright_start_index + t];
                        //cout << "Neg If: " << "temp modified value: " << temp_modified_value[t] << " and col: " << temp_modified_col[t] << endl;
                    }
                    else
                    {
                        //cout << sparse_mat->values[current_iright_start_index + t] << endl;
                        temp_modified_value[t] = sparse_mat->values[found_iright_start_index + t];
                        temp_modified_col[t] = sparse_mat->col_index[found_iright_start_index + t];
                        //cout << "Neg Else: " << "temp modified value: " << temp_modified_value[t] << " and col: " << temp_modified_col[t] << endl;
                    }
                }
            }
            //sparse_mat->printMatrix();
            /*for (int p = 0; p < nnzs_current_iright; p++)
            {
                cout << "temp modified value: " << temp_modified_value[p] << " and col: " << temp_modified_col[p] << endl;
            }*/

            // if num_add < 0: no nnzs change, only change in values for iright at found row
            if (num_add < 0)
            {
                for (int neg = 0; neg < nnzs_found_iright; neg++)
                {
                    sparse_mat->values[found_iright_start_index + neg] = temp_modified_value[neg];
                    sparse_mat->col_index[found_iright_start_index + neg] = temp_modified_col[neg];
                }
            }
            else
            {
                vector<double> new_value(sparse_mat->nnzs + num_add);
                vector<int> new_col(sparse_mat->nnzs + num_add);

                // assign unchanged values/col numbers
                for (int ni = 0; ni < found_iright_start_index; ni++)
                {
                    new_value[ni] = sparse_mat->values[ni];
                    new_col[ni] = sparse_mat->col_index[ni];
                }

                // insert temp modified values and col from the start of the found iright non-zeros
                for (int temp_ni = 0; temp_ni < nnzs_current_iright; temp_ni++)
                {
                    new_value[found_iright_start_index + temp_ni] = temp_modified_value[temp_ni];
                    new_col[found_iright_start_index + temp_ni] = temp_modified_col[temp_ni];
                }

                // push back the rest of non-zeros
                for (int ni = found_iright_start_index + nnzs_found_iright; ni < sparse_mat->nnzs; ni++)
                {
                    new_value[ni + num_add] = sparse_mat->values[ni];
                    new_col[ni + num_add] = sparse_mat->col_index[ni];
                }

                for (int irp = 0; irp < sparse_mat->rows; irp++)
                {
                    if ((found_iright_start_index >= sparse_mat->row_position[irp]) && (found_iright_start_index <= sparse_mat->row_position[irp + 1]))
                    {
                        for (int irp_ = irp + 1; irp_ < sparse_mat->rows + 1; irp_++)
                        {
                            sparse_mat->row_position[irp_] += num_add;
                        }
                    }
                }

                sparse_mat->nnzs += num_add;
                sparse_mat->values.reset(new double[sparse_mat->nnzs + num_add]);
                sparse_mat->col_index.reset(new int[sparse_mat->nnzs + num_add]);

                for (int pp = 0; pp < sparse_mat->nnzs; pp++)
                {
                    sparse_mat->values[pp] = new_value[pp];
                    sparse_mat->col_index[pp] = new_col[pp];
                    //cout << "New value: " << new_value[pp] << " and new col: " << new_col[pp] << endl;
                }
            }
            //sparse_mat->printMatrix();
            cout << endl;   
        }
        if (i == sparse_mat->rows - 2) break;
    }
    cout << "Result: " << endl;
    sparse_mat->printMatrix();



    delete[] row_position;
    delete[] col_index;
    delete[] values;
    delete sparse_mat;

    return 0;
}