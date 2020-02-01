#pragma once
#include "solver.h"
#include <memory>


template<class T>
void swap_rows(Matrix<T>& A, int& j, int& k)
{
	const int N = A.cols;
	for (int i = 0; i < A.cols; i++)
	{
		T* num = new T;
		*num = A.values[k * A.rows + i];
		A.values[k * A.rows + i] = A.values[j * A.rows + i];
		A.values[j * A.rows + i] = *num;
		delete num;
	}
}

// find the index of the largest value
// and return it as the index of row
template<class T>
int argmax(int& k, Matrix<T>& A)
{
	int index = k * A.cols + k;
	for (int i = k; i < A.cols; i++)
	{
		if (abs(A.values[i * A.cols + k]) > abs(A.values[index]))
		{
			index = i * A.cols + k;
		}
	}
	return index/A.cols;
}

template<class T>
void LU_decomposition_pp(Matrix<T>& A, Matrix<T>& L, Matrix<T>& P_)
{
	const int m = A.cols;

	// create an empty matrix L
	for (int i = 0; i < m * m; i++)
	{
		L.values[i] = 0;
	}

	// create eye matrix P_

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) { P_.values[i * m + j] = 1; }
			else P_.values[i * m + j] = 0;
		}
	}

	for (int k = 0; k < m - 1; k++)
	{
//		cout << "\nLoop over rows+1: i = " << k << endl;
//		cout << "---------------" << endl;
		int j = argmax(k, A);
		swap_rows(A, j, k);
		swap_rows(P_, j, k);
		swap_rows(L, j, k);
//		A.printMatrix();
//		L.printMatrix();
		for (int i = k + 1; i < m; i++)
		{
			const double s = A.values[i * m + k] / A.values[k * m + k];
			L.values[i * m + k] = s;
			for (int j = k; j < m; j++)
			{
				A.values[i * m + j] -= A.values[k * m + j] * s;
			}
		}
//		A.printMatrix();
//		L.printMatrix();
	}

	// add diagnal eye into L
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) { L.values[i * m + j] = 1; }
		}
	}

	// get the transpose of P_
	P_.transpose(P_);
}

template<class T>
void backward_substitution(Matrix<T>& U, T* b, T* output)
{
	const int m = U.rows;

	for (int k = m - 1; k >= 0; --k)
	{
		T s(0);
		for (int j = k + 1; j < m; j++)
		{
//			cout << s << " ";
			s += U.values[k * m + j] * output[j];
		}
		output[k] = (b[k] - s) / U.values[k * m + k];
	}
}

template<class T>
void forward_substitution(Matrix<T>& L, T* b, T* output)
{
	const int m = L.rows;
//	cout << "Printing pinvb: " << endl;
//	for (int i = 0; i < m; i++)
//	{
//		cout << b[i] << " ";
//	}
//	cout << endl;
	for (int k = 0; k < m; k++)
	{
		T s(0);
		for (int j = 0; j < k; j++)
		{
//			cout << s << " ";
			s += L.values[k * m + j] * output[j];
		}
//		cout << endl;
		output[k] = (b[k] - s) / L.values[k * m + k];
	}
//	cout << "Printing forward sub vector output: " << endl;
//	for (int i = 0; i < m; i++)
//	{
//		cout << " " << output[i];
//	}
//	cout << endl;
}

template<class T>
void LU_dense(Matrix<T>& A, T* b, T* x)
{
	if (A.rows != A.cols)
	{
		cerr << "Cannot decomposite non-square matrix into LU. \n";
		return;
	}

	const int m = A.cols;
	auto* L = new Matrix<T>(m, m, true);
	auto* P = new Matrix<T>(m, m, true);
	LU_decomposition_pp(A, *L, *P);
	// print out P, L, U
	//L->printMatrix();
	
	// A is now the upper triangle
	//A.printMatrix();
	// get the inverse i.e. the transpose of P
	P->transpose(*P);
	//P->printMatrix();

	auto* pinvb = new double[m * 1];
	P->matVecMult(b, pinvb);

	auto* y = new double[m * 1];
    // using cblas to calculate triangular system
    //cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, m, 1, 1, L, m, pinvb, m);
	forward_substitution(*L, pinvb, y);
	backward_substitution(A, y, x);

	cout << "LU Solution Dense: \n";
	for (int i = 0; i < m; i++)
	{
		cout << "x" << i << ": " << x[i] << endl;
	}

	delete L;
	delete P;
	delete[] pinvb;
	delete[] y;
	/*delete[] solution;*/
}

template<class T>
void gauss_elimination(Matrix<T>& A, T* x, T* b)
{
    if (A.rows != A.cols)
    {
        cerr << "Cannot apply Gaussian Elimination on non-square matrix. \n";
        return;
    }

    const int m = A.cols;
    int i, j, k;

    for (i = 0; i < m; i++) {                   //Pivotisation
        for (k = i + 1; k < m; k++) {
            if (abs(A.values[i * m + i]) < abs(A.values[k * m + i])) {
                for (j = 0; j < m; j++) {
                    auto* temp = new double;
                    *temp = A.values[i * m + j];
                    A.values[i * m + j] = A.values[k * m + j];
                    A.values[k * m + j] = *temp;
                    delete temp;
                }
                auto* temp2 = new double;
                *temp2 = b[i];
                *(b + i) = b[k];
                *(b + k) = *temp2;
                delete temp2;
            }
        }
    }

    // Perform Gauss-Elimination
    for (i = 0; i < m - 1; i++) {
        for (k = i + 1; k < m; k++) {
            auto* t = new double;
            *t = A.values[k * m + i] / A.values[i * m + i];
            for (j = 0; j < m; j++) {

                //make the elements below the pivot elements equal to zero or eliminate the variables
                A.values[k * m + j] = A.values[k * m + j] - *t * A.values[i * m + j];
            }
            *(b + k) = b[k] - *t * b[i];
            delete t;
        }
    }

    //back-substitution
    for (i = m - 1; i >= 0; i--) {

        // make the variable to be calculated equal to the rhs of the last equation
        *(x + i) = b[i];
        for (j = i + 1; j < m; j++) {

            // subtract all the lhs values except the coefficient of the variable whose value is being calculated
            if (j != i) {
                *(x + i) = x[i] - A.values[i * m + j] * x[j];
            }
        }
        //now finally divide the rhs by the coefficient of the variable to be calculated
        *(x + i) = x[i] / A.values[i * m + i];
    }
}

template<class T>
void gauss_seidel_dense(Matrix<T>& A, T* x, T* b, T er, T urf) {
    // gauss_idel(Matrix<T> & A, T * x, T * b, er, urf, nmax, rmx)

    // solves systems of linear eqs., with the Gauss-Seidel method

    /*
    // input
    n = number of eqs.
    a  = matrix of coefficients
    b  = matrix of constants
    x  = initial values for the unknown (nq values)
    er = termination criterion
    urf = relaxation factor
    nmax = maximum number of iterations
    // output
    x = solution
    rmx = maximum residual (verification)
    */

    if (A.rows != A.cols)
    {
        cerr << "Cannot apply Gauss-Seidel on non-square matrix. \n";
        return;
    }

    const int n = A.cols;
    const int m = A.cols;
    int i, j, k, niter;
    double sum, xold;
    int nmax = 1000;
    double rmx = 0.0;

    // initialisation of x
    for (i = 0; i < n; i++) {
        *(x + i) = 0.0;
    }

    //Pivotisation
    for (i = 0; i < m; i++) {
        for (k = i + 1; k < m; k++) {
            if (abs(A.values[i * m + i]) < abs(A.values[k * m + i])) {
                for (j = 0; j < m; j++) {
                    auto* temp = new double;
                    *temp = A.values[i * m + j];
                    A.values[i * m + j] = A.values[k * m + j];
                    A.values[k * m + j] = *temp;
                    delete temp;
                }
                auto* temp2 = new double;
                *temp2 = b[i];
                *(b + i) = b[k];
                *(b + k) = *temp2;
                delete temp2;
            }
        }
        cout << endl;
    }

    // start of iterations
    for (niter = 0; niter < nmax; niter++) {

        double ea = 0.0;

        // cout << "iter " << niter << endl;  // print iterations until convergence

        // new x's calculation
        for (i = 0; i < n; i++) {
            xold = x[i];
            sum = b[i];

            // cout << "i " << i << endl;

            for (j = 0; j <= i - 1; j++) {
                sum += -A.values[i * n + j] * x[j];
                // cout << "iter j1: " << j << ", A: " << A.values[i * n + j] << ", x: " << x[j] << ", sum: " << sum << endl;
            }

            for (j = i + 1; j < n; j++) {
                sum += -A.values[i * n + j] * x[j];
                // cout << "iter j2: " << j << ", A: " << A.values[i * n + j] << ", x: " << x[j] << ", sum: " << sum << endl;
            }

            *(x + i) = sum / A.values[i * n + i];

            // error and underelaxation
            double ern = abs(*(x + i) - xold);

            // cout << i << " " <<  ern << endl;
            // cout << endl;
            // ea = __max(ea, ern);
            ea = max(ea, ern);
            *(x + i) = *(x + i) * urf + xold * (1. - urf);
        }

        // Checks for exit
        if (ea < er) {
            break;
        }

    }

    if (niter >= nmax) {
        cout << "Warning! Iterations' limit";
    }

    //// Verification
    //for (i = 0; i < n; i++) {
    //    sum = - b[i];
    //    for (j = 0; j < n; j++) {
    //        sum += A.values[i * n + j] * x[j];
    //    }
    //    rmx = __max(abs(sum), rmx);
    //}
}

template<class T>
void jacobi_dense(Matrix<T>& A, T* x, T* b, int maxit, double tolerance)
//user can input Matrix A and array b and can also define iteration tolerance and number of iteration times
{
    //set the condition that matrix A must be a square matrix
    if (A.rows != A.cols)
    {
        cerr << "Input matrix A must be a sqaure matrix!" << endl;
    }
    //have a guess of the size of output matrix, which has the same size of the rows size of matrix A
    unique_ptr<T[]> x_new_array(new T[A.rows]);
    
    //initialize the x matrix and new x matrix as zero matrixs
    for (int i = 0; i < A.rows; i++)
    {
        x[i] = 0;
        x_new_array[i] = 0;
        //cerr << x_new_array[i]<< endl;
    }

    unique_ptr<T[]> total_sum(new T[A.rows]);
    for (int k = 0; k <maxit; k++) //record the number of iteration
    {
        for (int i = 0; i < A.rows; i++)
        {
            total_sum[i] = 0;
        }
        for (int i = 0; i < A.rows; i++)
        {
            //Define a multiple production by A[i,:i] and x[:i]
            vector<T> mul_Ax(A.rows); // Since that A[i,:i] and x[:i] are both 1d arrays,
            //so the multiple production should be dot multiple production, which returns one value
            //so we create a vector to store values
            //auto* mul_Ax = new double[i];
            if (i == 0) //when i =0; A[i,:i] and x[:i] are 0
            {
                mul_Ax[i] = 0;
            }
            else
            {
                for (int j = 0; j < i; j++)
                {
                    mul_Ax[i] += A.values[i * A.cols + j] * x[j];
                }
            }
            //Define a multiple production by A[i,i+1:] and x[i+1:]
            // Since that A[i,i+1:] and x[i+1:] are both 1d arrays,
            //so the multiple production should be dot multiple production, which returns one value
            //so we create a vector to store values
            vector<T> mul2_Ax(A.cols);
            if (i == A.rows - 1)//when i is the last index; A[i,i+1:] and x[i+1:] are 0
            {
                mul2_Ax[i] = 0;
            }
            else
            {
                for (int j = i + 1; j < A.cols; j++)
                {
                    mul2_Ax[i] += A.values[i * A.cols + j] * x[j];
                }
            }
            x_new_array[i] = (1. / A.values[i * A.cols + i]) * (b[i] - mul_Ax[i] - mul2_Ax[i]);
        }
        double pow_sum = 0;
        for (int i = 0; i < A.rows; i++)
        {
            for (int j = 0; j < A.cols; j++)
            {
                total_sum[i] += A.values[i * A.cols + j] * x_new_array[j];

            }
            pow_sum += pow(total_sum[i] - b[i], 2);
        }
        //calculate the error
        double residual = sqrt(pow_sum);
//        cout << "residual" << k << " = " << residual << endl;
        if (residual < tolerance)
        {
            break;
        }
        for (int i = 0; i < A.rows; i++)
        {
            x[i] = x_new_array[i];
            cout << x[i] << endl;
        }
    }
}

template<class T>
void jacobi_dense_blas(Matrix<T>& A, T* x, T* b, int maxit, double tolerance)
//user can input Matrix A and array b and can also define iteration tolerance and number of iteration times
{
    //set the condition that matrix A must be a square matrix
    if (A.rows != A.cols)
    {
        cerr << "Input matrix A must be a sqaure matrix!" << endl;
    }
    //have a guess of the size of output matrix, which has the same size of the rows size of matrix A
    unique_ptr<T[]> x_new_array(new T[A.rows]);
    
    //initialize the x matrix and new x matrix as zero matrixs
    for (int i = 0; i < A.rows; i++)
    {
        x[i] = 0;
        x_new_array[i] = 0;
        //cerr << x_new_array[i]<< endl;
    }

    unique_ptr<T[]> total_sum(new T[A.rows]);
    for (int k = 0; k <maxit; k++) //record the number of iteration
    {
        for (int i = 0; i < A.rows; i++)
        {
            total_sum[i] = 0;
        }
        for (int i = 0; i < A.rows; i++)
        {
            //Define a multiple production by A[i,:i] and x[:i]
            vector<T> mul_Ax(A.rows);
            // Since that A[i,:i] and x[:i] are both 1d arrays,
            //so the multiple production should be dot multiple production, which returns one value
            //so we create a vector to store values
            //auto* mul_Ax = new double[i];
            if (i == 0) //when i =0; A[i,:i] and x[:i] are 0
            {
                mul_Ax[i] = 0;
            }
            else
            {
                for (int j = 0; j < i; j++)
                {
                    mul_Ax[i] += A.values[i * A.cols + j] * x[j];
                }
            }
            //Define a multiple production by A[i,i+1:] and x[i+1:]
            // Since that A[i,i+1:] and x[i+1:] are both 1d arrays,
            //so the multiple production should be dot multiple production, which returns one value
            //so we create a vector to store values
            vector<T> mul2_Ax(A.cols);
            if (i == A.rows - 1)//when i is the last index; A[i,i+1:] and x[i+1:] are 0
            {
                mul2_Ax[i] = 0;
            }
            else
            {
                for (int j = i + 1; j < A.cols; j++)
                {
                    mul2_Ax[i] += A.values[i * A.cols + j] * x[j];
                }
            }
            x_new_array[i] = (1. / A.values[i * A.cols + i]) * (b[i] - mul_Ax[i] - mul2_Ax[i]);
        }
        double pow_sum = 0;
        for (int i = 0; i < A.rows; i++)
        {
            for (int j = 0; j < A.cols; j++)
            {
                total_sum[i] += A.values[i * A.cols + j] * x_new_array[j];

            }
            pow_sum += pow(total_sum[i] - b[i], 2);
        }
        //calculate the error
        double residual = sqrt(pow_sum);
//        cout << "residual" << k << " = " << residual << endl;
        if (residual < tolerance)
        {
            break;
        }
        for (int i = 0; i < A.rows; i++)
        {
            x[i] = x_new_array[i];
            cout << x[i] << endl;
        }
    }
}


template<class T>
void gauss_seidel_sparse(CSRMatrix<T>& A, T* x, T* b, int maxit, double tolerance)
{
    for (int i = 0; i < A.rows; i++) //initialize x array
    {
        x[i] = 0;
    }

    unique_ptr<int[]> row_diff(new int[A.cols]);
    vector<int> diag_values;
    //get the numbers of nonzero values for each row
    for (int i = 0; i < A.cols; i++)
    {
        row_diff[i] = A.row_position[i + 1] - A.row_position[i];
    }

    int m = A.rows;
    int n = A.cols;
    unique_ptr<T[]> sum(new T[m]);
    unique_ptr<T[]> total_sum(new double[A.rows]);
    for (int k = 0; k < maxit; k++) //the number of iterations
    {
        for (int i = 0; i < m; i++)
        {
            sum[i] = 0;
            total_sum[i] = 0;
        }
        int counter = 0;//use counter to track the index of each value
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < row_diff[i]; j++)
            {
                if (A.col_index[counter] == i) //find diagonal values(A[i,i])
                {
                    diag_values.push_back(A.values[counter]);
                }
                else
                {
                    //cout << "counter=" << counter;
                    sum[i] += A.values[counter] * x[A.col_index[counter]];//the result of ( A[i, :i] @ x[:i] ) + ( A[i, i+1:] @ x[i+1:] )
                    //which is actually a sum of dot products
                }
                counter++;
            }
            x[i] = 1. / diag_values[i] * (b[i] - sum[i]); //Gauss_seidle method is very similar with Jacobi method
            // The difference is that Gauss-Seidle method uses the latest updated values in the iteration while the
            //Jacobi method stores the values in one new array and uses the value obtained from the last step
            
            //cout << x[i] << endl;
        }
        double pow_sum = 0;
        for (int i = 0; i < A.rows; i++)
        {
            for (counter = A.row_position[i]; counter < A.row_position[i] + row_diff[i]; counter++)
            {
                total_sum[i] += A.values[counter] * x[A.col_index[counter]];
            }
            pow_sum += pow(total_sum[i] - b[i], 2);
        }
        //calculate the error
        double residual = sqrt(pow_sum);
        if (residual < tolerance)
        {
            break;
        }
    }
    
    cout << "Gauss-Seidel Solution Sparse: " << endl;
    for (int i = 0; i < m; i++)
    {
        cout << "x" << i << ": " << x[i] << endl;
        
    }
}

template<class T>
void jacobi_sparse(CSRMatrix<T>& A, T* x, T* b, int maxit, double tolerance)
{
    int m = A.rows;
    unique_ptr<T[]> x_new_array(new T[A.cols]);
    for (int i = 0; i < m; i++)
    {
        x[i] = 0;
        x_new_array[i] = 0;
    }
    unique_ptr<int[]> row_diff(new int[A.cols]);
    vector<int> diag_values;
    //get the numbers of nonzero values for each row
    for (int i = 0; i < A.cols; i++)
    {
        row_diff[i] = A.row_position[i + 1] - A.row_position[i];
    }

    shared_ptr<T[]> sum(new T[m]);
    shared_ptr<T[]> total_sum(new double[A.rows]);
    cout << "\nJacobi Solution Sparse: " << endl;
    for (int k = 0; k < 7; k++) //the number of iterations
    {
        for (int i = 0; i < m; i++)
        {
            sum[i] = 0;
            total_sum[i] = 0;
        }
        int counter = 0;//use counter to track the index of each value
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < row_diff[i]; j++)
            {
                if (A.col_index[counter] == i) //find diagonal values(A[i,i])
                {
                    diag_values.push_back(A.values[counter]);
                }
                else
                {
                    //cout << "counter=" << counter;
                    sum[i] += A.values[counter] * x[A.col_index[counter]];
                    //cout << "sum=" << sum[i] << endl;
                }
                counter++;
            }
            x_new_array[i] = (1. / diag_values[i]) * (b[i] - sum[i]); //Jacobi stores the array in one new array
        }
        double pow_sum = 0;
        for (int i = 0; i < A.rows; i++)
        {
            for (counter = A.row_position[i]; counter < A.row_position[i] + row_diff[i]; counter++)
            {
                total_sum[i] += A.values[counter] * x_new_array[A.col_index[counter]];
            }
//            cout << "total_sum_sparse="<<total_sum[i] << endl;
            pow_sum += pow(total_sum[i] - b[i], 2);
        }
        //calculate the error
        double residual = sqrt(pow_sum);
        if (residual < tolerance)
        {
            break;
        }
        
        for (int i = 0; i < m; i++)
        {
            x[i] = x_new_array[i];//copy the new array to the original one after each iteration, which is different to Gauss_seidle method
        }
    
    }
    for (int i = 0; i < m; i++)
    {
        cout << "x" << i << ": " << x[i] << endl;

    }
}

template<class T>
void thomas(Matrix<T>& A, T* x, T* r)
{
    int i;
    double term;
    const int n = A.cols;
    A.values[2 * n] = A.values[2 * n] / A.values[1 * n];
    r[0] = r[0] / A.values[1 * n];

    //-- forward elimination
    for (i = 1; i < n; i++) {
        term = A.values[1 * n + i] - A.values[0 * n + i] * A.values[2 * n + i - 1];
        A.values[2 * n + i] = A.values[2 * n + i] / term;
        r[i] = (r[i] - A.values[0 * n + i] * r[i - 1]) / term;
    }
    //-- back substitution
    x[n - 1] = r[n - 1];
    for (i = n - 2; i >= 0; i--) {
            x[i] = r[i] - A.values[2 * n + i] * x[i + 1];
        }
}
