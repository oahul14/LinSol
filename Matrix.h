#pragma once    
#include <vector>
#include <tuple>
#include <memory>

using namespace std;

template<class T>
class Matrix
{
public:
	/////// public methods
	
	// constructor where we want to preallocate memory: own our own memory
	Matrix(int rows, int cols, bool preallocate);
	
	// constructor where we want to preallocate memoryalready preallocated memory outside
	// dont own our own memory
	Matrix(int rows, int cols, T* values_ptr);

	// copy constructor
	Matrix(Matrix &B);

	// destructor
	// virtual ~Matrix() = 0; pure virtual function: the sub class MUST overwrite this func
	virtual ~Matrix(); //the sub class CAN overwrite this func

	// print matrix values
	void printValues();
	virtual void printMatrix();

	// some basic functions
	void transpose(Matrix<T>& itself);
    // generate random n x n matrix
    void genRanDense(bool dom);
    void genRanSparse(double sparsity, bool dom);
	void matMatMult(Matrix& mat_right, Matrix& output);
	void matVecMult(T* vec, T* output);
	void matMatMult_colMajor(Matrix& mat_right, Matrix& output);
	
	/////// public variables
	// matrix size;
	// explicitly using the c++11 nullptr;
	unique_ptr<T[]> values;
	int rows = -1;
	int cols = -1;

	// private variables: no need for other classes to know
	// values
protected:
	bool preallocated = false;
private:
	int size_of_values = -1;
	

};

