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
#include "testing.cpp"
#include "testing.h"

using namespace std;

int main()
{
    double* mat_array = new double[15*15]{98,   0,   0,   0,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,  91,   0,  13,   0,   0,   0,   0,   0,   5,   6,  6,   0,   0,   0,
    0,   0, 108,   0,   0,   0,   0,   0,   0,   0,   5,   0,  13,   0,   0,
    0,  13,   0, 145,   0,   0,   0,   0,   0,   0,   0,  10,  11,   0,   6,
    5,   0,   0,   0, 118,   0,   8,  11,   0,   9,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0, 133,   0,   0,  11,  12,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   8,   0, 142,  13,   0,   0,  11,   0,   6,   0,   0,
    0,   0,   0,   0,  11,   0,  13,  82,   0,   0,   0,   5,   0,   0,  10,
    0,   0,   0,   0,   0,  11,   0,   0,  93,   0,   0,   7,   0,   5,   0,
    0,   5,   0,   0,   9,  12,   0,   0,   0,  72,   0,   9,   0,   0,   0,
    0,   6,   5,   0,   0,   0,  11,   0,   0,   0, 114,   0,   0,   0,   0,
    0,   6,   0,  10,   0,   0,   0,   5,   7,   9,   0,  52,  11,   6,   0,
    0,   0,  13,  11,   0,   0,   6,   0,   0,   0,   0,  11, 106,   9,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   5,   0,   0,   6,   9, 116,   7,
        0,   0,   0,   6,   0,   0,   0,  10,   0,   0,   0,   0,   0,   7,  58};
    double* b = new double[15]{5, 2, 6, 5, 10, 8, 3, 8, 6, 4, 1, 1, 3, 10, 1};
    //time_all_given(15, mat_array, b, true);
    time_all_rand(10, "dense", true);
    delete[] mat_array;
    delete[] b;
    return 0;
}
