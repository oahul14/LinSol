#include <iostream>
#include <Accelerate/Accelerate.h>

using namespace std;
 
int main() {
 
    double A[9] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };
 
    double B[9] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };
 
    double C[9];
     
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, A, 3, B, 3, 0.0, C, 3);
    cout << "---- cblas_dgemm ----" << endl;
    cout << C[0] << " " << C[1] << " " << C[2] << endl;
    cout << C[3] << " " << C[4] << " " << C[5] << endl;
    cout << C[6] << " " << C[7] << " " << C[8] << endl;
 
    return 0;
}
