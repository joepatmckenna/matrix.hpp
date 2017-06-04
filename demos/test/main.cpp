#include <iostream> // cout
#include <typeinfo> // typeid
#include "matrix.hpp"
using namespace std;

typedef double real;

int main ()
{
    // TODO:
    // slicing
    // non-template declaration: matrix A(3,3);
    // ones, identity, pd, spd, orthogonal
    
    { // default constructor
        matrix<real> A;
        // matrix B;
    }
    { // parameterized constructor
        matrix<real> A(3,3);
    }
    { // copy constructor
        matrix<real> A;
        matrix<real> B(A), C=A;
        matrix<float> D;
        matrix<double> E(D), F=D;
    }
    { // copy assignment
        matrix<real> A(3,3), B(3,3);
        A=B;
        matrix<float> D;
        matrix<double> E;
        D=E;
    }
    { // move constructor
        matrix<real> A(matrix<real>(3,3)), B=matrix<real>(3,3);
        matrix<float> C(matrix<double>(3,3)), D=matrix<double>(3,3);
    }
    { // move assignment
        matrix<real> A=matrix<real>(3,3);
        matrix<float> B=matrix<double>(3,3);
    }
    { // access
        matrix<real> A(3,3), B(3,1);
        A(0,0); B(0); 
    }
    { // overloaded assignment
        matrix<real> A(3,3);
        A=0; 
        A(0,0)=A(1,1)=A(2,2)=1.;
    }
    { // matrix-matrix operations
        matrix<float> A(3,3), B(3,3);
        A+B; A*B; A-B; 
        A = 0.; A(0,0) = A(1,1) = A(2,2) = 2.;
        A=A^4;
        print(A);
        
    }
    { // elementwise operations
        matrix<double> A(3,3);
        matrix<float> B(3,3);
        B.times(A); A.times(B); 1.*A; A*1.; A.power(2);
        B*1.; 1.*B; 
        A='a'+true;
        B=(B=1)*true;
    }
    { // print
//         matrix<real> A(3,3);
//         print(A=0);
    }
    { // assignment type conversion
        matrix<float> A(3,3); matrix<double> B(3,3);
        float a=0; double b=0;
        A=B; B=A;
        A=b; B=a;
    }
    { // matrix-matrix operation type conversion
        matrix<float> A(3,3); matrix<double> B(3,3);
        A+B; A*B; A-B;
    }
    {
        matrix<double> A;
        matrix<float> B;
        static_cast<matrix<double>> (B);
//         static_cast<matrix<double>> (b);
        matrix<matrix<bool>> Q(1,10);
    }
    {
        matrix<real> A(3,3);
        A = 0.;
        A(0,0)=1.; A(0,1)=2.; A(0,2)=3.;
        print(A[0]);
        A[0]=1.;
        print(A[0]);
        
    }
    return 0;
}
