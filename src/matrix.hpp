#ifndef MATRIX_H
#define MATRIX_H
#include<math.h> // pow
#include<iostream>
#include<typeinfo> // typeid
using namespace std;

typedef double real;

//
// declaration
//
# define inf (1./0.)

// matrix slice class
class slice {
public:
    // 0,1: row,col; 0: lower bound
    // 0,1: row,col; 1: upper bound
    // 0,1: row,col; 2: stride
    int params[2][3];
    slice() {
        params[0][0] = params[1][0] = -inf;
        params[0][1] = params[1][1] = inf;
        params[0][2] = params[1][2] = 1;
    };
    slice(int i00, int i01, int i02, int i10, int i11, int i12) {
        params[0][0] = i00; params[1][0] = i10;
        params[0][1] = i01; params[1][1] = i11;
        params[0][2] = i02; params[1][2] = i12;
    }
};

// matrix class
template<typename T=float>
class matrix {
public:
    T **data;
    int shape[2];
    // default constructor
    matrix();
    // constructor(s) with arguments
    matrix(int,int);
    matrix(T*);
    // destructor
    ~matrix();
    // copy constructor
    matrix(const matrix&);
    // copy constructor (type conversion)
    template<typename U> matrix<T>(const matrix<U>&);
    // copy assignment
    matrix& operator=(const matrix&);
    // copy assignment (type conversion)
    template<typename U> matrix<T>& operator=(const matrix<U>&);
//     // type cast operator
//     template <typename U> operator matrix<U>() {return matrix<U>();} // matrix<U> = matrix<T>
//     // move constructor
//     matrix(matrix&&)=default;
//     // move assignment
//     matrix& operator=(matrix&&);
//     T& operator[](int[2]);
    // access A[i]
    matrix<T>& operator[](int);
    // access: A(0,0)
    T& operator()(int,int);
    // access: A(0) if A is vector
    T& operator()(int);
    // slicing
    matrix& operator[](slice);
    // leading plus sign
    matrix& operator+();
    // negation
    matrix& operator-();
    // scalar assignment: A=c
    matrix& operator=(T);
    // elementwise matrix-matrix product: A.times(B)
    matrix& times(const matrix&);
    // elementwise matrix power: A.power(2)
    matrix& power(const double);
//     real trace(const matrix&);
//     real det(const matrix&);
//     matrix& inv(const matrix&);
};

template <typename T>
matrix<T>::matrix(T* A) {
    shape[0] = 1;
    shape[1] = sizeof(A)/sizeof(*A);
    data = new T*[1];
    data[0] = A;
}

// // matrix scalar multiplication: c*A or A*c
template <typename T> matrix<T>& operator*(const T,matrix<T>&);
template <typename T> matrix<T>& operator*(matrix<T>&,const T);
template <typename T, typename U> auto operator*(matrix<T>& A,const U c) -> matrix<decltype(**A.data*c)>;
template <typename T, typename U> auto operator*(const U c, matrix<T>& A) -> matrix<decltype(**A.data*c)>;
// matrix-matrix addition: A+B
template <typename T> matrix<T>& operator+(const matrix<T>&,const matrix<T>&);
template <typename T, typename U> auto operator+(const matrix<T>& A,const matrix<U>& B) -> matrix<decltype(**A.data+**B.data)>;
// matrix-matrix subtraction: A-B
template <typename T> matrix<T>& operator-(const matrix<T>&,const matrix<T>&);
template <typename T, typename U> auto operator-(const matrix<T>& A,const matrix<U>& B) -> matrix<decltype(**A.data-**B.data)>;
// matrix-matrix multiplication: A*B
template <typename T> matrix<T>& operator*(const matrix<T>&,const matrix<T>&);
template <typename T, typename U> auto operator*(const matrix<T>& A,const matrix<U>& B) -> matrix<decltype(**A.data*(**B.data))>;
// matrix power: A^n
template <typename T> matrix<T>& operator^(const matrix<T>& A, int n);
// zero matrix: zeros(A)
template <typename T> matrix<T>& zeros(int m, int n);
// print matrix: print(A)
template <typename T> void print(matrix<T>& A);

//
// implementation
//

// default constructor
template <typename T> matrix<T>::matrix() {
    shape[0]=1; shape[1]=1;
    data = new T*[1];
    data[0]=new T[1];
}

// parameterized constructor(s)
template <typename T> matrix<T>::matrix(int m, int n) {
    shape[0] = m; shape[1] = n;
    data = new T*[m];
    for(int row=0; row<m; row++) {
        data[row] = new T[n];
    }
}

// destructor
template <typename T> matrix<T>::~matrix() {
    for(int row=0; row<shape[0]; row++) {
        delete[] data[row];
    }
    delete[] data;
}

// copy constructor
template <typename T> matrix<T>::matrix(const matrix& A) {
    shape[0] = A.shape[0]; shape[1] = A.shape[1];
    data = new T*[shape[0]];
    for(int row=0; row<shape[0]; row++) {
        data[row] = new T[shape[1]];
        for(int col=0; col<shape[1]; col++) {
            data[row][col] = A.data[row][col];
        } 
    } 
}

// copy constructor (type conversion)
template <typename T> template <typename U>
matrix<T>::matrix(const matrix<U>& A) {
    shape[0] = A.shape[0]; shape[1] = A.shape[1];
    data = new T*[shape[0]];
    for(int row=0; row<shape[0]; row++) {
        data[row] = new T[shape[1]];
        for(int col=0; col<shape[1]; col++) {
            data[row][col] = A.data[row][col];
        } 
    } 
}


// copy assignment
template <typename T>
matrix<T>& matrix<T>::operator=(const matrix<T>& A) {
    if (shape[0] == A.shape[0] && shape[1] == A.shape[1]) {
        for(int row=0; row<shape[0]; row++) {
            for(int col=0; col<shape[1]; col++) {
                data[row][col] = A.data[row][col];
            }
        } 
    } else {
        cout << "matrix::operator=(const matrix&) :";
        cout << " can't assign " << A.shape[0] << "x" << A.shape[1];
        cout << " to " << shape[0] << "x" << shape[1] << " matrix." << endl;
    } 
    return *this;
}
// copy assignment (type conversion)
template<typename T> template<typename U>
matrix<T>& matrix<T>::operator=(const matrix<U>& A) {
    if (shape[0] == A.shape[0] && shape[1] == A.shape[1]) {
        for(int row=0; row<shape[0]; row++) {
            for(int col=0; col<shape[1]; col++) {
                data[row][col] = A.data[row][col];
            }
        } 
    } else {
        cout << "matrix::operator=(const matrix&) :";
        cout << " can't assign " << A.shape[0] << "x" << A.shape[1];
        cout << " to " << shape[0] << "x" << shape[1] << " matrix." << endl;
    } 
    return *this;
}

// // move constructor
// template <typename T> matrix<T>::operator=(matrix<T>&& A) {}
// // move assignment
// template <typename T> matrix<T>& matrix<T>::operator=(matrix<T>&& A) {};

// negation
template <typename T>
matrix<T>& matrix<T>::operator-() {
    for(int row=0; row<shape[0]; row++) {
        for(int col=0; col<shape[1]; col++) {
            data[row][col] = -data[row][col];
        }
    } 
    return *this;
}

// leading plus sign
template <typename T>
matrix<T>& matrix<T>::operator+() {
    return *this;
}

// access
template <typename T> T&
matrix<T>::operator()(int i,int j) {
    T *datum = nullptr;
    if (i<shape[0] && j<shape[1]) {
        datum = &data[i][j];
    } else {
        cout << "matrix::operator()(int,int) :";
        cout << " element (" << i << "," << j << ") is out of bounds" << endl;
    }
    return *datum;
}
template <typename T>
T& matrix<T>::operator()(int i) {
    T *datum = nullptr;
    if (i<shape[0] && shape[1]==1) {
        datum = &data[i][0];
    } else if (i<shape[1] && shape[0]==1) {
        datum  = &data[0][i];
    } else {
        cout << "matrix::operator()(int) :";
        cout << " element " << i << " is ambiguous or out of bounds" << endl;
    }
    return *datum;
}

// access: A[i]
// need to find a better way to do this instead of construct new matrix every time i access row
template <typename T>
matrix<T>& matrix<T>::operator[](int i) {
    matrix<T> *R = new matrix<T>(1,shape[1]);
    R->data[0] = data[i];
    return *R;
}

// slicing
template <typename T>
matrix<T>& matrix<T>::operator[](slice s) {
    if (s.params[0][1]==inf) {
        s.params[0][1]=shape[0];
    }
    if (s.params[1][1]==inf) {
        s.params[1][1]=shape[1];
    }
    int nrow = (s.params[0][1]-s.params[0][0]+1)/s.params[0][2];
    int ncol = (s.params[1][1]-s.params[1][0]+1)/s.params[1][2];
    matrix<T> *A = new matrix<T>(nrow,ncol);
    int i=0, j=0;
    for (int row=s.params[0][0]; row<s.params[0][1]; row+=s.params[0][2]) {
        for (int col=s.params[1][0]; col<s.params[1][1]; col+=s.params[1][2]) {
            A->data[i][j] = &data[row][col];
        }
    }
    return *A;
}

// scalar assignment
template <typename T> matrix<T>& matrix<T>::operator=(T c) {
    for (int row=0; row<shape[0]; row++) {
        for (int col=0; col<shape[1]; col++) {
            data[row][col] = c;
        }
    }
    return *this;
}

// matrix-matrix addition: A+B
template <typename T> matrix<T>& operator+(const matrix<T>& A,const matrix<T>& B) {
    matrix<T> *C;
    if (A.shape[0]==B.shape[0] && A.shape[1]==B.shape[1]) {
        C = new matrix<T>(A.shape[0],A.shape[1]);
        for (int row=0; row<A.shape[0]; row++) {
            for (int col=0; col<A.shape[1]; col++) {
                C->data[row][col] = A.data[row][col]+B.data[row][col];
            }
        }
    } else {
        C = nullptr;
        cout << "matrix::operator+(const matrix&,const matrix&) :";
        cout << " can't add " << A.shape[0] << "x" << A.shape[1];
        cout << " and " << B.shape[0] << "x" << B.shape[1] << " matrices." << endl;
    }
    return *C;
}
template <typename T, typename U> auto operator+(const matrix<T>& A,const matrix<U>& B) -> matrix<decltype(**A.data+**B.data)> {
    typedef decltype(**A.data+**B.data) V;
    matrix<V> *C;
    if (A.shape[0]==B.shape[0] && A.shape[1]==B.shape[1]) {
        C = new matrix<V>(A.shape[0],A.shape[1]);
        for (int row=0; row<A.shape[0]; row++) {
            for (int col=0; col<A.shape[1]; col++) {
                C->data[row][col] = A.data[row][col]+B.data[row][col];
            }
        }
    } else {
        C = nullptr;
        cout << "matrix::operator+(const matrix&,const matrix&) :";
        cout << " can't add " << A.shape[0] << "x" << A.shape[1];
        cout << " and " << B.shape[0] << "x" << B.shape[1] << " matrices." << endl;
    }
    return *C;
}

// matrix-matrix subtraction: A-B
template <typename T> matrix<T>& operator-(const matrix<T>& A,const matrix<T>& B) {
    matrix<T> *C;
    if (A.shape[0]==B.shape[0] && A.shape[1]==B.shape[1]) {
        C = new matrix<T>(A.shape[0],A.shape[1]);
        for (int row=0; row<A.shape[0]; row++) {
            for (int col=0; col<A.shape[1]; col++) {
                C->data[row][col] = A.data[row][col]-B.data[row][col];
            }
        }
    } else {
        C = nullptr;
        cout << "matrix::operator+(const matrix&,const matrix&) :";
        cout << " can't subtract " << B.shape[0] << "x" << B.shape[1];
        cout << " from " << A.shape[0] << "x" << A.shape[1] << " matrix." << endl;
    }
    return *C;
}
template <typename T, typename U> auto operator-(const matrix<T>& A,const matrix<U>& B) -> matrix<decltype(**A.data-**B.data)> {
    typedef decltype(**A.data+**B.data) V;
    matrix<V> *C;
    if (A.shape[0]==B.shape[0] && A.shape[1]==B.shape[1]) {
        C = new matrix<V>(A.shape[0],A.shape[1]);
        for (int row=0; row<A.shape[0]; row++) {
            for (int col=0; col<A.shape[1]; col++) {
                C->data[row][col] = A.data[row][col]-B.data[row][col];
            }
        }
    } else {
        C = nullptr;
        cout << "matrix::operator+(const matrix&,const matrix&) :";
        cout << " can't subtract " << B.shape[0] << "x" << B.shape[1];
        cout << " from " << A.shape[0] << "x" << A.shape[1] << " matrix." << endl;
    }
    return *C;
}

// matrix matrix multiplication C=A*B
template <typename T>
matrix<T>& operator*(const matrix<T>& A,const matrix<T>& B) {
    matrix<T> *C;
    if (A.shape[1]==B.shape[0]) {
        C = &zeros<T>(A.shape[0],B.shape[1]);
        for (int row=0; row<C->shape[0]; row++) {
            for (int col=0; col<C->shape[1]; col++) {
                for(int i=0; i<A.shape[1]; i++) {
                    C->data[row][col] += A.data[row][i]*B.data[i][col];
                }
            }
        }
    } else {
        C = nullptr;
        cout << "matrix::operator*(const matrix& A,const matrix& B) :";
        cout << " can't multiply " << A.shape[0] << "x" << A.shape[1];
        cout << " and " << B.shape[0] << "x" << B.shape[1] << " matrices." << endl;
    }
    return *C;   
}
template <typename T, typename U>
auto operator*(const matrix<T>& A,const matrix<U>& B) -> matrix<decltype(**A.data*(**B.data))> {
    typedef decltype(**A.data*(**B.data)) V;
    matrix<V> *C;
    if (A.shape[1]==B.shape[0]) {
        C = &zeros<V>(A.shape[0],B.shape[1]);
        for (int row=0; row<C->shape[0]; row++) {
            for (int col=0; col<C->shape[1]; col++) {
                for(int i=0; i<A.shape[1]; i++) {
                    C->data[row][col] += A.data[row][i]*B.data[i][col];
                }
            }
        }
    } else {
        C = nullptr;
        cout << "matrix::operator*(const matrix& A,const matrix& B) :";
        cout << " can't multiply " << A.shape[0] << "x" << A.shape[1];
        cout << " and " << B.shape[0] << "x" << B.shape[1] << " matrices." << endl;
    }
    return *C;   
}

// matrix power: B=A^n
template <typename T>
matrix<T>& operator^(const matrix<T>& A, const int n) {
    matrix<T> *B, *C;
    if (A.shape[0]==A.shape[1] && n>=0) {
        B = new matrix<T>(A.shape[0],A.shape[1]); *B=A;
        C = new matrix<T>(A.shape[0],A.shape[1]); *C=0.;
        for(int power=1; power<n; power++) {
            for (int row=0; row<C->shape[0]; row++) {
                for (int col=0; col<C->shape[1]; col++) {
                    for(int i=0; i<A.shape[1]; i++) {
                        C->data[row][col] += A.data[row][i]*B->data[i][col];
                    }
                }
            }
            *B=*C;
            *C=0.;
        }
        C->~matrix();
    } else if (A.shape[0]!=A.shape[1]) {
        B = C = nullptr;
        cout << "matrix::operator^(const matrix&,const int) :";
        cout << " can't compute power of " << A.shape[0] << "x" << A.shape[1] << " matrix." << endl;
    } else {
        B = C = nullptr;
        cout << "matrix::operator^(const matrix&,const int) :";
        cout << " can't compute negative power of matrix right now." << endl;
    }
    return *B;   
}


// elementwise matrix matrix product C=A.times(B)
template <typename T> matrix<T>& matrix<T>::times(const matrix<T>& B) {
    matrix *C;
    if (shape[0]==B.shape[0] && shape[1]==B.shape[1]) {
        C = new matrix(shape[0],shape[1]);
        for (int row=0; row<shape[0]; row++) {
            for (int col=0; col<shape[1]; col++) {
                C->data[row][col] = data[row][col]*B.data[row][col];
            }
        }
    } else {
        C = nullptr;
        cout << "matrix::times(const matrix&) :";
        cout << " can't elementwise multiply " << shape[0] << "x" << shape[1];
        cout << " and " << B.shape[0] << "x" << B.shape[1] << " matrices." << endl;
    }

    return *C;
}

// matrix elementwise power B=A.pow(n)
template <typename T>
matrix<T>& matrix<T>::power(const double exp) {
    matrix<T> *B = new matrix<T>(shape[0],shape[1]);
    for (int row=0; row<B->shape[0]; row++) {
        for (int col=0; col<B->shape[1]; col++) {
            B->data[row][col] = pow(B->data[row][col],exp);
        }
    }
    return *B;
}

// matrix scalar multiplication B=c*A or B=A*c
template <typename T> matrix<T>& operator*(const T c,matrix<T>& A) {
    matrix<T> *B = new matrix<T>(A.shape[0],A.shape[1]);
    for (int row=0; row<B->shape[0]; row++) {
        for (int col=0; col<B->shape[1]; col++) {
            B->data[row][col] = c*A.data[row][col];
        }
    }
    return *B;
}
template <typename T> matrix<T>& operator*(matrix<T>& A,const T c) {
    matrix<T> *B = new matrix<T>(A.shape[0],A.shape[1]);
    for (int row=0; row<B->shape[0]; row++) {
        for (int col=0; col<B->shape[1]; col++) {
            B->data[row][col] = c*A.data[row][col];
        }
    }
    return *B;
}
template <typename T, typename U>
auto operator*(matrix<T>& A,const U c) -> matrix<decltype(**A.data*c)> {
    typedef decltype(**A.data*c) V;
    matrix<V> *B = new matrix<V>(A.shape[0],A.shape[1]);
    for (int row=0; row<B->shape[0]; row++) {
        for (int col=0; col<B->shape[1]; col++) {
            B->data[row][col] = c*A.data[row][col];
        }
    }
    return *B;    
}
template <typename T, typename U>
auto operator*(const U c,matrix<T>& A) -> matrix<decltype(**A.data*c)> {
    typedef decltype(**A.data*c) V;
    matrix<V> *B = new matrix<V>(A.shape[0],A.shape[1]);
    for (int row=0; row<B->shape[0]; row++) {
        for (int col=0; col<B->shape[1]; col++) {
            B->data[row][col] = c*A.data[row][col];
        }
    }
    return *B;    
}

// print
template <typename T> void print(matrix<T>& A) {
    for (int row=0; row<A.shape[0]; row++) {
        for (int col=0; col<A.shape[1]; col++) {
            cout << A(row,col) << " ";
        }
        cout << endl;
    }
}

// zero matrix
template <typename T> matrix<T>& zeros(int m, int n) {
    matrix<T> *A = new matrix<T>(m,n); *A = 0.; return *A;
}


template <typename T>
T norm(matrix<T> x) {
  T mag = 0;
  if (x.shape[0] == 1) {
    for (int i=0; i<x.shape[1]; i++) {
      mag += pow(x.data[0][i], 2);
    }
  } else if (x.shape[1] == 1) {
    for (int i=0; i<x.shape[0]; i++) {
      mag += pow(x.data[i][0], 2.);
    }
  } else {
    cout << "norm: only vectors supported right now" << endl;
  }
  mag = pow(mag, .5);
  return mag;
}

template <typename T>
T dot(matrix<T> x, matrix<T> y) {
  T prod = 0;
  if (x.shape[0]==y.shape[0] && x.shape[1]==y.shape[1]) {
    for (int row=0; row<x.shape[0]; row++) {
      for (int col=0; col<x.shape[1]; col++) {
        prod += x.data[row][col] * y.data[row][col];
      }
    }
  } else {
    cout << "dot: matrices not same size" << endl;
  }
  return prod;
}

template <typename T>
matrix<T>& transpose(const matrix<T>& A) {
  matrix<T> *At;
  At = new matrix<T>(A.shape[1],A.shape[0]);
  for (int row=0; row<A.shape[0]; row++) {
    for (int col=0; col<A.shape[1]; col++) {
      At->data[col][row] = A.data[row][col];
    }
  }
  return *At;
}

#endif
