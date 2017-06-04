#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include "../../src/matrix.hpp"
#include "cg.cpp"
using namespace std;

int main() {

  // initialize arrays
  int m=500, n=100;
  matrix<real> X(m,n+1), y(m,1);

  // read data to arrays
  ifstream file("conj_grad_data.txt");
  if(file.is_open()) {
    for(int i=0; i<m; i++) {
      X(i,0) = 1.;
      for(int j=1; j<=n; j++) {
        file >> X(i,j);
      }
      file >> y(i,0);
    }
  }

  matrix<real> A(n+1,n+1), b(n+1,1), c0(n+1,1), Xt = transpose(X);
  A = Xt*X;
  b = Xt*y;
  c0 = 0.;

  matrix<real> x = conjugate_gradient(A, b, c0, 1e-10);

  return 0;
}
