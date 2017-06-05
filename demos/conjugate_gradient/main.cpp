#include <fstream>
#include "../../src/matrix.hpp"
#include "cg.cpp"
using namespace std;

int main() {

  int m;
  char filename[64];

  for (int n=2; n<100; n++) {

    m = 10*n;
    matrix<real> X(m,n+1), Xt(n+1,m), y(m,1);
    matrix<real> A(n+1,n+1), b(n+1,1), c(n+1,1), c0(n+1,1);

    sprintf(filename, "data/c_%i.txt", n);
    ifstream file(filename);
    if(file.is_open()) {
      for(int i=0; i<=n; i++) {
        file >> c(i);
      }
    }

    for (int dataset=0; dataset<10; dataset++) {
      // read data
      sprintf(filename, "data/xy_%i_%i.txt", n, dataset);
      ifstream file(filename);
      if(file.is_open()) {
        for(int i=0; i<m; i++) {
          X(i,0) = 1.;
          for(int j=1; j<=n; j++) {
            file >> X(i,j);
          }
          file >> y(i,0);
        }
      }

      // conj grad
      Xt = transpose(X);
      A = Xt*X;
      b = Xt*y;
      c0 = 0.;
      c0 = conjugate_gradient(A, b, c0, 1e-10);
      cout << n << " " << dataset << " " <<  norm(c - c0) << endl;
    }
  }

  return 0;
}
