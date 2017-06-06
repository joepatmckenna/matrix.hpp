#include <fstream>
#include <ctime>
#include "../../src/matrix.hpp"
#include "cg.cpp"
using namespace std;

int main() {

  int m, start_s, stop_s;
  char filename[64];
  ifstream ifile;
  ofstream ofile;
  matrix<real> times(99,1);
  times=0;

  for (int n=1; n<100; n++) {

    start_s = clock();

    m = 10*n;
    matrix<real> X(m,n+1), Xt(n+1,m), y(m,1);
    matrix<real> A(n+1,n+1), b(n+1,1), c(n+1,1), c0(n+1,1);

    sprintf(filename, "data/c_%i.txt", n);
    ifile.open(filename, ifstream::in);
    if(ifile.is_open()) {
      for(int i=0; i<=n; i++) {
        ifile >> c(i);
      }
    }

    for (int dataset=0; dataset<10; dataset++) {

      // read data
      sprintf(filename, "data/xy_%i_%i.txt", n, dataset);
      ifile.open(filename, ifstream::in);
      if(ifile.is_open()) {
        for(int i=0; i<m; i++) {
          X(i,0) = 1.;
          for(int j=1; j<=n; j++) {
            ifile >> X(i,j);
          }
          ifile >> y(i,0);
        }
      }

      // conj grad
      Xt = transpose(X);
      A = Xt*X;
      b = Xt*y;
      c0 = 0.;
      c0 = conjugate_gradient(A, b, c0, 1e-10, n);

      if (n==1 && dataset == 0) {
        ofile.open("data/c_1_min.txt", ios::trunc);
        ofile << c0(0,0) << " " << c0(1,0) << endl;
        ofile.close();
      }

      stop_s = clock();
      times(n-1,0) = times(n-1,0) + (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }
  }

  times = 1000/10*times;

  ofile.open("data/times.txt", ios::trunc);
  for (int i=1; i<100; i++) {
    ofile << times(i-1,0) << endl;
  }
  ofile.close();

  return 0;
}
