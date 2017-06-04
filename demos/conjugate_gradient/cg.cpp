// conjugate gradient
template <typename T>
matrix<T> conjugate_gradient(matrix<T> A, matrix<T> b, matrix<T> c0, real tol) {
  int k=0, n=A.shape[0]-1;
  matrix<T> c(n+1,1), p(n+1,1), r(n+1,1), Ap(n+1,1);
  T alpha, beta, rtr;

  p = r = b - A*c0;
  rtr = dot(r,r);
  while (norm(r) > tol and k < n) {
    Ap = A*p;
    alpha = rtr / dot(p, Ap);
    c = c + alpha*p;
    beta = rtr;
    r = r - alpha*Ap;
    rtr = dot(r,r);
    beta = rtr / beta;
    p = r + beta*p;
    cout << k << " " << norm(r) << endl;
    k += 1;
  }
  return c;
}
