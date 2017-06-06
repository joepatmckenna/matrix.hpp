// conjugate gradient
template <typename T>
matrix<T> conjugate_gradient(matrix<T> A, matrix<T> b, matrix<T> c0, T tol, int max_iters) {
  int iters=0, n=A.shape[0]-1;
  matrix<T> c(n+1,1), p(n+1,1), r(n+1,1), Ap(n+1,1);
  T alpha, beta, rtr;

  p = r = b - A*c0;
  rtr = dot(r,r);
  while (norm(r) > tol and iters < max_iters) {
    cout << iters << " " << norm(r) << endl;
    Ap = A*p;
    alpha = rtr / dot(p, Ap);
    c = c + alpha*p;
    beta = rtr;
    r = r - alpha*Ap;
    rtr = dot(r,r);
    beta = rtr / beta;
    p = r + beta*p;
    iters += 1;
  }
  return c;
}
