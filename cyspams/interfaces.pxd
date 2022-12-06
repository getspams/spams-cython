# distutils: language = c++
# cython: language_level = 3

cdef extern from 'interfaces.h' nogil:
    void nnls(const double *A, const double *y, const int m, const int n, double *x, double &rnorm)
    void lasso(double *A, double *y, const int m, const int p, const int n, double *x)
    void lasso(double *A, double *y, const int m, const int p, const int n, double *x, const double lambda1, const double lambda2)
    void lasso(double *A, double *y, const int m, const int p, const int n, double *x, const double lambda1, const double lambda2, const int mode, const bint pos, const bint ols, const int max_length_path, const int L, const bint cholesky, const int n_threads, const bint verbose)
