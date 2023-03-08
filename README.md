# spams-cython
This package exposes convenient interfaces of some functions of the `SPArse Modeling Software (SPAMS)` C++ library, allowing them to be used from Cython code by releasing the Python Global Interpreter Lock (GIL)

## Installation from PyPI
```Shell
pip install spams-cython
```

## Installation from source
```Shell
git clone https://github.com/getspams/spams-cython.git
cd spams-cython
pip install .
```

## Use `spams-cython` into your project
### pyproject.toml
```TOML
[build-system]
requires = [
    "setuptools",
    "Cython",
    "spams-cython"
]
```
### setup.py
```python
from setuptools import setup, Extension
from Cython.Build import cythonize
import cyspams

extensions = [
    Extension(
        mymodule,
        sources['mymodule.pyx'],
        include_dirs=cyspams.get_include()
    )
]

setup(ext_modules=cythonize(extensions))
```
### mymodule.pyx
```cython
from cyspams.interfaces cimport nnls, lasso
```

> **Note**
>
> To build your project with `spams-cython` you need to have a BLAS/LAPACK library on your system (e.g. OpenBLAS, Intel MKL)

## Exposed functions
### `nnls`
```cython
void nnls(const double *A, const double *y, const int m, const int n, double *x, double &rnorm)
```
__Args in__:\
__A__ → Matrix '_A = (m, n)_' stored as 1D contiguous array (column-major order)\
__y__ → Vector '_y = (m)_'\
__m__ → Dimension '_m_'\
__n__ → Dimension '_n_'

__Args out__:\
__x__ → Solution vector '_x = (n)_'\
__rnorm__ → Squared Euclidean norm of the final residual vector

---
### `lasso`
```cython
void lasso(double *A, double *y, const int m, const int p, const int n, double *x)
void lasso(double *A, double *y, const int m, const int p, const int n, double *x, const double lambda1, const double lambda2)
void lasso(double *A, double *y, const int m, const int p, const int n, double *x, const double lambda1, const double lambda2, const int mode, const bint pos)
void lasso(double *A, double *y, const int m, const int p, const int n, double *x, const double lambda1, const double lambda2, const int mode, const bint pos, const bint ols, const int max_length_path, const int L, const bint cholesky, const int n_threads, const bint verbose)
```
__Args in__:\
__A__ → Matrix '_A = (m, p)_' stored as 1D contiguous array (column-major order)\
__y__ → Matrix '_y = (m, n)_' stored as 1D contiguous array (column-major order)\
__m__ → Dimension '_m_'\
__p__ → Dimension '_p_'\
__n__ → Dimension '_n_'\
__lambda1__ → Regularization parameter (default = 0.0)\
__lambda2__ → Parameter for solving the Elastic-Net (default = 0.0)\
__mode__ → 0 = L1COEFFS, 1 = L2ERROR, 2 = PENALTY, 3 = SPARSITY, 4 = L2ERROR2, 5 = PENALTY2, 6 = FISTAMODE (default = 2)\
__pos__ → Adds non-negativity constraints on the coefficients (default = true)\
__ols__ → Perform an orthogonal projection before returning the solution (default = false)\
__max_length_path__ → Maximum length of the path (default = -1)\
__L__ → Maximum number of steps of the homotopy algorithm. Can be used as a stopping criterion (default = -1)\
__cholesky__ → Choose between Cholesky implementation or one based on the matrix inversion Lemma (default = false)\
__n_threads__ → Number of threads to use (default = 1)\
__verbose__ → Verbose mode (default = false)

__Args out__:\
__x__ → Solution matrix '_x = (p, n)_' stored as 1D contiguous array (column-major order)

---
