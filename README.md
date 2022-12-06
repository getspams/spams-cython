# cyspams
This package exposes convenient interfaces of some functions of the `SPArse Modeling Software (SPAMS)` C++ library, allowing them to be used from Cython code by releasing the Python Global Interpreter Lock (GIL)

## Installation from PyPI
```Shell
pip install cyspams
```

## Installation from source
```Shell
git clone https://github.com/getspams/spams-cython.git
cd spams-cython
pip install .
```

## Use `cyspams` into your project
### pyproject.toml
```TOML
[build-system]
requires = [
    "setuptools",
    "Cython",
    "cyspams"
]
```
### your_module.pyx
```Cython
from cyspams.interfaces cimport nnls, lasso
```

> **Note**
>
> To build your project with `cyspams` you need to have a BLAS/LAPACK library on your system (e.g. OpenBLAS, Intel MKL)
