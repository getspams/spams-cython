/**
 * @file interfaces.h
 * Interfaces for NNLS and LASSO (SPAMS) solvers
 * 
 * Diffusion Imaging and Connectivity Estimation (DICE) lab
 * 
 * @date 29/11/2022
*/

#ifndef INTERFACES_H
#define INTERFACES_H

#include "nnls.h"
#include "spams.h"

/**
 * @brief Interface for the '_nnls()' function in 'nnls.h'
 * 
 * @param[in] A Matrix 'A = (m, n)' stored as 1D contiguous array (column-major order)
 * @param[in] y Vector 'y = (m)'
 * @param[in] m Dimension 'm'
 * @param[in] n Dimension 'n'
 * 
 * @param[out] x Solution vector 'x = (n)'
 * @param[out] rnorm Squared Euclidean norm of the final residual vector
*/
void nnls(const double *A, const double *y, const int m, const int n, double *x, double &rnorm)
{
    // NOTE work on a copy of 'A' and 'y' because the '_nnls()' call will modify them
    double *_A = new double[m * n];
    double *_y = new double[m];
    std::copy(A, A + m * n, _A);
    std::copy(y, y + m, _y);

    _nnls(_A, _y, m, n, x, rnorm);
    // NOTE throw exceptions
    // if (ret == 1)
    //     std::cout << "WARNING: NNLS max iterations exceeded" << std::endl;
    // if (ret == 2)
    //     throw std::runtime_error("NNLS failed");

    delete[] _A;
    delete[] _y;
}

/**
 * @brief Interface for the '_lassoD()' function in 'spams.h'
 * 
 * @param[in] A Matrix 'A = (m, p)' stored as 1D contiguous array (column-major order)
 * @param[in] y Matrix 'y = (m, n)' stored as 1D contiguous array (column-major order)
 * @param[in] m Dimension 'm'
 * @param[in] p Dimension 'p'
 * @param[in] n Dimension 'n'
 * @param[in] lambda1 Regularization parameter (default = 0.0)
 * @param[in] lambda2 Parameter for solving the Elastic-Net (default = 0.0)
 * @param[in] mode 0=L1COEFFS, 1=L2ERROR, 2=PENALTY, 3=SPARSITY, 4=L2ERROR2, 5=PENALTY2, 6=FISTAMODE (default = 2)
 * @param[in] pos Adds non-negativity constraints on the coefficients (default = true)
 * @param[in] ols Perform an orthogonal projection before returning the solution (default = false)
 * @param[in] max_length_path Maximum length of the path (default = -1)
 * @param[in] L Maximum number of steps of the homotopy algorithm. Can be used as a stopping criterion (default = -1)
 * @param[in] cholesky Choose between Cholesky implementation or one based on the matrix inversion Lemma (default = false)
 * @param[in] n_threads Number of threads to use (default = 1)
 * @param[in] verbose Verbose mode (default = false)
 * 
 * @param[out] x Solution matrix 'x = (p, n)' stored as 1D contiguous array (column-major order)
*/
void lasso(
    double *A,
    double *y,
    const int m,
    const int p,
    const int n,
    double *x,
    const double lambda1 = 0.0,
    const double lambda2 = 0.0,
    const int mode = 2,
    const bool pos = true,
    const bool ols = false,
    const int max_length_path = -1,
    const int L = -1,
    const bool cholesky = false,
    const int n_threads = 1,
    const bool verbose = false
    )
{
    Matrix<double> *_A = new Matrix<double>(A, m, p);
    Matrix<double> *_y = new Matrix<double>(y, m, n);
    Matrix<double> *path;
    const bool return_reg_path = false;
    Matrix<double> *_x = new Matrix<double>(x, p, n);
    _lassoD(_y, _A, &path, return_reg_path, L, lambda1, lambda2, (constraint_type)mode, pos, ols, n_threads, max_length_path, verbose, cholesky)->toFull(*_x);

    delete _A;
    delete _y;
    delete _x;
}

#endif
