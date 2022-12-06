/**
 * @file nnls.h
 * C++ implementation of the Non-Negative Least Squares (NNLS) algorithm
 * Based on the Fortran code in:
 * LAWSON, Charles L.; HANSON, Richard J. Solving least squares problems. Society for Industrial and Applied Mathematics, 1995.
 * 
 * Diffusion Imaging and Connectivity Estimation (DICE) lab
 * 
 * @date 29/11/2022
*/

#ifndef NNLS_H
#define NNLS_H

#include <cstdlib>
#include <cmath>
#include <algorithm>

void _solve_triangular_system(const int nsetp, double *zz, const double *a, int jj, const int m, const int *index);
int _householder_transform(const int mode, const int lpivot, const int l1, const int m, double *u, const int u_dim1, double &up, double *cm, const int ice, const int icv, const int ncv);

/**
 * @brief Non-Negative Least Squares (NNLS) algorithm
 * 
 * @param[in, out] A Matrix 'a = (m, n)' stored as 1D contiguous array (column-major order). On exit 'a' contains the product matrix 'q*a', where 'q' is an '(m, m)' orthogonal matrix generated implicitly by this function
 * @param[in, out] y Vector 'b = (m)'. On exit 'b' contains 'q*b', where 'q' is an '(m, m)' orthogonal matrix generated implicitly by this function
 * @param[in] m Dimension 'm'
 * @param[in] n Dimension 'n'
 * 
 * @param[out] x Solution vector 'x = (n)'
 * @param[out] rnorm Squared Euclidean norm of the final residual vector
*/
int _nnls(double *a, double *b, const int m, const int n, double *x, double &rnorm)
{
    // Check data
    if (m <= 0 || n <= 0 || a == NULL || b == NULL || x == NULL)
        return 2;

    // Allocate and zero-initialize the working space
    int *index = new int[n]{};
    double *w = new double[n]{};
    double *zz = new double[m]{};

    int ni;
    int iz1 = 0;
    int iz2 = n - 1;
    int nsetp = 0;
    int iz;
    double sm;
    int mi;
    int npp1 = 0;
    int izmax;
    double wmax;
    int j;
    double asave;
    double up;
    double unorm;
    double d;
    double ztest;
    int jz;
    int jj;
    int ip;
    int ii;
    int iter = 0;
    int itmax = n < 3 ? n * 3 : n * m;
    double alpha;
    double t;
    int pfeas;
    int k;
    double cc;
    double ss;
    double xr;
    double yr;
    int nj;
    double temp;

    // Initialize 'x' and 'index'
    for (ni = 0; ni < n; ni++)
    {
        x[ni] = 0.0;
        index[ni] = ni;
    }

    // Main loop; quit if all coefficients are already in the solution or if 'm' cols of 'a' have been triangulated
    while (iz1 <= iz2 && nsetp < m)
    {
        // Compute components of the dual (negative gradient) vector 'w'
        for (iz = iz1; iz <= iz2; iz++)
        {
            ni = index[iz];
            sm = 0.0;
            for (mi = npp1; mi < m; mi++)
                sm += a[ni * m + mi] * b[mi];
            w[ni] = sm;
        }

        izmax = 0;
        for (;;)
        {
            // Find largest positive 'w[j]'
            wmax = 0.0;
            for (iz = iz1; iz <= iz2; iz++)
            {
                if (w[index[iz]] > wmax)
                {
                    wmax = w[index[iz]];
                    izmax = iz;
                }
            }
            // Terminate if 'wmax <= 0.0' (it indicates satisfaction of the Kuhn-Tucker conditions)
            if (wmax <= 0.0)
                break;

            // The sign of 'w[j]' is ok for 'j' to be moved to set 'p'. Begin the transformation and check new diagonal element to avoid near linear dependence
            j = index[izmax];
            asave = a[j * m + npp1];
            up = 0.0;
            _householder_transform(1, npp1, npp1 + 1, m, &a[j * m], 1, up, NULL, 1, 1, 0);
            unorm = 0.0;
            if (nsetp != 0)
                for (mi = 0; mi < nsetp; mi++)
                    unorm += pow(a[j * m + mi], 2);
            unorm = sqrt(unorm);
            d = unorm + fabs(a[j * m + npp1]) * 0.01;
            if (d - unorm > 0.0)
            {
                // Col 'j' is sufficiently independent. Copy 'b' into 'zz', update 'zz' and solve for 'ztest' (= proposed new value for 'x[j]')
                std::copy(b, b + m, zz);
                _householder_transform(2, npp1, npp1 + 1, m, &a[j * m], 1, up, zz, 1, 1, 1);
                ztest = zz[npp1] / a[j * m + npp1];
                if (ztest > 0.0)
                    break;
            }
            // Reject 'j' as a candidate to be moved from 'zz' to 'p'. Restore 'a[npp1, j]', set 'w[j] = 0.0', and loop back to test dual coefficients again
            a[j * m + npp1] = asave;
            w[j] = 0.0;
        }
        // Terminate if 'wmax <= 0.0' (it indicates satisfaction of the Kuhn-Tucker conditions)
        if (wmax <= 0.0)
            break;

        // Index 'j = index[izmax]' has been selected to be moved from 'zz' to 'p'. Update 'b' and indices, apply householder transformations to cols in new set 'zz', zero sub-diagonal elements in col 'j', set 'w[j] = 0.0'
        std::copy(zz, zz + m, b);
        index[izmax] = index[iz1];
        index[iz1] = j;
        iz1++;
        nsetp = npp1 + 1;
        npp1++;
        if (iz1 <= iz2)
            for (jz = iz1; jz <= iz2; jz++)
            {
                jj = index[jz];
                _householder_transform(2, nsetp - 1, npp1, m, &a[j * m], 1, up, &a[jj * m], 1, m, 1);
            }
        if (nsetp != m)
            for (mi = npp1; mi < m; mi++)
                a[j * m + mi] = 0.0;
        w[j] = 0.0;
        // Solve the triangular system; store the solution temporarily in 'zz'
        _solve_triangular_system(nsetp, zz, a, jj, m, index);

        // Secondary loop
        while (++iter < itmax)
        {
            // See if all new constrained coefficients are feasible; if not, compute alpha
            alpha = 2.0;
            for (ip = 0; ip < nsetp; ip++)
            {
                ni = index[ip];
                if (zz[ip] <= 0.0)
                {
                    t = -x[ni] / (zz[ip] - x[ni]);
                    if (alpha > t)
                    {
                        alpha = t;
                        jj = ip - 1;
                    }
                }
            }
            // If all new constrained coefficients are feasible then still 'alpha == 2'. If so, then exit from the secondary loop to main loop
            if (alpha == 2.0)
                break;

            // Use '0.0 < alpha < 1.0' to interpolate between old 'x' and new 'zz'
            for (ip = 0; ip < nsetp; ip++)
            {
                ni = index[ip];
                x[ni] += alpha * (zz[ip] - x[ni]);
            }

            // Modify 'zz' and 'b' and 'index' to move coefficient 'i' from 'p' to 'zz'
            pfeas = 1;
            k = index[jj + 1];
            do
            {
                x[k] = 0.0;
                if (jj != nsetp - 1)
                {
                    jj++;
                    for (ni = jj + 1; ni < nsetp; ni++)
                    {
                        ii = index[ni];
                        index[ni - 1] = ii;
                        // Compute orthogonal rotation matrix
                        if (fabs(a[ii * m + ni - 1]) > fabs(a[ii * m + ni]))
                        {
                            xr = a[ii * m + ni] / a[ii * m + ni - 1];
                            yr = hypot(xr, 1.0);
                            cc = copysign(1.0 / yr, a[ii * m + ni - 1]);
                            ss = cc * xr;
                            a[ii * m + ni - 1] = fabs(a[ii * m + ni - 1]) * yr;
                        }
                        else if (a[ii * m + ni] != 0.0)
                        {
                            xr = a[ii * m + ni - 1] / a[ii * m + ni];
                            yr = hypot(xr, 1.0);
                            ss = copysign(1.0 / yr, a[ii * m + ni]);
                            cc = ss * xr;
                            a[ii * m + ni - 1] = fabs(a[ii * m + ni]) * yr;
                        }
                        else
                        {
                            a[ii * m + ni - 1] = 0.0;
                            cc = 0.0;
                            ss = 1.0;
                        }
                        a[ii * m + ni] = 0.0;
                        for (nj = 0; nj < n; nj++)
                            if (nj != ii)
                            {
                                temp = a[nj * m + ni - 1];
                                a[nj * m + ni - 1] = cc * temp + ss * a[nj * m + ni];
                                a[nj * m + ni] = -ss * temp + cc * a[nj * m + ni];
                            }
                        temp = b[ni - 1];
                        b[ni - 1] = cc * temp + ss * b[ni];
                        b[ni] = -ss * temp + cc * b[ni];
                    }
                }
                npp1 = nsetp - 1;
                nsetp--;
                iz1--;
                index[iz1] = k;
                // See if the remaining coefficients in 'p' are feasible; they should be because of the way alpha was determined. If any are infeasible it is due to round-off error. Any that are non-positive will be set to zero and moved from 'p' to 'zz'
                for (jj = 0, pfeas = 1; jj < nsetp; jj++)
                {
                    k = index[jj];
                    if (x[k] <= 0.0)
                    {
                        pfeas = 0;
                        break;
                    }
                }
            } while (pfeas == 0);

            // Copy 'b' into 'zz', then solve again and loop back
            std::copy(b, b + m, zz);
            _solve_triangular_system(nsetp, zz, a, jj, m, index);
        }
        if (iter >= itmax)
            break;

        for (ip = 0; ip < nsetp; ip++)
        {
            k = index[ip];
            x[k] = zz[ip];
        }
    }

    // Compute the norm of the final residual vector
    if (rnorm != NULL)
    {
        sm = 0.0;
        if (npp1 < m)
            for (mi = npp1; mi < m; mi++)
                sm += pow(b[mi], 2);
        else
            for (ni = 0; ni < n; ni++)
                w[ni] = 0.0;
        rnorm = sqrt(sm);
    }

    // Free working space
    delete[] index;
    delete[] w;
    delete[] zz;

    if (iter >= itmax)
        return 1;
    return 0;
}

void _solve_triangular_system(const int nsetp, double *zz, const double *a, int jj, const int m, const int *index)
{
    int i;
    int ii;
    int ip;
    for (i = 0; i < nsetp; i++)
    {
        ip = nsetp - i - 1;
        if (i != 0)
            for (ii = 0; ii <= ip; ii++)
                zz[ii] -= a[jj * m + ii] * zz[ip + 1];
        jj = index[ip];
        zz[ip] /= a[jj * m + ip];
    }
}

int _householder_transform(const int mode, const int lpivot, const int l1, const int m, double *u, const int u_dim1, double &up, double *cm, const int ice, const int icv, const int ncv)
{
    // Check params
    if (mode != 1 && mode != 2)
        return 1;
    if (m < 1 || u == NULL || u_dim1 < 1)
        return 1;
    if (lpivot < 0 || lpivot >= l1 || l1 > m)
        return 1;

    double cl = fabs(u[lpivot * u_dim1]);
    int i;
    double clinv;
    double sm;
    double b;
    int j;

    if (mode == 2)
    {
        if (cl <= 0.0)
            return 0;
    }
    else
    {
        // Construct the transformation
        for (i = l1; i < m; i++)
            cl = fmax(fabs(u[i * u_dim1]), cl);
        if (cl <= 0.0)
            return 0;
        clinv = 1.0 / cl;
        sm = pow(u[lpivot * u_dim1] * clinv, 2);
        for (i = l1; i < m; i++)
            sm += pow(u[i * u_dim1] * clinv, 2);
        cl *= sqrt(sm);
        if (u[lpivot * u_dim1] > 0.0)
            cl = -cl;
        up = u[lpivot * u_dim1] - cl;
        u[lpivot * u_dim1] = cl;
    }

    // if no vectors where to apply then change pivot vector
    b = up * u[lpivot * u_dim1];
    // 'b' must be non-positive here; if 'b >= 0.0', then return
    if (b >= 0.0)
        return 0;
    // ok, for all vectors we want to apply
    if (cm == NULL)
        return 2;

    for (i = 0; i < ncv; i++)
    {
        sm = cm[lpivot * ice + i * icv] * up;
        for (j = l1; j < m; j++)
            sm += cm[j * ice + i * icv] * u[j * u_dim1];
        if (sm != 0.0)
        {
            sm *= 1.0 / b;
            cm[lpivot * ice + i * icv] += sm * up;
            for (j = l1; j < m; j++)
                cm[j * ice + i * icv] += u[j * u_dim1] * sm;
        }
    }
    return 0;
}

#endif
