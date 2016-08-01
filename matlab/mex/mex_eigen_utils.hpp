/* MEX_EIGEN_UTILS

 * Originally developed for a C++ implementation to sample (2D) SLDS with Gaussian likelihoods
 * and (optional) spatial term. See
 *   J. F. P. Kooij, T-PAMI 2015
 *
 * Requires Eigen and GNU Scientific Library.
 * Compile with mex -IC:\Libraries\eigen mex_slds.cpp
 *
 * (c) Julian Kooij, University of Amsterdam
 *
 */

#pragma once

#include "mex.h"

#include <iostream>
#include <vector>
#include <limits>

/*The mex library*/
#include <math.h>
#include <stdlib.h>

#include <Eigen/Eigen>

extern void _main();

#if 0
    #define DEBUG_PRINTF mexPrintf
    #define DEBUG_SHOW_D(x) DEBUG_PRINTF(#x " = %d\n", (x));
    #define DEBUG_SHOW_F(x) DEBUG_PRINTF(#x " = %.4f\n", (x));
#else
    #define DEBUG_PRINTF
    #define DEBUG_SHOW_D(x)
    #define DEBUG_SHOW_F(x)
#endif

#define ASSERT_RUNTIME(cond, msg) if (!(cond)) { throw std::exception("ERROR: " msg " (failed (" #cond ") )"); }

/* ------------------------------ */
/* Setup typedefs                 */
/* ------------------------------ */

typedef Eigen::Vector2d Vec2;
typedef Eigen::Matrix2d Mat2x2;
typedef Eigen::Map<Eigen::MatrixXd> MapXd;
typedef Eigen::Map<Mat2x2> MapMat2x2;
typedef Eigen::Map<Vec2> MapVec2;
typedef Eigen::ArrayXd DistType;
typedef Eigen::Ref<Eigen::ArrayXi> LabelRef;

template <typename Derived>
void print_eigen(const Eigen::DenseBase<Derived> &m, const char *label)
{
    std::ostringstream ss;
    ss << m;
    mexPrintf(" --- %s (%d x %d) ---\n", label, m.rows(), m.cols());
    mexPrintf("%s\n", ss.str().c_str());
    mexPrintf(" --- %s (%d x %d) ---\n", label, m.rows(), m.cols());
}

/* ------------------------------ */
/* Mex conversion utilities       */
/* ------------------------------ */

template <typename Derived>
Eigen::Map<Derived> mex_to_eigen_map(const mxArray *arr)
{
    size_t rows, cols;

    if (Derived::RowsAtCompileTime == Eigen::Dynamic || Derived::ColsAtCompileTime == Eigen::Dynamic)
    {
        // determine size from mxArray
        const mwSize * dims = mxGetDimensions(arr);
        rows = dims[0];
        cols = dims[1];
    }
    else
    {
        rows = Derived::RowsAtCompileTime;
        cols = Derived::ColsAtCompileTime;
    }

    Eigen::Map<Derived> map(mxGetPr(arr), rows, cols);
    return map;
}

template <typename Derived>
Eigen::Map<Derived> mex_to_eigen_map(const mxArray *arr, size_t poffset)
{
    /* With pointer offset (only supported when size is known at compile time */

    ASSERT_RUNTIME(Derived::RowsAtCompileTime != Eigen::Dynamic, "number of rows must be known at runtime to use pointer offset parameter");
    ASSERT_RUNTIME(Derived::ColsAtCompileTime != Eigen::Dynamic, "number of columns must be known at runtime to use pointer offset parameter");

    // size is known at compile-time
    const size_t rows = Derived::RowsAtCompileTime;
    const size_t cols = Derived::ColsAtCompileTime;

    // NOTE: here, pointer offset is expressed in terms of output items
    Eigen::Map<Derived> map(mxGetPr(arr) + poffset*rows*cols, rows, cols);
    return map;
}

template <typename Derived>
void mex_to_eigen(const mxArray *arr, Eigen::DenseBase<Derived> &out)
{
    MapXd map( mex_to_eigen_map<Eigen::MatrixXd>(arr) );
    out = map.cast<Derived::Scalar>();
}

template <typename Derived>
mxArray * mex_from_eigen(const Eigen::MatrixBase<Derived> &m)
{
    const mwSize rows = m.rows();
    const mwSize cols = m.cols();
    mwSize dims[2] = {rows, cols};
    mxArray *arr = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

    MapXd map( mex_to_eigen_map<Eigen::MatrixXd>(arr) );
    map = m.cast<double>();

    return arr;
}

template <typename Derived>
mxArray * mex_from_eigen(const Eigen::ArrayBase<Derived> &m)
{
    return mex_from_eigen(m.matrix());
}

template <typename T>
void mex_to_vector(const mxArray *arr, std::vector<int>& vec, int valoffset = 0)
{
        size_t T = mxGetNumberOfElements(arr);
        vec.resize(T);

        double *pVals = mxGetPr(arr);
        for(size_t t = 0; t < T; ++t, ++pVals) {
            vec[t] = static_cast<int>(*pVals) + valoffset;
        }
}

template <typename T>
mxArray * mex_from_vector(const std::vector<T> &v, int valoffset = 0)
{
    mxArray *arr = mxCreateNumericMatrix(1, v.size(), mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(arr);

    std::vector<int>::const_iterator itr;
    for(itr = v.begin(); itr != v.end(); ++itr, ++p) {
        *p = (*itr) + valoffset;
    }

    return arr;
}

template <typename T>
mxArray * mex_from_scalar(const T& value)
{
    mxArray *arr = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(arr);
    *p = static_cast<double>(value);
    return arr;
}

/* ------------------------------ */
/* Utilities                      */
/* ------------------------------ */

template <typename T>
inline bool isNaN(const T &x) {
    return (x != x);
}
