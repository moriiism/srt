#ifndef MORIIISM_SRT_SRTMATHLIB_SRTMATHLIB_H_
#define MORIIISM_SRT_SRTMATHLIB_SRTMATHLIB_H_

#include "mib_blas.h"
#include "mi_sort.h"

namespace SrtMathlib
{
    double GetSum(
        long narr,
        const double* const val_arr);

    double GetRootMeanSquareError(
        const double* const det_arr,
        const double* const val_arr,
        int ndet);

} // namespace SrtMathlib

#endif // MORIIISM_SRT_SRTMATHLIB_SRTMATHLIB_H_
