#ifndef MORIIISM_SRT_SRTLIB_SMTH_ZAL_CUDA_H_
#define MORIIISM_SRT_SRTLIB_SMTH_ZAL_CUDA_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"

namespace SrtlibSmthZalCuda
{
    __global__
    void GetDerivUBetaArr(
        const double* const sky_dash_dev_arr,
        int nskyx, int nskyy,
        double* const beta_dev_arr);
    
} // namespace SrtlibSmthZalCuda

#endif // MORIIISM_SRT_SRTLIB_SMTH_ZAL_CUDA_H_
