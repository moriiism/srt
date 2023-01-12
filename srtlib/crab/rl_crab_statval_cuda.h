#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_STATVAL_CUDA_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_STATVAL_CUDA_H_

#include "mib_blas.h"
#include "mir_math.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"

namespace SrtlibCrabRlCrabStatvalCuda
{
    __global__
    void VecDiffSqrt(const double* const vec1_dev_arr,
                     const double* const vec2_dev_arr,
                     int nsize,
                     double* const vec3_dev_arr);
    
    double GetHellingerDist(
        cublasHandle_t handle,    
        const double* const sky_pre_dev_arr,
        const double* const flux_pre_dev_arr,
        const double* const sky_new_dev_arr,
        const double* const flux_new_dev_arr,
        int nsky, int nphase);
        

} // SrtlibCrabRlCrabStatvalCuda

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_STATVAL_CUDA_H_
