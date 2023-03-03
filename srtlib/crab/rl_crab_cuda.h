#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_CUDA_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_CUDA_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
//#include "mir_math.h"
#include "../smth_zal.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"


namespace SrtlibRlCrabCuda
{
    __global__
    void VecDiv(
        const double* const vec1_arr,
        const double* const vec2_arr,
        int nsize,
        double* const vec3_arr);

    __global__
    void VecMul(
        const double* const vec1_arr,
        const double* const vec2_arr,
        int nsize,
        double* const vec3_arr);
    
    void GetDetArr(
        cublasHandle_t handle,
        const double* const sky_dev_arr,
        const double* const resp_norm_mat_dev_arr,
        int ndet, int nsky,
        double* const det_dev_arr);

    void GetDenArr(
        cublasHandle_t handle,
        const double* const sky_dev_arr,
        const double* const flux_dev_arr,
        const double* const det_0_dev_arr,
        const double* const bg_dev_arr,
        const double* const resp_norm_mat_dev_arr,
        int ndet, int nsky, int nphase,
        double* const* const den_dev_arr);

    void GetYDashArr(
        const double* const* const data_dev_arr,
        const double* const* const den_dev_arr,
        int ndet, int nphase,
        double* const* const y_dash_dev_arr);

    void GetMvalArr(
        cublasHandle_t handle,
        const double* const* const y_dash_dev_arr,
        const double* const resp_norm_mat_dev_arr,
        const double* const sky_dev_arr,
        int ndet, int nsky, int nphase,
        double* const mval_dev_arr);

    void GetNvalArr(
        cublasHandle_t handle,
        const double* const* const y_dash_dev_arr,
        const double* const flux_dev_arr,
        const double* const det_0_dev_arr,
        int ndet, int nphase,
        double* const nval_dev_arr);

} // namespace SrtlibRlCrabCuda

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_CUDA_H_
