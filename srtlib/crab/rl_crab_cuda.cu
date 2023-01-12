#include<unistd.h>
#include "mir_math.h"
#include "rl_crab_cuda.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"

__global__
void SrtlibRlCrabCuda::VecDiv(
    const double* const vec1_arr,
    const double* const vec2_arr,
    int nsize,
    double* const vec3_arr)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < nsize){
        vec3_arr[index] = vec1_arr[index] / vec2_arr[index];
    }
    __syncthreads();
}

__global__
void SrtlibRlCrabCuda::VecMul(
    const double* const vec1_arr,
    const double* const vec2_arr,
    int nsize,
    double* const vec3_arr)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < nsize){    
        vec3_arr[index] = vec1_arr[index] * vec2_arr[index];
    }
    __syncthreads();
}

void SrtlibRlCrabCuda::GetDetArr(
    cublasHandle_t handle,
    const double* const sky_dev_arr,
    const double* const resp_norm_mat_dev_arr,
    int ndet, int nsky,
    double* const det_dev_arr) // ndet
{
    double alpha = 1.0;
    double beta = 0.0;
    cublasDgemv(handle, CUBLAS_OP_N,
		ndet, nsky, &alpha,
		resp_norm_mat_dev_arr, ndet,
		sky_dev_arr, 1,
		&beta, det_dev_arr, 1);
}

void SrtlibRlCrabCuda::GetDenArr(
    cublasHandle_t handle,
    const double* const sky_dev_arr,
    const double* const flux_dev_arr,
    const double* const det_0_dev_arr,
    const double* const bg_dev_arr,
    const double* const resp_norm_mat_dev_arr,
    int ndet, int nsky, int nphase,
    double* const* const den_dev_arr)
{
    double* det_dev_arr = NULL;
    size_t mem_size_ndet = ndet * sizeof(double);
    cudaMalloc((void **)&det_dev_arr, mem_size_ndet);
    SrtlibRlCrabCuda::GetDetArr(handle, sky_dev_arr,
                                resp_norm_mat_dev_arr,
                                ndet, nsky, det_dev_arr);
    double* flux_arr = new double[nphase];
    cublasGetVector(nphase, sizeof(double), flux_dev_arr, 1,
                    flux_arr, 1);
    for(int iphase = 0; iphase < nphase; iphase++){
        double flux_tmp = flux_arr[iphase];
        cublasDcopy(handle, ndet,
                    bg_dev_arr, 1,
                    den_dev_arr[iphase], 1);
        cublasDaxpy(handle, ndet,
                    &flux_tmp,
                    det_0_dev_arr, 1,
                    den_dev_arr[iphase], 1);
        double alpha = 1.0;
        cublasDaxpy(handle, ndet,
                    &alpha, det_dev_arr, 1,
                    den_dev_arr[iphase], 1);
    }
    delete [] flux_arr;
    cudaFree(det_dev_arr);
}


void SrtlibRlCrabCuda::GetYDashArr(
    const double* const* const data_dev_arr,
    const double* const* const den_dev_arr,
    int ndet, int nphase,
    double* const* const y_dash_dev_arr)
{
    int blocksize = 512;
    dim3 block (blocksize, 1, 1);
    dim3 grid (ndet / block.x + 1, 1, 1);
    for(int iphase = 0; iphase < nphase; iphase++){
        SrtlibRlCrabCuda::VecDiv<<<grid,block>>>(
            data_dev_arr[iphase], den_dev_arr[iphase],
            ndet, y_dash_dev_arr[iphase]);
    }
}

void SrtlibRlCrabCuda::GetMvalArr(
    cublasHandle_t handle,
    const double* const* const y_dash_dev_arr,
    const double* const resp_norm_mat_dev_arr,
    const double* const sky_dev_arr,
    int ndet, int nsky, int nphase,
    double* const mval_dev_arr)
{
    double* y_dash_sum_dev_arr = NULL;
    size_t mem_size_ndet = ndet * sizeof(double);
    cudaMalloc((void **)&y_dash_sum_dev_arr, mem_size_ndet);
    cublasDcopy(handle, ndet,
                y_dash_dev_arr[0], 1,
                y_dash_sum_dev_arr, 1);
    for(int iphase = 1; iphase < nphase; iphase++){
        double alpha = 1.0;
        cublasDaxpy(handle, ndet,
                    &alpha, y_dash_dev_arr[iphase], 1,
                    y_dash_sum_dev_arr, 1);
    }
    double* coeff_dev_arr = NULL;
    size_t mem_size_nsky = nsky * sizeof(double);
    cudaMalloc((void **)&coeff_dev_arr, mem_size_nsky);
    double alpha = 1.0;
    double beta = 0.0;
    cublasDgemv(handle, CUBLAS_OP_T,
                ndet, nsky, &alpha,
                resp_norm_mat_dev_arr, ndet,
                y_dash_sum_dev_arr, 1,
                &beta, coeff_dev_arr, 1);
    int blocksize = 512;
    dim3 block (blocksize, 1, 1);
    dim3 grid (nsky / block.x + 1, 1, 1);
    SrtlibRlCrabCuda::VecMul<<<grid,block>>>(
        coeff_dev_arr, sky_dev_arr, nsky, mval_dev_arr);
    cudaFree(y_dash_sum_dev_arr);
    cudaFree(coeff_dev_arr);
}

void SrtlibRlCrabCuda::GetNvalArr(
    cublasHandle_t handle,
    const double* const* const y_dash_dev_arr,
    const double* const flux_dev_arr,
    const double* const det_0_dev_arr,
    int ndet, int nphase,
    double* const nval_dev_arr)
{
    double* dot_arr = new double[nphase];
    for(int iphase = 0; iphase < nphase; iphase++){
        double dot = 0.0;
        cublasDdot(handle, ndet,
                   y_dash_dev_arr[iphase], 1,
                   det_0_dev_arr, 1, &dot);
        dot_arr[iphase] = dot;
    }
    double* dot_dev_arr = NULL;
    size_t mem_size_nphase = nphase * sizeof(double);
    cudaMalloc((void **)&dot_dev_arr, mem_size_nphase);
    cublasSetVector(nphase, sizeof(double), dot_arr, 1,
		    dot_dev_arr, 1);
    int blocksize = 512;
    dim3 block (blocksize, 1, 1);
    dim3 grid (nphase / block.x + 1, 1, 1);
    SrtlibRlCrabCuda::VecMul<<<grid,block>>>(flux_dev_arr,
                                             dot_dev_arr,
                                             nphase,
                                             nval_dev_arr);
    delete [] dot_arr;
    cudaFree(dot_dev_arr);
}
