#include "rl_crab_statval_cuda.h"

__global__
void SrtlibCrabRlCrabStatvalCuda::VecDiffSqrt(
    const double* const vec1_dev_arr,
    const double* const vec2_dev_arr,
    int nsize,
    double* const vec3_dev_arr)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < nsize){
        vec3_dev_arr[index] = sqrt(vec1_dev_arr[index])
            - sqrt(vec2_dev_arr[index]);
    }
    __syncthreads();
}

double SrtlibCrabRlCrabStatvalCuda::GetHellingerDist(
    cublasHandle_t handle,    
    const double* const sky_pre_dev_arr,
    const double* const flux_pre_dev_arr,
    const double* const sky_new_dev_arr,
    const double* const flux_new_dev_arr,
    int nsky, int nphase)
{
    double sum_sky_pre = 0.0;
    double sum_sky_new = 0.0;

    double* dummy_sky_dev_arr = NULL;
    double* dummy_flux_dev_arr = NULL;
    size_t mem_size_dummy_sky = nsky * sizeof(double);
    size_t mem_size_dummy_flux = nphase * sizeof(double);
    cudaMalloc((void **)&dummy_sky_dev_arr, mem_size_dummy_sky);
    cudaMalloc((void **)&dummy_flux_dev_arr, mem_size_dummy_flux);
    double* dummy_sky_arr = new double[nsky];
    double* dummy_flux_arr = new double[nphase];
    for(int isky = 0; isky < nsky; isky++){
        dummy_sky_arr[isky] = 1.0;
    }
    for(int iphase = 0; iphase < nphase; iphase++){
        dummy_flux_arr[iphase] = 1.0;
    }
    cublasSetVector(nsky, sizeof(double), dummy_sky_arr, 1,
		    dummy_sky_dev_arr, 1);
    cublasSetVector(nphase, sizeof(double), dummy_flux_arr, 1,
		    dummy_flux_dev_arr, 1);
    cublasDdot(handle, nsky, sky_pre_dev_arr, 1,
               dummy_sky_dev_arr, 1, &sum_sky_pre);
    cublasDdot(handle, nsky, sky_new_dev_arr, 1,
               dummy_sky_dev_arr, 1, &sum_sky_new);
    
    double sum_flux_pre = 0.0;
    double sum_flux_new = 0.0;
    cublasDdot(handle, nphase, flux_pre_dev_arr, 1,
               dummy_flux_dev_arr, 1, &sum_flux_pre);
    cublasDdot(handle, nphase, flux_new_dev_arr, 1,
               dummy_flux_dev_arr, 1, &sum_flux_new);
    double sum_pre = sum_sky_pre + sum_flux_pre;
    double sum_new = sum_sky_new + sum_flux_new;

    double* sky_pre_norm_dev_arr = NULL;
    double* sky_new_norm_dev_arr = NULL;
    double* flux_pre_norm_dev_arr = NULL;
    double* flux_new_norm_dev_arr = NULL;
    size_t mem_size_nsky = nsky * sizeof(double);
    size_t mem_size_nphase = nphase * sizeof(double);    
    cudaMalloc((void **)&sky_pre_norm_dev_arr, mem_size_nsky);
    cudaMalloc((void **)&sky_new_norm_dev_arr, mem_size_nsky);    
    cudaMalloc((void **)&flux_pre_norm_dev_arr, mem_size_nphase);
    cudaMalloc((void **)&flux_new_norm_dev_arr, mem_size_nphase);
    cublasDcopy(handle, nsky,
                sky_pre_dev_arr, 1,
                sky_pre_norm_dev_arr, 1);
    cublasDcopy(handle, nsky,
                sky_new_dev_arr, 1,
                sky_new_norm_dev_arr, 1);
    cublasDcopy(handle, nphase,
                flux_pre_dev_arr, 1,
                flux_pre_norm_dev_arr, 1);
    cublasDcopy(handle, nphase,
                flux_new_dev_arr, 1,
                flux_new_norm_dev_arr, 1);

    double alpha = 1.0 / sum_pre;
    cublasDscal(handle, nsky, &alpha,
                sky_pre_norm_dev_arr, 1);
    cublasDscal(handle, nphase, &alpha,
                flux_pre_norm_dev_arr, 1);    
    alpha = 1.0 / sum_new;
    cublasDscal(handle, nsky, &alpha,
                sky_new_norm_dev_arr, 1);
    cublasDscal(handle, nphase, &alpha,
                flux_new_norm_dev_arr, 1);

    double* diff_sqrt_sky_dev_arr = NULL;
    double* diff_sqrt_flux_dev_arr = NULL;
    cudaMalloc((void **)&diff_sqrt_sky_dev_arr, mem_size_nsky);
    cudaMalloc((void **)&diff_sqrt_flux_dev_arr, mem_size_nphase);

    int blocksize = 512;
    dim3 block (blocksize, 1, 1);
    dim3 grid_sky (nsky / block.x + 1, 1, 1);
    SrtlibCrabRlCrabStatvalCuda::VecDiffSqrt<<<grid_sky,block>>>(
        sky_pre_norm_dev_arr,
        sky_new_norm_dev_arr,
        nsky,
        diff_sqrt_sky_dev_arr);
    dim3 grid_flux (nphase / block.x + 1, 1, 1);
    SrtlibCrabRlCrabStatvalCuda::VecDiffSqrt<<<grid_flux,block>>>(
        flux_pre_norm_dev_arr,
        flux_new_norm_dev_arr,
        nphase,
        diff_sqrt_flux_dev_arr);

    double sum_diff_sqrt_sky = 0.0;
    cublasDdot(handle, nsky, diff_sqrt_sky_dev_arr, 1,
               diff_sqrt_sky_dev_arr, 1, &sum_diff_sqrt_sky);
    double sum_diff_sqrt_flux = 0.0;
    cublasDdot(handle, nphase, diff_sqrt_flux_dev_arr, 1,
               diff_sqrt_flux_dev_arr, 1, &sum_diff_sqrt_flux);
    double sum_diff_sqrt = sum_diff_sqrt_sky + sum_diff_sqrt_flux;
    double ans = sqrt(sum_diff_sqrt);

    delete [] dummy_sky_arr;
    delete [] dummy_flux_arr;
    cudaFree(dummy_sky_dev_arr);
    cudaFree(dummy_flux_dev_arr);
    cudaFree(sky_pre_norm_dev_arr); 
    cudaFree(sky_new_norm_dev_arr);
    cudaFree(flux_pre_norm_dev_arr);
    cudaFree(flux_new_norm_dev_arr);
    cudaFree(diff_sqrt_sky_dev_arr);
    cudaFree(diff_sqrt_flux_dev_arr);
    
    return (ans);
}
