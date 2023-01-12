#include "../smth_zal_cuda.h"
#include "rl_crab_cuda.h"
#include "rl_crab_statval_cuda.h"
#include "rl_crab_smth_pf_zal_cuda.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"

__global__
void SrtlibRlCrabSmthPfZalCuda::GetSkyNewCuda(
    const double* const alpha_dev_arr,
    const double* const beta_dev_arr,
    const double* const mval_dev_arr,
    int nsky, double mu,
    double* const sky_new_dev_arr)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < nsky){ 
        double bval = 1.0 - mu * beta_dev_arr[index];
        double num = -1.0 * bval
            + sqrt(bval * bval
                   + 4.0 * mu * alpha_dev_arr[index]
                   * mval_dev_arr[index]);
        double den = 2.0 * mu * alpha_dev_arr[index];
        sky_new_dev_arr[index] = num / den;
    }
    __syncthreads();    
}

void SrtlibRlCrabSmthPfZalCuda::GetSkyNewArr(
    cublasHandle_t handle,
    const double* const sky_dev_arr,
    const double* const mval_dev_arr,
    int nskyx, int nskyy, double mu,
    double* const sky_new_dev_arr)
{
    int nsky = nskyx * nskyy;
    double* alpha_arr = new double[nsky];
    SrtlibSmthZal::GetDerivUAlphaArr(
        nskyx, nskyy, alpha_arr);

    double* alpha_dev_arr = NULL;
    size_t mem_size_nsky = nsky * sizeof(double);
    cudaMalloc((void **)&alpha_dev_arr, mem_size_nsky);
    cublasSetVector(nsky, sizeof(double), alpha_arr, 1,
		    alpha_dev_arr, 1);
    
    double* beta_dev_arr = NULL;
    cudaMalloc((void **)&beta_dev_arr, mem_size_nsky);
    int blocksize = 512;
    dim3 block (blocksize, 1, 1);
    dim3 grid  (nsky / block.x + 1, 1, 1);
    SrtlibSmthZalCuda::GetDerivUBetaArr<<<grid,block>>>(
        sky_dev_arr, nskyx, nskyy, beta_dev_arr);
    SrtlibRlCrabSmthPfZalCuda::GetSkyNewCuda<<<grid,block>>>(
        alpha_dev_arr,
        beta_dev_arr,
        mval_dev_arr,
        nsky, mu,
        sky_new_dev_arr);

    cudaFree(alpha_dev_arr);
    cudaFree(beta_dev_arr);
    delete [] alpha_arr;
}

__global__
void SrtlibRlCrabSmthPfZalCuda::GetFluxNewArr(
    const double* const nval_dev_arr,
    const double* const flux_target_dev_arr,
    const double* const phase_dev_arr,
    int nphase, double gamma,
    double* const flux_new_dev_arr)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < nphase){ 
        double bval = phase_dev_arr[index]
            - 2.0 * gamma * flux_target_dev_arr[index];
        double num = -1.0 * bval
            + sqrt(bval * bval
                   + 8.0 * gamma * nval_dev_arr[index]);
        double den = 4.0 * gamma;
        flux_new_dev_arr[index] = num / den;
    }
    __syncthreads();
}

void SrtlibRlCrabSmthPfZalCuda::GetSkyFluxNewArr(
    cublasHandle_t handle,
    const double* const sky_pre_dev_arr,
    const double* const flux_pre_dev_arr,
    const double* const* const data_dev_arr,
    const double* const bg_dev_arr,
    const double* const flux_target_dev_arr,
    const double* const phase_dev_arr,
    const double* const det_0_dev_arr,
    const double* const resp_norm_mat_dev_arr,    
    int ndet, int nskyx, int nskyy, int nphase,
    double mu, double gamma,
    double* const sky_new_dev_arr,
    double* const flux_new_dev_arr)
{
    int nsky = nskyx * nskyy;
    double** den_dev_arr = new double*[nphase];
    double** y_dash_dev_arr = new double*[nphase];
    size_t mem_size_ndet = ndet * sizeof(double);
    size_t mem_size_nsky = nsky * sizeof(double);
    size_t mem_size_nphase = nphase * sizeof(double);
    for(int iphase = 0; iphase < nphase; iphase++){
        cudaMalloc((void **)&den_dev_arr[iphase], mem_size_ndet);
        cudaMalloc((void **)&y_dash_dev_arr[iphase], mem_size_ndet);
    }
    double* mval_dev_arr = NULL;
    double* nval_dev_arr = NULL;
    cudaMalloc((void **)&mval_dev_arr, mem_size_nsky);
    cudaMalloc((void **)&nval_dev_arr, mem_size_nphase);
        
    SrtlibRlCrabCuda::GetDenArr(handle,
                                sky_pre_dev_arr,
                                flux_pre_dev_arr,
                                det_0_dev_arr,
                                bg_dev_arr,
                                resp_norm_mat_dev_arr,
                                ndet, nsky, nphase,
                                den_dev_arr);
    SrtlibRlCrabCuda::GetYDashArr(data_dev_arr,
                                  den_dev_arr,
                                  ndet, nphase,
                                  y_dash_dev_arr);
    SrtlibRlCrabCuda::GetMvalArr(handle,
                                 y_dash_dev_arr,
                                 resp_norm_mat_dev_arr,
                                 sky_pre_dev_arr,
                                 ndet, nsky, nphase,
                                 mval_dev_arr);
    SrtlibRlCrabCuda::GetNvalArr(handle,
                                 y_dash_dev_arr,
                                 flux_pre_dev_arr,
                                 det_0_dev_arr,
                                 ndet, nphase,
                                 nval_dev_arr);
    SrtlibRlCrabSmthPfZalCuda::GetSkyNewArr(
        handle,
        sky_pre_dev_arr,
        mval_dev_arr,
        nskyx, nskyy, mu,
        sky_new_dev_arr);
    int blocksize = 512;
    dim3 block (blocksize, 1, 1);
    dim3 grid  (nphase / block.x + 1, 1, 1);
    SrtlibRlCrabSmthPfZalCuda::GetFluxNewArr<<<grid,block>>>(
        nval_dev_arr,
        flux_target_dev_arr,
        phase_dev_arr,
        nphase, gamma,
        flux_new_dev_arr);
    for(int iphase = 0; iphase < nphase; iphase++){
        cudaFree(den_dev_arr[iphase]);
        cudaFree(y_dash_dev_arr[iphase]);
    }
    cudaFree(den_dev_arr);
    cudaFree(y_dash_dev_arr);
    cudaFree(mval_dev_arr);
    cudaFree(nval_dev_arr);
}

void SrtlibRlCrabSmthPfZalCuda::RichlucyCrabSmthPfZal(
    FILE* const fp_log,
    const double* const sky_init_arr,
    const double* const flux_init_arr,
    const double* const* const data_arr,
    const double* const bg_arr,
    const double* const flux_target_arr,
    const double* const phase_arr,
    const double* const det_0_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nskyx, int nskyy, int nphase,
    double mu, double gamma,
    string outdir,
    string outfile_head,
    int nem, double tol_em,
    double* const sky_new_arr,
    double* const flux_new_arr)
{
    int nsky = nskyx * nskyy;

    double* sky_init_dev_arr = NULL;
    double* flux_init_dev_arr = NULL;
    double** data_dev_arr = new double*[nphase];
    double* bg_dev_arr = NULL;
    double* flux_target_dev_arr = NULL;
    double* phase_dev_arr = NULL;
    double* det_0_dev_arr = NULL;
    double* resp_norm_mat_dev_arr = NULL;
    double* sky_new_dev_arr = NULL;
    double* flux_new_dev_arr = NULL;
    
    size_t mem_size_nsky = nsky * sizeof(double);
    size_t mem_size_ndet = ndet * sizeof(double);
    size_t mem_size_nsky_ndet = nsky * ndet * sizeof(double);
    size_t mem_size_nphase = nphase * sizeof(double);
    cudaMalloc((void **)&sky_init_dev_arr, mem_size_nsky);
    cudaMalloc((void **)&flux_init_dev_arr, mem_size_nphase);
    for(int iphase = 0; iphase < nphase; iphase++){
        cudaMalloc((void **)&data_dev_arr[iphase], mem_size_ndet);
    }
    cudaMalloc((void **)&bg_dev_arr, mem_size_ndet);
    cudaMalloc((void **)&flux_target_dev_arr, mem_size_nphase);
    cudaMalloc((void **)&phase_dev_arr, mem_size_nphase);
    cudaMalloc((void **)&det_0_dev_arr, mem_size_ndet);
    cudaMalloc((void **)&resp_norm_mat_dev_arr, mem_size_nsky_ndet);
    cudaMalloc((void **)&sky_new_dev_arr, mem_size_nsky);
    cudaMalloc((void **)&flux_new_dev_arr, mem_size_nphase);

    cublasSetVector(nsky, sizeof(double), sky_init_arr, 1,
		    sky_init_dev_arr, 1);
    cublasSetVector(nphase, sizeof(double), flux_init_arr, 1,
		    flux_init_dev_arr, 1);
    for(int iphase = 0; iphase < nphase; iphase++){
        cublasSetVector(ndet, sizeof(double), data_arr[iphase], 1,
                        data_dev_arr[iphase], 1);
    }
    cublasSetVector(ndet, sizeof(double), bg_arr, 1,
		    bg_dev_arr, 1);
    cublasSetVector(nphase, sizeof(double), flux_target_arr, 1,
		    flux_target_dev_arr, 1);
    cublasSetVector(nphase, sizeof(double), phase_arr, 1,
		    phase_dev_arr, 1);
    cublasSetVector(ndet, sizeof(double), det_0_arr, 1,
		    det_0_dev_arr, 1);
    cublasSetVector(ndet * nsky, sizeof(double),
                    resp_norm_mat_arr, 1,
		    resp_norm_mat_dev_arr, 1);
    
    cublasHandle_t handle;
    cublasCreate(&handle);

    double* sky_pre_dev_arr = NULL;
    double* flux_pre_dev_arr = NULL;
    cudaMalloc((void **)&sky_pre_dev_arr, mem_size_nsky);
    cudaMalloc((void **)&flux_pre_dev_arr, mem_size_nphase);
    cublasDcopy(handle, nsky, sky_init_dev_arr, 1,
                sky_pre_dev_arr, 1);
    cublasDcopy(handle, nphase, flux_init_dev_arr, 1,
                flux_pre_dev_arr, 1);
    
    for(int iem = 0; iem < nem; iem ++){
        SrtlibRlCrabSmthPfZalCuda::GetSkyFluxNewArr(
            handle,
	    sky_pre_dev_arr,
            flux_pre_dev_arr,
            data_dev_arr,
            bg_dev_arr,
            flux_target_dev_arr,
            phase_dev_arr,
            det_0_dev_arr,
            resp_norm_mat_dev_arr,    
            ndet, nskyx, nskyy, nphase,
            mu, gamma,
            sky_new_dev_arr,
            flux_new_dev_arr);
        double helldist
            = SrtlibCrabRlCrabStatvalCuda::GetHellingerDist(
                handle,
                sky_pre_dev_arr, flux_pre_dev_arr,
                sky_new_dev_arr, flux_new_dev_arr,
                nsky, nphase);
        if (access( "/tmp/rl_crab_smth_pf_zal_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_crab_smth_pf_zal_stop file is found,"
                "then stop.\n");
            break;
        }
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n",
                             iem, helldist);
            break;
         }
        cublasDcopy(handle, nsky, sky_new_dev_arr, 1,
                    sky_pre_dev_arr, 1);
        cublasDcopy(handle, nphase, flux_new_dev_arr, 1,
                    flux_pre_dev_arr, 1);
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n",
                         iem, helldist);
    }
    cudaFree(sky_pre_dev_arr);
    cudaFree(flux_pre_dev_arr);

    int status = cublasDestroy(handle);

    cublasGetVector(nsky, sizeof(double), sky_new_dev_arr, 1,
                    sky_new_arr, 1);
    cublasGetVector(nphase, sizeof(double), flux_new_dev_arr, 1,
                    flux_new_arr, 1);
    cudaFree(sky_init_dev_arr);
    cudaFree(flux_init_dev_arr);
    for(int iphase = 0; iphase < nphase; iphase++){
        cudaFree(data_dev_arr);
    }
    cudaFree(bg_dev_arr);
    cudaFree(flux_target_dev_arr);
    cudaFree(phase_dev_arr);
    cudaFree(det_0_dev_arr);
    cudaFree(resp_norm_mat_dev_arr);
    cudaFree(sky_new_dev_arr);
    cudaFree(flux_new_dev_arr);
}

//// Accerelated by Zhou-Alexander-Lange,
//// which is H.Zhou, D.Alexander, K.Lange,
//// "A quasi-Newton acceleration for high-dimensional
//// optimization algorithms", Stat Comput (2011) 21, 261.
//// Case: q = 1
//void SrtlibRlCrabSmthPfZal::RichlucyCrabSmthPfZalQ1(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const flux_init_arr,
//    const double* const* const data_arr,
//    const double* const bg_arr,    
//    const double* const flux_target_arr,
//    const double* const phase_arr,
//    const double* const det_0_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nskyx, int nskyy, int nphase,
//    double mu, double gamma,
//    string outdir,
//    string outfile_head,
//    int nem, double tol_em,
//    double* const sky_new_arr,
//    double* const flux_new_arr)
//{
//    int nsky = nskyx * nskyy;
//    double* sky_0_arr  = new double[nsky];
//    double* sky_1_arr  = new double[nsky];
//    double* sky_2_arr  = new double[nsky];
//    double* u_sky_arr  = new double[nsky];
//    double* v_sky_arr  = new double[nsky];
//    double* diff_sky_arr  = new double[nsky];
//    double* sky_0_new_arr  = new double[nsky];
//
//    double* flux_0_arr  = new double[nphase];
//    double* flux_1_arr  = new double[nphase];
//    double* flux_2_arr  = new double[nphase];
//    double* u_flux_arr  = new double[nphase];
//    double* v_flux_arr  = new double[nphase];
//    double* diff_flux_arr  = new double[nphase];
//    double* flux_0_new_arr  = new double[nphase];
//    
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_0_arr, 1);
//    dcopy_(nphase, const_cast<double*>(flux_init_arr), 1, flux_0_arr, 1);
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_0_arr,
//            flux_0_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_1_arr,
//            flux_1_arr);
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_1_arr,
//            flux_1_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_2_arr,
//            flux_2_arr);
//
//        MibBlas::Sub(sky_1_arr, sky_0_arr, nsky, u_sky_arr);
//        MibBlas::Sub(flux_1_arr, flux_0_arr, nphase, u_flux_arr);        
//        MibBlas::Sub(sky_2_arr, sky_1_arr, nsky, v_sky_arr);
//        MibBlas::Sub(flux_2_arr, flux_1_arr, nphase, v_flux_arr);        
//        MibBlas::Sub(v_sky_arr, u_sky_arr, nsky, diff_sky_arr);
//        MibBlas::Sub(v_flux_arr, u_flux_arr, nphase, diff_flux_arr);
//
//        double num = ddot_(nsky, u_sky_arr, 1, u_sky_arr, 1)
//            + ddot_(nphase, u_flux_arr, 1, u_flux_arr, 1);
//        double den = ddot_(nsky, u_sky_arr, 1, diff_sky_arr, 1)
//            + ddot_(nphase, u_flux_arr, 1, diff_flux_arr, 1);
//        double cval = -1.0 * num / den;
//        MiIolib::Printf2(fp_log, "cval = %e\n", cval);
//
//        if (cval < 0.0){
//            // usual update
//            MiIolib::Printf2(fp_log, "cval < 0, then update by sky1 flux1\n");
//            dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//        } else {
//            int nk = 100;
//            double eta = 0.5;
//            int flag_find = 0;
//            for (int ik = 0; ik < nk; ik ++){
//                double cval0 = cval * pow(eta, ik);
//                dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//                dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//                dscal_(nsky, (1.0 - cval0), sky_0_new_arr, 1);
//                dscal_(nphase, (1.0 - cval0), flux_0_new_arr, 1);
//                daxpy_(nsky, cval0, sky_2_arr, 1, sky_0_new_arr, 1);
//                daxpy_(nphase, cval0, flux_2_arr, 1, flux_0_new_arr, 1);
//                
//                int nneg_tmp = 0;
//                for(int isky = 0; isky < nsky; isky ++){
//                    if(sky_0_new_arr[isky] < 0.0){
//                        nneg_tmp ++;
//                    }
//                }
//                for(int iphase = 0; iphase < nphase; iphase ++){
//                    if(flux_0_new_arr[iphase] < 0.0){
//                        nneg_tmp ++;
//                    }
//                }
//                if (nneg_tmp > 0){
//                    continue;
//                } else{
//                    MiIolib::Printf2(fp_log, "cval0 = %e, ik = %d\n",
//                                     cval0, ik);
//                    flag_find = 1;
//                    break;
//                }
//            }
//            if(flag_find == 0){
//                // usual update
//                MiIolib::Printf2(fp_log, "flag_find = 0, then update by sky2 flux2\n");
//                dcopy_(nsky, sky_2_arr, 1, sky_0_new_arr, 1);
//                dcopy_(nphase, flux_2_arr, 1, flux_0_new_arr, 1);
//            }
//        }
//        double helldist  = SrtlibCrabRlCrabStatval::GetHellingerDist(
//            sky_0_arr, flux_0_arr,
//            sky_0_new_arr, flux_0_new_arr,
//            nsky, nphase);
//        MiIolib::Printf2(
//            fp_log, "iem = %d, helldist = %e\n",
//            iem, helldist);
//        if (helldist < tol_em){
//            MiIolib::Printf2(
//                fp_log, "iem = %d, helldist = %e\n",
//                iem, helldist);
//            break;
//        }
//        dcopy_(nsky, sky_0_new_arr, 1, sky_0_arr, 1);
//        dcopy_(nphase, flux_0_new_arr, 1, flux_0_arr, 1);
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//    }
//    dcopy_(nsky, sky_0_new_arr, 1, sky_new_arr, 1);
//    dcopy_(nphase, flux_0_new_arr, 1, flux_new_arr, 1);
//
//    delete [] sky_0_arr;
//    delete [] sky_1_arr;
//    delete [] sky_2_arr;
//    delete [] u_sky_arr;
//    delete [] v_sky_arr;
//    delete [] diff_sky_arr;
//    delete [] sky_0_new_arr;
//
//    delete [] flux_0_arr;
//    delete [] flux_1_arr;
//    delete [] flux_2_arr;
//    delete [] u_flux_arr;
//    delete [] v_flux_arr;
//    delete [] diff_flux_arr;
//    delete [] flux_0_new_arr;
//}
//
//// Accerelated by Squarem S3
//void SrtlibRlCrabSmthPfZal::RichlucyCrabSmthPfSqS3(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const flux_init_arr,
//    const double* const* const data_arr,
//    const double* const bg_arr,    
//    const double* const flux_target_arr,
//    const double* const phase_arr,
//    const double* const det_0_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nskyx, int nskyy, int nphase,
//    double mu, double gamma,
//    string outdir,
//    string outfile_head,
//    int nem, double tol_em,
//    double* const sky_new_arr,
//    double* const flux_new_arr)
//{
//    int nsky = nskyx * nskyy;
//    double* sky_0_arr  = new double[nsky];
//    double* sky_1_arr  = new double[nsky];
//    double* sky_2_arr  = new double[nsky];
//    double* u_sky_arr  = new double[nsky];
//    double* v_sky_arr  = new double[nsky];
//    double* diff_sky_arr  = new double[nsky];
//    double* sky_0_new_arr  = new double[nsky];
//
//    double* flux_0_arr  = new double[nphase];
//    double* flux_1_arr  = new double[nphase];
//    double* flux_2_arr  = new double[nphase];
//    double* u_flux_arr  = new double[nphase];
//    double* v_flux_arr  = new double[nphase];
//    double* diff_flux_arr  = new double[nphase];
//    double* flux_0_new_arr  = new double[nphase];
//    
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_0_arr, 1);
//    dcopy_(nphase, const_cast<double*>(flux_init_arr), 1, flux_0_arr, 1);
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_0_arr,
//            flux_0_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_1_arr,
//            flux_1_arr);
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_1_arr,
//            flux_1_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_2_arr,
//            flux_2_arr);
//
//        MibBlas::Sub(sky_1_arr, sky_0_arr, nsky, u_sky_arr);
//        MibBlas::Sub(flux_1_arr, flux_0_arr, nphase, u_flux_arr);        
//        MibBlas::Sub(sky_2_arr, sky_1_arr, nsky, v_sky_arr);
//        MibBlas::Sub(flux_2_arr, flux_1_arr, nphase, v_flux_arr);        
//        MibBlas::Sub(v_sky_arr, u_sky_arr, nsky, diff_sky_arr);
//        MibBlas::Sub(v_flux_arr, u_flux_arr, nphase, diff_flux_arr);
//
//        double num = ddot_(nsky, u_sky_arr, 1, u_sky_arr, 1)
//            + ddot_(nphase, u_flux_arr, 1, u_flux_arr, 1);
//        double den = ddot_(nsky, diff_sky_arr, 1, diff_sky_arr, 1)
//            + ddot_(nphase, diff_flux_arr, 1, diff_flux_arr, 1);
//        double sval = -1.0 * sqrt(num / den);
//        MiIolib::Printf2(fp_log, "sval = %e\n", sval);
//
//        int nk = 100;
//        double eta = 0.5;
//        int flag_find = 0;
//        for (int ik = 0; ik < nk; ik ++){
//            double sval0 = sval * pow(eta, ik);
//            dcopy_(nsky, sky_0_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_0_arr, 1, flux_0_new_arr, 1);
//            daxpy_(nsky, -2.0 * sval0, u_sky_arr, 1, sky_0_new_arr, 1);
//            daxpy_(nphase, -2.0 * sval0, u_flux_arr, 1, flux_0_new_arr, 1);
//            daxpy_(nsky, sval0 * sval0, diff_sky_arr, 1, sky_0_new_arr, 1);
//            daxpy_(nphase, sval0 * sval0, diff_flux_arr, 1, flux_0_new_arr, 1);
//
//            int nneg_tmp = 0;
//            for(int isky = 0; isky < nsky; isky ++){
//                if(sky_0_new_arr[isky] < 0.0){
//                    nneg_tmp ++;
//                }
//            }
//            for(int iphase = 0; iphase < nphase; iphase ++){
//                if(flux_0_new_arr[iphase] < 0.0){
//                    nneg_tmp ++;
//                }
//            }
//            if (nneg_tmp > 0){
//                continue;
//            } else{
//                MiIolib::Printf2(fp_log, "sval0 = %e, ik = %d\n",
//                                 sval0, ik);
//                flag_find = 1;
//                break;
//            }
//        }
//        if(flag_find == 0){
//            // usual update
//            MiIolib::Printf2(fp_log,
//                             "flag_find = 0, then update by sky2 flux2\n");
//            dcopy_(nsky, sky_2_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_2_arr, 1, flux_0_new_arr, 1);
//        }
//        double helldist  = SrtlibCrabRlCrabStatval::GetHellingerDist(
//            sky_0_arr, flux_0_arr,
//            sky_0_new_arr, flux_0_new_arr,
//            nsky, nphase);
//        MiIolib::Printf2(
//            fp_log, "iem = %d, helldist = %e\n",
//            iem, helldist);
//        if (helldist < tol_em){
//            MiIolib::Printf2(
//                fp_log, "iem = %d, helldist = %e\n",
//                iem, helldist);
//            break;
//        }
//        dcopy_(nsky, sky_0_new_arr, 1, sky_0_arr, 1);
//        dcopy_(nphase, flux_0_new_arr, 1, flux_0_arr, 1);
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//    }
//    dcopy_(nsky, sky_0_new_arr, 1, sky_new_arr, 1);
//    dcopy_(nphase, flux_0_new_arr, 1, flux_new_arr, 1);
//
//    delete [] sky_0_arr;
//    delete [] sky_1_arr;
//    delete [] sky_2_arr;
//    delete [] u_sky_arr;
//    delete [] v_sky_arr;
//    delete [] diff_sky_arr;
//    delete [] sky_0_new_arr;
//
//    delete [] flux_0_arr;
//    delete [] flux_1_arr;
//    delete [] flux_2_arr;
//    delete [] u_flux_arr;
//    delete [] v_flux_arr;
//    delete [] diff_flux_arr;
//    delete [] flux_0_new_arr;
//}

