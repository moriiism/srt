#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_ZAL_CUDA_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_ZAL_CUDA_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
//#include "mir_math.h"

namespace SrtlibRlCrabSmthPfZalCuda
{
    __global__
    void GetSkyNewCuda(
        const double* const alpha_dev_arr,
        const double* const beta_dev_arr,
        const double* const mval_dev_arr,
        double live_time_ratio_ave,
        int nsky, double mu,
        double* const sky_new_dev_arr);

    void GetSkyNewArr(
        cublasHandle_t handle,
        const double* const sky_dev_arr,
        const double* const mval_dev_arr,
        double live_time_ratio_ave,
        int nskyx, int nskyy, double mu,
        double* const sky_new_dev_arr);
    
    __global__
    void GetFluxNewArr(
        const double* const nval_dev_arr,
        const double* const flux_target_dev_arr,
        const double* const phase_dev_arr,
        const double* const live_time_ratio_dev_arr,
        int nphase, double gamma,
        double* const flux_new_dev_arr);

    void GetSkyFluxNewArr(
        cublasHandle_t handle,
        const double* const sky_pre_dev_arr,
        const double* const flux_pre_dev_arr,
        const double* const* const data_dev_arr,
        const double* const bg_dev_arr,
        const double* const flux_target_dev_arr,
        const double* const phase_dev_arr,
        const double* const live_time_ratio_dev_arr,
        const double* const det_0_dev_arr,
        const double* const resp_norm_mat_dev_arr,    
        int ndet, int nskyx, int nskyy, int nphase,
        double mu, double gamma,
        double* const sky_new_dev_arr,
        double* const flux_new_dev_arr);
    
    void RichlucyCrabSmthPfZal(
        FILE* const fp_log,
        const double* const sky_init_arr,
        const double* const flux_init_arr,
        const double* const* const data_arr,
        const double* const bg_arr,
        const double* const flux_target_arr,
        const double* const phase_arr,
        const double* const live_time_ratio_arr,
        const double* const det_0_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nskyx, int nskyy, int nphase,
        double mu, double gamma,
        string outdir,
        string outfile_head,
        int nem, double tol_em,
        double* const sky_new_arr,
        double* const flux_new_arr);
    

//    void RichlucyCrabSmthPfZalQ1(
//        FILE* const fp_log,
//        const double* const sky_init_arr,
//        const double* const flux_init_arr,
//        const double* const* const data_arr,
//        const double* const bg_arr,    
//        const double* const flux_target_arr,
//        const double* const phase_arr,
//        const double* const det_0_arr,
//        const double* const resp_norm_mat_arr,
//        int ndet, int nskyx, int nskyy, int nphase,
//        double mu, double gamma,
//        string outdir,
//        string outfile_head,
//        int nem, double tol_em,
//        double* const sky_new_arr,
//        double* const flux_new_arr);
//
//    void RichlucyCrabSmthPfSqS3(
//        FILE* const fp_log,
//        const double* const sky_init_arr,
//        const double* const flux_init_arr,
//        const double* const* const data_arr,
//        const double* const bg_arr,    
//        const double* const flux_target_arr,
//        const double* const phase_arr,
//        const double* const det_0_arr,
//        const double* const resp_norm_mat_arr,
//        int ndet, int nskyx, int nskyy, int nphase,
//        double mu, double gamma,
//        string outdir,
//        string outfile_head,
//        int nem, double tol_em,
//        double* const sky_new_arr,
//        double* const flux_new_arr);

} // namespace SrtlibRlCrabSmthPfZalCuda

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_ZAL_CUDA_H_
