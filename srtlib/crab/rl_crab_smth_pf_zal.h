#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_ZAL_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_ZAL_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
//#include "mir_math.h"

namespace SrtlibRlCrabSmthPfZal
{
    void GetSkyNew(
        const double* const alpha_arr,
        const double* const beta_arr,
        const double* const mval_arr,
        double live_time_ratio_ave,
        int nsky, double mu,
        double* const sky_new_arr);
    void GetSkyNewArr(
        const double* const sky_arr,
        const double* const mval_arr,
        double live_time_ratio_ave,
        int nskyx, int nskyy, double mu,
        double* const sky_new_arr);
    void GetFluxNewArr(
        const double* const nval_arr,
        const double* const flux_target_arr,
        const double* const phase_arr,
        const double* const live_time_ratio_arr,        
        int nphase, double gamma,
        double* const flux_new_arr);

    void GetSkyFluxNewArr(
        const double* const sky_pre_arr,
        const double* const flux_pre_arr,
        const double* const* const data_arr,
        const double* const bg_arr,
        const double* const flux_target_arr,
        const double* const phase_arr,
        const double* const live_time_ratio_arr,        
        const double* const det_0_arr,
        const double* const resp_norm_mat_arr,    
        int ndet, int nskyx, int nskyy, int nphase,
        double mu, double gamma,
        double* const sky_new_arr,
        double* const flux_new_arr);
    
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

    void RichlucyCrabSmthPfZalQ1(
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
        double* const flux_new_arr);

    void RichlucyCrabSmthPfSqS3(
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
        double* const flux_new_arr);

} // namespace SrtlibRlCrabSmthPfZal

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_ZAL_H_
