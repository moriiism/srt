#ifndef MORIIISM_SRT_SRTLIB_RL_DET2_H_
#define MORIIISM_SRT_SRTLIB_RL_DET2_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"

namespace SrtlibRlDet2
{
    void GetSkyNewFactorArr(
        const double* const sky_arr,
        const double* const data_arr,
        const double* const bg_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nsky,
        double* const sky_new_factor_arr);

    void GetSkyNewArr(
        const double* const sky_arr,
        const double* const data1_arr,
        const double* const data2_arr,
        const double* const bg1_arr,
        const double* const bg2_arr,
        const double* const resp_norm_mat_det1_arr,
        const double* const resp_norm_mat_det2_arr,
        int ndet, int nsky,
        double* const sky_new_arr);

    void Richlucy(
        FILE* const fp_log,
        const double* const sky_init_arr,
        const double* const data1_arr,
        const double* const data2_arr,
        const double* const bg1_arr,
        const double* const bg2_arr,
        const double* const resp_norm_mat_det1_arr,
        const double* const resp_norm_mat_det2_arr,    
        int ndet, int nsky,
        string outdir, string outfile_head,
        int nem, double tol_em,
        double* const sky_new_arr);

} // namespace SrtlibRlDet2

#endif // MORIIISM_SRT_SRTLIB_RL_DET2_H_
