#ifndef MORIIISM_SRT_SRTLIB_RL_H_
#define MORIIISM_SRT_SRTLIB_RL_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"

namespace SrtlibRl
{
    void GetDetArr(const double* const sky_arr,
                   const double* const resp_norm_mat_arr,
                   int ndet, int nsky,
                   double* const det_arr);

    void GetDenArr(const double* const sky_arr,
                   const double* const bg_arr,
                   const double* const resp_norm_mat_arr,
                   int ndet, int nsky,
                   double* const den_arr);

    void GetYDashArr(const double* const data_arr,
                     const double* const den_arr,
                     int ndet,
                     double* const y_dash_arr);

    void GetSkyNewArr(const double* const sky_arr,
                      const double* const data_arr,
                      const double* const bg_arr,
                      const double* const resp_norm_mat_arr,
                      int ndet, int nsky,
                      double* const sky_new_arr);

    void Richlucy(
        FILE* const fp_log,
        const double* const sky_init_arr,
        const double* const data_arr,
        const double* const bg_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nsky,
        string outdir, string outfile_head,
        int nem, double tol_em,
        double* const sky_new_arr);

    // accerelated richardson lucy (SQUAREM)
    void RichlucyAccSquarem(
        FILE* const fp_log,
        const double* const sky_init_arr,
        const double* const data_arr,
        const double* const bg_arr,        
        const double* const resp_norm_mat_arr,
        int ndet, int nsky,
        string outdir, string outfile_head,
        int nem, double tol_em,
        double* const sky_new_arr);
    
    // accerelated richardson lucy (Ikeda)
    void RichlucyAccIkeda(
        FILE* const fp_log,
        const double* const sky_init_arr,
        const double* const data_arr,
        const double* const bg_arr,        
        const double* const resp_norm_mat_arr,
        int ndet, int nsky,
        string outdir, string outfile_head,
        int nem, double tol_em,
        int nph_data, int rand_seed,
        double* const sky_new_arr);

    void GetInvVec(const double* const vec_arr, int nelm,
                   double* const inv_arr);

    // accerelated richardson lucy
    void RichlucyAccKuroda(
        FILE* const fp_log,
        const double* const sky_init_arr,
        const double* const data_arr,
        const double* const bg_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nsky,
        string outdir, string outfile_head,
        int nem, double tol_em,
        int k_restart, double delta_restart,
        double* const sky_new_arr);

} // namespace SrtlibRl

#endif // MORIIISM_SRT_SRTLIB_RL_H_
