#ifndef MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_EM_H_
#define MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_EM_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

namespace SrtlibRlBg2SmthEm
{
    void RichlucyBg2Smth(
        FILE* const fp_log,
        const double* const rho_init_arr,
        double nu_init,
        const double* const data_arr,
        const double* const bg_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nskyx, int nskyy, double mu,
        string outdir,
        string outfile_head,
        int nem, double tol_em,
        int npm, double tol_pm,
        int nnewton, double tol_newton,
        double* const rho_new_arr,
        double* const nu_new_ptr);

    // accerelation by SQUAREM
    void RichlucyBg2Smth_Acc(
        FILE* const fp_log,
        const double* const rho_init_arr,
        double nu_init,
        const double* const data_arr,
        const double* const bg_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nskyx, int nskyy, double mu,
        string outdir,
        string outfile_head,
        int nem, double tol_em,
        int npm, double tol_pm,
        int nnewton, double tol_newton,
        double* const rho_new_arr,
        double* const nu_new_ptr);
    
    // get m_arr & n_val
    void GetMArrNval(
        const double* const rho_arr, double nu,
        const double* const data_arr,
        const double* const bg_arr,
        const double* const resp_norm_mat_arr,        
        int ndet, int nsky,
        double* const mval_arr,
        double* const nval_ptr);

} // namespace SrtlibRlBg2SmthEm

#endif // MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_EM_H_

