#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_MM_EM_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_MM_EM_H_

#include "mib_blas.h"
#include "mif_fits.h"
//#include "mir_math.h"

void RichlucyFpsrcSmthBgMM(FILE* const fp_log,
                           const double* const rho_init_arr,
                           const double* const nu_init_arr,
                           double phi_init,
                           const double* const data_arr,
                           const double* const bg_arr,
                           const double* const* const det_fpsrc_arr,
                           const double* const resp_norm_mat_arr,
                           int ndet, int nskyx, int nskyy, int nsrc,
                           double mu,
                           string outdir,
                           string outfile_head,
                           int nem, double tol_em,
                           int npm, double tol_pm,
                           int nnewton, double tol_newton,
                           double* const rho_new_arr,
                           double* const nu_new_arr,
                           double* const phi_new_ptr);

// accerelation by SQUAREM
void RichlucyFpsrcSmthBgMM_Acc(FILE* const fp_log,
                               const double* const rho_init_arr,
                               const double* const nu_init_arr,
                               double phi_init,
                               const double* const data_arr,
                               const double* const bg_arr,
                               const double* const* const det_fpsrc_arr,
                               const double* const resp_norm_mat_arr,
                               int ndet, int nskyx, int nskyy, int nsrc,
                               double mu,
                               string outdir,
                               string outfile_head,
                               int nem, double tol_em,
                               int npm, double tol_pm,
                               int nnewton, double tol_newton,
                               double* const rho_new_arr,
                               double* const nu_new_arr,
                               double* const phi_new_ptr);

#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_MM_EM_H_
