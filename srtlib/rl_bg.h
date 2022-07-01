#ifndef MORIIISM_SRT_SRTLIB_RL_BG_H_
#define MORIIISM_SRT_SRTLIB_RL_BG_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

double GetAlpha(const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                const double* const bg_arr,
                const double* const data_arr,
                int nsky, int ndet);

void GetRhoNu_New(const double* const rho_arr, double nu,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const bg_arr,
                  int ndet, int nsky,
                  double* const rho_new_arr,
                  double* const nu_new_ptr);

void RichlucyBg(FILE* const fp_log,
                const double* const rho_init_arr,
                double nu_init,
                const double* const data_arr,
                const double* const bg_arr,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky,
                string outdir,
                string outfile_head,
                int nem,                
                double tol_em,
                double* const rho_new_arr,
                double* const nu_new_ptr);

void RichlucyBgAccSQUAREM(FILE* const fp_log,
                          const double* const rho_init_arr,
                          double nu_init,
                          const double* const data_arr,
                          const double* const bg_arr,
                          const double* const resp_norm_mat_arr,
                          int ndet, int nsky,
                          string outdir,
                          string outfile_head,
                          int nem,                
                          double tol_em,
                          double* const rho_new_arr,
                          double* const nu_new_ptr);

#endif // MORIIISM_SRT_SRTLIB_RL_BG_H_
