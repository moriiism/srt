#ifndef MORIIISM_SRT_SRTLIB_RL_H_
#define MORIIISM_SRT_SRTLIB_RL_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

// accerelated richardson lucy (Ikeda)
void RichlucyAccIkeda(FILE* const fp_log,
                      const double* const rho_init_arr,
                      const double* const data_arr,
                      const double* const resp_norm_mat_arr,
                      int ndet, int nsky,
                      string outdir, string outfile_head,
                      int nem, double tol_em,
                      int nph_data, int rand_seed,
                      double* const rho_new_arr);


// accerelated richardson lucy (SQUAREM)
void RichlucyAccSQUAREM(FILE* const fp_log,
                        const double* const rho_init_arr,
                        const double* const data_arr,
                        const double* const resp_norm_mat_arr,
                        int ndet, int nsky,
                        string outdir, string outfile_head,
                        int nem, double tol_em,
                        double* const rho_new_arr);

void GetInvVec(const double* const vec_arr, int nelm,
               double* const inv_arr);

// accerelated richardson lucy
void RichlucyAcc(FILE* const fp_log,
                 const double* const rho_init_arr,
                 const double* const data_arr,
                 const double* const resp_norm_mat_arr,
                 int ndet, int nsky,
                 string outdir, string outfile_head,
                 int nem, double tol_em,
                 int k_restart, double delta_restart,
                 double* const rho_new_arr);

void Richlucy(FILE* const fp_log,
              const double* const rho_init_arr,
              const double* const data_arr,
              const double* const resp_norm_mat_arr,
              int ndet, int nsky,
              string outdir, string outfile_head,
              int nem, double tol_em,
              double* const rho_new_arr);

void GetRhoNewArr(const double* const rho_arr,
                  const double* const data_arr,
                  const double* const resp_norm_mat_arr,
                  int ndet, int nsky,
                  double* const rho_new_arr);

void GetDetArr(const double* const rho_arr,
               const double* const resp_norm_mat_arr,
               int ndet, int nsky,
               double* const det_arr);

#endif // MORIIISM_SRT_SRTLIB_RL_H_
