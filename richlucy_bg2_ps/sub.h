#ifndef MORIIISM_SRT_RICHLUCY_BG_SUB_BG_H_
#define MORIIISM_SRT_RICHLUCY_BG_SUB_BG_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

void GetDetArr(const double* const rho_arr,
               const double* const resp_norm_mat_arr,
               int ndet, int nsky,
               double* const out_arr);

void GetRhoNu_New(const double* const rho_arr, double nu,
                  const double* const data_arr,
                  const double* const resp_norm_mat_arr,
                  const double* const bg_arr,
                  int ndet, int nsky,
                  double* const rho_new_arr,
                  double* const nu_new_ptr);

void RichlucyBg(const double* const rho_init_arr,
                double nu_init,
                const double* const data_arr,
                const double* const bg_arr,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky,
                int niter,
                string outdir,
                string outfile_head,
                double tol,
                double* const rho_new_arr,
                double* const nu_new_ptr);

double GetHellingerDist(const double* const rho_arr, double nu, 
                        const double* const rho_new_arr, double nu_new,
                        int nsky);

double GetFuncL(const double* const data_arr,
                const double* const bg_arr,
                const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky);

#endif // MORIIISM_SRT_RICHLUCY_BG_SUB_BG_H_
