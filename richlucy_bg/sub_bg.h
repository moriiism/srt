#ifndef MORIIISM_SRT_RICHLUCY_BG_SUB_BG_H_
#define MORIIISM_SRT_RICHLUCY_BG_SUB_BG_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

double GetNextNb(const double* const rho_arr,
                 const double* const data_arr,
                 const double* const resp_mat_arr,
                 const double* const bg_arr,
                 int ndet, int nsky,
                 double N_B,
                 int niter_newton,
                 double tol_newton);

double GetB(const double* const bg_arr, int ndet);

double GetN(const double* const rho_arr, int nsky);

double GetDerivF_NB(const double* const rho_arr,
                    const double* const data_arr,
                    const double* const resp_mat_arr,
                    const double* const bg_arr,
                    int ndet, int nsky,
                    double N_B);

double GetDeriv2F_NB(const double* const rho_arr,
                     const double* const data_arr,
                     const double* const resp_mat_arr,
                     const double* const bg_arr,
                     int ndet, int nsky,
                     double N_B);

void GetDetArr(const double* const rho_arr,
               const double* const resp_mat_arr,
               int ndet, int nsky,
               double* const out_arr);

void GetNextRhoArr(const double* const rho_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   const double* const bg_arr,
                   int ndet, int nsky,
                   double N_B,
                   double* const out_arr);

void RichlucyBg(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                const double* const bg_arr,
                int niter_main, int niter_em, int niter_newton,
                string outdir, string outfile_head,
                int ndet, int nskyx, int nskyy,
                double tol_main, double tol_em, double tol_newton,
                double* const out_arr, double* const N_B_ptr);

void LoadResp(string respdir, int nskyx, int nskyy,
              double epsilon,
              double** const mat_arr_ptr,
              int* const ndetx_ptr,
              int* const ndety_ptr);

void GetNdet(string respdir, int* const ndetx_ptr, int* const ndety_ptr);

double GetHellingerDist(const double* const rho_arr,
                        const double* const rho_new_arr,
                        int nsky);

#endif // MORIIISM_SRT_RICHLUCY_BG_SUB_BG_H_
