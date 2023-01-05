#ifndef MORIIISM_SRT_SRTLIB_SUB_H_
#define MORIIISM_SRT_SRTLIB_SUB_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

void GetDetArr(const double* const rho_arr,
               const double* const resp_norm_mat_arr,
               int ndet, int nsky,
               double* const out_arr);

void GetMArrNval(const double* const rho_arr, double nu,
                 const double* const data_arr,
                 const double* const resp_norm_mat_arr,
                 const double* const bg_arr,
                 int ndet, int nsky,
                 double* const mval_arr,
                 double* const nval_ptr);

double GetFindLipConst(const double* const rho_arr, double nu,
                       const double* const mval_arr, double nval,
                       double mu,
                       int nskyx, int nskyy,
                       double lip_const, double lambda,
                       int nnewton, double tol_newton);

double GetQMinusF(const double* const rho_new_arr, double nu_new,
                  const double* const rho_arr, double nu,
                  double mu, double lip_const,
                  int nskyx, int nskyy);

double GetFuncF(const double* const rho_arr,
                double mu,
                int nskyx, int nskyy);

void GetDiffF(const double* const rho_arr,
              double mu,
              int nskyx, int nskyy,
              double* const out_arr);

void GetRhoNu_New(const double* const rho_arr, double nu,
                  const double* const data_arr,
                  const double* const resp_norm_mat_arr,
                  const double* const bg_arr,
                  int ndet, int nskyx, int nskyy,
                  int npm, double tol_pm,
                  int nnewton, double tol_newton,
                  double* const rho_new_arr,
                  double* const nu_new_ptr);

void RichlucyBg2Smooth(const double* const rho_init_arr,
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

double GetHellingerDist(const double* const rho_arr, double nu, 
                        const double* const rho_new_arr, double nu_new,
                        int nsky);

double GetFuncL(const double* const data_arr,
                const double* const bg_arr,
                const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky);

#endif // MORIIISM_SRT_SRTLIB_SUB_H_
