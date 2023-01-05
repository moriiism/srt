#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PM_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PM_H_

#include "mib_blas.h"
#include "mif_fits.h"
//#include "mir_math.h"

void GetVvalArr(const double* const rho_arr,
                int nskyx, int nskyy,
                double mu, double lip_const,
                double* const vval_arr);

void GetWvalArr(const double* const nu_arr,
                int nsrc,
                double* const wval_arr);

double GetZval(double phi,
               double lip_const,
               double B_val);

double GetFindLipConst(FILE* const fp_log,
                       const double* const rho_arr,
                       const double* const nu_arr,
                       double phi,
                       const double* const mval_arr,
                       const double* const nval_arr,
                       double pval,
                       double phi_val,
                       int nph, double B_val,
                       double mu,
                       int nskyx, int nskyy, int nsrc,
                       double lip_const, double lambda,
                       int nnewton, double tol_newton);

double GetQMinusF(const double* const rho_new_arr,
                  const double* const nu_new_arr,
                  double phi_new,
                  const double* const rho_arr,
                  const double* const nu_arr,
                  double phi,
                  double mu, double lip_const,
                  double B_val,
                  int nskyx, int nskyy, int nsrc);

double GetFuncF(const double* const rho_arr, double phi,
                double mu, double B_val, 
                int nskyx, int nskyy);

void GetDiffF(const double* const rho_arr,
              double mu,
              int nskyx, int nskyy,
              double* const out_arr);


void GetRhoNuPhi_ByPM(FILE* const fp_log,
                      const double* const rho_arr,
                      const double* const nu_arr,
                      double phi,
                      const double* const mval_arr,
                      const double* const nval_arr,
                      double pval,
                      int nph, double B_val,
                      int ndet, int nskyx, int nskyy, int nsrc,
                      double mu,
                      int npm, double tol_pm,
                      int nnewton, double tol_newton,
                      double* const rho_new_arr,
                      double* const nu_new_arr,
                      double* const phi_new_ptr);

// nesterov
void GetRhoNuPhi_ByPM_Nesterov(FILE* const fp_log,
                               const double* const rho_arr,
                               const double* const nu_arr,
                               double phi,
                               const double* const mval_arr,
                               const double* const nval_arr,
                               double pval,
                               int nph, double B_val,
                               int ndet, int nskyx, int nskyy, int nsrc,
                               double mu,
                               int npm, double tol_pm,
                               int nnewton, double tol_newton,
                               double* const rho_new_arr,
                               double* const nu_new_arr,
                               double* const phi_new_ptr,
                               double* const helldist_ptr,
                               int* const flag_converge_ptr);

#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PM_H_
