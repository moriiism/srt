#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PMPROJ_PM_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PMPROJ_PM_H_

#include "mib_blas.h"
#include "mif_fits.h"
//#include "mir_math.h"

void GetDerivFpm(const double* const rho_arr,
                 double phi,
                 int nskyx, int nskyy, int nsrc,
                 double mu, double B_val,
                 double* const deriv_fpm_sky_arr,
                 double* const deriv_fpm_src_arr,
                 double* const deriv_fpm_phi_ptr);

void GetProjectedDerivFpm(const double* const deriv_fpm_sky_arr,
                          const double* const deriv_fpm_src_arr,
                          double deriv_fpm_phi,
                          int nsky, int nsrc,
                          double* const proj_deriv_fpm_sky_arr,
                          double* const proj_deriv_fpm_src_arr,
                          double* const proj_deriv_fpm_phi_ptr);

void GetVvalArr(const double* const rho_arr,
                const double* const proj_deriv_fpm_rho_arr,
                int nsky, double lip_const,
                double* const vval_arr);

void GetWvalArr(const double* const nu_arr,
                const double* const proj_deriv_fpm_nu_arr,
                int nsrc, double lip_const,
                double* const wval_arr);

double GetZval(double phi,
               double proj_deriv_fpm_phi,
               double lip_const);

void GetVvalArrWvalArrZval(const double* const rho_arr,
                           const double* const nu_arr,
                           double phi,
                           int nskyx, int nskyy, int nsrc,
                           double mu, double B_val,
                           double lip_const,
                           double* const vval_arr,
                           double* const wval_arr,
                           double* const zval_ptr);

double GetFindLipConst(const double* const rho_arr,
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

void GetRhoNuPhi_ByPM(const double* const rho_arr,
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

#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PMPROJ_PM_H_
