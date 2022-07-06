#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_MM_PM_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_MM_PM_H_

#include "mib_blas.h"
#include "mif_fits.h"

namespace FpsrcSmthBgMmPm
{
    double GetFindLipConst_Mm(FILE* const fp_log,
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

    void GetRhoNuPhi_ByPm_Mm(FILE* const fp_log,
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

    // nesterov
    void GetRhoNuPhi_ByPm_Mm_Nesterov(FILE* const fp_log,
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

} // namespace FpsrcSmthBgMmPm
    
#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_MM_PM_H_
