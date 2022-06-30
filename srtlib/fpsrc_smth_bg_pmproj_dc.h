#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PMPROJ_DC_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PMPROJ_DC_H_

#include "mib_blas.h"
#include "mif_fits.h"
//#include "mir_math.h"

void GetRhoNuPhi_ByDC_Pmproj(FILE* const fp_log,
                             const double* const rho_arr,
                             const double* const nu_arr,
                             double phi,
                             const double* const mval_arr,
                             const double* const nval_arr,
                             double pval,
                             int nph, double B_val,
                             int ndet, int nskyx, int nskyy, int nsrc,
                             double mu,
                             int ndc, double tol_dc,
                             int npm, double tol_pm,
                             int nnewton, double tol_newton,
                             double* const rho_new_arr,
                             double* const nu_new_arr,
                             double* const phi_new_ptr,
                             double* const helldist_ptr,
                             int* const flag_converge_ptr);

#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_PMPROJ_DC_H_
