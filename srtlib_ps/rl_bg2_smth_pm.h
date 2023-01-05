#ifndef MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_PM_H_
#define MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_PM_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

namespace SrtlibRlBg2SmthPm
{
    void GetVvalArr(const double* const rho_arr,
                    int nskyx, int nskyy,
                    double mu, double lip_const,
                    double* const vval_arr);

    double GetWval(double nu);

    double GetFindLipConst(
        FILE* const fp_log,
        const double* const rho_arr, double nu,
        const double* const mval_arr, double nval,
        double mu,
        int nskyx, int nskyy,
        double lip_const, double lambda,
        int nnewton, double tol_newton);

    double GetQMinusF(const double* const rho_new_arr,
                      double nu_new,
                      const double* const rho_arr,
                      double nu,
                      double mu, double lip_const,
                      int nskyx, int nskyy);

    double GetFuncF(const double* const rho_arr,
                    double mu,
                    int nskyx, int nskyy);

    void GetDiffF(const double* const rho_arr,
                  double mu,
                  int nskyx, int nskyy,
                  double* const out_arr);
    
    void GetRhoNu_ByPm(
        FILE* const fp_log,
        const double* const rho_arr, double nu,
        const double* const mval_arr, double nval,
        int nskyx, int nskyy,
        double mu,
        int npm, double tol_pm,
        int nnewton, double tol_newton,
        double* const rho_new_arr,
        double* const nu_new_ptr,
        double* const helldist_ptr,
        int* const flag_converge_ptr);

    void GetRhoNu_ByPm_Nesterov(
        FILE* const fp_log,
        const double* const rho_arr, double nu,
        const double* const mval_arr, double nval,
        int nskyx, int nskyy,
        double mu,
        int npm, double tol_pm,
        int nnewton, double tol_newton,
        double* const rho_new_arr,
        double* const nu_new_ptr,
        double* const helldist_ptr,
        int* const flag_converge_ptr);
    
} // namespace SrtlibRlBg2SmthPm

#endif // MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_PM_H_
