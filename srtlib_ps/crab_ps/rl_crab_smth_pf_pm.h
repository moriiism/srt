#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_PM_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_PM_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

namespace SrtlibRlCrabSmthPfPm
{
    void GetVvalArr(const double* const rho_arr,
                    int nskyx, int nskyy,
                    double mu, double lip_const,
                    double* const vval_arr);

    void GetWvalArr(const double* const nu_arr,
                    const double* const nu_0_arr,
                    int nphase,
                    double gamma, double lip_const,
                    double* const wval_arr);

    double GetFindLipConst(
        FILE* const fp_log,
        const double* const rho_arr,
        const double* const nu_arr,
        const double* const mval_arr,
        const double* const nval_arr,
        const double* const nu_0_arr,
        double mu, double gamma,
        int nskyx, int nskyy, int nphase,
        double lip_const, double lambda,
        int nnewton, double tol_newton);

    double GetQMinusF(const double* const rho_new_arr,
                      const double* const nu_new_arr,
                      const double* const rho_arr,
                      const double* const nu_arr,
                      const double* const nu_0_arr,
                      double mu, double gamma,
                      double lip_const,
                      int nskyx, int nskyy, int nphase);

    double GetTermD(const double* const nu_arr,
                    const double* const nu_0_arr,
                    int nphase);
    void GetDiffTermD(const double* const nu_arr,
                      const double* const nu_0_arr,
                      int nphase,
                      double* const term_d_diff_arr);
    
    void GetRhoNu_ByPm(
        FILE* const fp_log,
        const double* const rho_arr,
        const double* const nu_arr,
        const double* const mval_arr,
        const double* const nval_arr,
        const double* const nu_0_arr,
        int nskyx, int nskyy, int nphase,
        double mu, double gamma,
        int npm, double tol_pm,
        int nnewton, double tol_newton,
        double* const rho_new_arr,
        double* const nu_new_arr,
        double* const helldist_ptr,
        int* const flag_converge_ptr);

//    void GetRhoNu_ByPm_Nesterov(
//        FILE* const fp_log,
//        const double* const rho_arr,
//        const double* const nu_arr,
//        const double* const mval_arr,
//        const double* const nval_arr,
//        int nskyx, int nskyy, int nphase,
//        double mu, double gamma,
//        int npm, double tol_pm,
//        int nnewton, double tol_newton,
//        double* const rho_new_arr,
//        double* const nu_new_arr,
//        double* const helldist_ptr,
//        int* const flag_converge_ptr);
    
} // namespace SrtlibRlCrabSmthPfPm

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_PM_H_
