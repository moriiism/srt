#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_EM_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_EM_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

namespace SrtlibRlCrabSmthPfEm
{
    void RichlucyCrabSmthPf(
        FILE* const fp_log,
        const double* const rho_init_arr,
        const double* const nu_init_arr,
        const double* const* const data_arr,
        const double* const nu_0_arr,        
        const double* const phase_arr,
        const double* const det_0_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nskyx, int nskyy, int nphase,
        double mu, double gamma,
        string outdir,
        string outfile_head,
        int nem, double tol_em,
        int npm, double tol_pm,
        int nnewton, double tol_newton,
        double* const rho_new_arr,
        double* const nu_new_arr);

} // namespace SrtlibRlCrabSmthPfEm

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_EM_H_
