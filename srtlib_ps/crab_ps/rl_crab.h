#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"
#include "../rl.h"

namespace SrtlibRlCrab
{
    void GetRhoNuNewNumArr(const double* const rho_arr,
                           const double* const nu_arr,
                           const double* const* const data_arr,
                           const double* const phase_arr,
                           const double* const det_0_arr,
                           const double* const resp_norm_mat_arr,
                           int ndet, int nsky, int nphase,
                           double* const rho_new_arr,
                           double* const nu_new_arr);
    
    void GetRhoNuNewArr(const double* const rho_arr,
                        const double* const nu_arr,
                        const double* const* const data_arr,
                        const double* const phase_arr,
                        const double* const det_0_arr,
                        const double* const resp_norm_mat_arr,
                        int ndet, int nsky, int nphase,
                        double* const rho_new_arr,
                        double* const nu_new_arr);

    void RichlucyCrab(FILE* const fp_log,
                      const double* const rho_init_arr,
                      const double* const nu_init_arr,
                      const double* const* const data_arr,
                      const double* const phase_arr,
                      const double* const det_0_arr,
                      const double* const resp_norm_mat_arr,
                      int ndet, int nsky, int nphase,
                      string outdir, string outfile_head,
                      int nem, double tol_em,
                      double* const rho_new_arr,
                      double* const nu_new_arr);

    void RichlucyCrabAccZALq1(
        FILE* const fp_log,
        const double* const rho_init_arr,
        const double* const nu_init_arr,
        const double* const* const data_arr,
        const double* const phase_arr,
        const double* const det_0_arr,
        const double* const resp_norm_mat_arr,
        int ndet, int nsky, int nphase,
        string outdir, string outfile_head,
        int nem, double tol_em,
        double* const rho_new_arr,
        double* const nu_new_arr);
    

} // namespace SrtlibRl

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_H_


