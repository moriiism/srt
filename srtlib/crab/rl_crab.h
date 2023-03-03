#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
//# #include "mir_math.h"
#include "../smth_zal.h"

namespace SrtlibRlCrab
{
    void GetDetArr(const double* const sky_arr,
                   const double* const resp_norm_mat_arr,
                   int ndet, int nsky,
                   double* const det_arr);
    
    void GetDenArr(const double* const sky_arr,
                   const double* const flux_arr,
                   const double* const det_0_arr,
                   const double* const bg_arr,
                   const double* const resp_norm_mat_arr,
                   int ndet, int nsky, int nphase,
                   double* const* const den_arr);

    void GetYDashArr(const double* const* const data_arr,
                     const double* const* const den_arr,
                     int ndet, int nphase,
                     double* const* const y_dash_arr);

    void GetMvalArr(const double* const* const y_dash_arr,
                    const double* const resp_norm_mat_arr,
                    const double* const sky_arr,
                    int ndet, int nsky, int nphase,
                    double* const mval_arr);

    void GetNvalArr(const double* const* const y_dash_arr,
                    const double* const flux_arr,
                    const double* const det_0_arr,
                    int ndet, int nphase,
                    double* const nval_arr);

} // namespace SrtlibRlCrab

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_H_


