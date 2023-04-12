#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_DET2_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_DET2_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "../smth_zal.h"

namespace SrtlibRlCrabDet2
{
    void GetMvalArr(
        const double* const* const y_dash_det1_arr,
        const double* const* const y_dash_det2_arr,    
        const double* const resp_norm_mat_det1_arr,
        const double* const resp_norm_mat_det2_arr,
        const double* const sky_arr,
        int ndet, int nsky, int nphase,
        double* const mval_arr);

    void GetNvalArr(
        const double* const* const y_dash_det1_arr,
        const double* const* const y_dash_det2_arr,    
        const double* const flux_arr,
        const double* const det_0_det1_arr,
        const double* const det_0_det2_arr,    
        int ndet, int nphase,
        double* const nval_arr);

} // namespace SrtlibRlCrabDet2

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_DET2_H_


