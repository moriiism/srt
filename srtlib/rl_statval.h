#ifndef MORIIISM_SRT_SRTLIB_RL_STATVAL_H_
#define MORIIISM_SRT_SRTLIB_RL_STATVAL_H_

#include "mi_base.h"
#include "mib_blas.h"

namespace SrtlibRlStatval
{
    double GetHellingerDist(const double* const sky_pre_arr,
                            const double* const sky_new_arr,
                            int nsky);

} // SrtlibRlStatval

#endif // MORIIISM_SRT_SRTLIB_RL_STATVAL_H_

