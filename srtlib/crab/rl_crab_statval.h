#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_STATVAL_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_STATVAL_H_

#include "mib_blas.h"
#include "mir_math.h"

namespace SrtlibCrabRlCrabStatval
{
    double GetHellingerDist(const double* const sky_pre_arr,
                            const double* const flux_pre_arr,
                            const double* const sky_new_arr,
                            const double* const flux_new_arr,
                            int nsky, int nphase);
} // SrtlibCrabRlCrabStatval

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_STATVAL_H_
