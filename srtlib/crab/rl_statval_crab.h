#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_STATVAL_CRAB_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_STATVAL_CRAB_H_

#include "mib_blas.h"
#include "mir_math.h"

namespace SrtlibRlStatvalCrab
{
    double GetHellingerDist(const double* const rho_arr,
                            const double* const nu_arr,
                            const double* const rho_new_arr,
                            const double* const nu_new_arr,
                            int nsky, int nphase);
} // SrtlibRlStatvalCrab

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_STATVAL_CRAB_H_
