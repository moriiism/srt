#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_STATVAL_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_STATVAL_H_

#include "mib_blas.h"
#include "mir_math.h"

double GetHellingerDist(const double* const rho_arr,
                        const double* const nu_arr,
                        double phi,
                        const double* const rho_new_arr,
                        const double* const nu_new_arr,
                        double phi_new,
                        int nsky, int nsrc);

#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_STATVAL_H_

