#ifndef MORIIISM_SRT_SRTLIB_RL_STATVAL_H_
#define MORIIISM_SRT_SRTLIB_RL_STATVAL_H_

#include "mib_blas.h"
#include "mir_math.h"

double GetHellingerDist(const double* const rho_arr,
                        const double* const rho_new_arr,
                        int nsky);

double GetNegLogLike(const double* const rho_arr,
                     const double* const data_arr,
                     const double* const resp_norm_mat_arr,
                     int ndet, int nsky, double epsilon);

#endif // MORIIISM_SRT_SRTLIB_RL_STATVAL_H_

