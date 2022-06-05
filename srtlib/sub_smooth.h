#ifndef MORIIISM_SRT_SRTLIB_SUB_SMOOTH_H_
#define MORIIISM_SRT_SRTLIB_SUB_SMOOTH_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

int GetIbin(int ibinx, int ibiny, int nbinx);

double GetTermV(const double* const rho_arr, int nskyx, int nskyy);

void GetDiffTermV(const double* const rho_arr, int nskyx, int nskyy,
                  double* const rho_diff_arr);

#endif // MORIIISM_SRT_SRTLIB_SUB_SMOOTH_H_
