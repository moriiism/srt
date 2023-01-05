#ifndef MORIIISM_SRT_SRTLIB_SMTH_H_
#define MORIIISM_SRT_SRTLIB_SMTH_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

namespace SrtlibSmth
{

    int GetIbin(int ibinx, int ibiny, int nbinx);

    double GetTermV(const double* const rho_arr, int nskyx, int nskyy);

    void GetDiffTermV(const double* const rho_arr, int nskyx, int nskyy,
                      double* const rho_diff_arr);
} // namespace SrtlibSmth

#endif // MORIIISM_SRT_SRTLIB_SMTH_H_
