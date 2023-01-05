#ifndef MORIIISM_SRT_SRTLIB_SMTH_ZAL_H_
#define MORIIISM_SRT_SRTLIB_SMTH_ZAL_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

namespace SrtlibSmthZal
{
    int GetIbin(int ibinx, int ibiny, int nbinx);
    double GetDerivUAlpha(int iskyx, int iskyy,
                          int nskyx, int nskyy);
    void GetDerivUAlphaArr(int nskyx, int nskyy,
                           double* const alpha_arr);
    void GetDerivUBetaArr(const double* const sky_dash_arr,
                          int nskyx, int nskyy,
                          double* const beta_arr);
    
} // namespace SrtlibSmthZal

#endif // MORIIISM_SRT_SRTLIB_SMTH_ZAL_H_
