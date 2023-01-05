#ifndef MORIIISM_SRT_SRTLIB_SUB_PM_H_
#define MORIIISM_SRT_SRTLIB_SUB_PM_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

void GetVvalArr(const double* const rho_arr,
                int nskyx, int nskyy,
                double mu, double lip_const,
                double* const vval_arr);

double GetWval(double nu);

#endif // MORIIISM_SRT_SRTLIB_SUB_PM_H_
