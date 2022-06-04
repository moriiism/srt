#ifndef MORIIISM_SRT_RICHLUCY_BG2_SMOOTH_SUB_PM_H_
#define MORIIISM_SRT_RICHLUCY_BG2_SMOOTH_SUB_PM_H_

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

double GetTermV(const double* const rho_arr, int nskyx, int nskyy);

void GetDiffTermV(const double* const rho_arr, int nskyx, int nskyy,
                  double* const rho_diff_arr);

int GetIbin(int ibinx, int ibiny, int nbinx);


#endif // MORIIISM_SRT_RICHLUCY_BG2_SMOOTH_SUB_PM_H_
