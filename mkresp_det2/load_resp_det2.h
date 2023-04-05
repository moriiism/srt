#ifndef MORIIISM_SRT_MKRESP_DET2_SUB_H_
#define MORIIISM_SRT_MKRESP_DET2_SUB_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

void LoadRespDet2(string respdir1, string respdir2,
                  int nskyx, int nskyy,
                  int nphoton_input,
                  double** const resp_norm_det1_arr_ptr,
                  double** const resp_norm_det2_arr_ptr,
                  double** const eff_arr_ptr,
                  int* const ndetx_ptr,
                  int* const ndety_ptr);

void GetNdet(string respdir, int* const ndetx_ptr,
             int* const ndety_ptr);

#endif // MORIIISM_SRT_MKRESP_DET2_SUB_H_
