#ifndef MORIIISM_SRT_EMGIST_SUB_EMGIST_H_
#define MORIIISM_SRT_EMGIST_SUB_EMGIST_H_

#include "mi_rand.h"
#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"
#include "TRandom3.h"

void GetNdet(string respdir, int* const ndetx_ptr, int* const ndety_ptr);
void LoadResp(string respdir, int nskyx, int nskyy,
              double** const mat_arr_ptr,
              int* const ndetx_ptr,
              int* const ndety_ptr);
void GenRandomEvtFromProbDist(const double* const prob_arr, int nbin,
                              int nevt, int rand_seed,
                              double* const out_arr);
void GenCVImage(const double* const prob_arr, int nbin,
                int nevt, int rand_seed, int nfold,
                double** const out_tr_arr,
                double** const out_vl_arr,
                double*  const out_arr);

#endif // MORIIISM_SRT_EMGIST_SUB_EMGIST_H_
