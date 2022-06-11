#ifndef MORIIISM_SRT_SRTLIB_SIM_H_
#define MORIIISM_SRT_SRTLIB_SIM_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"
#include "TRandom3.h"

void GenRandomEvtFromProbDist(const double* const prob_arr, int nbin,
                              int nevt, int rand_seed,
                              double* const out_bin_arr,
                              int*   const out_evt_arr);
void GenCVImageByPartition(int* const evt_arr, int nevt,
                           int nfold, int rand_seed_partition,
                           int nbin,
                           double** const out_tr_arr,
                           double** const out_vl_arr);
void GenCVImage(const double* const prob_arr, int nbin,
                int nevt, int rand_seed, int nfold,
                double** const out_tr_arr,
                double** const out_vl_arr,
                double*  const out_arr);


#endif // MORIIISM_SRT_SRTLIB_SIM_H_

