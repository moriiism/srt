#ifndef MORIIISM_SRT_RICHLUCY_SUB_RICHLUCY_H_
#define MORIIISM_SRT_RICHLUCY_SUB_RICHLUCY_H_

#include "mi_rand.h"
#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

void GetNdet(string respdir, int* const ndetx_ptr, int* const ndety_ptr);
void LoadResp(string respdir, int nskyx, int nskyy,
              double epsilon,
              double** const mat_arr_ptr,
              int* const ndetx_ptr,
              int* const ndety_ptr);
void Richlucy(const double* const rho_arr, int nph,
              const double* const data_arr,
              const double* const resp_mat_arr,
              int nem, double tol_em,
              double tol_diff_l_var,
              int flag_line_search,
              string outdir, string outfile_head,
              int ndet, int nskyx, int nskyy, double epsilon,
              int bitpix, double* const out_arr);

void GetLineSearch(const double* const xval_arr,
                   const double* const xval_new_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   vector<int> index_supp_vec,
                   double beta, double mu,
                   int ndet, int nskyx, int nskyy,
                   double epsilon,
                   double* const out_arr);

double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                int ndet, int nskyx, int nskyy, double epsilon);

double GetKLDiv(const double* const rho_arr,
                const double* const rho_new_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky);

void GetNextRhoArr(const double* const rho_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   int ndet, int nsky,
                   double* const out_arr);

double GetSumLogPlus(const double* const data_arr, int ndata,
                     const int* const supp_arr);

void GetInvArr(const double* const data_arr, int ndata,
               double* const out_arr);

void GetInvArrSupp(const double* const data_arr, int ndata,
                   const int* const supp_arr, 
                   double* const out_arr);

void GetSuppArr(const double* const data_arr, int ndata,
                double epsilon,
                int* const out_arr);

int GetNZeroArr(const double* const data_arr, int ndata,
                double epsilon);

int GetIbin(int ibinx, int ibiny, int nbinx);

void GetMinMax(const double* const data_arr, int ndata,
               double* min_ptr, double* max_ptr);

void GetMinMaxSupp(const double* const data_arr,
                   const int* const index_supp_arr,
                   int ndata,
                   double* min_ptr, double* max_ptr);

#endif // MORIIISM_SRT_RICHLUCY_SUB_RICHLUCY_H_
