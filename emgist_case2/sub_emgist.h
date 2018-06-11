#ifndef MORIIISM_SRT_EMGIST_SUB_EMGIST_H_
#define MORIIISM_SRT_EMGIST_SUB_EMGIST_H_

#include "mi_rand.h"
#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

void GetNdet(string respdir, int* const ndetx_ptr, int* const ndety_ptr);
void LoadResp(string respdir, int nskyx, int nskyy,
              double** const mat_arr_ptr,
              int* const ndetx_ptr,
              int* const ndety_ptr);
void SolveByEM(const double* const rho_arr, int nph,
               const double* const data_arr,
               const double* const resp_mat_arr,
               double beta, double mu, double lconst_init,
               int nem, double tol_em, int npm, double tol_pm,
               int flag_line_search,
               string outdir, string outfile_head,
               int ndet, int nskyx, int nskyy, double epsilon,
               int bitpix, double* const out_arr);

void GetLineSearch(const double* const xval_arr,
                   const double* const xval_new_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   double beta, double mu,
                   int ndet, int nskyx, int nskyy,
                   double epsilon,
                   double* const out_arr,
                   int* flag_saturate_ptr);

void IsRhoAgainPlus(const double* const rho_arr,
                    const double* const data_arr,
                    const double* const resp_mat_arr,
                    int ndet, int nsky,
                    int* const out_arr);

double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu,
                int ndet, int nskyx, int nskyy, double epsilon);

void GetDiffL(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              double beta, double mu,
              int ndet, int nskyx, int nskyy, double epsilon,
              double* const out_arr);

void GetSuppArrByTrunc(const double* const rho_arr,
                       const double* const mval_arr,
                       double beta, int nsky,
                       int* const out_arr,
                       vector<int>* const index_supp_vec_ptr);

double GetFuncLsub(const double* const rho_arr,
                   const double* const mval_arr,
                   double beta, double mu,
                   int nskyx, int nskyy);

double GetHellingerDist(const double* const rho_arr,
                        const double* const rho_new_arr,
                        int nsky);

double GetKLDiv(const double* const rho_arr,
                const double* const rho_new_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky);

double GetFindLconst(const double* const rho_arr,
                     const double* const mval_arr,
                     double beta, double mu,
                     int nskyx, int nskyy,
                     double lconst_init,
                     double tau_init,
                     double* const pLy_arr,
                     double* const tau_out_ptr);

void GetFuncM(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              int ndet, int nsky,
              double* const out_arr);

void GetFuncSigma(const double* const rho_arr,
                  double mu,
                  double lconst, int nskyx, int nskyy,
                  double* out_arr);

double GetTau(const double* const mval_arr,
              const double* const sigma_arr, int nsky,
              double lconst, double beta,
              double tau_init);

double GetFuncS(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double beta, double lconst);

double GetFuncDiffS(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double beta, double lconst);

void GetFuncRho(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double beta, double lconst,
                double* const out_arr);

void GetFuncDiffRho(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double beta, double lconst,
                    double* const out_arr);

double GetQMinusF(const double* const rho_new_arr,
                  const double* const rho_arr,
                  double mu, double lconst,
                  int nskyx, int nskyy);

double GetFuncF(const double* const rho_arr,
                double mu,
                int nskyx, int nskyy);

void GetDiffF(const double* const rho_arr,
              double mu,
              int nskyx, int nskyy,
              double* const out_arr);

double GetFuncG(const double* const rho_arr,
                const double* const mval_arr,
                int nsky, double beta);

void GetDiffG(const double* const rho_arr,
              const double* const mval_arr,
              int nsky, double beta,
              double* const out_arr);

double GetTermV(const double* const rho_arr, int nskyx, int nskyy);

void GetDiffTermV(const double* const rho_arr, int nskyx, int nskyy,
                  double* const rho_diff_arr);

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

#endif // MORIIISM_SRT_EMGIST_SUB_EMGIST_H_
