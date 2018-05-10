#ifndef MORIIISM_SRT_PMMN_SUB_PMMN_CASE2_H_
#define MORIIISM_SRT_PMMN_SUB_PMMN_CASE2_H_

#include "mi_rand.h"
#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"

void GetNdet(string respdir, int* const ndetx_ptr, int* const ndety_ptr);
void LoadResp(string respdir, int nskyx, int nskyy,
              double** const mat_arr_ptr,
              int* const ndetx_ptr,
              int* const ndety_ptr);
void SolveByProxMapMN(const double* const rho_arr, int nph,
                      const double* const data_arr,
                      const double* const resp_mat_arr,
                      double beta, double mu, double lconst,
                      double tol, double tol_em, int nstep,
                      string outdir, string outfile_head,
                      int ndet, int nskyx, int nskyy, double epsilon,
                      int bitpix, double* const out_arr);
int GetFindIk(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              double beta, double mu, double lconst, double eta,
              double tol_em,
              int ndet, int nskyx, int nskyy, double epsilon);

double GetFindIkBisect(const double* const rho_arr,
                       const double* const data_arr,
                       const double* const resp_mat_arr,
                       double beta, double mu, double lconst, double eta,
                       double tol_em,
                       int ndet, int nskyx, int nskyy, double epsilon);

void GetProxMap(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu, double lconst,
                int ndet, int nskyx, int nskyy, double epsilon,
                double tol_em,
                double* const out_arr, int* flag_good_ptr);
void GetLineSearch(const double* const xval_arr,
                   const double* const xval_new_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   const double* const sigma_arr,
                   int ndet, int nsky, double lconst,
                   double epsilon,
                   double* const out_arr,
                   int* flag_saturate_ptr);
double GetKLDiv(const double* const rho_arr,
                const double* const rho_new_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky);
void GetFuncM(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              int ndet, int nsky,
              double* const out_arr);
void GetFuncRho(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon,
                double* const out_arr);
void GetFuncDiffRho(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double lconst, double epsilon,
                    double* const out_arr);
void GetFuncTauThres(const double* const sigma_arr,
                     const double* const mval_arr,
                     int nsky, double lconst, double epsilon,
                     double* const out_arr);
double GetFuncS(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon);
double GetFuncDiffS(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double lconst, double epsilon);
double GetFuncLem(const double* const rho_arr,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const sigma_arr,
                  int ndet, int nsky,
                  double lconst);
double GetFuncLem(const double* const rho_arr,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const sigma_arr,
                  int ndet, int nsky,
                  double lconst);
double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu,
                int ndet, int nskyx, int nskyy);
double GetFuncLSupp(const double* const rho_arr,
                    const int* const rho_supp_arr,
                    const double* const data_arr,
                    const double* const resp_mat_arr,
                    double beta, double mu,
                    int ndet, int nskyx, int nskyy);
double GetQMinusF(const double* const rho_new_arr,
                  const double* const rho_arr,
                  double beta, double mu, double lconst,
                  int nskyx, int nskyy);
void GetFuncSigma(const double* const rho_arr,
                  const double* const data_arr,
                  double beta, double mu,
                  double lconst, int nskyx, int nskyy,
                  double* out_arr);
void GetFuncSigmaSupp(const double* const rho_arr,
                      const int* const rho_supp_arr,
                      const double* const data_arr,
                      double beta, double mu,
                      double lconst, int nskyx, int nskyy,
                      double* out_arr);
double GetFuncF(const double* const rho_arr,
                double beta, double mu,
                int nskyx, int nskyy);
double GetFuncFSupp(const double* const rho_arr,
                    const int* const rho_supp_arr,
                    double beta, double mu,
                    int nskyx, int nskyy);
void GetDiffF(const double* const rho_arr,
              double beta, double mu,
              int nskyx, int nskyy,
              double* const out_arr);
void GetDiffFSupp(const double* const rho_arr,
                  const int* const rho_supp_arr,
                  double beta, double mu,
                  int nskyx, int nskyy,
                  double* const out_arr);
double GetFuncG(const double* const rho_arr, 
                const double* const data_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky);
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

#endif // MORIIISM_SRT_PMMN_SUB_PMMN_CASE2_H_

