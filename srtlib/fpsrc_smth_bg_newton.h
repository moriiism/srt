#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_NEWTON_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_NEWTON_H_

#include "mib_blas.h"
#include "mif_fits.h"
//#include "mir_math.h"

void GetRhoArr_FromLambda(double lambda,
                          double lip_const,
                          const double* const vval_arr,
                          const double* const mval_arr,
                          int nsky,
                          double* const rho_arr);

void GetNuArr_FromLambda(double lambda,
                         double lip_const,
                         const double* const wval_arr,
                         const double* const nval_arr,
                         int nsrc,
                         double* const nu_arr);

double GetPhi_FromLambda(double lambda,
                         double nph,
                         double lip_const,
                         double zval,
                         double phi_val,
                         double pval);

// get derivative of rho from lambda
void GetDerivRhoArr_FromLambda(double lambda,
                               double lip_const,
                               const double* const vval_arr,
                               const double* const mval_arr,
                               int nsky,
                               double* const deriv_rho_arr);

// get derivative of nu from lambda
void GetDerivNuArr_FromLambda(double lambda,
                              double lip_const,
                              const double* const wval_arr,
                              const double* const nval_arr,
                              int nsrc,
                              double* const deriv_nu_arr);

// get derivative of phi from lambda
double GetDerivPhi_FromLambda(double lip_const);

double GetSval_FromLambda(double lambda,
                          double lip_const,
                          const double* const vval_arr,
                          const double* const wval_arr,
                          double zval,
                          const double* const mval_arr,
                          const double* const nval_arr,
                          double pval,
                          double phi_val,
                          int nsky, int nsrc, int nph);

double GetDerivSval_FromLambda(double lambda,
                               double lip_const,
                               const double* const vval_arr,
                               const double* const wval_arr,
                               const double* const mval_arr,
                               const double* const nval_arr,
                               int nsky, int nsrc);

double GetLambdaUpdate_ByNewton(double lambda,
                                double lip_const,
                                const double* const vval_arr,
                                const double* const wval_arr,
                                double zval,
                                const double* const mval_arr,
                                const double* const nval_arr,
                                double pval,
                                double phi_val,
                                int nsky, int nsrc, int nph);

double GetLambda_ByNewton(double lambda_init,
                          const double* const vval_arr,
                          const double* const wval_arr,
                          double zval,
                          const double* const mval_arr,
                          const double* const nval_arr,
                          double pval,
                          double phi_val,
                          int nsky, int nsrc, int nph,
                          double lip_const,
                          int nnewton, double tol_newton);

void GetRhoArrNuArrPhi_ByNewton(const double* const vval_arr,
                                const double* const wval_arr,
                                double zval,
                                const double* const mval_arr,
                                const double* const nval_arr,
                                double pval,
                                double phi_val,
                                int nsky, int nsrc, int nph,
                                double lip_const,
                                int nnewton, double tol_newton,
                                double lambda,
                                double* const rho_new_arr,
                                double* const nu_new_arr,
                                double* const phi_new_ptr,
                                double* const lambda_new_ptr);

#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_NEWTON_H_
