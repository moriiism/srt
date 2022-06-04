#ifndef MORIIISM_SRT_RICHLUCY_BG2_SMOOTH_SUB_NEWTON_H_
#define MORIIISM_SRT_RICHLUCY_BG2_SMOOTH_SUB_NEWTON_H_

#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"

void GetRhoArr_FromLambda(double lambda,
                          double lip_const,
                          const double* const vval_arr,
                          const double* const mval_arr,
                          int nsky,
                          double* const rho_arr);

double GetNu_FromLambda(double lambda,
                        double lip_const,
                        double wval,
                        double nval);

void GetDerivRhoArr_FromLambda(double lambda,
                               double lip_const,
                               const double* const vval_arr,
                               const double* const mval_arr,
                               int nsky,
                               double* const deriv_rho_arr);

double GetDerivNu_FromLambda(double lambda,
                             double lip_const,
                             double wval,
                             double nval);

double GetSval_FromLambda(double lambda,
                          double lip_const,
                          const double* const vval_arr, double wval,
                          const double* const mval_arr, double nval,
                          int nsky);

double GetDerivSval_FromLambda(double lambda,
                               double lip_const,
                               const double* const vval_arr, double wval,
                               const double* const mval_arr, double nval,
                               int nsky);

double GetLambdaNew_ByNewton(double lambda,
                             double lip_const,
                             const double* const vval_arr, double wval,
                             const double* const mval_arr, double nval,
                             int nsky);

double GetLambda_ByNewton(double lambda_init,
                          const double* const vval_arr, double wval,
                          const double* const mval_arr, double nval,
                          int nsky, double lip_const,
                          int nnewton, double tol_newton);

void GetRhoArrNu_ByNewton(const double* const vval_arr, double wval,
                          const double* const mval_arr, double nval,
                          int nsky,
                          double lip_const,
                          int nnewton, double tol_newton,
                          double lambda,
                          double* const rho_new_arr,
                          double* const nu_new_ptr,
                          double* const lambda_new_ptr);


#endif // MORIIISM_SRT_RICHLUCY_BG2_SMOOTH_SUB_NEWTON_H_
