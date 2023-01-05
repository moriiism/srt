#ifndef MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_NEWTON_H_
#define MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_NEWTON_H_

#include "mir_math.h"

namespace SrtlibRlBg2SmthNewton
{
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

    // get derivative of rho from lambda
    void GetDerivRhoArr_FromLambda(double lambda,
                                   double lip_const,
                                   const double* const vval_arr,
                                   const double* const mval_arr,
                                   int nsky,
                                   double* const deriv_rho_arr);

    // get derivative of nu from lambda
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

    double GetLambda_ByNewton(FILE* const fp_log,
                              double lambda_init,
                              const double* const vval_arr, double wval,
                              const double* const mval_arr, double nval,
                              int nsky, double lip_const,
                              int nnewton, double tol_newton);

    void GetRhoArrNu_ByNewton(FILE* const fp_log,
                              const double* const vval_arr, double wval,
                              const double* const mval_arr, double nval,
                              int nsky,
                              double lip_const,
                              int nnewton, double tol_newton,
                              double lambda,
                              double* const rho_new_arr,
                              double* const nu_new_ptr,
                              double* const lambda_new_ptr);
    

} // namespace SrtlibRlBg2SmthNewton

#endif // MORIIISM_SRT_SRTLIB_RL_BG2_SMTH_NEWTON_H_
