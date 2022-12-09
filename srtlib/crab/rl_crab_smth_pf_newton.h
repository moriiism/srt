#ifndef MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_NEWTON_H_
#define MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_NEWTON_H_

#include "mi_iolib.h"
#include "mir_math.h"

namespace SrtlibRlCrabSmthPfNewton
{
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
                             int nphase,
                             double* const nu_arr);

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
                                  int nphase,
                                  double* const deriv_nu_arr);

    double GetSval_FromLambda(double lambda,
                              double lip_const,
                              const double* const vval_arr,
                              const double* const wval_arr,
                              const double* const mval_arr,
                              const double* const nval_arr,
                              int nsky, int nphase);

    double GetDerivSval_FromLambda(double lambda,
                                   double lip_const,
                                   const double* const vval_arr,
                                   const double* const wval_arr,
                                   const double* const mval_arr,
                                   const double* const nval_arr,
                                   int nsky, int nphase);

    double GetLambdaNew_ByNewton(double lambda,
                                 double lip_const,
                                 const double* const vval_arr,
                                 const double* const wval_arr,
                                 const double* const mval_arr,
                                 const double* const nval_arr,
                                 int nsky, int nphase);

    double GetLambda_ByNewton(FILE* const fp_log,
                              double lambda_init,
                              const double* const vval_arr,
                              const double* const wval_arr,
                              const double* const mval_arr,
                              const double* const nval_arr,
                              int nsky, int nphase,
                              double lip_const,
                              int nnewton, double tol_newton);

    void GetRhoArrNuArr_ByNewton(FILE* const fp_log,
                                 const double* const vval_arr,
                                 const double* const wval_arr,
                                 const double* const mval_arr,
                                 const double* const nval_arr,
                                 int nsky, int nphase,
                                 double lip_const,
                                 int nnewton, double tol_newton,
                                 double lambda,
                                 double* const rho_new_arr,
                                 double* const nu_new_arr,
                                 double* const lambda_new_ptr);
    

} // namespace SrtlibRlCrabSmthPfNewton

#endif // MORIIISM_SRT_SRTLIB_CRAB_RL_CRAB_SMTH_PF_NEWTON_H_
