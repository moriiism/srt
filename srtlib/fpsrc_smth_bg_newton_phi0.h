#ifndef MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_NEWTON_PHI0_H_
#define MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_NEWTON_PHI0_H_

#include "mib_blas.h"
#include "mif_fits.h"
//#include "mir_math.h"

double GetSval_FromLambda_Phi0(double lambda,
                               double lip_const,
                               const double* const vval_arr,
                               const double* const wval_arr,
                               const double* const mval_arr,
                               const double* const nval_arr,
                               int nsky, int nsrc, int nph);

double GetDerivSval_FromLambda_Phi0(double lambda,
                                    double lip_const,
                                    const double* const vval_arr,
                                    const double* const wval_arr,
                                    const double* const mval_arr,
                                    const double* const nval_arr,
                                    int nsky, int nsrc);

double GetLambdaUpdate_ByNewton_Phi0(double lambda,
                                     double lip_const,
                                     const double* const vval_arr,
                                     const double* const wval_arr,
                                     const double* const mval_arr,
                                     const double* const nval_arr,
                                     int nsky, int nsrc, int nph);

double GetLambda_ByNewton_Phi0(double lambda_init,
                               const double* const vval_arr,
                               const double* const wval_arr,
                               const double* const mval_arr,
                               const double* const nval_arr,
                               int nsky, int nsrc, int nph,
                               double lip_const,
                               int nnewton, double tol_newton);

void GetRhoArrNuArr_ByNewton_Phi0(const double* const vval_arr,
                                  const double* const wval_arr,
                                  const double* const mval_arr,
                                  const double* const nval_arr,
                                  int nsky, int nsrc, int nph,
                                  double lip_const,
                                  int nnewton, double tol_newton,
                                  double lambda,
                                  double* const rho_new_arr,
                                  double* const nu_new_arr,
                                  double* const lambda_new_ptr);

#endif // MORIIISM_SRT_SRTLIB_FPSRC_SMTH_BG_NEWTON_PHI0_H_
