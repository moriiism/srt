#include "fpsrc_smth_bg_newton_phi0.h"
#include "fpsrc_smth_bg_newton.h"
//#include "fpsrc_smth_bg_statval.h"

double GetSval_FromLambda_Phi0(double lambda,
                               double lip_const,
                               const double* const vval_arr,
                               const double* const wval_arr,
                               const double* const mval_arr,
                               const double* const nval_arr,
                               int nsky, int nsrc, int nph)
{
    double* rho_arr = new double[nsky];
    GetRhoArr_FromLambda(lambda,
                         lip_const,
                         vval_arr,
                         mval_arr,
                         nsky,
                         rho_arr);
    double* nu_arr = new double[nsrc];
    GetNuArr_FromLambda(lambda,
                        lip_const,
                        wval_arr,
                        nval_arr,
                        nsrc,
                        nu_arr);
    double sval = -1.0;
    for(int isky = 0; isky < nsky; isky ++){
        sval += rho_arr[isky];
    }
    for(int isrc = 0; isrc < nsrc; isrc ++){
        sval += nu_arr[isrc];
    }
    delete [] rho_arr;
    delete [] nu_arr;    
    return sval;
}

double GetDerivSval_FromLambda_Phi0(double lambda,
                                    double lip_const,
                                    const double* const vval_arr,
                                    const double* const wval_arr,
                                    const double* const mval_arr,
                                    const double* const nval_arr,
                                    int nsky, int nsrc)
{
    double* deriv_rho_arr = new double[nsky];
    GetDerivRhoArr_FromLambda(lambda,
                              lip_const,
                              vval_arr,
                              mval_arr,
                              nsky,
                              deriv_rho_arr);
    double* deriv_nu_arr = new double[nsrc];
    GetDerivNuArr_FromLambda(lambda,
                             lip_const,
                             wval_arr,
                             nval_arr,
                             nsrc,
                             deriv_nu_arr);
    double deriv_sval = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        deriv_sval += deriv_rho_arr[isky];
    }
    for(int isrc = 0; isrc < nsrc; isrc ++){
        deriv_sval += deriv_nu_arr[isrc];
    }
    delete [] deriv_rho_arr;
    delete [] deriv_nu_arr;
    return deriv_sval;
}

double GetLambdaUpdate_ByNewton_Phi0(double lambda,
                                     double lip_const,
                                     const double* const vval_arr,
                                     const double* const wval_arr,
                                     const double* const mval_arr,
                                     const double* const nval_arr,
                                     int nsky, int nsrc, int nph)
{
    double sval = GetSval_FromLambda_Phi0(lambda,
                                          lip_const,
                                          vval_arr, wval_arr,
                                          mval_arr, nval_arr,
                                          nsky, nsrc, nph);
    double deriv_sval = GetDerivSval_FromLambda_Phi0(lambda,
                                                     lip_const,
                                                     vval_arr, wval_arr,
                                                     mval_arr, nval_arr,
                                                     nsky, nsrc);
    double lambda_new = lambda - sval / deriv_sval;
    return lambda_new;
}

double GetLambda_ByNewton_Phi0(double lambda_init,
                               const double* const vval_arr,
                               const double* const wval_arr,
                               const double* const mval_arr,
                               const double* const nval_arr,
                               int nsky, int nsrc, int nph,
                               double lip_const,
                               int nnewton, double tol_newton)
{
    double lambda = lambda_init;
    double lambda_new = lambda_init;
    int flag_converge = 0;
    double sval = 0.0;
    for(int inewton = 0; inewton < nnewton; inewton ++){
        lambda_new = GetLambdaUpdate_ByNewton_Phi0(lambda,
                                                   lip_const,
                                                   vval_arr, wval_arr,
                                                   mval_arr, nval_arr,
                                                   nsky, nsrc, nph);
        sval = GetSval_FromLambda_Phi0(lambda_new,
                                       lip_const,
                                       vval_arr, wval_arr,
                                       mval_arr, nval_arr,
                                       nsky, nsrc, nph);
        if(fabs(sval) < tol_newton){
            flag_converge = 1;
            break;
        }
        lambda = lambda_new;
    }
    if(flag_converge != 1){
        printf("newton: sval = %e\n", sval);
    }
    return lambda_new;
}


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
                                  double* const lambda_new_ptr)
{
    double lambda_new = GetLambda_ByNewton_Phi0(lambda,
                                                vval_arr, wval_arr,
                                                mval_arr, nval_arr,
                                                nsky, nsrc, nph,
                                                lip_const,
                                                nnewton, tol_newton);
    GetRhoArr_FromLambda(lambda_new,
                         lip_const,
                         vval_arr,
                         mval_arr,
                         nsky,
                         rho_new_arr);
    GetNuArr_FromLambda(lambda_new,
                        lip_const,
                        wval_arr,
                        nval_arr,
                        nsrc,
                        nu_new_arr);
    *lambda_new_ptr = lambda_new;
}
