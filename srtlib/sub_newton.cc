#include "sub_newton.h"

void GetRhoArr_FromLambda(double lambda,
                          double lip_const,
                          const double* const vval_arr,
                          const double* const mval_arr,
                          int nsky,
                          double* const rho_arr)
{
    for(int isky = 0; isky < nsky; isky++){
        double bval = lambda - lip_const * vval_arr[isky];
        double num = sqrt(bval * bval + 4 * lip_const * mval_arr[isky]) - bval;
        double den = 2 * lip_const;
        rho_arr[isky] = num / den;
    }
}

double GetNu_FromLambda(double lambda,
                        double lip_const,
                        double wval,
                        double nval)
{
    double bval = lambda - lip_const * wval;
    double num = sqrt(bval * bval + 4 * lip_const * nval) - bval;
    double den = 2 * lip_const;
    double nu = num / den;
    return nu;
}

// get derivative of rho from lambda
void GetDerivRhoArr_FromLambda(double lambda,
                               double lip_const,
                               const double* const vval_arr,
                               const double* const mval_arr,
                               int nsky,
                               double* const deriv_rho_arr)
{
    double* rho_arr = new double[nsky];
    GetRhoArr_FromLambda(lambda,
                         lip_const,
                         vval_arr,
                         mval_arr,
                         nsky,
                         rho_arr);
    for(int isky = 0; isky < nsky; isky++){
        double bval = lambda - lip_const * vval_arr[isky];
        double den = sqrt(bval * bval + 4 * lip_const * mval_arr[isky]);
        double num = -1 * rho_arr[isky];
        deriv_rho_arr[isky] = num / den;
    }
    delete [] rho_arr;
}

// get derivative of nu from lambda
double GetDerivNu_FromLambda(double lambda,
                             double lip_const,
                             double wval,
                             double nval)
{
    double nu = GetNu_FromLambda(lambda,
                                 lip_const,
                                 wval,
                                 nval);
    double bval = lambda - lip_const * wval;
    double den = sqrt(bval * bval + 4 * lip_const * nval);
    double num = -1 * nu;
    double deriv_nu = num / den;
    return deriv_nu;
}

double GetSval_FromLambda(double lambda,
                          double lip_const,
                          const double* const vval_arr, double wval,
                          const double* const mval_arr, double nval,
                          int nsky)
{
    double* rho_arr = new double[nsky];
    GetRhoArr_FromLambda(lambda,
                         lip_const,
                         vval_arr,
                         mval_arr,
                         nsky,
                         rho_arr);
    double nu = GetNu_FromLambda(lambda,
                                 lip_const,
                                 wval,
                                 nval);
    double sval = nu - 1.0;
    for(int isky = 0; isky < nsky; isky ++){
        sval += rho_arr[isky];
    }
    delete [] rho_arr;
    return sval;
}

double GetDerivSval_FromLambda(double lambda,
                               double lip_const,
                               const double* const vval_arr, double wval,
                               const double* const mval_arr, double nval,
                               int nsky)
{
    double* deriv_rho_arr = new double[nsky];
    GetDerivRhoArr_FromLambda(lambda,
                              lip_const,
                              vval_arr,
                              mval_arr,
                              nsky,
                              deriv_rho_arr);
    double deriv_nu = GetDerivNu_FromLambda(lambda,
                                            lip_const,
                                            wval,
                                            nval);
    double deriv_sval = deriv_nu;
    for(int isky = 0; isky < nsky; isky ++){
        deriv_sval += deriv_rho_arr[isky];
    }
    delete [] deriv_rho_arr;
    return deriv_sval;
}

double GetLambdaNew_ByNewton(double lambda,
                             double lip_const,
                             const double* const vval_arr, double wval,
                             const double* const mval_arr, double nval,
                             int nsky)
{
    double sval = GetSval_FromLambda(lambda,
                                     lip_const,
                                     vval_arr, wval,
                                     mval_arr, nval,
                                     nsky);
    double deriv_sval = GetDerivSval_FromLambda(lambda,
                                                lip_const,
                                                vval_arr, wval,
                                                mval_arr, nval,
                                                nsky);
    double lambda_new = lambda - sval / deriv_sval;
    return lambda_new;
}

double GetLambda_ByNewton(double lambda_init,
                          const double* const vval_arr, double wval,
                          const double* const mval_arr, double nval,
                          int nsky, double lip_const,
                          int nnewton, double tol_newton)
{
    double lambda = lambda_init;
    double lambda_new = lambda_init;
    for(int inewton = 0; inewton < nnewton; inewton ++){
        lambda_new = GetLambdaNew_ByNewton(lambda,
                                           lip_const,
                                           vval_arr, wval,
                                           mval_arr, nval,
                                           nsky);
        double sval = GetSval_FromLambda(lambda_new,
                                         lip_const,
                                         vval_arr, wval,
                                         mval_arr, nval,
                                         nsky);
        if(fabs(sval) < tol_newton){
            break;
        }
        lambda = lambda_new;
    }
    return lambda_new;
}


void GetRhoArrNu_ByNewton(const double* const vval_arr, double wval,
                          const double* const mval_arr, double nval,
                          int nsky,
                          double lip_const,
                          int nnewton, double tol_newton,
                          double lambda,
                          double* const rho_new_arr,
                          double* const nu_new_ptr,
                          double* const lambda_new_ptr)
{
    double lambda_new = GetLambda_ByNewton(lambda,
                                           vval_arr, wval,
                                           mval_arr, nval,
                                           nsky, lip_const,
                                           nnewton, tol_newton);
    GetRhoArr_FromLambda(lambda_new,
                         lip_const,
                         vval_arr,
                         mval_arr,
                         nsky,
                         rho_new_arr);
    double nu_new = GetNu_FromLambda(lambda_new,
                                     lip_const,
                                     wval,
                                     nval);
    *nu_new_ptr = nu_new;
    *lambda_new_ptr = lambda_new;
}
