#include "rl_crab_smth_pf_newton.h"

void SrtlibRlCrabSmthPfNewton::GetRhoArr_FromLambda(
    double lambda,
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

void SrtlibRlCrabSmthPfNewton::GetNuArr_FromLambda(
    double lambda,
    double lip_const,
    const double* const wval_arr,
    const double* const nval_arr,
    int nphase,
    double* const nu_arr)
{
    for(int iphase = 0; iphase < nphase; iphase++){
        double bval = lambda - lip_const * wval_arr[iphase];
        double num = sqrt(bval * bval + 4 * lip_const * nval_arr[iphase]) - bval;
        double den = 2 * lip_const;
        nu_arr[iphase] = num / den;
    }
}

// get derivative of rho from lambda
void SrtlibRlCrabSmthPfNewton::GetDerivRhoArr_FromLambda(
    double lambda,
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
double SrtlibRlCrabSmthPfNewton::GetDerivNuArr_FromLambda(
    double lambda,
    double lip_const,
    const double* const wval_arr,
    const double* const nval_arr,
    int nphase,
    double* const deriv_nu_arr)
{
    double* nu_arr = new double[nphase];
    GetNuArr_FromLambda(lambda,
                        lip_const,
                        wval_arr,
                        nval_arr,
                        nphase,
                        nu_arr);
    for(int iphase = 0; iphase < nphase; iphase++){
        double bval = lambda - lip_const * wval_arr[iphase];
        double den = sqrt(bval * bval + 4 * lip_const * nval_arr[iphase]);
        double num = -1 * nu_arr[iphase];
        deriv_nu_arr[iphase] = num / den;
    }
    delete [] nu_arr;
}

double SrtlibRlCrabSmthPfNewton::GetSval_FromLambda(
    double lambda,
    double lip_const,
    const double* const vval_arr,
    const double* const wval_arr,
    const double* const mval_arr,
    const double* const nval_arr,
    int nsky, int nphase)
{
    double* rho_arr = new double[nsky];
    double* nu_arr = new double[nsky];
    GetRhoArr_FromLambda(lambda,
                         lip_const,
                         vval_arr,
                         mval_arr,
                         nsky,
                         rho_arr);
    GetNuArr_FromLambda(lambda,
                        lip_const,
                        wval_arr,
                        nval_arr,
                        nphase,
                        nu_arr);
    double sval = -1.0;
    for(int isky = 0; isky < nsky; isky ++){
        sval += rho_arr[isky];
    }
    for(int iphase = 0; iphase < nphase; iphase ++){
        sval += nu_arr[iphase];
    }
    delete [] rho_arr;
    delete [] nu_arr;
    return sval;
}

double SrtlibRlCrabSmthPfNewton::GetDerivSval_FromLambda(
    double lambda,
    double lip_const,
    const double* const vval_arr,
    const double* const wval_arr,
    const double* const mval_arr,
    const double* const nval_arr,
    int nsky, int nphase)
{
    double* deriv_rho_arr = new double[nsky];
    double* deriv_nu_arr = new double[nphase];
    GetDerivRhoArr_FromLambda(lambda,
                              lip_const,
                              vval_arr,
                              mval_arr,
                              nsky,
                              deriv_rho_arr);
    GetDerivNuArr_FromLambda(lambda,
                             lip_const,
                             wval_arr,
                             nval_arr,
                             nphase,
                             deriv_nu_arr);
    double deriv_sval = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        deriv_sval += deriv_rho_arr[isky];
    }
    for(int iphase = 0; iphase < nphase; iphase ++){
        deriv_sval += deriv_nu_arr[iphase];
    }
    delete [] deriv_rho_arr;
    delete [] deriv_nu_arr;
    return deriv_sval;
}

double SrtlibRlCrabSmthPfNewton::GetLambdaNew_ByNewton(
    double lambda,
    double lip_const,
    const double* const vval_arr,
    const double* const wval_arr,
    const double* const mval_arr,
    const double* const nval_arr,
    int nsky, int nphase)
{
    double sval = GetSval_FromLambda(lambda,
                                     lip_const,
                                     vval_arr, wval_arr,
                                     mval_arr, nval_arr,
                                     nsky, nphase);
    double deriv_sval = GetDerivSval_FromLambda(lambda,
                                                lip_const,
                                                vval_arr, wval_arr,
                                                mval_arr, nval_arr,
                                                nsky, nphase);
    double lambda_new = lambda - sval / deriv_sval;

    return lambda_new;
}

double SrtlibRlCrabSmthPfNewton::GetLambda_ByNewton(
    FILE* const fp_log,
    double lambda_init,
    const double* const vval_arr,
    const double* const wval_arr,
    const double* const mval_arr,
    const double* const nval_arr,
    int nsky, int nphase,
    double lip_const,
    int nnewton, double tol_newton)
{
    double lambda = lambda_init;
    double lambda_new = lambda_init;
    int flag_converge = 0;
    double sval = 0.0;
    for(int inewton = 0; inewton < nnewton; inewton ++){
        lambda_new = GetLambdaNew_ByNewton(lambda,
                                           lip_const,
                                           vval_arr, wval_arr,
                                           mval_arr, nval_arr,
                                           nsky, nphase);
        sval = GetSval_FromLambda(lambda_new,
                                  lip_const,
                                  vval_arr, wval_arr,
                                  mval_arr, nval_arr,
                                  nsky, nphase);
        if(fabs(sval) < tol_newton){
            flag_converge = 1;
            break;
        }
        lambda = lambda_new;
    }
    if(flag_converge != 1){
        MiIolib::Printf2(fp_log,
                         "      newton: not converge: sval = %e\n", sval);
        MiIolib::Printf2(fp_log,
                         "            : lambda_new = %e\n", lambda_new);
        MiIolib::Printf2(fp_log,
                         "            : lambda_init = %e\n", lambda_init);
        MiIolib::Printf2(fp_log,
                         "            : lip_const = %e\n", lip_const);
        MiIolib::Printf2(fp_log,
                         "      abort.\n");        
        abort();
    }
    
    return lambda_new;
}

void SrtlibRlCrabSmthPfNewton::GetRhoArrNuArr_ByNewton(
    FILE* const fp_log,
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
    double* const lambda_new_ptr)
{
    double lambda_new = GetLambda_ByNewton(fp_log,
                                           lambda,
                                           vval_arr, wval_arr,
                                           mval_arr, nval_arr,
                                           nsky, nphase,
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
                        nphase,
                        nu_new_arr);
    *lambda_new_ptr = lambda_new;
}
