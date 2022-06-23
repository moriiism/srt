#include "fpsrc_smth_bg_newton.h"
//#include "fpsrc_smth_bg_statval.h"

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

void GetNuArr_FromLambda(double lambda,
                         double lip_const,
                         const double* const wval_arr,
                         const double* const nval_arr,
                         int nsrc,
                         double* const nu_arr)
{
    for(int isrc = 0; isrc < nsrc; isrc++){
        double bval = lambda - lip_const * wval_arr[isrc];
        double num = sqrt(bval * bval + 4 * lip_const * nval_arr[isrc]) - bval;
        double den = 2 * lip_const;
        nu_arr[isrc] = num / den;
    }
}

double GetPhi_FromLambda(double lambda,
                         double nph,
                         double lip_const,
                         double zval,
                         double phi_val,
                         double pval)
{
    double bval = lambda - lip_const * zval;
    double num = (pval - nph) / phi_val - bval;
    double den = lip_const;
    double phi = num / den;
    return phi;
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
void GetDerivNuArr_FromLambda(double lambda,
                              double lip_const,
                              const double* const wval_arr,
                              const double* const nval_arr,
                              int nsrc,
                              double* const deriv_nu_arr)
{
    double* nu_arr = new double[nsrc];
    GetNuArr_FromLambda(lambda,
                        lip_const,
                        wval_arr,
                        nval_arr,
                        nsrc,
                        nu_arr);
    for(int isrc = 0; isrc < nsrc; isrc++){
        double bval = lambda - lip_const * wval_arr[isrc];
        double den = sqrt(bval * bval + 4 * lip_const * nval_arr[isrc]);
        double num = -1 * nu_arr[isrc];
        deriv_nu_arr[isrc] = num / den;
    }
    delete [] nu_arr;
}

// get derivative of phi from lambda
double GetDerivPhi_FromLambda(double lip_const)
{
    double deriv_phi = -1.0 / lip_const; 
    return deriv_phi;
}

double GetSval_FromLambda(double lambda,
                          double lip_const,
                          const double* const vval_arr,
                          const double* const wval_arr,
                          double zval,
                          const double* const mval_arr,
                          const double* const nval_arr,
                          double pval,
                          double phi_val,
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
    double phi = GetPhi_FromLambda(lambda,
                                   nph,
                                   lip_const,
                                   zval,
                                   phi_val,
                                   pval);
    double sval = phi - 1.0;
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

double GetDerivSval_FromLambda(double lambda,
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
    double deriv_phi = GetDerivPhi_FromLambda(lip_const);
    double deriv_sval = deriv_phi;
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

double GetLambdaUpdate_ByNewton(double lambda,
                                double lip_const,
                                const double* const vval_arr,
                                const double* const wval_arr,
                                double zval,
                                const double* const mval_arr,
                                const double* const nval_arr,
                                double pval,
                                double phi_val,
                                int nsky, int nsrc, int nph)
{
    double sval = GetSval_FromLambda(lambda,
                                     lip_const,
                                     vval_arr, wval_arr, zval,
                                     mval_arr, nval_arr, pval,
                                     phi_val,
                                     nsky, nsrc, nph);
    double deriv_sval = GetDerivSval_FromLambda(lambda,
                                                lip_const,
                                                vval_arr, wval_arr,
                                                mval_arr, nval_arr,
                                                nsky, nsrc);
    double lambda_new = lambda - sval / deriv_sval;
    return lambda_new;
}

double GetLambda_ByNewton(FILE* const fp_log,
                          double lambda_init,
                          const double* const vval_arr,
                          const double* const wval_arr,
                          double zval,
                          const double* const mval_arr,
                          const double* const nval_arr,
                          double pval,
                          double phi_val,
                          int nsky, int nsrc, int nph,
                          double lip_const,
                          int nnewton, double tol_newton)
{
    double lambda = lambda_init;
    double lambda_new = lambda_init;
    int flag_converge = 0;
    double sval = 0.0;
    for(int inewton = 0; inewton < nnewton; inewton ++){
        lambda_new = GetLambdaUpdate_ByNewton(lambda,
                                              lip_const,
                                              vval_arr, wval_arr, zval,
                                              mval_arr, nval_arr, pval,
                                              phi_val,
                                              nsky, nsrc, nph);
        sval = GetSval_FromLambda(lambda_new,
                                  lip_const,
                                  vval_arr, wval_arr, zval,
                                  mval_arr, nval_arr, pval,
                                  phi_val,
                                  nsky, nsrc, nph);
        if(fabs(sval) < tol_newton){
            flag_converge = 1;
            break;
        }
        lambda = lambda_new;
    }
    if(flag_converge != 1){
        MiIolib::Printf2(fp_log,
                         "      newton: not converge: sval = %e\n", sval);
    }
    
    return lambda_new;
}


void GetRhoArrNuArrPhi_ByNewton(FILE* const fp_log,
                                const double* const vval_arr,
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
                                double* const lambda_new_ptr)
{
    double lambda_new = GetLambda_ByNewton(fp_log,
                                           lambda,
                                           vval_arr, wval_arr, zval,
                                           mval_arr, nval_arr, pval,
                                           phi_val,
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
    double phi_new = GetPhi_FromLambda(lambda_new,
                                       nph,
                                       lip_const,
                                       zval,
                                       phi_val,
                                       pval);
    *phi_new_ptr = phi_new;
    *lambda_new_ptr = lambda_new;
}
