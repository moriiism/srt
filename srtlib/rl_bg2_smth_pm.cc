#include "smth.h"
#include "rl_statval.h"
#include "rl_bg2_smth_pm.h"
#include "rl_bg2_smth_newton.h"

void SrtlibRlBg2SmthPm::GetVvalArr(const double* const rho_arr,
                                   int nskyx, int nskyy,
                                   double mu, double lip_const,
                                   double* const vval_arr)
{
    int nsky = nskyx * nskyy;
    double* deriv_rho_arr = new double[nsky];
    SrtlibSmth::GetDiffTermV(rho_arr, nskyx, nskyy, deriv_rho_arr);
    for(int isky = 0; isky < nsky; isky++){
        vval_arr[isky] = rho_arr[isky] - mu / lip_const * deriv_rho_arr[isky];
    }
    delete [] deriv_rho_arr;
}

double SrtlibRlBg2SmthPm::GetWval(double nu)
{
    double wval = nu;
    return wval;
}


double SrtlibRlBg2SmthPm::GetFindLipConst(
    FILE* const fp_log,
    const double* const rho_arr, double nu,
    const double* const mval_arr, double nval,
    double mu,
    int nskyx, int nskyy,
    double lip_const, double lambda,
    int nnewton, double tol_newton)
{
    int nsky = nskyx * nskyy;
    int ik_max = 1000;
    double eta = 1.2;

    double* rho_new_arr = new double[nsky];
    double nu_new = 0.0;
    double lip_const_new = 1.0;
    int flag_find_lip_const = 0;
    for(int ik = 0; ik < ik_max; ik ++){
        // printf("loop in findlipconst: ik = %d\n", ik);
        lip_const_new = lip_const * pow(eta, ik);

        double* vval_arr = new double[nsky];
        GetVvalArr(rho_arr,
                   nskyx, nskyy,
                   mu, lip_const_new,
                   vval_arr);
        double wval = GetWval(nu);
        double lambda_new = 0.0;
        SrtlibRlBg2SmthNewton::GetRhoArrNu_ByNewton(
            fp_log,
            vval_arr, wval,
            mval_arr, nval,
            nsky,
            lip_const_new,
            nnewton, tol_newton,
            lambda,
            rho_new_arr,
            &nu_new,
            &lambda_new);
        double qminusf = GetQMinusF(rho_new_arr, nu_new,
                                    rho_arr, nu,
                                    mu, lip_const_new,
                                    nskyx, nskyy);
        delete [] vval_arr;
        if(qminusf >= 0.0){
            flag_find_lip_const = 1;
            break;
        }
    }
    if (flag_find_lip_const == 0){
        printf("warning: GetFindLipConst: no lip_const is found.\n");
        printf("warning: lip_const_new = %e\n", lip_const_new);
        abort();
    }
    
    delete [] rho_new_arr;
    return lip_const_new;
}

double SrtlibRlBg2SmthPm::GetQMinusF(const double* const rho_new_arr,
                                     double nu_new,
                                     const double* const rho_arr,
                                     double nu,
                                     double mu, double lip_const,
                                     int nskyx, int nskyy)
{
    double term1 = GetFuncF(rho_arr, mu, nskyx, nskyy);
    double term2 = -1 * GetFuncF(rho_new_arr, mu, nskyx, nskyy);
    int nsky = nskyx * nskyy;
    double* diff_rho_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, diff_rho_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(rho_arr), 1, diff_rho_arr, 1);

    double* diff_f_arr = new double[nsky];
    GetDiffF(rho_arr, mu, nskyx, nskyy, diff_f_arr);
    double term3 = ddot_(nsky, diff_rho_arr, 1, diff_f_arr, 1);

    double term4 = lip_const / 2.0 *
        (ddot_(nsky, diff_rho_arr, 1, diff_rho_arr, 1)
         + pow(nu_new - nu, 2));
    double ans = term1 + term2 + term3 + term4;
    delete [] diff_rho_arr;
    delete [] diff_f_arr;
    return ans;
}

double SrtlibRlBg2SmthPm::GetFuncF(const double* const rho_arr,
                                   double mu,
                                   int nskyx, int nskyy)
{
    double ans = mu * SrtlibSmth::GetTermV(rho_arr, nskyx, nskyy);
    return ans;
}

void SrtlibRlBg2SmthPm::GetDiffF(const double* const rho_arr,
                                 double mu,
                                 int nskyx, int nskyy,
                                 double* const out_arr)
{
    int nsky = nskyx * nskyy;
    SrtlibSmth::GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);
}


void SrtlibRlBg2SmthPm::GetRhoNu_ByPm(
    FILE* const fp_log,
    const double* const rho_arr, double nu,
    const double* const mval_arr, double nval,
    int nskyx, int nskyy,
    double mu,
    int npm, double tol_pm,
    int nnewton, double tol_newton,
    double* const rho_new_arr,
    double* const nu_new_ptr,
    double* const helldist_ptr,
    int* const flag_converge_ptr)
{
    int nsky = nskyx * nskyy;
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    double nu_pre = nu;
    double nu_new = 0.0;

    double lip_const = 1.0;
    double lip_const_new = 1.0;
    // lambda must be > 0, because bval must be > 0
    // at the functions, GetDerivRhoArr_FromLambda.
    double lambda = -1.0;
    double lambda_new = -1.0;
    int flag_converge = 0;
    double helldist = 0.0;    
    for(int ipm = 0; ipm < npm; ipm++){
        lip_const_new = GetFindLipConst(
            fp_log, rho_pre_arr, nu_pre,
            mval_arr, nval, mu,
            nskyx, nskyy, lip_const,
            lambda, nnewton, tol_newton);
        // printf("lip_const_new = %e\n", lip_const_new);
        double* vval_arr = new double[nsky];
        GetVvalArr(rho_pre_arr,
                   nskyx, nskyy,
                   mu, lip_const_new,
                   vval_arr);
        double wval = GetWval(nu_pre);

        SrtlibRlBg2SmthNewton::GetRhoArrNu_ByNewton(
            fp_log, vval_arr, wval,
            mval_arr, nval,
            nsky,
            lip_const_new,
            nnewton, tol_newton,
            lambda,
            rho_new_arr,
            &nu_new,
            &lambda_new);
        delete [] vval_arr;

        // set neg val to zero
        for(int isky = 0; isky < nsky; isky ++){
            if(rho_new_arr[isky] < 0.0){
                rho_new_arr[isky] = 0.0;
            }
        }
        if(nu_new < 0.0){
            nu_new = 0.0;
        }
        
        helldist = SrtlibRlStatval::GetHellingerDist(
            rho_pre_arr, nu_pre,
            rho_new_arr, nu_new, nsky);
        if (helldist < tol_pm){
            flag_converge = 1;
            MiIolib::Printf2(
                fp_log,
                "    ipm = %d, helldist = %.2e, lip_const_new = %.2e\n",
                ipm, helldist, lip_const_new);
            break;
        }
        dcopy_(nsky, rho_new_arr, 1, rho_pre_arr, 1);
        nu_pre = nu_new;
        lambda = lambda_new;
        lip_const = lip_const_new;
    }
    delete [] rho_pre_arr;
    *nu_new_ptr = nu_new;
    *helldist_ptr = helldist;
    *flag_converge_ptr = flag_converge;
}
