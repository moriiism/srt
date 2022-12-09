#include "smth.h"
#include "rl_statval.h"
#include "rl_crab_smth_pf_pm.h"
#include "rl_crab_smth_pf_newton.h"

void SrtlibRlCrabSmthPfPm::GetVvalArr(
    const double* const rho_arr,
    int nskyx, int nskyy,
    double mu, double lip_const,
    double* const vval_arr)
{
    int nsky = nskyx * nskyy;
    double* deriv_rho_arr = new double[nsky];
    SrtlibSmth::GetDiffTermV(rho_arr, nskyx, nskyy, deriv_rho_arr);
    for(int isky = 0; isky < nsky; isky++){
        vval_arr[isky] = rho_arr[isky]
            - mu / lip_const * deriv_rho_arr[isky];
    }
    delete [] deriv_rho_arr;
}

void SrtlibRlCrabSmthPfPm::GetWvalArr(
    const double* const nu_arr,
    const double* const nu_0_arr,
    int nphase,
    double gamma, double lip_const,
    double* const wval_arr)
{
    double* deriv_nu_arr = new double[nphase];
    for(int iphase = 0; iphase < nphase; iphase++){
        deriv_nu_arr[iphase] = 2.0 * (nu_arr[iphase] - nu_0_arr[iphase]);
    }
    for(int iphase = 0; iphase < nphase; iphase++){
        wval_arr[iphase] = nu_arr[iphase]
            - gamma / lip_const * deriv_nu_arr[iphase];
    }
    delete [] deriv_nu_arr;
}


double SrtlibRlCrabSmthPfPm::GetFindLipConst(
    FILE* const fp_log,
    const double* const rho_arr,
    const double* const nu_arr,
    const double* const mval_arr,
    const double* const nval_arr,
    const double* const nu_0_arr,
    double mu, double gamma,
    int nskyx, int nskyy, int nphase,
    double lip_const, double lambda,
    int nnewton, double tol_newton)
{
    int nsky = nskyx * nskyy;
    int ik_max = 1000;
    double eta = 1.2;

    double* rho_new_arr = new double[nsky];
    double* nu_new_arr = new double[nphase];
    double lip_const_new = 1.0;
    int flag_find_lip_const = 0;
    for(int ik = 0; ik < ik_max; ik ++){
        // printf("loop in findlipconst: ik = %d\n", ik);
        lip_const_new = lip_const * pow(eta, ik);

        double* vval_arr = new double[nsky];
        double* wval_arr = new double[nphase];
        GetVvalArr(rho_arr,
                   nskyx, nskyy,
                   mu, lip_const_new,
                   vval_arr);
        GetWvalArr(nu_arr,
                   nu_0_arr,
                   nphase,
                   gamma, lip_const_new,
                   wval_arr);        
        double lambda_new = 0.0;
        SrtlibRlCrabSmthPfNewton::GetRhoArrNuArr_ByNewton(
            fp_log,
            vval_arr, wval_arr,
            mval_arr, nval_arr,
            nsky, nphase,
            lip_const_new,
            nnewton, tol_newton,
            lambda,
            rho_new_arr,
            nu_new_arr,
            &lambda_new);
        double qminusf = GetQMinusF(rho_new_arr, nu_new_arr,
                                    rho_arr, nu_arr,
                                    mu, gamma, lip_const_new,
                                    nskyx, nskyy, nphase);
        delete [] vval_arr;
        delete [] wval_arr;
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
    delete [] nu_new_arr;
    return lip_const_new;
}

double SrtlibRlCrabSmthPfPm::GetQMinusF(
    const double* const rho_new_arr,
    const double* const nu_new_arr,
    const double* const rho_arr,
    const double* const nu_arr,
    double mu, double gamma, double lip_const,
    int nskyx, int nskyy, int nphase)
{
    double term1 = GetFuncF1(rho_arr, mu, nskyx, nskyy)
        + GetFuncF2(nu_arr, nu_0_arr, gamma, nphase);
    double term2 = -1 * (GetFuncF1(rho_new_arr, mu, nskyx, nskyy)
                         +
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

double SrtlibRlCrabSmthPfPm::GetFuncF1(
    const double* const rho_arr,
    double mu, int nskyx, int nskyy)
{
    double ans = mu * SrtlibSmth::GetTermV(rho_arr, nskyx, nskyy);
    return ans;
}

double SrtlibRlCrabSmthPfPm::GetFuncF2(
    const double* const nu_arr,
    const double* const nu_0_arr,
    double gamma, int nphase)
{
    double ans = gamma * SrtlibRlCrabSmthPfPm::GetTermD(
        nu_arr, nu_0_arr, nphase);
    return ans;
}

void SrtlibRlCrabSmthPfPm::GetDiffF1(
    const double* const rho_arr,
    double mu,
    int nskyx, int nskyy,
    double* const out_arr)
{
    int nsky = nskyx * nskyy;
    SrtlibSmth::GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);
}

void SrtlibRlCrabSmthPfPm::GetDiffF2(
    const double* const nu_arr,
    const double* const nu_0_arr,
    double gamma, int nphase,
    double* const out_arr)
{
    SrtlibRlCrabSmthPfPm::GetDiffTermD(
        nu_arr, nu_0_arr, nphase, out_arr);
    dscal_(nphase, gamma, out_arr, 1);
}

double SrtlibRlCrabSmthPfPm::GetTermD(
    const double* const nu_arr,
    const double* const nu_0_arr,
    int nphase)
{
    double termd = 0.0;
    for(int iphase = 0; iphase < nphase; iphase++){
        double diff = nu_arr[iphase] - nu_0_arr[iphase];
        termd += diff * diff;
    }
    return termd;
}

void SrtlibRlCrabSmthPfPm::GetDiffTermD(
    const double* const nu_arr,
    const double* const nu_0_arr,
    int nphase,
    double* const termd_diff_arr)
{
    for(int iphase = 0; iphase < nphase; iphase++){
        termd_diff_arr[iphase] = 2.0 * (nu_arr[iphase] - nu_0_arr[iphase]);
    }
}


void SrtlibRlCrabSmthPfPm::GetRhoNu_ByPm(
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
    // lambda must be < 0, because bval must be < 0
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

        // output of newton method sometimes become
        // very slightly negative value. So,
        // set negative value to zero.
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


void SrtlibRlCrabSmthPfPm::GetRhoNu_ByPm_Nesterov(
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

    // xval_pre
    double* rho_x_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_x_pre_arr, 1);
    double nu_x_pre = nu;
    // xval
    double* rho_x_arr = new double[nsky];
    double nu_x = 0.0;
    // yval
    double* rho_y_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_y_arr, 1);
    double nu_y = nu;
    // yval_new
    double* rho_y_new_arr = new double[nsky];
    double nu_y_new = 0.0;

    double tval = 1.0;

    double* rho_diff_arr = new double[nsky];
    double lip_const = 1.0;
    double lip_const_new = 1.0;

    // lambda must be < 0, because bval must be < 0
    // at the functions, GetDerivRhoArr_FromLambda.
    double lambda = -1.0;
    double lambda_new = -1.0;
    int flag_converge = 0;
    double helldist = 0.0;    
    for(int ipm = 0; ipm < npm; ipm++){
        lip_const_new = GetFindLipConst(
            fp_log, rho_y_arr, nu_y,
            mval_arr, nval, mu,
            nskyx, nskyy, lip_const,
            lambda, nnewton, tol_newton);
        // printf("lip_const_new = %e\n", lip_const_new);
        double* vval_arr = new double[nsky];
        GetVvalArr(rho_y_arr,
                   nskyx, nskyy,
                   mu, lip_const_new,
                   vval_arr);
        double wval = GetWval(nu_y);

        SrtlibRlBg2SmthNewton::GetRhoArrNu_ByNewton(
            fp_log, vval_arr, wval,
            mval_arr, nval,
            nsky,
            lip_const_new,
            nnewton, tol_newton,
            lambda,
            rho_x_arr,
            &nu_x,
            &lambda_new);
        delete [] vval_arr;

        double tval_new = ( 1.0 + sqrt(1.0 + 4.0 * tval * tval) ) / 2.0;
        // printf("ipm = %d, tval_new = %e, phi_x = %e\n", ipm, tval_new, phi_x);
        // something is wrong.

        double coeff = (tval - 1.0) / tval_new;
        MibBlas::Sub(rho_x_arr, rho_x_pre_arr, nsky, rho_diff_arr);
        double nu_diff = nu_x - nu_x_pre;

        int nk = 100;
        double eta = 0.8;
        int ifind_nonneg = 0;
        int nneg = 0;
        for(int ik = 0; ik < nk; ik ++){
            double coeff0 = coeff * pow(eta, ik);
            dcopy_(nsky, rho_x_arr, 1, rho_y_new_arr, 1);
            daxpy_(nsky, coeff0, rho_diff_arr, 1, rho_y_new_arr, 1);
            nu_y_new = nu_x + coeff0 * nu_diff;
            nneg = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_y_new_arr[isky] < 0.0){
                    nneg ++;
                }
            }
            if(nu_y_new < 0.0){
                nneg ++;
            }
            if (nneg == 0){
                ifind_nonneg = 1;
                break;
            }
        }
        if(ifind_nonneg == 0){
            // set neg val to zero
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_y_new_arr[isky] < 0.0){
                    rho_y_new_arr[isky] = 0.0;
                }
            }
            if(nu_y_new < 0.0){
                nu_y_new = 0.0;
            }
        }
        
        helldist = SrtlibRlStatval::GetHellingerDist(
            rho_y_arr, nu_y,
            rho_y_new_arr, nu_y_new, nsky);
        if (helldist < tol_pm){
            flag_converge = 1;
            MiIolib::Printf2(
                fp_log,
                "    ipm = %d, helldist = %.2e, lip_const_new = %.2e\n",
                ipm, helldist, lip_const_new);
            break;
        }
        dcopy_(nsky, rho_y_new_arr, 1, rho_y_arr, 1);
        nu_y = nu_y_new;

        dcopy_(nsky, rho_x_arr, 1, rho_x_pre_arr, 1);
        nu_x_pre = nu_x;

        tval = tval_new;
        
        lambda = lambda_new;
        lip_const = lip_const_new;
    }

    dcopy_(nsky, rho_y_new_arr, 1, rho_new_arr, 1);
    double nu_new = nu_y_new;

    delete [] rho_x_pre_arr;
    delete [] rho_x_arr;
    delete [] rho_y_arr;
    delete [] rho_y_new_arr;
    delete [] rho_diff_arr;
    
    *nu_new_ptr = nu_new;
    *helldist_ptr = helldist;
    *flag_converge_ptr = flag_converge;
}
