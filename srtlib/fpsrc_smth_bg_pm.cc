#include "fpsrc_smth_bg_pm.h"
#include "fpsrc_smth_bg_newton.h"
#include "sub_smooth.h"
#include "fpsrc_smth_bg_statval.h"

void GetVvalArr(const double* const rho_arr,
                int nskyx, int nskyy,
                double mu, double lip_const,
                double* const vval_arr)
{
    int nsky = nskyx * nskyy;
    double* deriv_rho_arr = new double[nsky];
    GetDiffTermV(rho_arr, nskyx, nskyy, deriv_rho_arr);
    for(int isky = 0; isky < nsky; isky++){
        vval_arr[isky] = rho_arr[isky] - mu / lip_const * deriv_rho_arr[isky];
    }
    delete [] deriv_rho_arr;
}

void GetWvalArr(const double* const nu_arr,
                  int nsrc,
                  double* const wval_arr)
{
    for(int isrc = 0; isrc < nsrc; isrc++){
        wval_arr[isrc] = nu_arr[isrc];
    }
}

double GetZval(double phi,
               double lip_const,
               double B_val)

{
    double zval = phi + B_val / (lip_const * phi * phi);
    return zval;
}


double GetFindLipConst(FILE* const fp_log,
                       const double* const rho_arr,
                       const double* const nu_arr,
                       double phi,
                       const double* const mval_arr,
                       const double* const nval_arr,
                       double pval,
                       double phi_val,
                       int nph, double B_val,
                       double mu,
                       int nskyx, int nskyy, int nsrc,
                       double lip_const, double lambda,
                       int nnewton, double tol_newton)
{
    int nsky = nskyx * nskyy;
    int ik_max = 100;
    double eta = 1.2;

    double* rho_new_arr = new double[nsky];
    double* nu_new_arr = new double[nsrc];
    double phi_new = 0.0;
    for(int ik = 0; ik < ik_max; ik ++){
        lip_const = lip_const * pow(eta, ik);

        double* vval_arr = new double[nsky];
        GetVvalArr(rho_arr,
                   nskyx, nskyy,
                   mu, lip_const,
                   vval_arr);
        double* wval_arr = new double[nsrc];
        GetWvalArr(nu_arr,
                   nsrc,
                   wval_arr);
        double zval = GetZval(phi,
                              lip_const,
                              B_val);
        double lambda_new = 0.0;
        GetRhoArrNuArrPhi_ByNewton(fp_log,
                                   vval_arr, wval_arr, zval,
                                   mval_arr, nval_arr, pval,
                                   phi_val,
                                   nsky, nsrc, nph,
                                   lip_const,
                                   nnewton, tol_newton,
                                   lambda,
                                   rho_new_arr,
                                   nu_new_arr,
                                   &phi_new,
                                   &lambda_new);
        double qminusf = GetQMinusF(rho_new_arr, nu_new_arr, phi_new,
                                    rho_arr, nu_arr, phi,
                                    mu, lip_const, B_val,
                                    nskyx, nskyy, nsrc);
        delete [] vval_arr;
        delete [] wval_arr;        
        if(qminusf >= 0.0 && phi_new > 0.0){
            break;
        }
    }
    delete [] rho_new_arr;
    delete [] nu_new_arr;
    return lip_const;
}


double GetQMinusF(const double* const rho_new_arr,
                  const double* const nu_new_arr,
                  double phi_new,
                  const double* const rho_arr,
                  const double* const nu_arr,
                  double phi,
                  double mu, double lip_const,
                  double B_val,
                  int nskyx, int nskyy, int nsrc)
{
    double term1 = GetFuncF(rho_arr, phi, mu, B_val, nskyx, nskyy);
    double term2 = -1 * GetFuncF(rho_new_arr, phi_new, mu,
                                 B_val, nskyx, nskyy);
    int nsky = nskyx * nskyy;
    double* diff_rho_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, diff_rho_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(rho_arr), 1, diff_rho_arr, 1);

    double* diff_f_arr = new double[nsky];
    GetDiffF(rho_arr, mu, nskyx, nskyy, diff_f_arr);
    double term3 = ddot_(nsky, const_cast<double*>(diff_rho_arr), 1,
                         const_cast<double*>(diff_f_arr), 1);
    term3 += -1 * (phi_new - phi) * B_val / (phi * phi);

    double* diff_nu_arr = new double[nsrc];
    dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, diff_nu_arr, 1);
    daxpy_(nsrc, -1.0, const_cast<double*>(nu_arr), 1, diff_nu_arr, 1);
    
    double term4 = lip_const / 2.0 *
        (ddot_(nsky, const_cast<double*>(diff_rho_arr), 1,
               const_cast<double*>(diff_rho_arr), 1)
         + ddot_(nsrc, const_cast<double*>(diff_nu_arr), 1,
               const_cast<double*>(diff_nu_arr), 1)
         + pow(phi_new - phi, 2));
    double ans = term1 + term2 + term3 + term4;
    delete [] diff_rho_arr;
    delete [] diff_f_arr;
    delete [] diff_nu_arr;
    return(ans);
}


double GetFuncF(const double* const rho_arr, double phi,
                double mu, double B_val, 
                int nskyx, int nskyy)
{
    double ans = mu * GetTermV(rho_arr, nskyx, nskyy) + B_val / phi;
    return(ans);
}

void GetDiffF(const double* const rho_arr,
              double mu,
              int nskyx, int nskyy,
              double* const out_arr)
{
    int nsky = nskyx * nskyy;
    GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);
}

void GetRhoNuPhi_ByPM(FILE* const fp_log,
                      const double* const rho_arr,
                      const double* const nu_arr,
                      double phi,
                      const double* const mval_arr,
                      const double* const nval_arr,
                      double pval,
                      int nph, double B_val,
                      int ndet, int nskyx, int nskyy, int nsrc,
                      double mu,
                      int npm, double tol_pm,
                      int nnewton, double tol_newton,
                      double* const rho_new_arr,
                      double* const nu_new_arr,
                      double* const phi_new_ptr)
{
    double phi_val = phi;
    
    int nsky = nskyx * nskyy;
    double* rho_pre_arr = new double[nsky];
    double* nu_pre_arr = new double[nsrc];    
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_arr), 1, nu_pre_arr, 1);
    double phi_pre = phi;
    double phi_new = 0.0;

    double lip_const = 1.0;
    double lambda = 0.0;
    double lambda_new = 0.0;
    double lip_const_new = 1.0;
    int flag_converge = 0;
    double helldist = 0.0;
    for(int ipm = 0; ipm < npm; ipm++){
        lip_const_new = GetFindLipConst(fp_log,
                                        rho_pre_arr, nu_pre_arr, phi_pre,
                                        mval_arr, nval_arr, pval,
                                        phi_val, nph, B_val, mu,
                                        nskyx, nskyy, nsrc, lip_const,
                                        lambda, nnewton, tol_newton);
        // printf("lip_const_new = %e\n", lip_const_new);
        double* vval_arr = new double[nsky];
        GetVvalArr(rho_pre_arr,
                   nskyx, nskyy,
                   mu, lip_const_new,
                   vval_arr);
        double* wval_arr = new double[nsrc];
        GetWvalArr(nu_pre_arr,
                   nsrc,
                   wval_arr);
        double zval = GetZval(phi_pre,
                              lip_const_new,
                              B_val);
        
        GetRhoArrNuArrPhi_ByNewton(fp_log,
                                   vval_arr, wval_arr, zval,
                                   mval_arr, nval_arr, pval,
                                   phi_val,
                                   nsky, nsrc, nph,
                                   lip_const_new,
                                   nnewton, tol_newton,
                                   lambda,
                                   rho_new_arr,
                                   nu_new_arr,
                                   &phi_new,
                                   &lambda_new);
        delete [] vval_arr;
        delete [] wval_arr;
        // printf("pm out: ipm = %d, phi_new = %e\n", ipm, phi_new);
        helldist = GetHellingerDist(rho_pre_arr,
                                    nu_pre_arr,
                                    phi_pre,
                                    rho_new_arr,
                                    nu_new_arr,
                                    phi_new,
                                    nsky, nsrc);
        // printf("ipm = %d, helldist = %e\n", ipm, helldist);
        if (helldist < tol_pm){
            flag_converge = 1;
            MiIolib::Printf2(
                fp_log,
                "    ipm = %d, helldist = %.2e, lip_const_new = %.2e\n",
                ipm, helldist, lip_const_new);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
        phi_pre = phi_new;
        lambda = lambda_new;
        lip_const = lip_const_new;
    }
    if (flag_converge == 0){
        MiIolib::Printf2(
            fp_log, "    pm: not converged: helldist = %.2e\n",
            helldist);        
    }
    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
    *phi_new_ptr = phi_new;
}


// nesterov
void GetRhoNuPhi_ByPM_Nesterov(FILE* const fp_log,
                               const double* const rho_arr,
                               const double* const nu_arr,
                               double phi,
                               const double* const mval_arr,
                               const double* const nval_arr,
                               double pval,
                               int nph, double B_val,
                               int ndet, int nskyx, int nskyy, int nsrc,
                               double mu,
                               int npm, double tol_pm,
                               int nnewton, double tol_newton,
                               double* const rho_new_arr,
                               double* const nu_new_arr,
                               double* const phi_new_ptr)
{
    double phi_val = phi;
    
    int nsky = nskyx * nskyy;
    // xval_pre
    double* rho_x_pre_arr = new double[nsky];
    double* nu_x_pre_arr = new double[nsrc];    
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_x_pre_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_arr), 1, nu_x_pre_arr, 1);
    double phi_x_pre = phi;

    // xval
    double* rho_x_arr = new double[nsky];
    double* nu_x_arr = new double[nsrc];
    double phi_x = 0.0;

    // yval
    double* rho_y_arr = new double[nsky];
    double* nu_y_arr = new double[nsrc];    
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_y_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_arr), 1, nu_y_arr, 1);
    double phi_y = phi;

    // yval_new
    double* rho_y_new_arr = new double[nsky];
    double* nu_y_new_arr = new double[nsrc];    
    double phi_y_new = 0.0;

    double tval = 1.0;

    double* rho_diff_arr = new double[nsky];
    double* nu_diff_arr = new double[nsrc];
    double lip_const = 1.0;
    double lambda = 0.0;
    double lambda_new = 0.0;
    double lip_const_new = 1.0;
    int flag_converge = 0;
    double helldist = 0.0;
    for(int ipm = 0; ipm < npm; ipm++){
        lip_const_new = GetFindLipConst(fp_log,
                                        rho_y_arr, nu_y_arr, phi_y,
                                        mval_arr, nval_arr, pval,
                                        phi_val, nph, B_val, mu,
                                        nskyx, nskyy, nsrc, lip_const,
                                        lambda, nnewton, tol_newton);
        // printf("lip_const_new = %e\n", lip_const_new);
        if(lip_const_new > 1.0e+30){
            break;
        }
        double* vval_arr = new double[nsky];
        GetVvalArr(rho_y_arr,
                   nskyx, nskyy,
                   mu, lip_const_new,
                   vval_arr);
        double* wval_arr = new double[nsrc];
        GetWvalArr(nu_y_arr,
                   nsrc,
                   wval_arr);
        double zval = GetZval(phi_y,
                              lip_const_new,
                              B_val);
        
        GetRhoArrNuArrPhi_ByNewton(fp_log,
                                   vval_arr, wval_arr, zval,
                                   mval_arr, nval_arr, pval,
                                   phi_val,
                                   nsky, nsrc, nph,
                                   lip_const_new,
                                   nnewton, tol_newton,
                                   lambda,
                                   rho_x_arr,
                                   nu_x_arr,
                                   &phi_x,
                                   &lambda_new);
        delete [] vval_arr;
        delete [] wval_arr;

        double tval_new = ( 1.0 + sqrt(1.0 + 4.0 * tval * tval) ) / 2.0;

        double coeff = (tval - 1.0) / tval_new;
        MibBlas::Sub(rho_x_arr, rho_x_pre_arr, nsky, rho_diff_arr);
        MibBlas::Sub(nu_x_arr, nu_x_pre_arr, nsrc, nu_diff_arr);
        double phi_diff = phi_x - phi_x_pre;

        int nk = 100;
        double eta = 0.8;
        int ifind_nonneg = 0;        
        for(int ik = 0; ik < nk; ik ++){
            double coeff0 = coeff * pow(eta, ik);
            dcopy_(nsky, rho_x_arr, 1, rho_y_new_arr, 1);
            dcopy_(nsrc, nu_x_arr, 1, nu_y_new_arr, 1);
            daxpy_(nsky, coeff0, rho_diff_arr, 1, rho_y_new_arr, 1);
            daxpy_(nsrc, coeff0, nu_diff_arr, 1, nu_y_new_arr, 1);
            phi_y_new = phi_x + coeff0 * phi_diff;
            int nneg = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_y_new_arr[isky] < 0.0){
                    nneg ++;
                }
            }
            for(int isrc = 0; isrc < nsrc; isrc ++){
                if(nu_y_new_arr[isrc] < 0.0){
                    nneg ++;
                }
            }
            if(phi_y_new < 0.0){
                nneg ++;
            }
            if (nneg == 0){
                ifind_nonneg = 1;
                break;
            }
        }
        if(ifind_nonneg == 0){
            MiIolib::Printf2(fp_log, "warning: ipm = %d, ifind_nonneg == 0\n", ipm);
            dcopy_(nsky, rho_x_arr, 1, rho_y_new_arr, 1);
            dcopy_(nsrc, nu_x_arr, 1, nu_y_new_arr, 1);
            phi_y_new = phi_x;
            tval_new = tval;
        }
        
        // printf("pm out: ipm = %d, phi_new = %e\n", ipm, phi_new);
        helldist = GetHellingerDist(rho_y_arr,
                                    nu_y_arr,
                                    phi_y,
                                    rho_y_new_arr,
                                    nu_y_new_arr,
                                    phi_y_new,
                                    nsky, nsrc);
        // printf("ipm = %d, helldist = %e\n", ipm, helldist);
        if (helldist < tol_pm){
            flag_converge = 1;
            //MiIolib::Printf2(
            //    fp_log,
            //    "    ipm = %d, helldist = %.2e, lip_const_new = %.2e\n",
            //    ipm, helldist, lip_const_new);
            break;
        }
        dcopy_(nsky, rho_y_new_arr, 1, rho_y_arr, 1);
        dcopy_(nsrc, nu_y_new_arr, 1, nu_y_arr, 1);
        phi_y = phi_y_new;

        dcopy_(nsky, rho_x_arr, 1, rho_x_pre_arr, 1);
        dcopy_(nsrc, nu_x_arr, 1, nu_x_pre_arr, 1);
        phi_x_pre = phi_x;

        tval = tval_new;
        
        lambda = lambda_new;
        lip_const = lip_const_new;
    }
    if (flag_converge == 0){
        MiIolib::Printf2(
            fp_log, "    pm: not converged: helldist = %.2e\n",
            helldist);        
    }

    dcopy_(nsky, rho_y_new_arr, 1, rho_new_arr, 1);
    dcopy_(nsrc, nu_y_new_arr, 1, nu_new_arr, 1);
    double phi_new = phi_y_new;

    delete [] rho_x_pre_arr;
    delete [] nu_x_pre_arr;
    delete [] rho_x_arr;
    delete [] nu_x_arr;
    delete [] rho_y_arr;
    delete [] nu_y_arr;
    delete [] rho_y_new_arr;
    delete [] nu_y_new_arr;
    delete [] rho_diff_arr;
    delete [] nu_diff_arr;
    
    *phi_new_ptr = phi_new;
}
