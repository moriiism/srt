#include "fpsrc_smth_bg_dc.h"
#include "fpsrc_smth_bg_pm.h"
#include "fpsrc_smth_bg_statval.h"

// accerelation by SQUAREM
void GetRhoNuPhi_ByDC_Acc(FILE* const fp_log,
                          const double* const rho_arr,
                          const double* const nu_arr,
                          double phi,
                          const double* const mval_arr,
                          const double* const nval_arr,
                          double pval,
                          int nph, double B_val,
                          int ndet, int nskyx, int nskyy, int nsrc,
                          double mu,
                          int ndc, double tol_dc,
                          int npm, double tol_pm,
                          int nnewton, double tol_newton,
                          double* const rho_new_arr,
                          double* const nu_new_arr,
                          double* const phi_new_ptr)
{
    int nsky = nskyx * nskyy;

    double* rho_0_arr  = new double[nsky];
    double* rho_1_arr  = new double[nsky];
    double* rho_2_arr  = new double[nsky];
    double* rho_dash_arr  = new double[nsky];
    double* r_rho_arr  = new double[nsky];
    double* r2_rho_arr  = new double[nsky];
    double* v_rho_arr  = new double[nsky];
    double* rho_0_new_arr  = new double[nsky];

    double* nu_0_arr  = new double[nsrc];
    double* nu_1_arr  = new double[nsrc];
    double* nu_2_arr  = new double[nsrc];
    double* nu_dash_arr  = new double[nsrc];
    double* r_nu_arr  = new double[nsrc];
    double* r2_nu_arr  = new double[nsrc];
    double* v_nu_arr  = new double[nsrc];
    double* nu_0_new_arr  = new double[nsrc];

    double phi_0 = 0.0;
    double phi_1 = 0.0;
    double phi_2 = 0.0;
    double phi_dash = 0.0;
    double r_phi = 0.0;
    double r2_phi = 0.0;
    double v_phi = 0.0;
    double phi_0_new = 0.0;
    
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_0_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_arr), 1, nu_0_arr, 1);
    phi_0 = phi;
    
    int flag_converge = 0;
    double helldist  = 0.0;
    for(int idc = 0; idc < ndc; idc++){
        GetRhoNuPhi_ByPM_Nesterov(fp_log,
                                  rho_0_arr, nu_0_arr, phi_0,
                                  mval_arr, nval_arr, pval,
                                  nph, B_val,
                                  ndet, nskyx, nskyy, nsrc,
                                  mu,
                                  npm, tol_pm,
                                  nnewton, tol_newton,
                                  rho_1_arr,
                                  nu_1_arr,
                                  &phi_1);
        GetRhoNuPhi_ByPM_Nesterov(fp_log,
                                  rho_1_arr, nu_1_arr, phi_1,
                                  mval_arr, nval_arr, pval,
                                  nph, B_val,
                                  ndet, nskyx, nskyy, nsrc,
                                  mu,
                                  npm, tol_pm,
                                  nnewton, tol_newton,
                                  rho_2_arr,
                                  nu_2_arr,
                                  &phi_2);
        MibBlas::Sub(rho_1_arr, rho_0_arr, nsky, r_rho_arr);
        MibBlas::Sub(rho_2_arr, rho_1_arr, nsky, r2_rho_arr);
        MibBlas::Sub(r2_rho_arr, r_rho_arr, nsky, v_rho_arr);

        MibBlas::Sub(nu_1_arr, nu_0_arr, nsrc, r_nu_arr);
        MibBlas::Sub(nu_2_arr, nu_1_arr, nsrc, r2_nu_arr);
        MibBlas::Sub(r2_nu_arr, r_nu_arr, nsrc, v_nu_arr);

        r_phi = phi_1 - phi_0;
        r2_phi = phi_2 - phi_1;
        v_phi = r2_phi - r_phi;

        double r_norm2 = ddot_(nsky, r_rho_arr, 1, r_rho_arr, 1)
            + ddot_(nsrc, r_nu_arr, 1, r_nu_arr, 1)
            + r_phi * r_phi;
        double v_norm2 = ddot_(nsky, v_rho_arr, 1, v_rho_arr, 1)
            + ddot_(nsrc, v_nu_arr, 1, v_nu_arr, 1)
            + v_phi * v_phi;
        double alpha = -1.0 * sqrt(r_norm2 / v_norm2);

        int nk = 100;
        double eta = 0.8;
        int ifind_nonneg = 0;
        for (int ik = 0; ik < nk; ik ++){
            double alpha0 = alpha * pow(eta, ik);
            dcopy_(nsky, rho_0_arr, 1, rho_dash_arr, 1);
            dcopy_(nsrc, nu_0_arr, 1, nu_dash_arr, 1);
            phi_dash = phi_0;
            daxpy_(nsky, -2 * alpha0, r_rho_arr, 1, rho_dash_arr, 1);
            daxpy_(nsrc, -2 * alpha0, r_nu_arr, 1, nu_dash_arr, 1);
            phi_dash += -2 * alpha0 * r_phi;
            daxpy_(nsky, alpha0 * alpha0, v_rho_arr, 1, rho_dash_arr, 1);
            daxpy_(nsrc, alpha0 * alpha0, v_nu_arr, 1, nu_dash_arr, 1);
            phi_dash += alpha0 * alpha0 * v_phi;

            int nneg_tmp = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_dash_arr[isky] < 0.0){
                    nneg_tmp ++;
                }
            }
            for(int isrc = 0; isrc < nsrc; isrc ++){
                if(nu_dash_arr[isrc] < 0.0){
                    nneg_tmp ++;
                }
            }
            if(phi_dash < 0.0){
                nneg_tmp ++;
            }
            if (nneg_tmp > 0){
                continue;
            }

            GetRhoNuPhi_ByPM_Nesterov(fp_log,
                                      rho_dash_arr, nu_dash_arr, phi_dash,
                                      mval_arr, nval_arr, pval,
                                      nph, B_val,
                                      ndet, nskyx, nskyy, nsrc,
                                      mu,
                                      npm, tol_pm,
                                      nnewton, tol_newton,
                                      rho_0_new_arr,
                                      nu_0_new_arr,
                                      &phi_0_new);
            int nneg = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_0_new_arr[isky] < 0.0){
                    nneg ++;
                }
            }
            for(int isrc = 0; isrc < nsrc; isrc ++){
                if(nu_0_new_arr[isrc] < 0.0){
                    nneg ++;
                }
            }
            if(phi_0_new < 0.0){
                nneg ++;
            }
            if (nneg == 0){
                ifind_nonneg = 1;
                break;
            }
        }
        if(ifind_nonneg == 0){
            MiIolib::Printf2(fp_log, "warning: idc = %d, ifind_nonneg == 0\n",
                             idc);
        }
        
        
        helldist  = GetHellingerDist(rho_0_arr, nu_0_arr, phi_0,
                                     rho_0_new_arr, nu_0_new_arr, phi_0_new,
                                     nsky, nsrc);
        if (helldist < tol_dc){
            flag_converge = 1;
            MiIolib::Printf2(fp_log, "  idc = %d, helldist = %.2e\n",
                             idc, helldist);
            break;
        }
        dcopy_(nsky, rho_0_new_arr, 1, rho_0_arr, 1);
        dcopy_(nsrc, nu_0_new_arr, 1, nu_0_arr, 1);
        phi_0 = phi_0_new;
    }
    if (flag_converge == 0){
        MiIolib::Printf2(fp_log, "  dc: not converged: helldist = %.2e\n",
                         helldist);
    }

    dcopy_(nsky, rho_0_new_arr, 1, rho_new_arr, 1);
    dcopy_(nsrc, nu_0_new_arr, 1, nu_new_arr, 1);
    double phi_new = phi_0_new;
    
    delete [] rho_0_arr;
    delete [] rho_1_arr;
    delete [] rho_2_arr;
    delete [] rho_dash_arr;
    delete [] r_rho_arr;
    delete [] r2_rho_arr;
    delete [] v_rho_arr;
    delete [] rho_0_new_arr;

    delete [] nu_0_arr;
    delete [] nu_1_arr;
    delete [] nu_2_arr;
    delete [] nu_dash_arr;
    delete [] r_nu_arr;
    delete [] r2_nu_arr;
    delete [] v_nu_arr;
    delete [] nu_0_new_arr;


    *phi_new_ptr = phi_new;
}


void GetRhoNuPhi_ByDC(FILE* const fp_log,
                      const double* const rho_arr,
                      const double* const nu_arr,
                      double phi,
                      const double* const mval_arr,
                      const double* const nval_arr,
                      double pval,
                      int nph, double B_val,
                      int ndet, int nskyx, int nskyy, int nsrc,
                      double mu,
                      int ndc, double tol_dc,
                      int npm, double tol_pm,
                      int nnewton, double tol_newton,
                      double* const rho_new_arr,
                      double* const nu_new_arr,
                      double* const phi_new_ptr)
{
    int nsky = nskyx * nskyy;    
    double* rho_pre_arr = new double[nsky];
    double* nu_pre_arr = new double[nsrc];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_arr), 1, nu_pre_arr, 1);
    double phi_pre = phi;
    double phi_new = 0.0;
    int flag_converge = 0;
    double helldist  = 0.0;
    for(int idc = 0; idc < ndc; idc++){
        GetRhoNuPhi_ByPM(fp_log,
                         rho_pre_arr, nu_pre_arr, phi_pre,
                         mval_arr, nval_arr, pval,
                         nph, B_val,
                         ndet, nskyx, nskyy, nsrc,
                         mu,
                         npm, tol_pm,
                         nnewton, tol_newton,
                         rho_new_arr,
                         nu_new_arr,
                         &phi_new);
        helldist  = GetHellingerDist(rho_pre_arr, nu_pre_arr, phi_pre,
                                     rho_new_arr, nu_new_arr, phi_new,
                                     nsky, nsrc);
        if (helldist < tol_dc){
            flag_converge = 1;
            MiIolib::Printf2(fp_log, "  idc = %d, helldist = %.2e\n",
                             idc, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
        phi_pre = phi_new;
    }
    if (flag_converge == 0){
        MiIolib::Printf2(fp_log, "  dc: not converged: helldist = %.2e\n",
                         helldist);
    }
    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
    *phi_new_ptr = phi_new;
}
