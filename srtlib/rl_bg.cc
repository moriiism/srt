#include<unistd.h>
#include "rl.h"
#include "rl_bg.h"
#include "rl_statval.h"
#include "sim.h"

double GetAlpha(const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                const double* const bg_arr,
                const double* const data_arr,
                int nsky, int ndet)
{
    double* tmp_arr = new double[ndet];
    GetDetArr(rho_arr, resp_norm_mat_arr, ndet, nsky, tmp_arr);
    double B_val = MibBlas::Sum(bg_arr, ndet);
    double* num_arr = new double[ndet];
    MibBlas::ElmWiseMul(ndet, 1.0, data_arr, tmp_arr, num_arr);
    daxpy_(ndet, nu / B_val, const_cast<double*>(bg_arr), 1, tmp_arr, 1);
    double alpha = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        alpha += num_arr[idet] / tmp_arr[idet];
    }
    delete [] tmp_arr;
    delete [] num_arr;
    return alpha;
}

void GetRhoNu_New(const double* const rho_arr, double nu,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const bg_arr,
                  int ndet, int nsky,
                  double* const rho_new_arr,
                  double* const nu_new_ptr)
{

    double B_val = MibBlas::Sum(bg_arr, ndet);
    
    // rho_new_arr[isky] = 1./(alpha + B) * rho_arr[isky]
    //                    * [t(R_mat) %*% (data_arr / den_arr)][isky]
    // new_new = B / (alpha + B)
    double* den_arr = new double[ndet];
    GetDetArr(rho_arr, resp_mat_arr, ndet, nsky, den_arr);
    daxpy_(ndet, nu / B_val, const_cast<double*>(bg_arr), 1, den_arr, 1);
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / den_arr[idet];
    }
    double* tmp_arr = new double[nsky];
    char transa[1];
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, tmp_arr, 1);
    double alpha = GetAlpha(rho_arr, nu,
                            resp_mat_arr,
                            bg_arr,
                            data_arr,
                            nsky, ndet);
    MibBlas::ElmWiseMul(nsky, 1.0 / (alpha + B_val),
                        tmp_arr, rho_arr, rho_new_arr);
    double nu_new = B_val / (alpha + B_val);

    delete [] den_arr;
    delete [] div_arr;
    delete [] tmp_arr;
    *nu_new_ptr = nu_new;
}


void RichlucyBg(FILE* const fp_log,
                const double* const rho_init_arr,
                double nu_init,
                const double* const data_arr,
                const double* const bg_arr,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky,
                string outdir,
                string outfile_head,
                int nem,                
                double tol_em,
                double* const rho_new_arr,
                double* const nu_new_ptr)
{
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    double nu_pre = nu_init;
    double nu_new = 0.0;
    for(int iem = 0; iem < nem; iem ++){
        GetRhoNu_New(rho_pre_arr, nu_pre,
                     data_arr,
                     resp_norm_mat_arr,
                     bg_arr,
                     ndet, nsky,
                     rho_new_arr,
                     &nu_new);
        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre,
                                            rho_new_arr, nu_new, nsky);
        if (access( "/tmp/rl_bg_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_bg_stop file is found, then stop.\n");
            break;
        }
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
            iem, helldist);
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
            break;
        }
        dcopy_(nsky, rho_new_arr, 1, rho_pre_arr, 1);
        nu_pre = nu_new;
    }
    delete [] rho_pre_arr;
    *nu_new_ptr = nu_new;
}


void RichlucyBgAccSQUAREM(FILE* const fp_log,
                          const double* const rho_init_arr,
                          double nu_init,
                          const double* const data_arr,
                          const double* const bg_arr,
                          const double* const resp_norm_mat_arr,
                          int ndet, int nsky,
                          string outdir,
                          string outfile_head,
                          int nem,                
                          double tol_em,
                          double* const rho_new_arr,
                          double* const nu_new_ptr)
{
    double* rho_0_arr  = new double[nsky];
    double* rho_1_arr  = new double[nsky];
    double* rho_2_arr  = new double[nsky];
    double* rho_dash_arr  = new double[nsky];
    double* r_rho_arr  = new double[nsky];
    double* r2_rho_arr  = new double[nsky];
    double* v_rho_arr  = new double[nsky];
    double* rho_0_new_arr  = new double[nsky];

    double nu_0 = 0.0;
    double nu_1 = 0.0;
    double nu_2 = 0.0;
    double nu_dash = 0.0;
    double r_nu = 0.0;
    double r2_nu = 0.0;
    double v_nu = 0.0;
    double nu_0_new = 0.0;

    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_0_arr, 1);
    nu_0 = nu_init;
    for(int iem = 0; iem < nem; iem ++){
        GetRhoNu_New(rho_0_arr, nu_0,
                     data_arr,
                     resp_norm_mat_arr,
                     bg_arr,
                     ndet, nsky,
                     rho_1_arr,
                     &nu_1);
        GetRhoNu_New(rho_1_arr, nu_1,
                     data_arr,
                     resp_norm_mat_arr,
                     bg_arr,
                     ndet, nsky,
                     rho_2_arr,
                     &nu_2);
        MibBlas::Sub(rho_1_arr, rho_0_arr, nsky, r_rho_arr);
        MibBlas::Sub(rho_2_arr, rho_1_arr, nsky, r2_rho_arr);
        MibBlas::Sub(r2_rho_arr, r_rho_arr, nsky, v_rho_arr);
        r_nu = nu_1 - nu_0;
        r2_nu = nu_2 - nu_1;
        v_nu = r2_nu - r_nu;

        double r_norm2 = ddot_(nsky, r_rho_arr, 1, r_rho_arr, 1)
            + r_nu * r_nu;
        double v_norm2 = ddot_(nsky, v_rho_arr, 1, v_rho_arr, 1)
            + v_nu * v_nu;
        double alpha = -1.0 * sqrt(r_norm2 / v_norm2);

        int nk = 100;
        double eta = 0.8;
        int ifind_nonneg = 0;
        for (int ik = 0; ik < nk; ik ++){
            double alpha0 = alpha * pow(eta, ik);
            dcopy_(nsky, rho_0_arr, 1, rho_dash_arr, 1);
            daxpy_(nsky, -2 * alpha0, r_rho_arr, 1, rho_dash_arr, 1);
            daxpy_(nsky, alpha0 * alpha0, v_rho_arr, 1, rho_dash_arr, 1);
            nu_dash = nu_0;
            nu_dash += -2 * alpha0 * r_nu;
            nu_dash += alpha0 * alpha0 * v_nu;
            
            int nneg_tmp = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_dash_arr[isky] < 0.0){
                    nneg_tmp ++;
                }
            }
            if(nu_dash < 0.0){
                nneg_tmp ++;
            }
            if (nneg_tmp > 0){
                continue;
            }
            GetRhoNu_New(rho_dash_arr, nu_dash,
                         data_arr,
                         resp_norm_mat_arr,
                         bg_arr,
                         ndet, nsky,
                         rho_0_new_arr,
                         &nu_0_new);
            int nneg = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_0_new_arr[isky] < 0.0){
                    nneg ++;
                }
            }
            if(nu_0_new < 0.0){
                nneg ++;
            }
            if (nneg == 0){
                ifind_nonneg = 1;
                break;
            }
        }
        if(ifind_nonneg == 0){
            MiIolib::Printf2(fp_log, "warning: iem = %d, ifind_nonneg == 0\n", iem);
        }

        double helldist  = GetHellingerDist(rho_0_arr, nu_0,
                                            rho_0_new_arr, nu_0_new, nsky);
        if (access( "/tmp/rl_bg_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_bg_stop file is found, then stop.\n");
            break;
        }
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
            iem, helldist);
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
            break;
        }
        dcopy_(nsky, rho_0_new_arr, 1, rho_0_arr, 1);
        nu_0 = nu_0_new;
    }
    dcopy_(nsky, rho_0_new_arr, 1, rho_new_arr, 1);
    double nu_new = nu_0_new;

    delete [] rho_0_arr;
    delete [] rho_1_arr;
    delete [] rho_2_arr;
    delete [] rho_dash_arr;
    delete [] r_rho_arr;
    delete [] r2_rho_arr;
    delete [] v_rho_arr;
    delete [] rho_0_new_arr;
    
    *nu_new_ptr = nu_new;
}

