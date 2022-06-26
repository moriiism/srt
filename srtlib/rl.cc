#include<unistd.h>
#include "rl.h"
#include "rl_statval.h"

void GetInvVec(const double* const vec_arr, int nelm,
               double* const inv_arr)
{
    dcopy_(nelm, const_cast<double*>(vec_arr), 1, inv_arr, 1);
    double norm2 = ddot_(nelm, const_cast<double*>(vec_arr), 1,
                         const_cast<double*>(vec_arr), 1);
    dscal_(nelm, 1.0 / norm2, inv_arr, 1);
}


// accerelated richardson lucy
void RichlucyAcc(FILE* const fp_log,
                 const double* const rho_init_arr,
                 const double* const data_arr,
                 const double* const resp_norm_mat_arr,
                 int ndet, int nsky,
                 string outdir, string outfile_head,
                 int nem, double tol_em,
                 int k_restart, double delta_restart,
                 double* const rho_new_arr)
{
    if (delta_restart <= tol_em){
        printf("delta_restart(=%e) <= tol_em(=%e), then abort.\n",
               delta_restart, tol_em);
        abort();
    }
    double* rho_0_arr  = new double[nsky];
    double* rho_1_arr  = new double[nsky];
    double* rho_2_arr  = new double[nsky];
    double* rho_dot_old_arr = new double[nsky];
    double* rho_dot_new_arr = new double[nsky];

    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_0_arr, 1);
    GetRhoNewArr(rho_0_arr, data_arr, resp_norm_mat_arr,
                 ndet, nsky, rho_1_arr);
    dcopy_(nsky, rho_1_arr, 1, rho_dot_old_arr, 1);
    
    for(int iem = 0; iem < nem; iem ++){
        GetRhoNewArr(rho_1_arr, data_arr, resp_norm_mat_arr,
                     ndet, nsky, rho_2_arr);

        // epsilon-acceleration step
        double* delta_rho_0_arr = new double[nsky];
        double* delta_rho_1_arr = new double[nsky];
        MibBlas::Sub(rho_1_arr, rho_0_arr, nsky, delta_rho_0_arr);
        MibBlas::Sub(rho_2_arr, rho_1_arr, nsky, delta_rho_1_arr);
        double* delta_rho_0_inv_arr = new double[nsky];
        double* delta_rho_1_inv_arr = new double[nsky];        
        GetInvVec(delta_rho_0_arr, nsky, delta_rho_0_inv_arr);
        GetInvVec(delta_rho_1_arr, nsky, delta_rho_1_inv_arr);

        double* delta_delta_rho_01_inv_arr = new double[nsky];
        MibBlas::Sub(delta_rho_1_inv_arr, delta_rho_0_inv_arr, nsky,
                     delta_delta_rho_01_inv_arr);
        double* delta_delta_rho_01_inv_inv_arr = new double[nsky];
        GetInvVec(delta_delta_rho_01_inv_arr, nsky,
                  delta_delta_rho_01_inv_inv_arr);
        MibBlas::Add(rho_1_arr, delta_delta_rho_01_inv_inv_arr, nsky,
                     rho_dot_new_arr);
        delete [] delta_delta_rho_01_inv_arr;
        delete [] delta_delta_rho_01_inv_inv_arr;
        delete [] delta_rho_0_inv_arr;
        delete [] delta_rho_1_inv_arr;
        delete [] delta_rho_0_arr;
        delete [] delta_rho_1_arr;

        // re-starting procedure
        double* delta_rho_dot_arr = new double[nsky];
        MibBlas::Sub(rho_dot_new_arr, rho_dot_old_arr, nsky,
                     delta_rho_dot_arr);
        double diff2 = ddot_(nsky, delta_rho_dot_arr, 1,
                             delta_rho_dot_arr, 1);
        delete [] delta_rho_dot_arr;
        MiIolib::Printf2(fp_log, "iem = %d, diff2 = %.2e\n",
                         iem, diff2);
        if (diff2 < delta_restart){
            if (diff2 < tol_em){
                MiIolib::Printf2(fp_log, "iem = %d, diff2 = %.2e\n",
                                 iem, diff2);
                break;
            }
            double* rho_tmp_arr = new double[nsky];
            GetRhoNewArr(rho_dot_new_arr, data_arr, resp_norm_mat_arr,
                         ndet, nsky, rho_tmp_arr);
            double lval_tmp = -1.0 * GetNegLogLike(rho_tmp_arr, data_arr,
                                                   resp_norm_mat_arr,
                                                   ndet, nsky, 1.0e-10);
            double lval_2 = -1.0 * GetNegLogLike(rho_2_arr, data_arr,
                                                 resp_norm_mat_arr,
                                                 ndet, nsky, 1.0e-10);
            if (lval_tmp > lval_2){
                dcopy_(nsky, rho_tmp_arr, 1, rho_2_arr, 1);
                dcopy_(nsky, rho_dot_new_arr, 1, rho_1_arr, 1);
                delta_restart *= pow(10.0, -1.0 * k_restart);
                printf("delta_restart = %e\n", delta_restart);
            }
            delete [] rho_tmp_arr;
        }
        dcopy_(nsky, rho_dot_new_arr, 1, rho_dot_old_arr, 1);
        dcopy_(nsky, rho_1_arr, 1, rho_0_arr, 1);
        dcopy_(nsky, rho_2_arr, 1, rho_1_arr, 1);
       
        if (access( "/tmp/rl_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_stop file is found, then stop.\n");
            break;
        }
       
    }
    dcopy_(nsky, rho_dot_new_arr, 1, rho_new_arr, 1);
    delete [] rho_0_arr;
    delete [] rho_1_arr;
    delete [] rho_2_arr;
    delete [] rho_dot_old_arr;
    delete [] rho_dot_new_arr;
}


void Richlucy(FILE* const fp_log,
              const double* const rho_init_arr,
              const double* const data_arr,
              const double* const resp_norm_mat_arr,
              int ndet, int nsky,
              string outdir, string outfile_head,
              int nem, double tol_em,
              double* const rho_new_arr)
{
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    for(int iem = 0; iem < nem; iem ++){
        GetRhoNewArr(rho_pre_arr, data_arr, resp_norm_mat_arr,
                     ndet, nsky, rho_new_arr);
        double helldist  = GetHellingerDist(rho_pre_arr, rho_new_arr, nsky);
        if (access( "/tmp/rl_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_stop file is found, then stop.\n");
            break;
        }
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
            iem, helldist);
        
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
    }
    delete [] rho_pre_arr;
}

void GetRhoNewArr(const double* const rho_arr,
                  const double* const data_arr,
                  const double* const resp_norm_mat_arr,
                  int ndet, int nsky,
                  double* const rho_new_arr)
{
    double* den_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet ++){
        den_arr[idet] = 0.0;
    }
    GetDetArr(rho_arr, resp_norm_mat_arr, ndet, nsky, den_arr);

    double* div_arr = new double[ndet];
    double sum = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / den_arr[idet];
        sum += data_arr[idet];
    }

    double* tmp_arr = new double[nsky];
    char transa[1];
    strcpy(transa, "T");
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, tmp_arr, 1);
    for(int isky = 0; isky < nsky; isky ++){
        rho_new_arr[isky] = tmp_arr[isky] * rho_arr[isky];
        rho_new_arr[isky] /= sum;
    }
    delete [] den_arr;
    delete [] div_arr;
    delete [] tmp_arr;
}

void GetDetArr(const double* const rho_arr,
               const double* const resp_norm_mat_arr,
               int ndet, int nsky,
               double* const det_arr) // ndet
{
    // det_arr = R_mat %*% rho_arr
    char transa[1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);
}


