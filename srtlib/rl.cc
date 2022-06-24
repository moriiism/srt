#include<unistd.h>
#include "rl.h"
#include "rl_statval.h"

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


