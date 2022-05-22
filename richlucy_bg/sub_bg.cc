#include "sub_bg.h"


double GetB(const double* const bg_arr, int ndet)
{
    double B = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        B += bg_arr[idet];
    }
    return(B);
}

void GetDetArr(const double* const rho_arr,
               const double* const resp_mat_arr,
               int ndet, int nsky,
               double* const out_arr) // ndet
{
    // det_arr = R_mat %*% rho_arr
    char* transa = new char [1];
    strcpy(transa, "N");    
    // double* det_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, out_arr, 1);
}

double GetAlpha(const double* const rho_arr,
                double nu,
                const double* const resp_mat_arr,
                const double* const bg_arr, double B,
                const double* const data_arr,
                int nsky, int ndet){

    double* det_arr = new double[ndet];
    GetDetArr(rho_arr, resp_mat_arr, ndet, nsky, det_arr);

    double alpha = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        double den = det_arr[idet] + bg_arr[idet] / B * nu;
        alpha += data_arr[idet] * det_arr[idet] / den;
    }
    return alpha;
}

void GetRhoNu_New(const double* const rho_arr, double nu,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const bg_arr, double B,
                  int ndet, int nsky,
                  double* const rho_new_arr,
                  double* const nu_new_ptr){

    // rho_new_arr[isky] = 1./(alpha + B) * rho_arr[isky]
    //                    * [t(R_mat) %*% (data_arr / den_arr)][isky]
    // new_new = B / (alpha + B)
    double* det_arr = new double[ndet];
    GetDetArr(rho_arr, resp_mat_arr, ndet, nsky, det_arr);
    // den_arr
    double* den_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        den_arr[idet] = det_arr[idet] + bg_arr[idet] / B * nu;
    }
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / den_arr[idet];
    }
    double* tmp_arr = new double[nsky];
    char* transa = new char [1];
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, tmp_arr, 1);
    double alpha = GetAlpha(rho_arr, nu,
                            resp_mat_arr,
                            bg_arr, B,
                            data_arr,
                            nsky, ndet);
    for(int isky = 0; isky < nsky; isky ++){
        rho_new_arr[isky] = tmp_arr[isky] * rho_arr[isky] / (alpha + B);
    }
    double nu_new = B / (alpha + B);

    delete [] det_arr;
    delete [] den_arr;
    delete [] div_arr;
    delete [] tmp_arr;
    delete [] transa;
    *nu_new_ptr = nu_new;
}


void RichlucyBg(const double* const rho_init_arr,
                double nu_init,
                const double* const data_arr,
                const double* const bg_arr,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky,
                int niter,
                string outdir,
                string outfile_head,
                double tol,
                double* const rho_new_arr,
                double* const nu_new_ptr)
{
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    double nu_pre = nu_init;
    double B_val = GetB(bg_arr, ndet);
    double nu_new = nu_init;
    for(int iiter = 0; iiter < niter; iiter ++){
        GetRhoNu_New(rho_pre_arr, nu_pre,
                     data_arr,
                     resp_norm_mat_arr,
                     bg_arr, B_val,
                     ndet, nsky,
                     rho_new_arr,
                     &nu_new);
        double helldist = GetHellingerDist(rho_pre_arr, rho_new_arr, nsky);
        if (helldist < tol){
            printf("iiter = %d, helldist = %e\n", iiter, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        nu_pre = nu_new;
        printf("iiter = %d, helldist = %e\n", iiter, helldist);
    }
    delete [] rho_pre_arr;
    *nu_new_ptr = nu_new;
}

double GetHellingerDist(const double* const rho_arr,
                        const double* const rho_new_arr,
                        int nsky)
{
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        double diff = sqrt(rho_arr[isky]) - sqrt(rho_new_arr[isky]);
        sum += diff * diff;
    }
    double ans = sqrt(sum);
    return (ans);
}
