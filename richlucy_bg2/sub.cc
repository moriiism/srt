#include "sub.h"

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


void GetRhoNu_New(const double* const rho_arr, double nu,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const bg_arr,
                  int ndet, int nsky,
                  double* const rho_new_arr,
                  double* const nu_new_ptr){

    double nph = MibBlas::Sum(data_arr, ndet);
    
    // rho_new_arr[isky] = 1/nph * rho_arr[isky]
    //                    * [t(R_mat) %*% (data_arr / den_arr)][isky]
    // nu_new = 1/nph * nu
    //         * [bg_arr * (data_arr / den_arr)]
    double* den_arr = new double[ndet];
    GetDetArr(rho_arr, resp_mat_arr, ndet, nsky, den_arr);
    daxpy_(ndet, nu, const_cast<double*>(bg_arr), 1, den_arr, 1);
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
    MibBlas::ElmWiseMul(nsky, 1.0 / nph,
                        tmp_arr, rho_arr, rho_new_arr);

    double nu_new = ddot_(ndet, div_arr, 1, const_cast<double*>(bg_arr), 1)
        * nu / nph;
    
    delete [] den_arr;
    delete [] div_arr;
    delete [] tmp_arr;
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
    double nu_new = nu_init;
    for(int iiter = 0; iiter < niter; iiter ++){
        GetRhoNu_New(rho_pre_arr, nu_pre,
                     data_arr,
                     resp_norm_mat_arr,
                     bg_arr,
                     ndet, nsky,
                     rho_new_arr,
                     &nu_new);

        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre,
                                            rho_new_arr, nu_new, nsky);
        if (helldist < tol){
            printf("iiter = %d, helldist = %e\n",
                   iiter, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        nu_pre = nu_new;

        double lval = 0.0;        
        if (iiter % 100 == 0){
            lval = GetFuncL(data_arr, bg_arr,
                            rho_new_arr, nu_new,
                            resp_norm_mat_arr,
                            ndet, nsky);
            printf("iiter = %d, helldist = %e, lval = %e\n",
                   iiter, helldist, lval);
        } else {
            printf("iiter = %d, helldist = %e\n",
                   iiter, helldist);
        }
    }
    delete [] rho_pre_arr;
    *nu_new_ptr = nu_new;
}

double GetHellingerDist(const double* const rho_arr, double nu, 
                        const double* const rho_new_arr, double nu_new,
                        int nsky)
{
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        double diff = sqrt(rho_arr[isky]) - sqrt(rho_new_arr[isky]);
        sum += diff * diff;
    }
    double diff = sqrt(nu) - sqrt(nu_new);
    sum += diff * diff;
    double ans = sqrt(sum);
    return (ans);
}


double GetFuncL(const double* const data_arr,
                const double* const bg_arr,
                const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky)
{
    // ans = - sum_v [ Y(v) log( sum_u t(v,u) rho_u + b(v) * nu ) ]
    double* tmp_arr = new double[ndet];
    GetDetArr(rho_arr, resp_norm_mat_arr,
              ndet, nsky, tmp_arr);
    daxpy_(ndet, nu, const_cast<double*>(bg_arr), 1, tmp_arr, 1);
    for(int idet = 0; idet < ndet; idet ++){
        tmp_arr[idet] = log(tmp_arr[idet]);
    }
    double ans = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1.0,
                              tmp_arr, 1.0);
    delete [] tmp_arr;
    return(ans);
}
