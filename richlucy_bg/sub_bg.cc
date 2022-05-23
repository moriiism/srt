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
    int flag_line_search = 1;
    
    double* rho_pre_arr = new double[nsky];
    double* rho_this_arr = new double[nsky];
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

        // line search
        if (flag_line_search == 1){
            dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_this_arr, 1);
            double nu_this = nu_new;
            GetLineSearch(rho_pre_arr,
                          rho_this_arr,
                          nu_pre,
                          nu_this,
                          data_arr,
                          resp_norm_mat_arr,
                          bg_arr,
                          B_val,
                          ndet, nsky,
                          rho_new_arr,
                          &nu_new);
        }
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


double GetFuncL(const double* const data_arr,
                const double* const bg_arr,
                const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky)
{
    double B_val = GetB(bg_arr, ndet);
    
    // term1 = - sum_v [ Y(v) log( sum_u t(v,u) rho_u + b(v) / B * nu ) ]
    double* det_arr = new double[ndet];
    GetDetArr(rho_arr, resp_norm_mat_arr,
              ndet, nsky, det_arr);
    double term1 = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        term1 += data_arr[idet] *
            log(det_arr[idet] + bg_arr[idet] / B_val * nu);
    }
    term1 *= -1.0;
    delete [] det_arr;

    // term2 = - sum_v Y(v) log(B/nu) + B/nu
    double nph = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        nph += data_arr[idet];
    }
    double term2 = -1.0 * nph * log(B_val / nu) + B_val / nu;

    double ans = term1 + term2;
    return(ans);
}


void GetLineSearch(const double* const rho_1_arr,
                   const double* const rho_2_arr,
                   double nu_1,
                   double nu_2,
                   const double* const data_arr,
                   const double* const resp_norm_mat_arr,
                   const double* const bg_arr,
                   double B_val,
                   int ndet, int nsky,
                   double* const rho_new_arr,
                   double* nu_new_ptr)
{
    double lval_2 = GetFuncL(data_arr, bg_arr,
                             rho_2_arr, nu_2,
                             resp_norm_mat_arr,
                             ndet, nsky);
    // back tracking method
    double* theta_1_arr = new double[nsky];
    double* theta_2_arr = new double[nsky];
    double* theta_this_arr = new double[nsky];
    double* rho_this_arr = new double[nsky];
    double nu_this = 0.0;
    for(int isky = 0; isky < nsky; isky++){
        theta_1_arr[isky] = log( rho_1_arr[isky] / (1.0 - nu_1) );
        theta_2_arr[isky] = log( rho_2_arr[isky] / (1.0 - nu_2) );
    }
    int nstep = 10;
    double eta = 0.8;
    double kval = 100.0;
    int ifind = 0;
    for(int istep = 0; istep < nstep; istep ++){
        if (kval < 1.0){
            // not find better kval
            ifind = 0;
            break;
        }
        kval *= eta;
        double log_nu = kval * (log(nu_2) - log(nu_1)) + log(nu_1);
        nu_this = exp(log_nu);
        for(int isky = 0; isky < nsky; isky++){
            theta_this_arr[isky] = kval
                * (theta_2_arr[isky] - theta_1_arr[isky])
                + theta_1_arr[isky];
        }
        for(int isky = 0; isky < nsky; isky++){
            rho_this_arr[isky] = (1.0 - nu_this) * exp(theta_this_arr[isky]);
        }
        double lval_this = GetFuncL(data_arr, bg_arr,
                                    rho_this_arr, nu_this,
                                    resp_norm_mat_arr,
                                    ndet, nsky);
        if (lval_this < lval_2){
            printf("find: lval_this = %e\n", lval_this);
            // find better kval
            ifind = 1;
            break;
        }
    }
    printf("ifind = %d,", ifind);
    
    double nu_new = 0.0;
    if (ifind == 0){
        for(int isky = 0; isky < nsky; isky++){
            rho_new_arr[isky] = rho_2_arr[isky];
        }
        nu_new = nu_2;
    } else if (ifind == 1){
        for(int isky = 0; isky < nsky; isky++){
            rho_new_arr[isky] = rho_this_arr[isky];
        }
        nu_new = nu_this;
    }

    *nu_new_ptr = nu_new;
}
