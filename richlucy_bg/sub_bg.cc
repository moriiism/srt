#include "sub_bg.h"

double GetB(const double* const bg_arr, int ndet)
{
    double B_val = MibBlas::Sum(bg_arr, ndet);
    return B_val;
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

double GetAlpha(const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                const double* const bg_arr,
                const double* const data_arr,
                int nsky, int ndet){
    double* tmp_arr = new double[ndet];
    GetDetArr(rho_arr, resp_norm_mat_arr, ndet, nsky, tmp_arr);
    double B_val = GetB(bg_arr, ndet);
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

//double GetAlpha(const double* const rho_arr,
//                double nu,
//                const double* const resp_mat_arr,
//                const double* const bg_arr, double B,
//                const double* const data_arr,
//                int nsky, int ndet){
//
//    double* det_arr = new double[ndet];
//    GetDetArr(rho_arr, resp_mat_arr, ndet, nsky, det_arr);
//
//    double alpha = 0.0;
//    for(int idet = 0; idet < ndet; idet++){
//        double den = det_arr[idet] + bg_arr[idet] / B * nu;
//        alpha += data_arr[idet] * det_arr[idet] / den;
//    }
//    return alpha;
//}


void GetRhoNu_New(const double* const rho_arr, double nu,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const bg_arr,
                  int ndet, int nsky,
                  double* const rho_new_arr,
                  double* const nu_new_ptr){

    double B_val = GetB(bg_arr, ndet);
    
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
    // int flag_line_search = 1;
    
    double* rho_pre_arr = new double[nsky];
    // double* rho_this_arr = new double[nsky];
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
        
        //// line search
        //if (flag_line_search == 1){
        //    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_this_arr, 1);
        //  double nu_this = nu_new;
        //    GetLineSearch(rho_pre_arr,
        //                 rho_this_arr,
        //                 nu_pre,
        //                 nu_this,
        //                 data_arr,
        //                 resp_norm_mat_arr,
        //                bg_arr,
        //                B_val,
        //                ndet, nsky,
        //                rho_new_arr,
        //                &nu_new);
        //
        //}

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
    double B_val = GetB(bg_arr, ndet);
    
    // term1 = - sum_v [ Y(v) log( sum_u t(v,u) rho_u + b(v) / B * nu ) ]
    double* tmp_arr = new double[ndet];
    GetDetArr(rho_arr, resp_norm_mat_arr,
              ndet, nsky, tmp_arr);
    daxpy_(ndet, nu / B_val, const_cast<double*>(bg_arr), 1, tmp_arr, 1);
    for(int idet = 0; idet < ndet; idet ++){
        tmp_arr[idet] = log(tmp_arr[idet]);
    }
    double term1 = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1.0,
                                tmp_arr, 1.0);
    
    //for(int idet = 0; idet < ndet; idet ++){
    //    term1 += data_arr[idet] *
    //        log(det_arr[idet] + bg_arr[idet] / B_val * nu);
    //}
    delete [] tmp_arr;
    
    // term2 = - sum_v Y(v) log(B/nu) + B/nu
    double nph = MibBlas::Sum(data_arr, ndet);
    //for(int idet = 0; idet < ndet; idet ++){
    //    nph += data_arr[idet];
    //}
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
    double nu_new = 0.0;
    double lval_2 = GetFuncL(data_arr, bg_arr,
                             rho_2_arr, nu_2,
                             resp_norm_mat_arr,
                             ndet, nsky);
    double* theta_1_arr = new double[nsky];
    double* theta_2_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky++){
        theta_1_arr[isky] = log( rho_1_arr[isky] / (1.0 - nu_1) );
        theta_2_arr[isky] = log( rho_2_arr[isky] / (1.0 - nu_2) );
    }
    // adjcent    
    double kval = 1.0;
    double eta = 0.8;
    // large
    double* theta_large_arr = new double[nsky];
    double* rho_large_arr = new double[nsky];
    double kval_large = kval / eta;
    double log_nu_large = kval_large * (log(nu_2) - log(nu_1)) + log(nu_1);
    
    double nu_large = exp(log_nu_large);
    for(int isky = 0; isky < nsky; isky++){
        theta_large_arr[isky] =
            kval_large * (theta_2_arr[isky] - theta_1_arr[isky])
            + theta_1_arr[isky];
    }
    for(int isky = 0; isky < nsky; isky++){
        rho_large_arr[isky] = (1.0 - nu_large) * exp(theta_large_arr[isky]);
    }

    double lval_large = GetFuncL(data_arr, bg_arr,
                                 rho_large_arr, nu_large,
                                 resp_norm_mat_arr,
                                 ndet, nsky);
    // small
    double* theta_small_arr = new double[nsky];
    double* rho_small_arr = new double[nsky];
    double kval_small = kval * eta;
    double log_nu_small = kval_small * (log(nu_2) - log(nu_1)) + log(nu_1);
    double nu_small = exp(log_nu_small);
    for(int isky = 0; isky < nsky; isky++){
        theta_small_arr[isky] =
            kval_small * (theta_2_arr[isky] - theta_1_arr[isky])
            + theta_1_arr[isky];
    }
    for(int isky = 0; isky < nsky; isky++){
        rho_small_arr[isky] = (1.0 - nu_small) * exp(theta_small_arr[isky]);
    }
    double lval_small = GetFuncL(data_arr, bg_arr,
                                 rho_small_arr, nu_small,
                                 resp_norm_mat_arr,
                                 ndet, nsky);
    if ( (lval_2 <= lval_small) &&
         (lval_2 <= lval_large) ){
        for(int isky = 0; isky < nsky; isky++){
            rho_new_arr[isky] = rho_2_arr[isky];
        }
        nu_new = nu_2;
        *nu_new_ptr = nu_new;
        printf("kval = 1.0\n");
        return;
    }
    if (lval_small < lval_2){
        eta = 0.8;
    } else if (lval_large < lval_2){
        eta = 1.25;
    }
    // search
    double* theta_this_arr = new double[nsky];
    double* rho_this_arr = new double[nsky];
    double* rho_pre_arr = new double[nsky];
    double nu_this = 0.0;
    double nu_pre = 0.0;
    kval = 1.0;
    int nstep = 10;
    double lval_pre = lval_2;
    for(int istep = 0; istep < nstep; istep ++){
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
        printf("lval_pre, lval_this = %e, %e\n", 
               lval_pre, lval_this);
        
        if (lval_pre < lval_this){
            break;
        }
        lval_pre = lval_this;
        for(int isky = 0; isky < nsky; isky++){
            rho_pre_arr[isky] = rho_this_arr[isky];
        }
        nu_pre = nu_this;
    }

    for(int isky = 0; isky < nsky; isky++){
        rho_new_arr[isky] = rho_pre_arr[isky];
    }
    nu_new = nu_pre;
    *nu_new_ptr = nu_new;

    
    return;
}
