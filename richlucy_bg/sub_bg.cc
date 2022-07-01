//#include "sub_bg.h"

//double GetFuncL(const double* const data_arr,
//                const double* const bg_arr,
//                const double* const rho_arr,
//                double nu,
//                const double* const resp_norm_mat_arr,
//                int ndet, int nsky)
//{
//    double B_val = MibBlas::Sum(bg_arr, ndet);
//    
//    // term1 = - sum_v [ Y(v) log( sum_u t(v,u) rho_u + b(v) / B * nu ) ]
//    double* tmp_arr = new double[ndet];
//    GetDetArr(rho_arr, resp_norm_mat_arr,
//              ndet, nsky, tmp_arr);
//    daxpy_(ndet, nu / B_val, const_cast<double*>(bg_arr), 1, tmp_arr, 1);
//    for(int idet = 0; idet < ndet; idet ++){
//        tmp_arr[idet] = log(tmp_arr[idet]);
//    }
//    double term1 = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1.0,
//                                tmp_arr, 1.0);
//    
//    //for(int idet = 0; idet < ndet; idet ++){
//    //    term1 += data_arr[idet] *
//    //        log(det_arr[idet] + bg_arr[idet] / B_val * nu);
//    //}
//    delete [] tmp_arr;
//    
//    // term2 = - sum_v Y(v) log(B/nu) + B/nu
//    double nph = MibBlas::Sum(data_arr, ndet);
//    //for(int idet = 0; idet < ndet; idet ++){
//    //    nph += data_arr[idet];
//    //}
//    double term2 = -1.0 * nph * log(B_val / nu) + B_val / nu;
//
//    double ans = term1 + term2;
//    return(ans);
//}
//

//void GetLineSearch(const double* const rho_1_arr,
//                   const double* const rho_2_arr,
//                   double nu_1,
//                   double nu_2,
//                   const double* const data_arr,
//                   const double* const resp_norm_mat_arr,
//                   const double* const bg_arr,
//                   double B_val,
//                   int ndet, int nsky,
//                   double* const rho_new_arr,
//                   double* nu_new_ptr)
//{
//    double nu_new = 0.0;
//    double lval_2 = GetFuncL(data_arr, bg_arr,
//                             rho_2_arr, nu_2,
//                             resp_norm_mat_arr,
//                             ndet, nsky);
//    double* theta_1_arr = new double[nsky];
//    double* theta_2_arr = new double[nsky];
//    for(int isky = 0; isky < nsky; isky++){
//        theta_1_arr[isky] = log( rho_1_arr[isky] / (1.0 - nu_1) );
//        theta_2_arr[isky] = log( rho_2_arr[isky] / (1.0 - nu_2) );
//    }
//    // adjcent    
//    double kval = 1.0;
//    double eta = 0.8;
//    // large
//    double* theta_large_arr = new double[nsky];
//    double* rho_large_arr = new double[nsky];
//    double kval_large = kval / eta;
//    double log_nu_large = kval_large * (log(nu_2) - log(nu_1)) + log(nu_1);
//    
//    double nu_large = exp(log_nu_large);
//    for(int isky = 0; isky < nsky; isky++){
//        theta_large_arr[isky] =
//            kval_large * (theta_2_arr[isky] - theta_1_arr[isky])
//            + theta_1_arr[isky];
//    }
//    for(int isky = 0; isky < nsky; isky++){
//        rho_large_arr[isky] = (1.0 - nu_large) * exp(theta_large_arr[isky]);
//    }
//
//    double lval_large = GetFuncL(data_arr, bg_arr,
//                                 rho_large_arr, nu_large,
//                                 resp_norm_mat_arr,
//                                 ndet, nsky);
//    // small
//    double* theta_small_arr = new double[nsky];
//    double* rho_small_arr = new double[nsky];
//    double kval_small = kval * eta;
//    double log_nu_small = kval_small * (log(nu_2) - log(nu_1)) + log(nu_1);
//    double nu_small = exp(log_nu_small);
//    for(int isky = 0; isky < nsky; isky++){
//        theta_small_arr[isky] =
//            kval_small * (theta_2_arr[isky] - theta_1_arr[isky])
//            + theta_1_arr[isky];
//    }
//    for(int isky = 0; isky < nsky; isky++){
//        rho_small_arr[isky] = (1.0 - nu_small) * exp(theta_small_arr[isky]);
//    }
//    double lval_small = GetFuncL(data_arr, bg_arr,
//                                 rho_small_arr, nu_small,
//                                 resp_norm_mat_arr,
//                                 ndet, nsky);
//    if ( (lval_2 <= lval_small) &&
//         (lval_2 <= lval_large) ){
//        for(int isky = 0; isky < nsky; isky++){
//            rho_new_arr[isky] = rho_2_arr[isky];
//        }
//        nu_new = nu_2;
//        *nu_new_ptr = nu_new;
//        printf("kval = 1.0\n");
//        return;
//    }
//    if (lval_small < lval_2){
//        eta = 0.8;
//    } else if (lval_large < lval_2){
//        eta = 1.25;
//    }
//    // search
//    double* theta_this_arr = new double[nsky];
//    double* rho_this_arr = new double[nsky];
//    double* rho_pre_arr = new double[nsky];
//    double nu_this = 0.0;
//    double nu_pre = 0.0;
//    kval = 1.0;
//    int nstep = 10;
//    double lval_pre = lval_2;
//    for(int istep = 0; istep < nstep; istep ++){
//        kval *= eta;
//        double log_nu = kval * (log(nu_2) - log(nu_1)) + log(nu_1);
//        nu_this = exp(log_nu);
//        
//        for(int isky = 0; isky < nsky; isky++){
//            theta_this_arr[isky] = kval
//                * (theta_2_arr[isky] - theta_1_arr[isky])
//                + theta_1_arr[isky];
//        }
//        for(int isky = 0; isky < nsky; isky++){
//            rho_this_arr[isky] = (1.0 - nu_this) * exp(theta_this_arr[isky]);
//        }
//        double lval_this = GetFuncL(data_arr, bg_arr,
//                                    rho_this_arr, nu_this,
//                                    resp_norm_mat_arr,
//                                    ndet, nsky);
//        printf("lval_pre, lval_this = %e, %e\n", 
//               lval_pre, lval_this);
//        
//        if (lval_pre < lval_this){
//            break;
//        }
//        lval_pre = lval_this;
//        for(int isky = 0; isky < nsky; isky++){
//            rho_pre_arr[isky] = rho_this_arr[isky];
//        }
//        nu_pre = nu_this;
//    }
//
//    for(int isky = 0; isky < nsky; isky++){
//        rho_new_arr[isky] = rho_pre_arr[isky];
//    }
//    nu_new = nu_pre;
//    *nu_new_ptr = nu_new;
//
//    
//    return;
//}
//
