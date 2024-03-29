#include "rl_crab_statval.h"

double SrtlibCrabRlCrabStatval::GetHellingerDist(
    const double* const sky_pre_arr,
    const double* const flux_pre_arr,
    const double* const sky_new_arr,
    const double* const flux_new_arr,
    int nsky, int nphase)
{
    double sum_pre = 0.0;
    double sum_new = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum_pre += sky_pre_arr[isky];
        sum_new += sky_new_arr[isky];
    }
    for(int iphase = 0; iphase < nphase; iphase++){
        sum_pre += flux_pre_arr[iphase];
        sum_new += flux_new_arr[iphase];
    }
    double* sky_pre_norm_arr = new double[nsky];
    double* flux_pre_norm_arr = new double[nphase];
    double* sky_new_norm_arr = new double[nsky];
    double* flux_new_norm_arr = new double[nphase];
    for(int isky = 0; isky < nsky; isky ++){
        sky_pre_norm_arr[isky] = sky_pre_arr[isky] / sum_pre;
        sky_new_norm_arr[isky] = sky_new_arr[isky] / sum_new;
    }
    for(int iphase = 0; iphase < nphase; iphase++){
        flux_pre_norm_arr[iphase] = flux_pre_arr[iphase] / sum_pre;
        flux_new_norm_arr[iphase] = flux_new_arr[iphase] / sum_new;
    }
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        double diff = sqrt(sky_pre_norm_arr[isky])
            - sqrt(sky_new_norm_arr[isky]);
        sum += diff * diff;
    }
    for(int iphase = 0; iphase < nphase; iphase++){
        double diff = sqrt(flux_pre_norm_arr[iphase])
            - sqrt(flux_new_norm_arr[iphase]);
        sum += diff * diff;
    }
    double ans = sqrt(sum);
    delete [] sky_pre_norm_arr;
    delete [] flux_pre_norm_arr;
    delete [] sky_new_norm_arr;
    delete [] flux_new_norm_arr;
    return (ans);
}

//double SrtlibRlStatval::GetNegLogLike(
//    const double* const rho_arr,
//    const double* const data_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nsky, double epsilon)
//{
//    double* rho_epsilon_arr = new double[nsky];
//    for(int isky = 0; isky < nsky; isky ++){
//        if(rho_arr[isky] < epsilon){
//            rho_epsilon_arr[isky] = epsilon;
//        } else {
//            rho_epsilon_arr[isky] = rho_arr[isky];
//        }
//    }
//    
//    // - log likelihood = - sum_v [ Y(v) log( sum_u t(v,u) rho_u ) ]
//    double* det_arr = new double[ndet];
//    char transa[1];
//    strcpy(transa, "N");
//    dgemv_(transa, ndet, nsky, 1.0,
//           const_cast<double*>(resp_norm_mat_arr), ndet,
//           const_cast<double*>(rho_epsilon_arr), 1,
//           0.0, det_arr, 1);
//    for(int idet = 0; idet < ndet; idet ++){
//        det_arr[idet] = log(det_arr[idet]);
//    }
//    double neg_log_like = -1.0 * ddot_(ndet,
//                                       const_cast<double*>(data_arr), 1,
//                                       const_cast<double*>(det_arr), 1);
//    delete [] det_arr;
//    delete [] rho_epsilon_arr;
//    return neg_log_like;
//}


//double SrtlibRlStatval::GetKLDiv(const double* const rho_arr,
//                           const double* const rho_new_arr,
//                           const double* const resp_mat_arr,
//                           int ndet, int nsky)
//{
//    char* transa = new char [1];
//    strcpy(transa, "N");
//    // q.vec = R.mat %*% y.vec
//    double* q_arr = new double[ndet];
//    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
//           const_cast<double*>(rho_arr), 1,
//           0.0, q_arr, 1);
//    // q.new.vec = R.mat %*% y.new.vec
//    double* q_new_arr = new double[ndet];
//    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
//           const_cast<double*>(rho_new_arr), 1,
//           0.0, q_new_arr, 1);
//
//    delete [] transa;
//
//    // q.vec = q.vec / sum(q.vec)
//    // q.new.vec = q.new.vec / sum(q.new.vec)
//    double sum_q = 0.0;
//    double sum_q_new = 0.0;
//    for(int idet = 0; idet < ndet; idet ++){
//        sum_q += q_arr[idet];
//        sum_q_new += q_new_arr[idet];
//    }
//    dscal_(ndet, 1.0/sum_q, q_arr, 1);
//    dscal_(ndet, 1.0/sum_q_new, q_new_arr, 1);
//    
//    double ans = 0.0;
//    for(int idet = 0; idet < ndet; idet ++){
//        if(q_new_arr[idet] > 0.0){
//            ans = ans + q_new_arr[idet] * log( q_new_arr[idet] / q_arr[idet] );
//        }
//    }
//    delete [] q_arr;
//    delete [] q_new_arr;
//    return (ans);
//}
