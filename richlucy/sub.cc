#include "sub.h"

double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                int ndet, int nskyx, int nskyy, double epsilon)
{
    int nsky = nskyx * nskyy;

    double* rho_epsilon_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        if(rho_arr[isky] < epsilon){
            rho_epsilon_arr[isky] = epsilon;
        } else {
            rho_epsilon_arr[isky] = rho_arr[isky];
        }
    }
    
    // term1 = - sum_v [ Y(v) log( sum_u t(v,u) rho_u ) ]
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_epsilon_arr), 1,
           0.0, det_arr, 1);
    for(int idet = 0; idet < ndet; idet ++){
        det_arr[idet] = log(det_arr[idet]);
    }
    double term1 = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1,
                                const_cast<double*>(det_arr), 1);
    delete [] det_arr;
    delete [] transa;

    double ans = term1;

    delete [] rho_epsilon_arr;
    return(ans);
}


double GetKLDiv(const double* const rho_arr,
                const double* const rho_new_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky)
{
    char* transa = new char [1];
    strcpy(transa, "N");
    // q.vec = R.mat %*% y.vec
    double* q_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, q_arr, 1);
    // q.new.vec = R.mat %*% y.new.vec
    double* q_new_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_new_arr), 1,
           0.0, q_new_arr, 1);

    delete [] transa;

    // q.vec = q.vec / sum(q.vec)
    // q.new.vec = q.new.vec / sum(q.new.vec)
    double sum_q = 0.0;
    double sum_q_new = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        sum_q += q_arr[idet];
        sum_q_new += q_new_arr[idet];
    }
    dscal_(ndet, 1.0/sum_q, q_arr, 1);
    dscal_(ndet, 1.0/sum_q_new, q_new_arr, 1);
    
    double ans = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        if(q_new_arr[idet] > 0.0){
            ans = ans + q_new_arr[idet] * log( q_new_arr[idet] / q_arr[idet] );
        }
    }
    delete [] q_arr;
    delete [] q_new_arr;
    return (ans);
}




