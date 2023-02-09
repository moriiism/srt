#include "rl_statval.h"

double SrtlibRlStatval::GetHellingerDist(
    const double* const sky_pre_arr,
    const double* const sky_new_arr,
    int nsky)
{
    double sum_pre = 0.0;
    double sum_new = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum_pre += sky_pre_arr[isky];
        sum_new += sky_new_arr[isky];
    }
    double* sky_pre_norm_arr = new double[nsky];
    double* sky_new_norm_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_pre_norm_arr[isky] = sky_pre_arr[isky] / sum_pre;
        sky_new_norm_arr[isky] = sky_new_arr[isky] / sum_new;
    }
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        double diff = sqrt(sky_pre_norm_arr[isky])
            - sqrt(sky_new_norm_arr[isky]);
        sum += diff * diff;
    }
    double ans = sqrt(sum);
    delete [] sky_pre_norm_arr;
    delete [] sky_new_norm_arr;
    return (ans);
}
