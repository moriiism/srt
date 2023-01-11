#include "rl_crab_statval_cuda.h"

double SrtlibCrabRlCrabStatvalCuda::GetHellingerDist(
    const double* const sky_pre_dev_arr,
    const double* const flux_pre_dev_arr,
    const double* const sky_new_dev_arr,
    const double* const flux_new_dev_arr,
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
