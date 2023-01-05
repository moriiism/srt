#include "rl_statval_crab.h"

double SrtlibRlStatvalCrab::GetHellingerDist(
    const double* const rho_arr,
    const double* const nu_arr,
    const double* const rho_new_arr,
    const double* const nu_new_arr,
    int nsky, int nphase)
{
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        double diff = sqrt(rho_arr[isky]) - sqrt(rho_new_arr[isky]);
        sum += diff * diff;
    }
    for(int iphase = 0; iphase < nphase; iphase ++){
        double diff = sqrt(nu_arr[iphase]) - sqrt(nu_new_arr[iphase]);
        sum += diff * diff;
    }
    double ans = sqrt(sum);
    return (ans);
}

