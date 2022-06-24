#include "rl_statval.h"

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

