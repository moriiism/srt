#include "fpsrc_smth_bg_statval.h"

double GetHellingerDist(const double* const rho_arr,
                        const double* const nu_arr,
                        double phi,
                        const double* const rho_new_arr,
                        const double* const nu_new_arr,
                        double phi_new,
                        int nsky, int nsrc)
{
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        double diff = sqrt(rho_arr[isky]) - sqrt(rho_new_arr[isky]);
        sum += diff * diff;
    }
    for(int isrc = 0; isrc < nsrc; isrc ++){
        double diff = sqrt(nu_arr[isrc]) - sqrt(nu_new_arr[isrc]);
        sum += diff * diff;
    }
    double diff = sqrt(phi) - sqrt(phi_new);
    sum += diff * diff;
    double ans = sqrt(sum);
    return (ans);
}

