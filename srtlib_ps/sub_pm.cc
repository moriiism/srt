#include "sub_pm.h"
#include "sub_smooth.h"

void GetVvalArr(const double* const rho_arr,
                int nskyx, int nskyy,
                double mu, double lip_const,
                double* const vval_arr)
{
    int nsky = nskyx * nskyy;
    double* deriv_rho_arr = new double[nsky];
    GetDiffTermV(rho_arr, nskyx, nskyy, deriv_rho_arr);
    for(int isky = 0; isky < nsky; isky++){
        vval_arr[isky] = rho_arr[isky] - mu / lip_const * deriv_rho_arr[isky];
    }
    delete [] deriv_rho_arr;
}

double GetWval(double nu)
{
    double wval = nu;
    return wval;
}

