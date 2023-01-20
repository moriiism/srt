#include "srtmathlib.h"

double SrtMathlib::GetSum(
    long narr,
    const double* const val_arr)
{
    MiBase::IsValidArray(narr, val_arr);
    double ans = 0.0;
    for(long idata = 0; idata < narr; idata++){
        ans += val_arr[idata];
    }
    return ans;
}

double SrtMathlib::GetRootMeanSquareError(
    const double* const det_arr,
    const double* const val_arr,
    int ndet)
{
    double sum = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        double diff = det_arr[idet] - val_arr[idet];
        sum += diff * diff;
    }
    double ave = sum / ndet;
    double ans = sqrt(ave);
    return ans;
}
