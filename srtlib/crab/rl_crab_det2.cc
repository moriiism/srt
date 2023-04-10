#include<unistd.h>
#include "rl_crab_det2.h"

void SrtlibRlCrabDet2::GetMvalArr(
    const double* const* const y_dash_det1_arr,
    const double* const* const y_dash_det2_arr,    
    const double* const resp_norm_mat_det1_arr,
    const double* const resp_norm_mat_det2_arr,
    const double* const sky_arr,
    int ndet, int nsky, int nphase,
    double* const mval_arr)
{
    double* y_dash_sum_det1_arr = new double[ndet];
    double* y_dash_sum_det2_arr = new double[ndet];
    dcopy_(ndet, const_cast<double*>(y_dash_det1_arr[0]), 1,
           y_dash_sum_det1_arr, 1);
    dcopy_(ndet, const_cast<double*>(y_dash_det2_arr[0]), 1,
           y_dash_sum_det2_arr, 1);    
    for(int iphase = 1; iphase < nphase; iphase++){
        daxpy_(ndet, 1.0,
               const_cast<double*>(y_dash_det1_arr[iphase]), 1,
               y_dash_sum_det1_arr, 1);
        daxpy_(ndet, 1.0,
               const_cast<double*>(y_dash_det2_arr[iphase]), 1,
               y_dash_sum_det2_arr, 1);
    }

    double* coeff_det1_arr = new double[nsky];
    double* coeff_det2_arr = new double[nsky];    
    char transa[2];
    strcpy(transa, "T");
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_det1_arr), ndet,
           y_dash_sum_det1_arr, 1,
           0.0, coeff_det1_arr, 1);
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_det2_arr), ndet,
           y_dash_sum_det2_arr, 1,
           0.0, coeff_det2_arr, 1);

    double* coeff_arr = new double[nsky];
    MibBlas::Add(coeff_det1_arr,
                 coeff_det2_arr,
                 nsky,
                 coeff_arr);
    MibBlas::ElmWiseMul(nsky, 1.0,
                        coeff_arr, sky_arr,
                        mval_arr);
    delete [] y_dash_sum_det1_arr;
    delete [] y_dash_sum_det2_arr;    
    delete [] coeff_det1_arr;
    delete [] coeff_det2_arr;
    delete [] coeff_arr;
}

void SrtlibRlCrabDet2::GetNvalArr(
    const double* const* const y_dash_det1_arr,
    const double* const* const y_dash_det2_arr,    
    const double* const flux_arr,
    const double* const det_0_det1_arr,
    const double* const det_0_det2_arr,    
    int ndet, int nphase,
    double* const nval_arr)
{
    for(int iphase = 0; iphase < nphase; iphase++){
        double dot_det1 = ddot_(
            ndet,
            const_cast<double*>(y_dash_det1_arr[iphase]), 1,
            const_cast<double*>(det_0_det1_arr), 1);
        double dot_det2 = ddot_(
            ndet,
            const_cast<double*>(y_dash_det2_arr[iphase]), 1,
            const_cast<double*>(det_0_det2_arr), 1);
        nval_arr[iphase] = flux_arr[iphase]
            * (dot_det1 + dot_det2);
    }
}


