#include "smth_zal.h"
#include "smth_zal_cuda.h"

__global__
void SrtlibSmthZalCuda::GetDerivUBetaArr(
    const double* const sky_dash_dev_arr,
    int nskyx, int nskyy,
    double* const beta_dev_arr)
{
    int isky = blockIdx.x * blockDim.x + threadIdx.x;
    int iskyx = isky % nskyx;
    int iskyy = isky / nskyx;
    beta_dev_arr[isky] = 0.0;

    double beta = 0.0;
    int isky_plus_x = (iskyx + 1) + iskyy * nskyx;
    int isky_minus_x = (iskyx - 1) + iskyy * nskyx;
    int isky_plus_y = iskyx + (iskyy + 1) * nskyx;
    int isky_minus_y = iskyx + (iskyy - 1) * nskyx;
            
    // term1
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        beta += sky_dash_dev_arr[isky] + sky_dash_dev_arr[isky_plus_x];
    }
    // term2
    if((1 <= iskyx && iskyx < nskyx) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        beta += sky_dash_dev_arr[isky_minus_x] + sky_dash_dev_arr[isky];
    }
    // term3
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        beta += sky_dash_dev_arr[isky] + sky_dash_dev_arr[isky_plus_y];
    }
    // term4
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (1 <= iskyy && iskyy < nskyy)){
        beta += sky_dash_dev_arr[isky_minus_x] + sky_dash_dev_arr[isky];
    }
    // term5
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (nskyy - 1 == iskyy)){
        beta += sky_dash_dev_arr[isky] + sky_dash_dev_arr[isky_plus_x];
    }
    // term6
    if((1 <= iskyx && iskyx < nskyx) &&
       (nskyy - 1 == iskyy)){
        beta += sky_dash_dev_arr[isky_minus_x] + sky_dash_dev_arr[isky];
    }
    // term7
    if((nskyx - 1 == iskyx) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        beta += sky_dash_dev_arr[isky] + sky_dash_dev_arr[isky_plus_y];
    }
    // term8
    if((nskyx - 1 == iskyx) &&
       (1 <= iskyy && iskyy < nskyy)){
        beta += sky_dash_dev_arr[isky_minus_y] + sky_dash_dev_arr[isky];
    }
    beta *= 2.0;
    beta_dev_arr[isky] = beta;
}
