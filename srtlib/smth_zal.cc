#include "smth_zal.h"

int SrtlibSmthZal::GetIbin(int ibinx, int ibiny, int nbinx)
{
    int ibin = ibinx + ibiny * nbinx;
    return(ibin);
}

double SrtlibSmthZal::GetDerivUAlpha(int iskyx, int iskyy,
                                     int nskyx, int nskyy)
{
    double alpha = 0.0;
    // term1
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        alpha += 4.0;
    }
    // term2
    if((1 <= iskyx && iskyx < nskyx) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        alpha += 4.0;
    }
    // term3
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        alpha += 4.0;
    }
    // term4
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (1 <= iskyy && iskyy < nskyy)){
        alpha += 4.0;
    }
    // term5
    if((0 <= iskyx && iskyx < nskyx - 1) &&
       (nskyy - 1 == iskyy)){
        alpha += 4.0;
    }
    // term6
    if((1 <= iskyx && iskyx < nskyx) &&
       (nskyy - 1 == iskyy)){
        alpha += 4.0;
    }
    // term7
    if((nskyx - 1 == iskyx) &&
       (0 <= iskyy && iskyy < nskyy - 1)){
        alpha += 4.0;
    }
    // term8
    if((nskyx - 1 == iskyx) &&
       (1 <= iskyy && iskyy < nskyy)){
        alpha += 4.0;
    }

    return alpha;
}

void SrtlibSmthZal::GetDerivUAlphaArr(int nskyx, int nskyy,
                                      double* const alpha_arr)
{
    for(int iskyx = 0; iskyx < nskyx; iskyx ++){
        for(int iskyy = 0; iskyy < nskyy; iskyy ++){
            int isky = GetIbin(iskyx, iskyy, nskyx);
            alpha_arr[isky] = SrtlibSmthZal::GetDerivUAlpha(iskyx, iskyy,
                                                            nskyx, nskyy);
        }
    }
}

void SrtlibSmthZal::GetDerivUBetaArr(const double* const sky_dash_arr,
                                     int nskyx, int nskyy,
                                     double* const beta_arr)
{
    for(int iskyx = 0; iskyx < nskyx; iskyx ++){
        for(int iskyy = 0; iskyy < nskyy; iskyy ++){
            int isky = GetIbin(iskyx, iskyy, nskyx);
            beta_arr[isky] = 0.0;
        }
    }
    for(int iskyx = 0; iskyx < nskyx; iskyx ++){
        for(int iskyy = 0; iskyy < nskyy; iskyy ++){
            double beta = 0.0;
            int isky = GetIbin(iskyx, iskyy, nskyx);
            int isky_plus_x = GetIbin(iskyx + 1, iskyy, nskyx);
            int isky_minus_x = GetIbin(iskyx - 1, iskyy, nskyx);
            int isky_plus_y = GetIbin(iskyx, iskyy + 1, nskyx);
            int isky_minus_y = GetIbin(iskyx, iskyy - 1, nskyx);
            
            // term1
            if((0 <= iskyx && iskyx < nskyx - 1) &&
               (0 <= iskyy && iskyy < nskyy - 1)){
                beta += sky_dash_arr[isky] + sky_dash_arr[isky_plus_x];
            }
            // term2
            if((1 <= iskyx && iskyx < nskyx) &&
               (0 <= iskyy && iskyy < nskyy - 1)){
                beta += sky_dash_arr[isky_minus_x] + sky_dash_arr[isky];
            }
            // term3
            if((0 <= iskyx && iskyx < nskyx - 1) &&
               (0 <= iskyy && iskyy < nskyy - 1)){
                beta += sky_dash_arr[isky] + sky_dash_arr[isky_plus_y];
            }
            // term4
            if((0 <= iskyx && iskyx < nskyx - 1) &&
               (1 <= iskyy && iskyy < nskyy)){
                beta += sky_dash_arr[isky_minus_x] + sky_dash_arr[isky];
            }
            // term5
            if((0 <= iskyx && iskyx < nskyx - 1) &&
               (nskyy - 1 == iskyy)){
                beta += sky_dash_arr[isky] + sky_dash_arr[isky_plus_x];
            }
            // term6
            if((1 <= iskyx && iskyx < nskyx) &&
               (nskyy - 1 == iskyy)){
                beta += sky_dash_arr[isky_minus_x] + sky_dash_arr[isky];
            }
            // term7
            if((nskyx - 1 == iskyx) &&
               (0 <= iskyy && iskyy < nskyy - 1)){
                beta += sky_dash_arr[isky] + sky_dash_arr[isky_plus_y];
            }
            // term8
            if((nskyx - 1 == iskyx) &&
               (1 <= iskyy && iskyy < nskyy)){
                beta += sky_dash_arr[isky_minus_y] + sky_dash_arr[isky];
            }
            beta *= 2.0;
            beta_arr[isky] = beta;
        }
    }
}

