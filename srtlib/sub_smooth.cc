#include "sub_smooth.h"

int GetIbin(int ibinx, int ibiny, int nbinx)
{
    int ibin = ibinx + ibiny * nbinx;
    return(ibin);
}

double GetTermV(const double* const rho_arr, int nskyx, int nskyy)
{
    double termv = 0.0;
    for(int iskyx = 0; iskyx < nskyx - 1; iskyx ++){
        for(int iskyy = 0; iskyy < nskyy - 1; iskyy ++){
            int isky = iskyx + iskyy * nskyx;
            int isky_plus_x = (iskyx + 1) + iskyy * nskyx;
            int isky_plus_y = iskyx + (iskyy + 1) * nskyx;
            double diff1 = rho_arr[isky] - rho_arr[isky_plus_x];
            double diff2 = rho_arr[isky] - rho_arr[isky_plus_y];
            termv += diff1 * diff1 + diff2 * diff2;
        }
    }
    for(int iskyx = 0; iskyx < nskyx - 1; iskyx ++){
        int isky = iskyx + (nskyy - 1) * nskyx;
        int isky_plus_x = (iskyx + 1) + (nskyy - 1) * nskyx;
        double diff1 = rho_arr[isky] - rho_arr[isky_plus_x];
        termv += diff1 * diff1;
    }
    for(int iskyy = 0; iskyy < nskyy - 1; iskyy ++){
        int isky = (nskyx - 1) + iskyy * nskyx;
        int isky_plus_y = (nskyx - 1) + (iskyy + 1) * nskyx;
        double diff2 = rho_arr[isky] - rho_arr[isky_plus_y];
        termv += diff2 * diff2;
    }
    return (termv);
}


void GetDiffTermV(const double* const rho_arr, int nskyx, int nskyy,
                  double* const rho_diff_arr)
{
    // iskyx = 0, iskyy = 0
    // isky_plus_x, isky_plus_y
    {
        int iskyx = 0;
        int iskyy = 0;
        int isky = GetIbin(iskyx, iskyy, nskyx);
        int isky_plus_x = GetIbin(iskyx + 1, iskyy, nskyx);
        int isky_plus_y = GetIbin(iskyx, iskyy + 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_plus_x])
            + (rho_arr[isky] - rho_arr[isky_plus_y]);
    }

    // iskyx = 0, 1 <= iskyy <= nskyy - 2
    // isky_plus_x, isky_plus_y, isky_minus_y
    {
        int iskyx = 0;
        for(int iskyy = 1; iskyy < nskyy - 1; iskyy ++){
            int isky = GetIbin(iskyx, iskyy, nskyx);
            int isky_plus_x = GetIbin(iskyx + 1, iskyy, nskyx);
            int isky_plus_y = GetIbin(iskyx, iskyy + 1, nskyx);
            int isky_minus_y =  GetIbin(iskyx, iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }

    // iskyx = 0, iskyy = nskyy - 1
    // isky_plus_x, isky_minus_y
    {
        int iskyx = 0;
        int iskyy = nskyy - 1;
        int isky = GetIbin(iskyx, iskyy, nskyx);
        int isky_plus_x = GetIbin(iskyx + 1, iskyy, nskyx);
        int isky_minus_y = GetIbin(iskyx, iskyy - 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_plus_x])
            + (rho_arr[isky] - rho_arr[isky_minus_y]);
    }

    // 1 <= iskyx <= nskyx - 2, iskyy = 0
    // isky_plus_x, isky_minus_x, isky_plus_y
    {
        int iskyy = 0;
        for(int iskyx = 1; iskyx < nskyx - 1; iskyx ++){
            int isky = GetIbin(iskyx, iskyy, nskyx);
            int isky_plus_x  = GetIbin(iskyx + 1, iskyy    , nskyx);
            int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_plus_y  = GetIbin(iskyx    , iskyy + 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y]);
        }
    }

    // 1 <= iskyx <= nskyx - 2, 1 <= iskyy <= nskyy - 2
    // isky_plus_x, isky_minus_x, isky_plus_y, isky_minus_y
    for(int iskyx = 1; iskyx < nskyx - 1; iskyx ++){
        for(int iskyy = 1; iskyy < nskyy - 1; iskyy ++){
            int isky         = GetIbin(iskyx    , iskyy    , nskyx);
            int isky_plus_x  = GetIbin(iskyx + 1, iskyy    , nskyx);
            int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_plus_y  = GetIbin(iskyx    , iskyy + 1, nskyx);
            int isky_minus_y = GetIbin(iskyx    , iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }

    // 1 <= iskyx <= nskyx - 2, iskyy = nskyy - 1
    // isky_plus_x, isky_minus_x, isky_minus_y
    {
        int iskyy = nskyy - 1;
        for(int iskyx = 1; iskyx < nskyx - 1; iskyx ++){
            int isky          = GetIbin(iskyx    , iskyy    , nskyx);
            int isky_plus_x   = GetIbin(iskyx + 1, iskyy    , nskyx);
            int isky_minus_x  = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_minus_y  = GetIbin(iskyx    , iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }

    // iskyx = nskyx - 1, iskyy = 0
    // isky_minus_x, isky_plus_y
    {
        int iskyx = nskyx - 1;
        int iskyy = 0;
        int isky         = GetIbin(iskyx    , iskyy    , nskyx);
        int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
        int isky_plus_y  = GetIbin(iskyx    , iskyy + 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_minus_x])
            + (rho_arr[isky] - rho_arr[isky_plus_y]);
    }

    // iskyx = nskyx - 1, 1 <= iskyy <= nskyy - 2
    // isky_minus_x, isky_plus_y, isky_minus_y
    {
        int iskyx = nskyx - 1;
        for(int iskyy = 1; iskyy < nskyy - 1; iskyy ++){
            int isky          = GetIbin(iskyx    , iskyy    , nskyx);
            int isky_minus_x  = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_plus_y   = GetIbin(iskyx    , iskyy + 1, nskyx);
            int isky_minus_y  = GetIbin(iskyx    , iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }
    
    // iskyx = nskyx - 1, iskyy = nskyy - 1
    // isky_minus_x, isky_minus_y
    {
        int iskyx = nskyx - 1;
        int iskyy = nskyy - 1;
        int isky         = GetIbin(iskyx    , iskyy    , nskyx);
        int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
        int isky_minus_y = GetIbin(iskyx    , iskyy - 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_minus_x])
            + (rho_arr[isky] - rho_arr[isky_minus_y]);
    }
}

