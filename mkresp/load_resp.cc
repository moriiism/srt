#include "load_resp.h"

void LoadResp(string respdir, int nskyx, int nskyy,
              int nphoton_input,
              double** const resp_arr_ptr,
              double** const resp_norm_arr_ptr,
              double** const eff_arr_ptr,
              int* const ndetx_ptr,
              int* const ndety_ptr)
{
    printf("--- LoadResp --- \n");
    //    row: detector
    //    col: sky
    int ndetx = 0;
    int ndety = 0;
    GetNdet(respdir, &ndetx, &ndety);
    printf("LoadResp: ndetx = %d\n", ndetx);
    printf("LoadResp: ndety = %d\n", ndety);
    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, ndetx, ndety);

    
    int ndet = ndetx * ndety;
    int nsky = nskyx * nskyy;
    printf("LoadResp: ndet = %d, nsky = %d\n", ndet, nsky);
    double* resp_arr = new double [ndet * nsky];
    double* resp_norm_arr = new double [ndet * nsky];
    double* eff_arr = new double [nsky];
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            char infile[kLineSize];
            sprintf(infile, "%s/gimage_%3.3d_%3.3d.img",
                    respdir.c_str(), iskyx, iskyy);
            int bitpix = 0;
            double* data_arr = NULL;
            MifFits::InFitsImageD(infile, img_info,
                                  &bitpix, &data_arr);
            // force to be non-negative
            for (int idet = 0; idet < ndet; idet ++){
                if(data_arr[idet] < 0.0){
                    data_arr[idet] = 0.0;
                }
            }
            // div by nphoton_input
            for(int idet = 0; idet < ndet; idet ++){
                data_arr[idet] /= nphoton_input;
            }
            double eff = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                eff += data_arr[idet];
            }
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;
            for(int idet = 0; idet < ndet; idet ++){
                resp_arr[imat + idet] = data_arr[idet];
                resp_norm_arr[imat + idet] = data_arr[idet] / eff;
            }
            delete [] data_arr;
            // fill efficiency matrix
            eff_arr[isky] = eff;
        }
    }
    delete img_info;

    // check
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;
            double resp_norm_sum = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                resp_norm_sum += resp_norm_arr[imat + idet];
            }
            // printf("resp_norm_sum = %e\n", resp_norm_sum);
            if ( fabs(resp_norm_sum - 1.0) > 1.e-10){
                printf("warning: resp_norm_sum = %e\n", resp_norm_sum);
            }
        }
    }

    // check
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;
            double resp_sum = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                resp_sum += resp_arr[imat + idet];
            }
            //printf("resp_sum = %e, eff_arr[isky] = %e\n",
            //    resp_sum, eff_arr[isky]);
            if ( fabs(resp_sum - eff_arr[isky]) > 1.e-10){
                printf("warning:\n");
                printf("resp_sum = %e, eff_arr[isky] = %e\n",
                       resp_sum, eff_arr[isky]);
            }
        }
    }

    // average efficiency
    double ave_eff = 0.0;
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            int isky = nskyx * iskyy + iskyx;
            ave_eff += eff_arr[isky];
        }
    }
    ave_eff /= nskyx * nskyy;
    printf("ave_eff = %e\n", ave_eff);
    
    
    *resp_arr_ptr = resp_arr;
    *resp_norm_arr_ptr = resp_norm_arr;
    *eff_arr_ptr = eff_arr;
    *ndetx_ptr = ndetx;
    *ndety_ptr = ndety;
}


void GetNdet(string respdir, int* const ndetx_ptr, int* const ndety_ptr)
{
    char infile[kLineSize];
    sprintf(infile, "%s/gimage_000_000.img", respdir.c_str());
    int naxis = MifFits::GetNaxis(infile);
    printf("GetNdet: naxis = %d\n", naxis);
    if(2 != naxis){
        printf("GetNdet: bad naxis = %d\n", naxis);
        abort();
    }
    int* ndet_arr = new int[naxis];
    for(int iaxis = 0; iaxis < naxis; iaxis ++){
        ndet_arr[iaxis] = MifFits::GetAxisSize(infile, iaxis);
    }
    *ndetx_ptr = ndet_arr[0];
    *ndety_ptr = ndet_arr[1];
    delete [] ndet_arr;
}

