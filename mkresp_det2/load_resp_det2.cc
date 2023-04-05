#include "load_resp_det2.h"

void LoadRespDet2(string respdir1, string respdir2,
                  int nskyx, int nskyy,
                  int nphoton_input,
                  double** const resp_norm_det1_arr_ptr,
                  double** const resp_norm_det2_arr_ptr,
                  double** const eff_arr_ptr,
                  int* const ndetx_ptr,
                  int* const ndety_ptr)
{
    printf("--- LoadRespDet2 --- \n");
    //    row: detector
    //    col: sky
    int ndetx = 0;
    int ndety = 0;
    GetNdet(respdir1, &ndetx, &ndety);
    printf("LoadResp: ndetx = %d\n", ndetx);
    printf("LoadResp: ndety = %d\n", ndety);
    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, ndetx, ndety);

    int ndet = ndetx * ndety;
    int nsky = nskyx * nskyy;
    printf("LoadResp: ndet = %d, nsky = %d\n", ndet, nsky);

    printf("ndet * nsky = %d\n", ndet * nsky);
    double* resp_norm_det1_arr = new double [ndet * nsky];
    double* resp_norm_det2_arr = new double [ndet * nsky];
    double* eff_arr = new double [nsky];
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            int bitpix = 0;
            char infile1[kLineSize];
            sprintf(infile1, "%s/gimage_%3.3d_%3.3d.img",
                    respdir1.c_str(), iskyx, iskyy);
            char infile2[kLineSize];
            sprintf(infile2, "%s/gimage_%3.3d_%3.3d.img",
                    respdir2.c_str(), iskyx, iskyy);
            double* data1_arr = NULL;
            MifFits::InFitsImageD(infile1, img_info,
                                  &bitpix, &data1_arr);
            double* data2_arr = NULL;
            MifFits::InFitsImageD(infile2, img_info,
                                  &bitpix, &data2_arr);            
            
            // force to be non-negative
            for (int idet = 0; idet < ndet; idet ++){
                if(data1_arr[idet] < 0.0){
                    data1_arr[idet] = 0.0;
                }
                if(data2_arr[idet] < 0.0){
                    data2_arr[idet] = 0.0;
                }
            }
            // div by nphoton_input
            for(int idet = 0; idet < ndet; idet ++){
                data1_arr[idet] /= nphoton_input;
                data2_arr[idet] /= nphoton_input;
            }
            double eff = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                eff += data1_arr[idet];
                eff += data2_arr[idet];
            }
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;
            for(int idet = 0; idet < ndet; idet ++){
                resp_norm_det1_arr[imat + idet] = data1_arr[idet] / eff;
                resp_norm_det2_arr[imat + idet] = data2_arr[idet] / eff;
            }
            delete [] data1_arr;
            delete [] data2_arr;
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
                resp_norm_sum += resp_norm_det1_arr[imat + idet];
                resp_norm_sum += resp_norm_det2_arr[imat + idet];
            }
            // printf("resp_norm_sum = %e\n", resp_norm_sum);
            if ( fabs(resp_norm_sum - 1.0) > 1.e-10){
                printf("warning: resp_norm_sum = %e\n", resp_norm_sum);
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
    
    *resp_norm_det1_arr_ptr = resp_norm_det1_arr;
    *resp_norm_det2_arr_ptr = resp_norm_det2_arr;
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

