#include "sub_bg.h"

double GetNextNb(const double* const rho_arr,
                 const double* const data_arr,
                 const double* const resp_mat_arr,
                 const double* const bg_arr,
                 int ndet, int nsky,
                 double N_B)
{
    double N_B_new = N_B;
    double B = GetB(bg_arr, ndet);
    printf("B = %e\n", B);
    
    double deriv_f_at_B = GetDerivF_NB(rho_arr, data_arr, resp_mat_arr,
                                       bg_arr, ndet, nsky, B);
    if(deriv_f_at_B >= 0.0){
        N_B_new = B;
    } else{
        // get next N_B by Newton Method
        for(int iter = 0; iter < 10; iter ++){
            double deriv_f_at_NB = GetDerivF_NB(rho_arr, data_arr,
                                                resp_mat_arr, bg_arr,
                                                ndet, nsky, N_B_new);
            double deriv2_f_at_NB = GetDeriv2F_NB(rho_arr, data_arr,
                                                  resp_mat_arr, bg_arr,
                                                  ndet, nsky, N_B_new);
            N_B_new -= deriv_f_at_NB / deriv2_f_at_NB;
        }
    }
    return(N_B_new);
}

double GetB(const double* const bg_arr, int ndet)
{
    double B = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        B += bg_arr[idet];
    }
    return(B);
}

double GetN(const double* const rho_arr, int nsky)
{
    double N = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        N += rho_arr[isky];
    }
    return(N);
}

double GetDerivF_NB(const double* const rho_arr,
                    const double* const data_arr,
                    const double* const resp_mat_arr,
                    const double* const bg_arr,
                    int ndet, int nsky,
                    double N_B)
{
    // det_arr = R_mat %*% rho_arr
    char* transa = new char [1];
    strcpy(transa, "N");    
    double* det_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);

    double sum = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        sum += data_arr[idet] * det_arr[idet] /
            (N_B * det_arr[idet] + bg_arr[idet]);
        //printf("data = %e\n ", data_arr[idet]);
        //printf("det = %e\n ", det_arr[idet]);
        //printf("N_B = %e\n ", N_B);
        //printf("bg  = %e\n ", bg_arr[idet]);
        //printf("sum = %e\n", sum);
    }
    double deriv_f = 1.0 - sum;
    delete [] transa;
    delete [] det_arr;
    return(deriv_f);
}

double GetDeriv2F_NB(const double* const rho_arr,
                     const double* const data_arr,
                     const double* const resp_mat_arr,
                     const double* const bg_arr,
                     int ndet, int nsky,
                     double N_B)
{
    // det_arr = R_mat %*% rho_arr
    char* transa = new char [1];
    strcpy(transa, "N");    
    double* det_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);

    double sum = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        sum += data_arr[idet] * pow(det_arr[idet], 2)
            / pow(N_B * det_arr[idet] + bg_arr[idet], 2);
    }
    double deriv2_f = sum;
    return(deriv2_f);
}

void GetDetArr(const double* const rho_arr,
               const double* const resp_mat_arr,
               int ndet, int nsky,
               double* const out_arr) // ndet
{
    // det_arr = R_mat %*% rho_arr
    char* transa = new char [1];
    strcpy(transa, "N");    
    // double* det_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, out_arr, 1);
}

void GetNextRhoArr(const double* const rho_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   const double* const bg_arr,
                   int ndet, int nsky,
                   double N_B,
                   double* const out_arr)
{
    // coeff
    double B = GetB(bg_arr, ndet);
    double num = 1.0 - B / N_B;
    double* det_arr = new double[ndet];
    GetDetArr(rho_arr, resp_mat_arr, ndet, nsky, det_arr);
    double den = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        den += data_arr[idet] * det_arr[idet] /
            (det_arr[idet] + bg_arr[idet] / N_B);
    }
    double coeff = num / den;

    // mval_arr = t(R_mat) %*% (data_arr / det_arr) * rho_arr
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / (det_arr[idet] + bg_arr[idet] / N_B);
    }
    char* transa = new char [1];    
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, out_arr, 1);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = out_arr[isky] * rho_arr[isky];
        out_arr[isky] *= coeff;
    }
    
    delete [] transa;
    delete [] det_arr;
    delete [] div_arr;
}

void RichlucyBg(const double* const rho_arr, int nph,
                const double* const data_arr,
                const double* const resp_mat_arr,
                const double* const bg_arr,
                int niter,
                string outdir, string outfile_head,
                int ndet, int nskyx, int nskyy,
                double* const out_arr)
{
    int nsky = nskyx * nskyy;
    double* rho_new_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_new_arr, 1);
    double N_B = GetN(rho_new_arr, nsky) + GetB(bg_arr, ndet);
    printf("N_B = %e\n", N_B);
    for(int iiter = 0; iiter < niter; iiter ++){
        printf("iiter = %d\n", iiter);
        N_B = GetNextNb(rho_new_arr, data_arr, resp_mat_arr,
                        bg_arr, ndet, nsky, N_B);
        printf("N_B = %e\n", N_B);
        int nem = 50;
        for(int iem = 0; iem < nem; iem ++){
            double* rho_next_arr = new double[nsky];
            GetNextRhoArr(rho_new_arr, data_arr, resp_mat_arr, bg_arr,
                          ndet, nsky, N_B, rho_next_arr);
            dcopy_(nsky, const_cast<double*>(rho_next_arr), 1, rho_new_arr, 1);
            delete [] rho_next_arr;
        }
    }
    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, out_arr, 1);
    delete [] rho_new_arr;
}


void LoadResp(string respdir, int nskyx, int nskyy,
              double epsilon,
              double** const mat_arr_ptr,
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
    double* mat_arr = new double [ndet * nsky];
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            char infile[kLineSize];
            sprintf(infile, "%s/gimage_%3.3d_%3.3d.img", respdir.c_str(), iskyx, iskyy);

            int bitpix = 0;
            double* data_arr = NULL;
            MifFits::InFitsImageD(infile, img_info,
                                  &bitpix, &data_arr);
            double sum = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                if(data_arr[idet] < epsilon){
                    data_arr[idet] = epsilon;
                }
                sum += data_arr[idet];
            }
            
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;
            for(int idet = 0; idet < ndet; idet ++){
                mat_arr[imat + idet] = data_arr[idet] / sum;
            }
            delete [] data_arr;
        }
    }
    delete img_info;
    
    *mat_arr_ptr = mat_arr;
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
