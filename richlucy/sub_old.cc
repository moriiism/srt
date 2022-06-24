#include "sub.h"

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

void Richlucy(const double* const rho_arr, int nph,
              const double* const data_arr,
              const double* const resp_mat_arr,
              int nem, double tol_em,
              double tol_diff_l_var,
              int flag_line_search,
              string outdir, string outfile_head,
              int ndet, int nskyx, int nskyy, double epsilon,
              int bitpix, double* const out_arr)
{
    char outfile_moni[kLineSize];
    char qdpfile_timelog[kLineSize];
    char qdpfile_delta_logl[kLineSize];
    char qdpfile_kldiv[kLineSize];
    char qdpfile_diff_l_var[kLineSize];
    
    sprintf(outfile_moni, "%s/%s_moni.dat",
            outdir.c_str(), outfile_head.c_str());
    sprintf(qdpfile_timelog, "%s/%s_time_logl.qdp",
            outdir.c_str(), outfile_head.c_str());
    sprintf(qdpfile_delta_logl, "%s/%s_delta_logl.qdp",
            outdir.c_str(), outfile_head.c_str());
    sprintf(qdpfile_kldiv, "%s/%s_kldiv.qdp",
            outdir.c_str(), outfile_head.c_str());
    sprintf(qdpfile_diff_l_var, "%s/%s_diff_l_var.qdp",
            outdir.c_str(), outfile_head.c_str());        
    
    FILE* fp_moni = fopen(outfile_moni, "w");
    FILE* fp_timelog = fopen(qdpfile_timelog, "w");
    FILE* fp_delta_logl = fopen(qdpfile_delta_logl, "w");
    FILE* fp_kldiv = fopen(qdpfile_kldiv, "w");
    FILE* fp_diff_l_var = fopen(qdpfile_diff_l_var, "w");
    setbuf(fp_moni, NULL);
    setbuf(fp_timelog, NULL);
    setbuf(fp_delta_logl, NULL);
    setbuf(fp_kldiv, NULL);
    setbuf(fp_diff_l_var, NULL);    
    
    fprintf(fp_moni, "# iem, nzero, kldiv, helldist, logl, logl.inc, delta.logl, tdiff, lconst, diff_l_var\n");
    fprintf(fp_timelog, "# tdiff logl.inc\n");
    fprintf(fp_delta_logl, "skip sing\n");
    fprintf(fp_kldiv, "skip sing\n");
    fprintf(fp_diff_l_var, "skip sing\n");        
    
    double time_st = MiTime::GetTimeSec();
    double logl_init = GetFuncL(rho_arr, data_arr, resp_mat_arr,
                                ndet, nskyx, nskyy, epsilon);
    
    int nsky = nskyx * nskyy;
    double* rho_new_arr = new double[nsky];
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_new_arr, 1);
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);


    for(int iem = 0; iem < nem; iem ++){
        double* mval_arr = new double[nsky];
        GetNextRhoArr(rho_new_arr, data_arr, resp_mat_arr, ndet, nsky, mval_arr);
        dcopy_(nsky, mval_arr, 1, rho_new_arr, 1);

        double kldiv    = GetKLDiv(rho_pre_arr, rho_new_arr, resp_mat_arr, ndet, nsky);
        double logl_pre = GetFuncL(rho_pre_arr, data_arr, resp_mat_arr,
                                   ndet, nskyx, nskyy, epsilon);
        double logl     = GetFuncL(rho_new_arr, data_arr, resp_mat_arr,
                                   ndet, nskyx, nskyy, epsilon);
        
        double delta_logl = logl - logl_pre;
        double logl_inc   = logl - logl_init;
        double time = MiTime::GetTimeSec();
        double tdiff = time - time_st;
        printf("iem = %d  kldiv = %e  logl = %.10e  logl - logl_init = %e "
               "delta_logl = %e  tdiff = %e \n",
               iem, kldiv, logl, logl_inc, delta_logl, tdiff);
        fprintf(fp_moni, "%d  %e  %.10e  %e  %e  %e  \n",
                iem, kldiv, logl, logl_inc, delta_logl, tdiff);
        fprintf(fp_timelog, "%e  %e\n", tdiff, logl_inc);
        fprintf(fp_delta_logl, "%d  %e\n", iem, delta_logl);
        fprintf(fp_kldiv, "%d  %e\n", iem, kldiv);

        if(kldiv < tol_em){
            printf("kldiv (%e) < tol_em (%e)\n", kldiv, tol_em);
            dcopy_(nsky, rho_new_arr, 1, out_arr, 1);

            char outfile_last_head[kLineSize];
            sprintf(outfile_last_head, "%s/%s_last_head.dat",
                    outdir.c_str(), outfile_head.c_str());
            FILE* fp_last_head = fopen(outfile_last_head, "w");
            fprintf(fp_last_head, "# mu  beta  iem  nzero  kldiv  helldist logl  logl_inc  delta_logl  tdiff  lconst  diff_l_var\n");
            fclose(fp_last_head);
            
            // last
            char outfile_last[kLineSize];
            sprintf(outfile_last, "%s/%s_last.dat",
                    outdir.c_str(), outfile_head.c_str());
            FILE* fp_last = fopen(outfile_last, "w");
            fprintf(fp_last, "%d  %e  %.10e  %e  %e  %e\n",
                    iem, kldiv, logl, logl_inc, delta_logl, tdiff);
            fclose(fp_last);
            break;
        }
        
        // for next
        dcopy_(nsky, rho_new_arr, 1, rho_pre_arr, 1);
        dcopy_(nsky, rho_new_arr, 1, out_arr, 1);

        delete [] mval_arr;
    }

    fclose(fp_moni);
    fclose(fp_timelog);
    fclose(fp_delta_logl);
    fclose(fp_kldiv);
    fclose(fp_diff_l_var);
    delete [] rho_new_arr;
    delete [] rho_pre_arr;
}


double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                int ndet, int nskyx, int nskyy, double epsilon)
{
    int nsky = nskyx * nskyy;

    double* rho_epsilon_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        if(rho_arr[isky] < epsilon){
            rho_epsilon_arr[isky] = epsilon;
        } else {
            rho_epsilon_arr[isky] = rho_arr[isky];
        }
    }
    
    // term1 = - sum_v [ Y(v) log( sum_u t(v,u) rho_u ) ]
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_epsilon_arr), 1,
           0.0, det_arr, 1);
    for(int idet = 0; idet < ndet; idet ++){
        det_arr[idet] = log(det_arr[idet]);
    }
    double term1 = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1,
                                const_cast<double*>(det_arr), 1);
    delete [] det_arr;
    delete [] transa;

    double ans = term1;

    delete [] rho_epsilon_arr;
    return(ans);
}


double GetKLDiv(const double* const rho_arr,
                const double* const rho_new_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky)
{
    char* transa = new char [1];
    strcpy(transa, "N");
    // q.vec = R.mat %*% y.vec
    double* q_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, q_arr, 1);
    // q.new.vec = R.mat %*% y.new.vec
    double* q_new_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_new_arr), 1,
           0.0, q_new_arr, 1);

    delete [] transa;

    // q.vec = q.vec / sum(q.vec)
    // q.new.vec = q.new.vec / sum(q.new.vec)
    double sum_q = 0.0;
    double sum_q_new = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        sum_q += q_arr[idet];
        sum_q_new += q_new_arr[idet];
    }
    dscal_(ndet, 1.0/sum_q, q_arr, 1);
    dscal_(ndet, 1.0/sum_q_new, q_new_arr, 1);
    
    double ans = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        if(q_new_arr[idet] > 0.0){
            ans = ans + q_new_arr[idet] * log( q_new_arr[idet] / q_arr[idet] );
        }
    }
    delete [] q_arr;
    delete [] q_new_arr;
    return (ans);
}

void GetNextRhoArr(const double* const rho_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   int ndet, int nsky,
                   double* const out_arr)
{
    // det_arr = R_mat %*% rho_arr
    char* transa = new char [1];
    strcpy(transa, "N");    
    double* det_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);

    // mval_arr = t(R_mat) %*% (data_arr / det_arr) * rho_arr
    double* div_arr = new double[ndet];
    double sum = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / det_arr[idet];
        sum += data_arr[idet];
    }
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, out_arr, 1);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = out_arr[isky] * rho_arr[isky];
        out_arr[isky] /= sum;
    }
    
    delete [] transa;
    delete [] det_arr;
    delete [] div_arr;
}


////////////////////////////////////////////////

double GetSumLogPlus(const double* const data_arr, int ndata,
                     const int* const supp_arr)
{
    double ans = 0.0;
    for(int idata = 0; idata < ndata; idata ++){
        if(supp_arr[idata] > 0){
            ans = ans + log(data_arr[idata]);
        }
    }
    return(ans);
}

void GetInvArr(const double* const data_arr, int ndata,
               double* const out_arr)
{
    for(int idata = 0; idata < ndata; idata ++){
        out_arr[idata] = 1.0 / data_arr[idata];
    }
}

void GetInvArrSupp(const double* const data_arr, int ndata,
                   const int* const supp_arr, 
                   double* const out_arr)
{
    for(int idata = 0; idata < ndata; idata ++){
        if(supp_arr[idata] > 0){
            out_arr[idata] = 1.0 / data_arr[idata];
        } else {
            out_arr[idata] = 0.0;
        }
    }
}
    
void GetSuppArr(const double* const data_arr, int ndata,
                double epsilon,
                int* const out_arr)
{
    for(int idata = 0; idata < ndata; idata ++){
        out_arr[idata] = 0;
        if(data_arr[idata] > epsilon){
            out_arr[idata] = 1;
        }
    }
}

int GetNZeroArr(const double* const data_arr, int ndata,
                double epsilon)
{
    int nzero = 0;
    for(int idata = 0; idata < ndata; idata ++){
        if(fabs(data_arr[idata]) <= epsilon){
            nzero = nzero + 1;
        }
    }
    return(nzero);
}

int GetIbin(int ibinx, int ibiny, int nbinx)
{
    int ibin = ibinx + ibiny * nbinx;
    return(ibin);
}

void GetMinMax(const double* const data_arr, int ndata,
               double* min_ptr, double* max_ptr)
{
    double min = data_arr[0];
    double max = data_arr[0];
    for(int idata = 0; idata < ndata; idata ++){
        if(min > data_arr[idata]){
            min = data_arr[idata];
        }
        if(max < data_arr[idata]){
            max = data_arr[idata];
        }
    }
    *min_ptr = min;
    *max_ptr = max;
}

void GetMinMaxSupp(const double* const data_arr,
                   const int* const index_supp_arr,
                   int ndata,
                   double* min_ptr, double* max_ptr)
{
    double min = 0.0;
    double max = 0.0;
    for(int idata = 0; idata < ndata; idata ++){
        if(1 == index_supp_arr[idata]){
            min = data_arr[idata];
            max = data_arr[idata];
            break;
        }
    }
    for(int idata = 0; idata < ndata; idata ++){
        if(1 == index_supp_arr[idata]){
            if(min > data_arr[idata]){
                min = data_arr[idata];
            }
            if(max < data_arr[idata]){
                max = data_arr[idata];
            }
        }
    }
    *min_ptr = min;
    *max_ptr = max;
}

