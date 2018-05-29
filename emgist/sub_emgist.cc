#include "sub_emgist.h"

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

////

void SolveByEM(const double* const rho_arr, int nph,
               const double* const data_arr,
               const double* const resp_mat_arr,
               double beta, double mu, double lconst,
               double tol, double tol_em, int nstep, int flag_line_search,
               string outdir, string outfile_head,
               int ndet, int nskyx, int nskyy, double epsilon,
               int bitpix, double* const out_arr)
{
    char outfile_moni[kLineSize];
    char outfile_timelog[kLineSize];
    sprintf(outfile_moni, "%s/%s_moni.dat",
            outdir.c_str(), outfile_head.c_str());
    sprintf(outfile_timelog, "%s/%s_time_logl.dat",
            outdir.c_str(), outfile_head.c_str());
    FILE* fp_moni = fopen(outfile_moni, "w");
    FILE* fp_timelog = fopen(outfile_timelog, "w");
    setbuf(fp_moni, NULL);
    setbuf(fp_timelog, NULL);
    
    fprintf(fp_moni, "# istep, kldiv, logl, logl.inc, delta.logl, tdiff, lconst, diff_l_var\n");
    fprintf(fp_timelog, "# tdiff logl.inc\n");
    
    double time_st = MiTime::GetTimeSec();
    double logl_init = GetFuncL(rho_arr, data_arr, resp_mat_arr, beta, mu,
                                ndet, nskyx, nskyy);
    fprintf(fp_moni, "0  0  %.10e  0.0  0.0  0.0  0.0  0.0\n", logl_init);
    
    int nsky = nskyx * nskyy;
    double* rho_new_arr = new double[nsky];
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_new_arr, 1);
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);

    double* pLy_arr = new double[nsky];
    int nem = 1000;
    for(int iem = 0; iem < nem; iem ++){
        double* mval_arr = new double[nsky];
        GetFuncM(rho_new_arr, data_arr, resp_mat_arr, ndet, nsky, mval_arr);

        double lconst = 1.0e-5;
        int npm = 1000;
        for(int ipm = 0; ipm < npm; ipm ++){
            lconst /= 1.2;
            lconst = GetFindLconst(rho_new_arr, mval_arr, beta, mu,
                                   nskyx, nskyy, epsilon, lconst, pLy_arr);

            double logl_sub_pre = GetFuncLsub(rho_new_arr, mval_arr, beta, mu, nskyx, nskyy);
            dcopy_(nsky, const_cast<double*>(pLy_arr), 1, rho_new_arr, 1);
            double logl_sub = GetFuncLsub(rho_new_arr, mval_arr, beta, mu, nskyx, nskyy);

            double diff_logl_sub = logl_sub - logl_sub_pre;
            // printf("diff_logl_sub = %e\n", logl_sub - logl_sub_pre);
            if(fabs(diff_logl_sub) < 3.0e-1){
                break;
            }
        }

        //
        // image of differential
        //
        double* diff_l_arr = new double[nsky];
        GetDiffL(rho_new_arr, data_arr, resp_mat_arr,
                 beta, mu, ndet, nskyx, nskyy, diff_l_arr);
        if(iem % 10 == 0){
            int naxis = 2;
            long* naxes = new long[naxis];
            naxes[0] = nskyx;
            naxes[1] = nskyy;
            char tag[kLineSize];
            sprintf(tag, "%4.4d", iem);
            MifFits::OutFitsImageD(outdir, outfile_head + "_diff",
                                   tag, 2,
                                   bitpix,
                                   naxes, diff_l_arr);
            delete [] naxes;
        }
        // delete for zero
        vector<double> diff_l_vec;
        for(int isky = 0; isky < nsky; isky++){
            if(rho_new_arr[isky] > epsilon){
                diff_l_vec.push_back(diff_l_arr[isky]);
            }
        }
        delete [] diff_l_arr;
        double diff_l_stddev = MirMath::GetStddev(diff_l_vec);
        double diff_l_mean   = MirMath::GetAMean(diff_l_vec);
        double diff_l_var = diff_l_stddev / diff_l_mean;
//        printf("diff_l_variability: stddev / mean = %e, (stddev, mean) = (%e, %e) \n",
//               diff_l_var, diff_l_stddev, diff_l_mean);

        
        double kldiv    = GetKLDiv(rho_pre_arr, rho_new_arr, resp_mat_arr, ndet, nsky);
        double logl_pre = GetFuncL(rho_pre_arr, data_arr, resp_mat_arr, beta, mu, ndet, nskyx, nskyy);
        double logl     = GetFuncL(rho_new_arr, data_arr, resp_mat_arr, beta, mu, ndet, nskyx, nskyy);
        
        double delta_logl = logl - logl_pre;
        double logl_inc   = logl - logl_init;
        double time = MiTime::GetTimeSec();
        double tdiff = time - time_st;
        printf("iem = %d  kldiv = %e  logl = %.10e  logl - logl_init = %e "
               "delta_logl = %e  tdiff = %e  lconst = %e, diff_l_var = %e\n",
               iem, kldiv, logl, logl_inc, delta_logl, tdiff, lconst, diff_l_var);
        fprintf(fp_moni, "%d  %e  %.10e  %e  %e  %e  %e  %e\n",
                iem, kldiv, logl, logl_inc, delta_logl, tdiff, lconst, diff_l_var);
        fprintf(fp_timelog, "%e  %e\n", tdiff, logl_inc);

        if(iem % 10 == 0){
            int naxis = 2;
            long* naxes = new long[naxis];
            naxes[0] = nskyx;
            naxes[1] = nskyy;
            char tag[kLineSize];
            sprintf(tag, "%4.4d", iem);
            MifFits::OutFitsImageD(outdir, outfile_head,
                                   tag, 2,
                                   bitpix,
                                   naxes, rho_new_arr);
            delete [] naxes;
        }

        if( delta_logl > 10.0){
            printf("delta_logl (%e) > 10.0, then break.\n", delta_logl);
            break;
        }
        if(kldiv < tol && fabs(diff_l_var) < 1.0e-1){
            printf("kldiv (%e) < tol (%e)\n", kldiv, tol);
            dcopy_(nsky, rho_new_arr, 1, out_arr, 1);
            break;
        }

        // for next
        dcopy_(nsky, rho_new_arr, 1, rho_pre_arr, 1);
        dcopy_(nsky, rho_new_arr, 1, out_arr, 1);
    }

    fclose(fp_moni);
    fclose(fp_timelog);
    delete [] rho_new_arr;
    delete [] rho_pre_arr;
}

double GetTau(const double* const mval_arr,
              const double* const sigma_arr, int nsky,
              double lconst, double epsilon,
              double tau_init)
{
    double* tau_thres_arr = new double[nsky];
    GetFuncTauThres(sigma_arr, mval_arr, nsky, lconst, epsilon, tau_thres_arr);
    double tau_thres_min = 0.0;
    double tau_thres_max = 0.0;
    GetMinMax(tau_thres_arr, nsky, &tau_thres_min, &tau_thres_max);

    double tau = 0.0;    
    if(tau_init > tau_thres_min){
        tau = tau_init;
    } else {
        tau = tau_thres_max;
    }
    
    int nnewton = 100;
    double tol_newton = 1.0e-3;
    for(int inewton = 0; inewton < nnewton; inewton ++){
        tau = tau - GetFuncS(tau, sigma_arr, mval_arr, nsky, lconst, epsilon)
            / GetFuncDiffS(tau, sigma_arr, mval_arr, tau_thres_arr, nsky, lconst, epsilon);
        if( fabs(GetFuncS(tau, sigma_arr, mval_arr, nsky, lconst, epsilon) ) < tol_newton){
            break;
        }
    }
    return(tau);
}



double GetFindLconst(const double* const rho_arr,
                     const double* const mval_arr,
                     double beta, double mu,
                     int nskyx, int nskyy, double epsilon,
                     double lconst_init,
                     double* const pLy_arr)
{
    int nsky = nskyx * nskyy;
    double eta = 1.2;
    int ik_max = 1000;
    double lconst = lconst_init;
    double tau_pre = 0.0;
    for(int ik = 0; ik < ik_max; ik ++){
        lconst *= eta;

        double* sigma_arr = new double[nsky];
        GetFuncSigma(rho_arr, beta, mu, lconst, nskyx, nskyy, sigma_arr);
        double tau = GetTau(mval_arr, sigma_arr, nsky, lconst, epsilon, tau_pre);
        GetFuncRho(tau, sigma_arr, mval_arr, nsky, lconst, epsilon, pLy_arr);
        double qminusf = GetQMinusF(pLy_arr, rho_arr,
                                    beta, mu, lconst,
                                    nskyx, nskyy);
        //printf("GetFindLconst: ik: (L, qminusf) = %d (%e, %e)\n",
        //       ik, lconst, qminusf);
        if(qminusf >= 0.0){
            break;
        }
        tau_pre = tau;
    }
    return(lconst);
}


//
//
//

double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu,
                int ndet, int nskyx, int nskyy)
{
    int nsky = nskyx * nskyy;
    
    // term1 = - sum_v [ Y(v) log( sum_u t(v,u) rho_u ) ]
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);
    for(int idet = 0; idet < ndet; idet ++){
        det_arr[idet] = log(det_arr[idet]);
    }
    double term1 = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1,
                                const_cast<double*>(det_arr), 1);
    delete [] det_arr;
    delete [] transa;

    // term2 = (1 - beta) sum_u log rho_u
    double term2 = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        term2 += log(rho_arr[isky]);
    }
    term2 *= (1.0 - beta);

    // term3 = mu V(rho)
    double term3 = mu * GetTermV(rho_arr, nskyx, nskyy);

    double ans = term1 + term2 + term3;
    return(ans);
}

void GetDiffL(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              double beta, double mu,
              int ndet, int nskyx, int nskyy,
              double* const out_arr)
{
    int nsky = nskyx * nskyy;

    // det_arr = R_mat %*% rho_arr    
    // term1_diff = -1 * t(R_mat) %*% (data_arr / det_arr)
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / det_arr[idet];
    }

    double* term1_diff_arr = new double[nsky];
    strcpy(transa, "T");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, term1_diff_arr, 1);
    for(int isky = 0; isky < nsky; isky++){
        term1_diff_arr[isky] = -1 * term1_diff_arr[isky];
    }

    delete [] det_arr;
    delete [] transa;
    delete [] div_arr;

    // term2_diff = (1 - beta) / rho_u
    double* term2_diff_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky++){
        term2_diff_arr[isky] = (1.0 - beta) / rho_arr[isky];
    }

    // term3_diff = mu * diff_v
    double* term3_diff_arr = new double[nsky];
    GetDiffTermV(rho_arr, nskyx, nskyy, term3_diff_arr);

    for(int isky = 0; isky < nsky; isky++){
        out_arr[isky] = term1_diff_arr[isky] + term2_diff_arr[isky] + term3_diff_arr[isky];
    }
    delete [] term1_diff_arr;
    delete [] term2_diff_arr;
    delete [] term3_diff_arr;
}


double GetFuncLsub(const double* const rho_arr,
                   const double* const mval_arr,
                   double beta, double mu,
                   int nskyx, int nskyy)
{
    int nsky = nskyx * nskyy;
    
    // term1 = - sum_u [ mval_u * log( rho_u ) ]
    double* log_rho_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        log_rho_arr[isky] = log(rho_arr[isky]);
    }
    double term1 = -1.0 * ddot_(nsky, const_cast<double*>(mval_arr), 1,
                                const_cast<double*>(log_rho_arr), 1);
    // term2 = (1 - beta) sum_u log(rho_u)
    double term2 = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        term2 += log_rho_arr[isky];
    }
    term2 *= (1.0 - beta);
    delete [] log_rho_arr;
    
    // term3 = mu V(rho)
    double term3 = mu * GetTermV(rho_arr, nskyx, nskyy);

    double ans = term1 + term2 + term3;
    return(ans);
}


void GetFuncSigma(const double* const rho_arr,
                  double beta, double mu,
                  double lconst, int nskyx, int nskyy,
                  double* out_arr)
{
    int nsky = nskyx * nskyy;
    double* rho_diff_sparse_arr = new double[nsky];
    GetDiffSparse(rho_arr, nsky, beta, rho_diff_sparse_arr);
    double* rho_diff_v_arr = new double[nsky];
    GetDiffTermV(rho_arr, nskyx, nskyy, rho_diff_v_arr);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = rho_arr[isky] 
            - (rho_diff_sparse_arr[isky] + mu * rho_diff_v_arr[isky]) / lconst;
    }
    delete [] rho_diff_sparse_arr;
    delete [] rho_diff_v_arr;
}

void GetDiffSparse(const double* const rho_arr, int nsky,
                   double beta, double* const out_arr)
{
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = (1.0 - beta) / rho_arr[isky];
    }
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


void GetFuncM(const double* const rho_arr,
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
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / det_arr[idet];
    }
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, out_arr, 1);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = out_arr[isky] * rho_arr[isky];
    }
    delete [] transa;
    delete [] det_arr;
    delete [] div_arr;
}



double GetFuncS(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon)
{
    double* rho_arr = new double[nsky];
    GetFuncRho(tau, sigma_arr, mval_arr, nsky, lconst, epsilon, rho_arr);
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum += rho_arr[isky];
    }
    double ans = sum - 1.0;
    delete [] rho_arr;
    return(ans);
}

void GetFuncRho(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon,
                double* const out_arr)
{
    for(int isky = 0; isky < nsky; isky ++){
        double termb = sigma_arr[isky] + tau / lconst;
        out_arr[isky] = ( termb + sqrt( termb * termb + 4.0 * mval_arr[isky] / lconst ) ) / 2.0;
        if(out_arr[isky] <= epsilon){
            out_arr[isky] = epsilon;
        }
    }
}

double GetFuncDiffS(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    const double* const tau_thres_arr,
                    int nsky, double lconst, double epsilon)
{
    double* diff_rho_arr = new double[nsky];
    GetFuncDiffRho(tau, sigma_arr, mval_arr, tau_thres_arr,
                   nsky, lconst, epsilon,
                   diff_rho_arr);
    double ans = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        ans += diff_rho_arr[isky];
    }
    delete [] diff_rho_arr;
    return(ans);
}

void GetFuncDiffRho(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    const double* const tau_thres_arr,
                    int nsky, double lconst, double epsilon,
                    double* const out_arr)
{
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = 0.0;
        if(tau <= tau_thres_arr[isky]){
            out_arr[isky] = 0.0;
        } else {
            double termb = sigma_arr[isky] + tau / lconst;
            double root = sqrt( termb * termb + 4.0 * mval_arr[isky] / lconst );
            out_arr[isky] = (termb + root) / (2.0 * lconst * root);
        }
    }
}

void GetFuncTauThres(const double* const sigma_arr,
                     const double* const mval_arr,
                     int nsky, double lconst, double epsilon,
                     double* const out_arr)
{
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = lconst * (epsilon - sigma_arr[isky]) - mval_arr[isky] / epsilon;
    }
}

/////////////////////////////////////////////////////////////////




double GetQMinusF(const double* const rho_new_arr,
                  const double* const rho_arr,
                  double beta, double mu, double lconst,
                  int nskyx, int nskyy)
{
    //term1 =      FuncF(rho.vec, beta, mu, nrow, ncol)
    //term2 = -1 * FuncF(rho.new.vec, beta, mu, nrow, ncol)
    
    double term1 = GetFuncF(rho_arr, beta, mu, nskyx, nskyy);
    double term2 = -1 * GetFuncF(rho_new_arr, beta, mu, nskyx, nskyy);

    // term3 = sum( (rho.new.vec - rho.vec) * DiffF(rho.vec, beta, mu, nrow, ncol) )
    // term4 = L / 2.0 * sum( (rho.new.vec - rho.vec) * (rho.new.vec - rho.vec) )
    int nsky = nskyx * nskyy;
    double* diff_rho_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, diff_rho_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(rho_arr), 1, diff_rho_arr, 1);
    double* diff_f_arr = new double[nsky];
    GetDiffF(rho_arr, beta, mu, nskyx, nskyy, diff_f_arr);
    double term3 = ddot_(nsky, const_cast<double*>(diff_rho_arr), 1, const_cast<double*>(diff_f_arr), 1);
    double term4 = lconst *
        ddot_(nsky, const_cast<double*>(diff_rho_arr), 1, const_cast<double*>(diff_rho_arr), 1) / 2.0;
    double ans = term1 + term2 + term3 + term4;
    delete [] diff_rho_arr;
    delete [] diff_f_arr;
    return(ans);
}

double GetFuncF(const double* const rho_arr,
                double beta, double mu,
                int nskyx, int nskyy)
{
    int nsky = nskyx * nskyy;
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum += log(rho_arr[isky]);
    }
    double ans = (1.0 - beta) * sum + mu * GetTermV(rho_arr, nskyx, nskyy);
    return(ans);
}

void GetDiffF(const double* const rho_arr,
              double beta, double mu,
              int nskyx, int nskyy,
              double* const out_arr)
{
    // term1 = (1 - beta) / rho.vec
    int nsky = nskyx * nskyy;
    double* rho_inv_arr = new double[nsky];
    GetInvArr(rho_arr, nsky, rho_inv_arr);

    // term2 = mu * DiffTermV(rho.vec, nrow, ncol)
    GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);

    // ans = term1 + term2
    daxpy_(nsky, 1.0 - beta, rho_inv_arr, 1, out_arr, 1);
    delete [] rho_inv_arr;
}

double GetFuncG(const double* const rho_arr,
                const double* const mval_arr,
                int nsky)
{
    double ans = 0.0;
    for(int isky = 0; isky < nsky; isky++){
        ans += mval_arr[isky] * log(rho_arr[isky]);
    }
    ans *= -1;
    return(ans);
}

void GetDiffG(const double* const rho_arr,
              const double* const mval_arr,
              int nsky,
              double* const out_arr)
{
    for(int isky = 0; isky < nsky; isky++){
        out_arr[isky] = -1 * mval_arr[isky] / rho_arr[isky];
    }
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
