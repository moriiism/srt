#include "sub_dca.h"

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

void SolveByDca(const double* const rho_arr, int nph,
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
    double logl_pre = logl_init;
    fprintf(fp_moni, "0  0  %.10e  0.0  0.0  0.0  0.0  0.0\n", logl_init);
    
    double eta = 1.2;
    int nsky = nskyx * nskyy;
    double* rho_new_arr = new double[nsky];
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_new_arr, 1);
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    for(int istep = 0; istep < nstep; istep ++){

        if(istep == 0){
            int nfind = 3;
            for(int ifind = 0; ifind < nfind; ifind ++){
                int flag_find = 0;
                lconst = GetFindIkBisect(rho_new_arr, data_arr, resp_mat_arr,
                                         beta, mu, lconst, eta,
                                         tol_em, flag_line_search,
                                         ndet, nskyx, nskyy, epsilon,
                                         &flag_find);
                if(1 == flag_find){
                    break;
                }
            }
        } else {
            
            double* pLy_arr = new double[nsky];
            int flag_good = 0;
            GetProxMap(rho_new_arr, data_arr, resp_mat_arr,
                       beta, mu, lconst,
                       ndet, nskyx, nskyy, epsilon,
                       tol_em, flag_line_search,
                       pLy_arr, &flag_good);
            double qminusf = GetQMinusF(pLy_arr, rho_new_arr,
                                        beta, mu, lconst,
                                        nskyx, nskyy);
            delete [] pLy_arr;
            printf("(L, qminusf, flag_good) = (%e, %e, %d)\n",
                   lconst, qminusf, flag_good);
            if(qminusf < 0){
                int nfind = 3;
                for(int ifind = 0; ifind < nfind; ifind ++){
                    int flag_find = 0;
                    lconst = GetFindIkBisect(rho_new_arr, data_arr, resp_mat_arr,
                                             beta, mu, lconst, eta,
                                             tol_em, flag_line_search,
                                             ndet, nskyx, nskyy, epsilon,
                                             &flag_find);
                    if(1 == flag_find){
                        break;
                    }
                }
            }
            
            lconst /= pow(eta, 1);            
            int ik = GetFindIk(rho_new_arr, data_arr, resp_mat_arr,
                               beta, mu, lconst, eta,
                               tol_em, flag_line_search,
                               ndet, nskyx, nskyy, epsilon);
            lconst = pow(eta, ik) * lconst;
        }

        
        printf("SolveByProxMap: L = %e\n", lconst);
        int flag_good = 0;

        double* rho_pmout_arr = new double[nsky];
        GetProxMap(rho_new_arr, data_arr, resp_mat_arr,
                   beta, mu, lconst,
                   ndet, nskyx, nskyy, epsilon,
                   tol_em, flag_line_search,
                   rho_pmout_arr, &flag_good);

        printf("SolveByProxMap: output of GetProxMap : flag_good = %d\n", flag_good);
        
        dcopy_(nsky, const_cast<double*>(rho_pmout_arr), 1, rho_new_arr, 1);
        delete [] rho_pmout_arr;


        //
        // image of differential
        //
        double* diff_l_arr = new double[nsky];
        GetDiffL(rho_new_arr, data_arr, resp_mat_arr,
                 beta, mu, ndet, nskyx, nskyy, diff_l_arr);
        {
            int naxis = 2;
            long* naxes = new long[naxis];
            naxes[0] = nskyx;
            naxes[1] = nskyy;
            char tag[kLineSize];
            sprintf(tag, "%4.4d", istep);
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
        printf("diff_l_variability: stddev /mean = %e, (stddev, mean) = (%e, %e) \n",
               diff_l_var, diff_l_stddev, diff_l_mean);

    
        
        double kldiv = GetKLDiv(rho_pre_arr, rho_new_arr, resp_mat_arr, ndet, nsky);
        logl_pre    = GetFuncL(rho_pre_arr, data_arr, resp_mat_arr, beta, mu, ndet, nskyx, nskyy);
        double logl = GetFuncL(rho_new_arr, data_arr, resp_mat_arr, beta, mu, ndet, nskyx, nskyy);
        
        double delta_logl = logl - logl_pre;
        double logl_inc   = logl - logl_init;
        double time = MiTime::GetTimeSec();
        double tdiff = time - time_st;
        printf("%d  kldiv = %e  logl = %.10e  logl - logl_init = %e "
               "delta_logl = %e  tdiff = %e  lconst = %e, diff_l_var = %e\n",
               istep, kldiv, logl, logl_inc, delta_logl, tdiff, lconst, diff_l_var);
        fprintf(fp_moni, "%d  %e  %.10e  %e  %e  %e  %e  %e\n",
                istep, kldiv, logl, logl_inc, delta_logl, tdiff, lconst, diff_l_var);
        fprintf(fp_timelog, "%e  %e\n", tdiff, logl_inc);


        
        if(istep % 10 == 0){
            int naxis = 2;
            long* naxes = new long[naxis];
            naxes[0] = nskyx;
            naxes[1] = nskyy;
            char tag[kLineSize];
            sprintf(tag, "%4.4d", istep);
            MifFits::OutFitsImageD(outdir, outfile_head,
                                   tag, 2,
                                   bitpix,
                                   naxes, rho_new_arr);
            delete [] naxes;
        }
        
        if( delta_logl > 1.0e-1){
            printf("delta_logl (%e) > 0.0, then break.\n", delta_logl);
            break;
        }
        if(kldiv < tol){
            if(fabs(diff_l_var) > 1.0e-2){
                double* rho_new2_arr = new double[nsky];
                MinToDiffL(rho_new_arr, data_arr, resp_mat_arr,
                           beta, mu, ndet, nskyx, nskyy, epsilon, rho_new2_arr);
                dcopy_(nsky, rho_new2_arr, 1, rho_new_arr, 1);
                
                lconst = 1.0;
                tol_em /= 2.0;
                tol    /= 2.0;
                if(tol_em < 1e-15){
                    tol_em = 1e-15;
                    tol    = 1e-15;
                }
                printf("L search again. tol_em = %e\n", tol_em);
            } else {
                printf("kldiv (%e) < tol (%e)\n", kldiv, tol);
                dcopy_(nsky, rho_new_arr, 1, out_arr, 1);
                break;
            }
        }

//        if(lconst > 1.0e+10){
//            lconst = 1.0;
//            printf("L > 1.0e+10, then L search again.\n");
//        }

        // for next
        dcopy_(nsky, rho_new_arr, 1, rho_pre_arr, 1);
        dcopy_(nsky, rho_new_arr, 1, out_arr, 1);
        logl_pre = logl;

    }
    
    fclose(fp_moni);
    fclose(fp_timelog);
    delete [] rho_new_arr;
    delete [] rho_pre_arr;
}


int GetFindIk(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              double beta, double mu, double lconst, double eta,
              double tol_em, int flag_line_search,
              int ndet, int nskyx, int nskyy, double epsilon)
{
    int ik_max = 100;
    double lconst_org = lconst;
    int nsky = nskyx * nskyy;
    int ik_new = 0;
    for(int ik = 0; ik < ik_max; ik ++){
        lconst = pow(eta, ik) * lconst_org;
        double* pLy_arr = new double[nsky];
        int flag_good = 0;
        GetProxMap(rho_arr, data_arr, resp_mat_arr,
                   beta, mu, lconst,
                   ndet, nskyx, nskyy, epsilon,
                   tol_em, flag_line_search,
                   pLy_arr, &flag_good);
        double qminusf = GetQMinusF(pLy_arr, rho_arr,
                                    beta, mu, lconst,
                                    nskyx, nskyy);
        delete [] pLy_arr;
        printf("FindIk: (ik, L, qminusf, flag_good) = (%d, %e, %e, %d)\n",
               ik, lconst, qminusf, flag_good);
        if(qminusf >= 0.0 && flag_good == 1){
            ik_new = ik;
            break;
        }
    }
    return(ik_new);
}


double GetFindIkBisect(const double* const rho_arr,
                       const double* const data_arr,
                       const double* const resp_mat_arr,
                       double beta, double mu, double lconst, double eta,
                       double tol_em, int flag_line_search,
                       int ndet, int nskyx, int nskyy, double epsilon,
                       int* const flag_find_ptr)
{
    int nsky = nskyx * nskyy;

    double log10_lconst_0 = 0.0;
    double log10_lconst = log10(lconst);
    double log10_lconst_1 = log10_lconst + 5;
    double lconst_good = pow(10, log10_lconst_1);
    int flag_find = 0;

    // check small side
    {
        double lconst_0 = pow(10, log10_lconst_0);
        double* pLy_arr = new double[nsky];
        int flag_good = 0;
        GetProxMap(rho_arr, data_arr, resp_mat_arr,
                   beta, mu, lconst_0,
                   ndet, nskyx, nskyy, epsilon,
                   tol_em, flag_line_search,
                   pLy_arr, &flag_good);
        double qminusf = GetQMinusF(pLy_arr, rho_arr,
                                    beta, mu, lconst,
                                    nskyx, nskyy);
        delete [] pLy_arr;
        if(qminusf > 0.0){
            flag_find = 1;
            *flag_find_ptr = flag_find;
            return(lconst_0);
        }
    }
    
    {
        double lconst_1 = pow(10, log10_lconst_1);
        double* pLy_arr = new double[nsky];
        int flag_good = 0;
        GetProxMap(rho_arr, data_arr, resp_mat_arr,
                   beta, mu, lconst_1,
                   ndet, nskyx, nskyy, epsilon,
                   tol_em, flag_line_search,
                   pLy_arr, &flag_good);
        double qminusf = GetQMinusF(pLy_arr, rho_arr,
                                    beta, mu, lconst,
                                    nskyx, nskyy);
        delete [] pLy_arr;
        if(qminusf < 0.0){
            flag_find = 0;
            *flag_find_ptr = flag_find;
            return(lconst_1);
        }
    }
    
    int nstep = 20;
    for(int istep = 0; istep < nstep; istep ++){
        lconst = pow(10, log10_lconst);
        double* pLy_arr = new double[nsky];
        int flag_good = 0;
        GetProxMap(rho_arr, data_arr, resp_mat_arr,
                   beta, mu, lconst,
                   ndet, nskyx, nskyy, epsilon,
                   tol_em, flag_line_search,
                   pLy_arr, &flag_good);
        double qminusf = GetQMinusF(pLy_arr, rho_arr,
                                    beta, mu, lconst,
                                    nskyx, nskyy);
        delete [] pLy_arr;
        printf("FindIk: (L, log10_lconst_0 (lconst0), log10_lconst_1 (lconst1), qminusf, flag_good) "
               " = (%e, %e(%e), %e(%e), %e, %d)\n",
               lconst,
               log10_lconst_0, pow(10, log10_lconst_0),
               log10_lconst_1, pow(10, log10_lconst_1),
               qminusf, flag_good);
        
        if(qminusf > 1.0){
            if(0 == flag_good){
                log10_lconst += (log10_lconst_1 - log10_lconst_0) * 0.1;
            } else {
                lconst_good = pow(10, log10_lconst);
                log10_lconst_1 = log10_lconst;
                log10_lconst -= (log10_lconst_1 - log10_lconst_0) * 0.5;
            }
        } else if(qminusf <= 0.0){
            if(0 == flag_good){
                log10_lconst += (log10_lconst_1 - log10_lconst_0) * 0.1;
            } else {
                log10_lconst_0 = log10_lconst;
                log10_lconst += (log10_lconst_1 - log10_lconst_0) * 0.5;
            }
        } else {
            if(0 == flag_good){
                log10_lconst += (log10_lconst_1 - log10_lconst_0) * 0.1;
            } else {
                lconst_good = pow(10, log10_lconst);
                flag_find = 1;
                break;
            }
        }
    }

    *flag_find_ptr = flag_find;
    return(lconst_good);
}





void GetProxMap(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu, double lconst,
                int ndet, int nskyx, int nskyy, double epsilon,
                double tol_em, int flag_line_search,
                double* const out_arr, int* flag_good_ptr)
{
    int nsky = nskyx * nskyy;
    double* sigma_arr = new double[nsky];
    GetFuncSigma(rho_arr, data_arr,
                 beta, mu, lconst, nskyx, nskyy,
                 sigma_arr);

    dcopy_(nsky, const_cast<double*>(rho_arr), 1, out_arr, 1);
    int nem = 1000;
    int flag_good = 1;
    double tau_pre = 1.0e-10;
    int flag_saturate = 0;
    for(int iem = 0; iem < nem; iem ++){
        double lem = GetFuncLem(out_arr, data_arr, resp_mat_arr, sigma_arr,
                                ndet, nsky, lconst);
        double* mval_arr = new double[nsky];
        GetFuncM(out_arr, data_arr, resp_mat_arr, ndet, nsky, mval_arr);
        double* tau_thres_arr = new double[nsky];
        GetFuncTauThres(sigma_arr, mval_arr, nsky, lconst, epsilon, tau_thres_arr);
        double tau_thres_min = 0.0;
        double tau_thres_max = 0.0;
        GetMinMax(tau_thres_arr, nsky, &tau_thres_min, &tau_thres_max);

    }

    delete [] sigma_arr;
    
    *flag_good_ptr = flag_good;
}

void GetLineSearch(const double* const xval_arr,
                   const double* const xval_new_arr,
                   const double* const data_arr,
                   const double* const resp_mat_arr,
                   const double* const sigma_arr,
                   int ndet, int nsky, double lconst,
                   double epsilon,
                   double* const out_arr,
                   int* flag_saturate_ptr)
{
    double xval0     = xval_arr[0];
    double xval0_new = xval_new_arr[0];
    double* xval2_arr = new double[nsky - 1];
    double* xval2_new_arr = new double[nsky - 1];
    double* theta_arr = new double[nsky - 1];
    double* theta_new_arr = new double[nsky - 1];
    
    for(int isky = 0; isky < nsky - 1; isky ++){
        xval2_arr[isky]     = xval_arr[isky + 1];
        xval2_new_arr[isky] = xval_new_arr[isky + 1];
    }
    for(int isky = 0; isky < nsky - 1; isky ++){        
        theta_arr[isky]     = log( xval2_arr[isky] / (1.0 - xval0) );
        theta_new_arr[isky] = log( xval2_new_arr[isky] / (1.0 - xval0_new) );
    }
    int nstep = 100;
    double lem_init = GetFuncLem(xval_arr, data_arr, resp_mat_arr, sigma_arr,
                                 ndet, nsky, lconst);
    double lem_pre = lem_init;
    double* xval_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(xval_arr), 1, xval_pre_arr, 1);

    char outfile[kLineSize];
    sprintf(outfile, "temp.dat"); 
    FILE* fp_out = fopen(outfile, "w");
    setbuf(fp_out, NULL);
    fprintf(fp_out, "skip sing\n");
    fprintf(fp_out, "read\n");

    int flag_saturate = 0;
    double eta = 3.0;
    for(int istep = 1; istep < nstep; istep ++){
        double factor = pow(eta, istep);
        double lxval0_this = factor * (log(xval0_new) - log(xval0)) + log(xval0);
        double xval0_this = exp(lxval0_this);

        double* theta_this_arr = new double[nsky - 1];
        for(int isky = 0; isky < nsky - 1; isky ++){
            theta_this_arr[isky] = factor * (theta_new_arr[isky] - theta_arr[isky]) + theta_arr[isky];
        }
        double alpha = 0.0;
        for(int isky = 0; isky < nsky - 1; isky ++){
            alpha += exp(theta_this_arr[isky]);
        }
        double* xval2_this_arr = new double[nsky - 1];
        for(int isky = 0; isky < nsky - 1; isky ++){
            xval2_this_arr[isky] = (1 - xval0_this) * exp(theta_this_arr[isky]) / alpha;
        }
        delete [] theta_this_arr;
        
        double* xval_this_arr = new double[nsky];
        xval_this_arr[0] = xval0_this;
        for(int isky = 0; isky < nsky - 1; isky ++){
            xval_this_arr[isky + 1] = xval2_this_arr[isky];
        }
        delete [] xval2_this_arr;

        double xval_this_min = 0.0;
        double xval_this_max = 0.0;
        GetMinMax(xval_this_arr, nsky, &xval_this_min, &xval_this_max);
        
        if(xval_this_min < epsilon){
            printf("xval_this_min < epsilon: factor(istep) = %e (%d)\n", factor, istep);
            if(istep != 1){
                dcopy_(nsky, const_cast<double*>(xval_pre_arr), 1, out_arr, 1);
            }
            flag_saturate = 1;
            break;
        }

        double lem = GetFuncLem(xval_this_arr, data_arr, resp_mat_arr, sigma_arr, ndet, nsky, lconst);
        fprintf(fp_out, "%d  %e\n", istep, lem - lem_init);

        if(lem_pre < lem){
            if(istep != 1){
                dcopy_(nsky, const_cast<double*>(xval_pre_arr), 1, out_arr, 1);
                // printf("linesearch step = %d\n", istep);
            }
            break;
        }
        lem_pre = lem;
        dcopy_(nsky, const_cast<double*>(xval_this_arr), 1, xval_pre_arr, 1);
        delete [] xval_this_arr;
    }

    delete [] xval2_arr;
    delete [] xval2_new_arr;
    delete [] theta_arr;
    delete [] theta_new_arr;
    delete [] xval_pre_arr;
    fclose(fp_out);

    *flag_saturate_ptr = flag_saturate;
}


void MinToDiffL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu,
                int ndet, int nskyx, int nskyy,
                double epsilon,
                double* const out_arr)
{
    int nsky = nskyx * nskyy;    
    double* diff_l_arr = new double[nsky];
    GetDiffL(rho_arr, data_arr, resp_mat_arr,
             beta, mu, ndet, nskyx, nskyy, diff_l_arr);
    double logl = GetFuncL(rho_arr, data_arr, resp_mat_arr,
                           beta, mu, ndet, nskyx, nskyy);
    double logl_min = logl;
    double eta = 1.0e-10;
    int nstep = 10;

    double* rho_new_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky++){
        out_arr[isky] = rho_arr[isky];
    }

    
    for(int istep = 0; istep < nstep; istep++){
        double sum_active = 0.0;
        for(int isky = 0; isky < nsky; isky++){
            if(rho_arr[isky] > epsilon){
                sum_active += rho_arr[isky];
            }
        }
        double sum_active_new = 0.0;
        for(int isky = 0; isky < nsky; isky++){
            double factor = istep * eta;
            if(rho_arr[isky] > epsilon){
                rho_new_arr[isky] = rho_arr[isky] - factor * diff_l_arr[isky];
                sum_active_new += rho_new_arr[isky];
            } else {
                rho_new_arr[isky] = epsilon;
            }
        }
        for(int isky = 0; isky < nsky; isky++){
            if(rho_arr[isky] > epsilon){
                rho_new_arr[isky] *= (sum_active / sum_active_new);
            }
        }
        int nzero = GetNZeroArr(rho_new_arr, nsky, epsilon);
        double logl_new = GetFuncL(rho_new_arr, data_arr, resp_mat_arr,
                                   beta, mu, ndet, nskyx, nskyy);
        printf("MinToDiffL: nzero: %d, logl_new, logl, logl_new - logl = (%e, %e, %e)\n",
               nzero, logl_new, logl, logl_new - logl);

        if(logl_new < logl_min){
            logl_min = logl_new;
            for(int isky = 0; isky < nsky; isky++){
                out_arr[isky] = rho_new_arr[isky];
            }
        }
    }
    delete [] diff_l_arr;
    delete [] rho_new_arr;
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

//
//
//

void GetFuncM(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              int ndet, int nsky,
              double* const out_arr)
{
    char* transa = new char [1];
    strcpy(transa, "N");    
    // num.vec = R.mat %*% rho.vec
    double* det_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);

    // ans.vec = t(R.mat) %*% (D.vec / num.vec) * rho.vec
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / det_arr[idet];
    }
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, out_arr, 1);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = out_arr[isky] * rho_arr[isky];
    }
    delete [] transa;
    delete [] det_arr;
    delete [] div_arr;
    
}

void GetFuncRho(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon,
                double* const out_arr)
{
    // termb.vec = sigma.vec + tau / L
    // ans.vec = ( termb.vec + sqrt( termb.vec * termb.vec + 4 * m.vec / L ) ) / 2.0
    // ans.vec = mapply(ThresEpsilon, ans.vec, epsilon)
    
    for(int isky = 0; isky < nsky; isky ++){
        double termb = sigma_arr[isky] + tau / lconst;
        out_arr[isky] = ( termb + sqrt( termb * termb + 4 * mval_arr[isky] / lconst ) ) / 2.0;
        if(out_arr[isky] < epsilon){
            out_arr[isky] = epsilon;
        }
    }
}


void GetFuncDiffRho(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double lconst, double epsilon,
                    double* const out_arr)
{
    double* tau_thres_arr = new double[nsky];
    GetFuncTauThres(sigma_arr, mval_arr, nsky, lconst, epsilon, tau_thres_arr);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = 0.0;
        if(tau <= tau_thres_arr[isky]){
            out_arr[isky] = 0.0;
        } else {
            double termb = sigma_arr[isky] + tau / lconst;
            double root = sqrt( termb * termb + 4 * mval_arr[isky] / lconst );
            out_arr[isky] = (termb + root) / (2 * lconst * root);
        }
    }
    delete [] tau_thres_arr;
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

double GetFuncS(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon)
{
    // ans = sum( FuncRho(tau, sigma.vec, m.vec, L, epsilon) ) - 1
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

double GetFuncDiffS(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double lconst, double epsilon)
{
    // ans = sum(FuncDiffRhoVec(tau, sigma.vec, m.vec, L, epsilon))
    
    double* diff_rho_arr = new double[nsky];
    GetFuncDiffRho(tau, sigma_arr, mval_arr,
                   nsky, lconst, epsilon,
                   diff_rho_arr);
    double ans = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        ans += diff_rho_arr[isky];
    }
    delete [] diff_rho_arr;
    return(ans);
}


//
//
//
 
double GetFuncLem(const double* const rho_arr,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const sigma_arr,
                  int ndet, int nsky,
                  double lconst)
{
    double* diff_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, diff_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(sigma_arr), 1, diff_arr, 1);
    double term1 = lconst
        * ddot_(nsky, const_cast<double*>(diff_arr), 1, const_cast<double*>(diff_arr), 1) / 2.0;
    double term2 = GetFuncG(rho_arr, data_arr, resp_mat_arr, ndet, nsky);
    double ans = term1 + term2;
    delete [] diff_arr;
    return(ans);
}

    
double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu,
                int ndet, int nskyx, int nskyy)
{
    int nsky = nskyx * nskyy;
    double term1 = GetFuncF(rho_arr, beta, mu, nskyx, nskyy);
    double term2 = GetFuncG(rho_arr, data_arr, resp_mat_arr, ndet, nsky);
    double ans = term1 + term2;
    return(ans);
}


double GetFuncLSupp(const double* const rho_arr,
                    const int* const rho_supp_arr,
                    const double* const data_arr,
                    const double* const resp_mat_arr,
                    double beta, double mu,
                    int ndet, int nskyx, int nskyy)
{
    int nsky = nskyx * nskyy;
    double term1 = GetFuncFSupp(rho_arr, rho_supp_arr, beta, mu, nskyx, nskyy);
    double term2 = GetFuncG(rho_arr, data_arr, resp_mat_arr, ndet, nsky);
    double ans = term1 + term2;
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
    double* diff_f_arr = new double[nsky];
    double* diff_g_arr = new double[nsky];
    GetDiffF(rho_arr, beta, mu,
             nskyx, nskyy,
             diff_f_arr);
    GetDiffG(rho_arr, data_arr, resp_mat_arr,
             ndet, nsky, diff_g_arr);
    for(int isky = 0; isky < nsky; isky++){
        out_arr[isky] = diff_f_arr[isky] + diff_g_arr[isky];
    }
    delete [] diff_f_arr;
    delete [] diff_g_arr;
}


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


void GetFuncSigma(const double* const rho_arr,
                  const double* const data_arr,
                  double beta, double mu,
                  double lconst, int nskyx, int nskyy,
                  double* out_arr)
{
//    term1 = rho.vec
//    term2 = -1.0 / L * DiffF(rho.vec, beta, mu, nrow, ncol)
//    ans.vec = term1 + term2

    int nsky = nskyx * nskyy;
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, out_arr, 1);
    double* difff_arr = new double[nsky];
    GetDiffF(rho_arr, beta, mu, nskyx, nskyy, difff_arr);
    daxpy_(nsky, -1.0 / lconst, const_cast<double*>(difff_arr), 1, out_arr, 1);
    delete [] difff_arr;
}

void GetFuncSigmaSupp(const double* const rho_arr,
                      const int* const rho_supp_arr,
                      const double* const data_arr,
                      double beta, double mu,
                      double lconst, int nskyx, int nskyy,
                      double* out_arr)
{
    int nsky = nskyx * nskyy;
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, out_arr, 1);
    double* difff_arr = new double[nsky];
    GetDiffFSupp(rho_arr, rho_supp_arr, beta, mu, nskyx, nskyy, difff_arr);
    daxpy_(nsky, -1.0 / lconst, const_cast<double*>(difff_arr), 1, out_arr, 1);
    delete [] difff_arr;
}

double GetFuncF(const double* const rho_arr,
                double beta, double mu,
                int nskyx, int nskyy)
{
    //    term1 = (1.0 - beta) * sum(log(rho.vec))
    //    term2 = mu * TermV(rho.vec, nrow, ncol)
    //    ans = term1 + term2
    
    int nsky = nskyx * nskyy;
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum += log(rho_arr[isky]);
    }
    double ans = (1.0 - beta) * sum + mu * GetTermV(rho_arr, nskyx, nskyy);
    return(ans);
}

double GetFuncFSupp(const double* const rho_arr,
                    const int* const rho_supp_arr,
                    double beta, double mu,
                    int nskyx, int nskyy)
{
    // term1 = (1.0 - beta) * SumLogPlus(rho.vec, rho.supp.vec)
    // term2 = mu * TermV(rho.vec, nrow, ncol)
    // ans = term1 + term2

    int nsky = nskyx * nskyy;
    double ans
        = (1.0 - beta) * GetSumLogPlus(rho_arr, nsky, rho_supp_arr)
        + mu * GetTermV(rho_arr, nskyx, nskyy);
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

void GetDiffFSupp(const double* const rho_arr,
                  const int* const rho_supp_arr,
                  double beta, double mu,
                  int nskyx, int nskyy,
                  double* const out_arr)
{
    // term1 = (1 - beta) / rho.vec
    int nsky = nskyx * nskyy;
    double* rho_inv_arr = new double[nsky];
    GetInvArrSupp(rho_arr, nsky, rho_supp_arr, rho_inv_arr);

    // term2 = mu * DiffTermV(rho.vec, nrow, ncol)
    GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);

    // ans = term1 + term2
    daxpy_(nsky, 1.0 - beta, rho_inv_arr, 1, out_arr, 1);
    delete [] rho_inv_arr;
}


double GetFuncG(const double* const rho_arr, 
                const double* const data_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky)
{
    //    num.vec = R.mat %*% (rho.vec)
    double* det_arr = new double[ndet];

    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);

    //ans = -1 * sum( D.vec * log( num.vec ) )
    for(int idet = 0; idet < ndet; idet ++){
        det_arr[idet] = log(det_arr[idet]);
    }
    double ans = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1, const_cast<double*>(det_arr), 1);

    delete [] det_arr;
    delete [] transa;
    return(ans);
}


void GetDiffG(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              int ndet, int nsky,
              double* const out_arr)
{
    //    num.vec = R.mat %*% (rho.vec)
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);
    
    // ans.vec = t(R.mat) %*% (D.vec / num.vec)
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / det_arr[idet];
    }
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, out_arr, 1);
    for(int isky = 0; isky < nsky; isky++){
        out_arr[isky] = -1 * out_arr[isky];
    }

    delete [] det_arr;
    delete [] transa;
    delete [] div_arr;
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
