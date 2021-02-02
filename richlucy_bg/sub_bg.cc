#include "sub_bg.h"

void GetNextNb(const double* const rho_arr,
               const double* const data_arr,
               const double* const resp_mat_arr,
               const duoble* const bg_arr,
               int ndet, int nsky,
               double N_B)
{
    double N_B_new = N_B;
    
    double B = GetB(bg_arr, ndet);
    double deriv_f_at_B = GetDerivF_NB(rho_arr, data_arr, resp_mat_arr,
                                       bg_arr, ndet, nsky, B);
    if(deriv_f_at_B >= 0.0){
        N_B_new = B;
    } else{
        // get next N_B by Newton Method
        for(int iter = 0; iter < 10; iter ++){
            double deriv_f_at_NB = GetDerivF_NB(rho_arr, data_arr, resp_mat_arr, bg_arr,
                                                ndet, nsky, N_B_new);
            double deriv2_f_at_NB = GetDeriv2F_NB(rho_arr, data_arr, resp_mat_arr, bg_arr,
                                                  ndet, nsky, N_B_new);
            N_B_new -= deriv_f_at_NB / deriv2_f_at_NB;
        }
    }
    return(N_B_new);
}

double GetB(const duoble* const bg_arr, int ndet)
{
    double B = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        B += bg_arr[idet];
    }
    return(B);
}

double GetDerivF_NB(const double* const rho_arr,
                    const double* const data_arr,
                    const double* const resp_mat_arr,
                    const duoble* const bg_arr,
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
        sum += data_arr[idet] * det_arr[idet] / (N_B * det_arr[idet] + bg_arr[idet]);
    }
    deriv_f = 1.0 - sum;
    return(deriv_f)
}

double GetDeriv2F_NB(const double* const rho_arr,
                     const double* const data_arr,
                     const double* const resp_mat_arr,
                     const duoble* const bg_arr,
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
    deriv2_f = sum;
    return(deriv2_f)
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
                   const duoble* const bg_arr,
                   int ndet, int nsky,
                   double N_B,
                   double* const out_arr)
{
    // coeff
    double num = 1.0 - B / N_B;
    double* det_arr = new double[ndet];
    GetDetArr(rho_arr, resp_mat_arr, ndet, nsky, det_arr);
    double den = 0.0;
    for(int idet = 0; idet < ndet; idet++){
        den += data_arr[idet] * det_arr[idet] / (det_arr[idet] + bg_arr[idet] / N_B);
    }
    double coeff = num / den;

    // mval_arr = t(R_mat) %*% (data_arr / det_arr) * rho_arr
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / (det_arr[idet] + bg_arr[idet] / N_B);
    }
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


        GetNextNb(rho_arr, data_arr, resp_mat_arr,
               const duoble* const bg_arr,
               int ndet, int nsky,
               double N_B)
        


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
