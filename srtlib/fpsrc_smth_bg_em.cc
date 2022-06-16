
// get mval_arr, nval_arr, pval
void GetMvalArrNvalArrPval(const double* const rho_arr,
                           const double* const nu_arr,
                           double phi,
                           const double* const data_arr,
                           const double* const resp_norm_mat_arr,
                           const double* const* const bval_arr,
                           const double* const bg_arr,
                           int ndet, int nsky, int nsrc,
                           double* const mval_arr,
                           double* const nval_arr,
                           double* const pval_ptr)
{
    double* den_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet ++){
        den_arr[idet] = 0.0;
    }
    GetDetArr(rho_arr, resp_norm_mat_arr, ndet, nsky, den_arr);
    daxpy_(ndet, nu, const_cast<double*>(bg_arr), 1, den_arr, 1);
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / den_arr[idet];
    }
    double* tmp_arr = new double[nsky];
    char transa[1];
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_arr), ndet,
           div_arr, 1,
           0.0, tmp_arr, 1);
    MibBlas::ElmWiseMul(nsky, 1.0,
                        tmp_arr, rho_arr, mval_arr);
    double nval = ddot_(ndet, div_arr, 1, const_cast<double*>(bg_arr), 1) * nu;

    delete [] den_arr;
    delete [] div_arr;
    delete [] tmp_arr;
    *nval_ptr = nval;
}


void RichlucyFpsrcSmthBg(const double* const rho_init_arr,
                         const double* const nu_init_arr,
                         double phi_init,
                         const double* const data_arr,
                         const double* const bg_arr,
                         const double* const* const det_fpsrc_arr,
                         const double* const resp_norm_mat_arr,
                         int ndet, int nskyx, int nskyy, int nsrc,
                         double mu,
                         string outdir,
                         string outfile_head,
                         int nem, double tol_em,
                         int ndc, double tol_dc,
                         int npm, double tol_pm,
                         int nnewton, double tol_newton,
                         double* const rho_new_arr,
                         double* const nu_new_arr,
                         double* const phi_new_ptr)
{
    int nsky = nskyx * nskyy;
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    double nu_pre = nu_init;
    double nu_new = nu_init;
    for(int iem = 0; iem < nem; iem ++){
        double* mval_arr = new double[nsky];
        double* nval_arr = new double[nsrc];
        GetMArrNval(rho_arr, nu, data_arr, resp_norm_mat_arr, bg_arr,
                    ndet, nsky, mval_arr, &nval);

        double pval = GetPval();
        
        GetRhoNuPhi_ByDC(rho_pre_arr, nu_pre_arr, phi_pre,
                        data_arr,
                        resp_norm_mat_arr,
                        bg_arr,
                        ndet, nskyx, nskyy,
                        mu,
                        ndc, tol_dc,
                        npm, tol_pm,
                        nnewton, tol_newton,
                        rho_new_arr,
                        nu_new_arr,
                        &phi_new);
        
        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre,
                                            rho_new_arr, nu_new, nsky);
        if (helldist < tol_em){
            printf("iem = %d, helldist = %e\n",
                   iem, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
        nu_pre = nu_new;

        double lval = 0.0;        
        if (iem % 100 == 0){
            lval = GetFuncL(data_arr, bg_arr,
                            rho_new_arr, nu_new,
                            resp_norm_mat_arr,
                            ndet, nsky);
            printf("iem = %d, helldist = %e, lval = %e\n",
                   iem, helldist, lval);
        } else {
            printf("iem = %d, helldist = %e\n",
                   iem, helldist);
        }
o    }
    delete [] rho_pre_arr;
    *nu_new_ptr = nu_new;
}

