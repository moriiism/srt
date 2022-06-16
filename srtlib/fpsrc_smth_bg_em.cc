

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


void GetRhoNuPhi_ByDC(rho_pre_arr, nu_pre_arr, phi_pre,
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
                      phi_new)
{
    for(int idc = 0; idc < ndc; idc++){
        GetRhoNuPhi_ByPM(rho_arr, nu_arr, phi,
                         data_arr,
                         resp_norm_mat_arr,
                         bg_arr,
                         ndet, nskyx, nskyy,
                         mu,
                         npm, tol_pm,
                         nnewton, tol_newton,
                         rho_new_arr,
                         nu_new_arr,
                         &phi_new);
        
        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre_arr, phi_pre,
                                            rho_new_arr, nu_new_arr, phi_new,
                                            nsky, nsrc);
        if (helldist < tol_dc){
            printf("idc = %d, helldist = %e\n",
                   idc, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
        nu_pre = nu_new;
    }

}



void GetRhoNuPhi_ByPM(const double* const rho_arr, double nu,
                     const double* const data_arr,
                     const double* const resp_norm_mat_arr,
                     const double* const bg_arr,
                     int ndet, int nskyx, int nskyy,
                     double mu,
                     int npm, double tol_pm,
                     int nnewton, double tol_newton,
                     double* const rho_new_arr,
                     double* const nu_new_ptr)
{
    int nsky = nskyx * nskyy;
    double* mval_arr = new double[nsky];
    double nval = 0.0;
    GetMArrNval(rho_arr, nu, data_arr, resp_norm_mat_arr, bg_arr,
                ndet, nsky, mval_arr, &nval);

    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    double nu_pre = nu;
    double nu_new = 0.0;
    double lambda = 0.0;

    double lip_const = 1.0;
    for(int ipm = 0; ipm < npm; ipm++){
        double lip_const_new = GetFindLipConst(rho_pre_arr, nu_pre,
                                               mval_arr, nval, mu,
                                               nskyx, nskyy, lip_const,
                                               lambda, nnewton, tol_newton);
        // printf("lip_const_new = %e\n", lip_const_new);
        double* vval_arr = new double[nsky];        
        GetVvalArr(rho_pre_arr,
                   nskyx, nskyy,
                   mu, lip_const_new,
                   vval_arr);
        double wval = GetWval(nu_pre);

        double lambda_new = 0.0;
        GetRhoArrNu_ByNewton(vval_arr, wval,
                             mval_arr, nval,
                             nsky,
                             lip_const_new,
                             nnewton, tol_newton,
                             lambda,
                             rho_new_arr,
                             &nu_new,
                             &lambda_new);

        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre,
                                            rho_new_arr, nu_new, nsky);
        delete [] vval_arr;
        if (helldist < tol_pm){
            printf("ipm = %d, helldist = %e\n",
                   ipm, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        nu_pre = nu_new;
        lambda = lambda_new;
        lip_const = lip_const_new;
    }
    *nu_new_ptr = nu_new;
    delete [] mval_arr;
    delete [] rho_pre_arr;
}
