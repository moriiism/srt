
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
        GetRhoArrNuArrPhi_ByNewton(vval_arr, wval_arr, zval,
                                   mval_arr, nval_arr, pval,
                                   phi_val,
                                   nsky, nsrc, nph,
                                   lip_const_new,
                                   nnewton, tol_newton,
                                   lambda,
                                   rho_new_arr,
                                   nu_new_arr,
                                   &phi_new,
                                   &lambda_new);

        double helldist  = GetHellingerDist(rho_pre_arr,
                                            nu_pre_arr,
                                            phi_pre,
                                            rho_new_arr,
                                            nu_new_arr,
                                            phi_new,
                                            nsky, nsrc);
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
