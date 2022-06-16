
void GetRhoNuPhi_ByDC(const double* const rho_arr,
                      const double* const nu_arr,
                      double phi,
                      const double* const mval_arr,
                      const double* const nval_arr,
                      double pval,
                      int nph, double B_val,
                      int ndet, int nskyx, int nskyy, int nsrc,
                      double mu,
                      int ndc, double tol_dc,
                      int npm, double tol_pm,
                      int nnewton, double tol_newton,
                      double* const rho_new_arr,
                      double* const nu_new_arr,
                      double* const phi_new_ptr)
{
    int nsky = nskyx * nskyy;    
    double* rho_pre_arr = new double[nsky];
    double* nu_pre_arr = new double[nsrc];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_arr), 1, nu_pre_arr, 1);
    double phi_pre = phi;
    double phi_new = 0.0;
    for(int idc = 0; idc < ndc; idc++){
        GetRhoNuPhi_ByPM(rho_pre_arr, nu_pre_arr, phi_pre,
                         mval_arr, nval_arr, pval,
                         phi_pre,
                         nph, B_val,
                         ndet, nskyx, nskyy, nsrc,
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
        phi_pre = phi_new;
    }
    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
    *phi_new_ptr = phi_new;
}
