
void GetVvalArr(const double* const rho_arr,
                int nskyx, int nskyy,
                double mu, double lip_const,
                double* const vval_arr)
{
    int nsky = nskyx * nskyy;
    double* deriv_rho_arr = new double[nsky];
    GetDiffTermV(rho_arr, nskyx, nskyy, deriv_rho_arr);
    for(int isky = 0; isky < nsky; isky++){
        vval_arr[isky] = rho_arr[isky] - mu / lip_const * deriv_rho_arr[isky];
    }
    delete [] deriv_rho_arr;
}

double GetWvalArr(const double* const nu_arr,
                  int nsrc,
                  double* const wval_arr)
{
    for(int isrc = 0; isrc < nsrc; isrc++){
        wval_arr[isrc] = nu_arr[isrc];
    }
}

double GetZval(double phi,
               double lip_const,
               double B_val)

{
    double zval = phi + B_val / (lip_const * phi * phi);
    return zval;
}


double GetFindLipConst(const double* const rho_arr, double nu,
                       const double* const mval_arr, double nval,
                       double mu,
                       int nskyx, int nskyy,
                       double lip_const, double lambda,
                       int nnewton, double tol_newton)
{
    int nsky = nskyx * nskyy;
    int ik_max = 100;
    double eta = 1.2;

    double* rho_new_arr = new double[nsky];
    double nu_new = 0.0;
    for(int ik = 0; ik < ik_max; ik ++){
        lip_const = lip_const * pow(eta, ik);

        double* vval_arr = new double[nsky];
        GetVvalArr(rho_arr,
                   nskyx, nskyy,
                   mu, lip_const,
                   vval_arr);
        double wval = GetWval(nu);
        double lambda_new = 0.0;
        GetRhoArrNu_ByNewton(vval_arr, wval,
                             mval_arr, nval,
                             nsky,
                             lip_const,
                             nnewton, tol_newton,
                             lambda,
                             rho_new_arr,
                             &nu_new,
                             &lambda_new);
        double qminusf = GetQMinusF(rho_new_arr, nu_new,
                                    rho_arr, nu, 
                                    mu, lip_const,
                                    nskyx, nskyy);
        delete [] vval_arr;
        if(qminusf >= 0.0){
            break;
        }
    }
    delete [] rho_new_arr;
    return lip_const;
}


double GetQMinusF(const double* const rho_new_arr, double nu_new,
                  const double* const rho_arr, double nu,
                  double mu, double lip_const,
                  int nskyx, int nskyy)
{
    double term1 = GetFuncF(rho_arr, mu, nskyx, nskyy);
    double term2 = -1 * GetFuncF(rho_new_arr, mu, nskyx, nskyy);
    int nsky = nskyx * nskyy;
    double* diff_rho_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, diff_rho_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(rho_arr), 1, diff_rho_arr, 1);
    double* diff_f_arr = new double[nsky];
    GetDiffF(rho_arr, mu, nskyx, nskyy, diff_f_arr);
    double term3 = ddot_(nsky, const_cast<double*>(diff_rho_arr), 1,
                         const_cast<double*>(diff_f_arr), 1);
    double term4 = lip_const / 2.0 *
        (ddot_(nsky, const_cast<double*>(diff_rho_arr), 1,
               const_cast<double*>(diff_rho_arr), 1) + pow(nu_new - nu, 2));
    double ans = term1 + term2 + term3 + term4;
    delete [] diff_rho_arr;
    delete [] diff_f_arr;
    return(ans);
}


double GetFuncF(const double* const rho_arr,
                double mu,
                int nskyx, int nskyy)
{
    double ans = mu * GetTermV(rho_arr, nskyx, nskyy);
    return(ans);
}

void GetDiffF(const double* const rho_arr,
              double mu,
              int nskyx, int nskyy,
              double* const out_arr)
{
    int nsky = nskyx * nskyy;
    GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);
}

void GetRhoNuPhi_ByPM(const double* const rho_arr,
                      const double* const nu_arr,
                      double phi,
                      const double* const mval_arr,
                      const double* const nval_arr,
                      double pval,
                      double phi_val,
                      int nph, double B_val,
                      int ndet, int nskyx, int nskyy, int nsrc,
                      double mu,
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
        double* wval_arr = new double[nsrc];
        GetWvalArr(nu_pre_arr,
                   nsrc,
                   wval_arr);
        double zval = GetZval(phi_pre,
                              lip_const_new,
                              B_val);
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
        delete [] wval_arr;
        if (helldist < tol_pm){
            printf("ipm = %d, helldist = %e\n",
                   ipm, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
        phi_pre = phi_new;
        lambda = lambda_new;
        lip_const = lip_const_new;
    }
    *nu_new_ptr = nu_new;
    delete [] mval_arr;
    delete [] rho_pre_arr;
}
