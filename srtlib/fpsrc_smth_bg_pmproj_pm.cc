#include "fpsrc_smth_bg_pmproj_pm.h"
#include "fpsrc_smth_bg_newton.h"
#include "sub_smooth.h"
#include "fpsrc_smth_bg_statval.h"

// f_pm(rho, nu, phi) = mu V(rho) + B / phi
void GetDerivFpm(const double* const rho_arr,
                 double phi,
                 int nskyx, int nskyy, int nsrc,
                 double mu, double B_val,
                 double* const deriv_fpm_sky_arr,
                 double* const deriv_fpm_src_arr,
                 double* const deriv_fpm_phi_ptr)
{
    int nsky = nskyx * nskyy;
    double* deriv_rho_arr = new double[nsky];
    GetDiffTermV(rho_arr, nskyx, nskyy, deriv_rho_arr);
    for(int isky = 0; isky < nsky; isky ++){
        deriv_fpm_sky_arr[isky] = mu * deriv_rho_arr[isky];
    }
    delete [] deriv_rho_arr;
    for(int isrc = 0; isrc < nsrc; isrc ++){
        deriv_fpm_src_arr[isrc] = 0.0;
    }
    double deriv_fpm_phi = -1 * B_val / (phi * phi);
    *deriv_fpm_phi_ptr = deriv_fpm_phi;
}

// proj
void GetProjectedDerivFpm(const double* const deriv_fpm_sky_arr,
                          const double* const deriv_fpm_src_arr,
                          double deriv_fpm_phi,
                          int nsky, int nsrc,
                          double* const proj_deriv_fpm_sky_arr,
                          double* const proj_deriv_fpm_src_arr,
                          double* const proj_deriv_fpm_phi_ptr)
{
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky++){
        sum += deriv_fpm_sky_arr[isky];
    }
    for(int isrc = 0; isrc < nsrc; isrc++){
        sum += deriv_fpm_src_arr[isrc];
    }
    sum += deriv_fpm_phi;

    for(int isky = 0; isky < nsky; isky++){
        proj_deriv_fpm_sky_arr[isky] = deriv_fpm_sky_arr[isky] 
            - sum / (nsky + nsrc + 1);
    }
    for(int isrc = 0; isrc < nsrc; isrc++){
        proj_deriv_fpm_src_arr[isrc] = deriv_fpm_src_arr[isrc]
            - sum / (nsky + nsrc + 1);
    }
    double proj_deriv_fpm_phi = deriv_fpm_phi 
        - sum / (nsky + nsrc + 1);

    *proj_deriv_fpm_phi_ptr = proj_deriv_fpm_phi;
}

void GetVvalArr(const double* const rho_arr,
                const double* const proj_deriv_fpm_rho_arr,
                int nsky, double lip_const,
                double* const vval_arr)
{
    for(int isky = 0; isky < nsky; isky++){
        vval_arr[isky] = rho_arr[isky]
            - proj_deriv_fpm_rho_arr[isky] / lip_const;
    }
}

void GetWvalArr(const double* const nu_arr,
                const double* const proj_deriv_fpm_nu_arr,
                int nsrc, double lip_const,
                double* const wval_arr)
{
    for(int isrc = 0; isrc < nsrc; isrc++){
        wval_arr[isrc] = nu_arr[isrc]
            - proj_deriv_fpm_nu_arr[isrc] / lip_const;
    }
}

double GetZval(double phi,
               double proj_deriv_fpm_phi,
               double lip_const)
{
    double zval = phi - proj_deriv_fpm_phi / lip_const;
    return zval;
}

void GetVvalArrWvalArrZval(const double* const rho_arr,
                           const double* const nu_arr,
                           double phi,
                           int nskyx, int nskyy, int nsrc,
                           double mu, double B_val,
                           double lip_const,
                           double* const vval_arr,
                           double* const wval_arr,
                           double* const zval_ptr)
{
    int nsky = nskyx * nskyy;    
    double* deriv_fpm_sky_arr = new double[nsky];
    double* deriv_fpm_src_arr = new double[nsrc];
    double deriv_fpm_phi = 0.0;
    GetDerivFpm(rho_arr,
                phi,
                nskyx, nskyy, nsrc,
                mu, B_val,
                deriv_fpm_sky_arr,
                deriv_fpm_src_arr,
                &deriv_fpm_phi);
    double* proj_deriv_fpm_sky_arr = new double[nsky];
    double* proj_deriv_fpm_src_arr = new double[nsrc];
    double proj_deriv_fpm_phi = 0.0;
    GetProjectedDerivFpm(deriv_fpm_sky_arr,
                         deriv_fpm_src_arr,
                         deriv_fpm_phi,
                         nsky, nsrc,
                         proj_deriv_fpm_sky_arr,
                         proj_deriv_fpm_src_arr,
                         &proj_deriv_fpm_phi);
    GetVvalArr(rho_arr,
               proj_deriv_fpm_sky_arr,
               nsky, lip_const,
               vval_arr);
    GetWvalArr(nu_arr,
               proj_deriv_fpm_src_arr,
               nsrc, lip_const,
               wval_arr);
    double zval = GetZval(phi,
                          proj_deriv_fpm_phi,
                          lip_const);
    delete [] deriv_fpm_sky_arr;
    delete [] deriv_fpm_src_arr;
    delete [] proj_deriv_fpm_sky_arr;
    delete [] proj_deriv_fpm_src_arr;
    *zval_ptr = zval;
}

double GetFindLipConst(const double* const rho_arr,
                       const double* const nu_arr,
                       double phi,
                       const double* const mval_arr,
                       const double* const nval_arr,
                       double pval,
                       double phi_val,
                       int nph, double B_val,
                       double mu,
                       int nskyx, int nskyy, int nsrc,
                       double lip_const, double lambda,
                       int nnewton, double tol_newton)
{
    int nsky = nskyx * nskyy;
    int ik_max = 100;
    double eta = 1.2;

    double* rho_new_arr = new double[nsky];
    double* nu_new_arr = new double[nsrc];
    double phi_new = 0.0;
    for(int ik = 0; ik < ik_max; ik ++){
        lip_const = lip_const * pow(eta, ik);

        double* vval_arr = new double[nsky];
        double* wval_arr = new double[nsrc];
        double zval = 0.0;
        GetVvalArrWvalArrZval(rho_arr, nu_arr, phi,
                              nskyx, nskyy, nsrc,
                              mu, B_val,
                              lip_const,
                              vval_arr,
                              wval_arr,
                              &zval);
        double lambda_new = 0.0;
        GetRhoArrNuArrPhi_ByNewton(vval_arr, wval_arr, zval,
                                   mval_arr, nval_arr, pval,
                                   phi_val,
                                   nsky, nsrc, nph,
                                   lip_const,
                                   nnewton, tol_newton,
                                   lambda,
                                   rho_new_arr,
                                   nu_new_arr,
                                   &phi_new,
                                   &lambda_new);
        double qminusf = GetQMinusF(rho_new_arr, nu_new_arr, phi_new,
                                    rho_arr, nu_arr, phi,
                                    mu, lip_const, B_val,
                                    nskyx, nskyy, nsrc);
        delete [] vval_arr;
        delete [] wval_arr;
        // debug
        if(qminusf >= 0.0 && phi_new > 0.0){
        // if(qminusf >= 0.0){
            break;
        }
    }
    delete [] rho_new_arr;
    delete [] nu_new_arr;
    return lip_const;
}


double GetQMinusF(const double* const rho_new_arr,
                  const double* const nu_new_arr,
                  double phi_new,
                  const double* const rho_arr,
                  const double* const nu_arr,
                  double phi,
                  double mu, double lip_const,
                  double B_val,
                  int nskyx, int nskyy, int nsrc)
{
    double term1 = GetFuncF(rho_arr, phi, mu, B_val, nskyx, nskyy);
    double term2 = -1 * GetFuncF(rho_new_arr, phi_new, mu,
                                 B_val, nskyx, nskyy);
    int nsky = nskyx * nskyy;
    double* diff_rho_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, diff_rho_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(rho_arr), 1, diff_rho_arr, 1);

    double* diff_nu_arr = new double[nsrc];
    dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, diff_nu_arr, 1);
    daxpy_(nsrc, -1.0, const_cast<double*>(nu_arr), 1, diff_nu_arr, 1);

    double* deriv_fpm_sky_arr = new double[nsky];
    double* deriv_fpm_src_arr = new double[nsrc];
    double deriv_fpm_phi = 0.0;
    GetDerivFpm(rho_arr, phi,
                nskyx, nskyy, nsrc,
                mu, B_val,
                deriv_fpm_sky_arr,
                deriv_fpm_src_arr,
                &deriv_fpm_phi);
    double* proj_deriv_fpm_sky_arr = new double[nsky];
    double* proj_deriv_fpm_src_arr = new double[nsrc];
    double proj_deriv_fpm_phi = 0.0;
    GetProjectedDerivFpm(deriv_fpm_sky_arr,
                         deriv_fpm_src_arr,
                         deriv_fpm_phi,
                         nsky, nsrc,
                         proj_deriv_fpm_sky_arr,
                         proj_deriv_fpm_src_arr,
                         &proj_deriv_fpm_phi);
    double term3 = ddot_(nsky, const_cast<double*>(diff_rho_arr), 1,
                         const_cast<double*>(proj_deriv_fpm_sky_arr), 1);
    term3 += ddot_(nsrc, const_cast<double*>(diff_nu_arr), 1,
                   const_cast<double*>(proj_deriv_fpm_src_arr), 1);
    term3 += (phi_new - phi) * proj_deriv_fpm_phi;

    delete [] deriv_fpm_sky_arr;
    delete [] deriv_fpm_src_arr;
    delete [] proj_deriv_fpm_sky_arr;
    delete [] proj_deriv_fpm_src_arr;
    
    double term4 = lip_const / 2.0 *
        (ddot_(nsky, const_cast<double*>(diff_rho_arr), 1,
               const_cast<double*>(diff_rho_arr), 1)
         + ddot_(nsrc, const_cast<double*>(diff_nu_arr), 1,
               const_cast<double*>(diff_nu_arr), 1)
         + pow(phi_new - phi, 2));
    double ans = term1 + term2 + term3 + term4;
    delete [] diff_rho_arr;
    delete [] diff_nu_arr;    
    return ans;
}


double GetFuncF(const double* const rho_arr, double phi,
                double mu, double B_val, 
                int nskyx, int nskyy)
{
    double ans = mu * GetTermV(rho_arr, nskyx, nskyy) + B_val / phi;
    return(ans);
}

void GetRhoNuPhi_ByPM(const double* const rho_arr,
                      const double* const nu_arr,
                      double phi,
                      const double* const mval_arr,
                      const double* const nval_arr,
                      double pval,
                      int nph, double B_val,
                      int ndet, int nskyx, int nskyy, int nsrc,
                      double mu,
                      int npm, double tol_pm,
                      int nnewton, double tol_newton,
                      double* const rho_new_arr,
                      double* const nu_new_arr,
                      double* const phi_new_ptr)
{
    double phi_val = phi;
    
    int nsky = nskyx * nskyy;
    double* rho_pre_arr = new double[nsky];
    double* nu_pre_arr = new double[nsrc];    
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_arr), 1, nu_pre_arr, 1);
    double phi_pre = phi;
    double phi_new = 0.0;

    double lip_const = 1.0;
    double lambda = 0.0;
    double lambda_new = 0.0;
    double lip_const_new = 1.0;
    for(int ipm = 0; ipm < npm; ipm++){
        lip_const_new = GetFindLipConst(rho_pre_arr, nu_pre_arr, phi_pre,
                                        mval_arr, nval_arr, pval,
                                        phi_val, nph, B_val, mu,
                                        nskyx, nskyy, nsrc, lip_const,
                                        lambda, nnewton, tol_newton);
        // printf("lip_const_new = %e\n", lip_const_new);
        double* vval_arr = new double[nsky];
        double* wval_arr = new double[nsrc];
        double zval = 0.0;
        GetVvalArrWvalArrZval(rho_arr, nu_arr, phi,
                              nskyx, nskyy, nsrc,
                              mu, B_val,
                              lip_const,
                              vval_arr,
                              wval_arr,
                              &zval);
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
        delete [] vval_arr;
        delete [] wval_arr;
        // printf("pm out: ipm = %d, phi_new = %e\n", ipm, phi_new);
        
        double helldist  = GetHellingerDist(rho_pre_arr,
                                            nu_pre_arr,
                                            phi_pre,
                                            rho_new_arr,
                                            nu_new_arr,
                                            phi_new,
                                            nsky, nsrc);
        // printf("ipm = %d, helldist = %e\n", ipm, helldist);
        if (helldist < tol_pm){
            printf("    ipm = %d, helldist = %e, lip_const_new = %e\n",
                   ipm, helldist, lip_const_new);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
        phi_pre = phi_new;
        lambda = lambda_new;
        lip_const = lip_const_new;
    }

    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
    *phi_new_ptr = phi_new;
}
