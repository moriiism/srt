#include "sub.h"
#include "sub_pm.h"
#include "sub_newton.h"
#include "sub_smooth.h"

void GetDetArr(const double* const rho_arr,
               const double* const resp_norm_mat_arr,
               int ndet, int nsky,
               double* const det_arr) // ndet
{
    // det_arr = R_mat %*% rho_arr
    char transa[1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);
}

// get m_arr & n_val
void GetMArrNval(const double* const rho_arr, double nu,
                 const double* const data_arr,
                 const double* const resp_norm_mat_arr,
                 const double* const bg_arr,
                 int ndet, int nsky,
                 double* const mval_arr,
                 double* const nval_ptr)
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


void GetRhoNu_New(const double* const rho_arr, double nu,
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


void RichlucyBg2Smooth(const double* const rho_init_arr,
                       double nu_init,
                       const double* const data_arr,
                       const double* const bg_arr,
                       const double* const resp_norm_mat_arr,
                       int ndet, int nskyx, int nskyy, double mu,
                       string outdir,
                       string outfile_head,
                       int nem, double tol_em,
                       int npm, double tol_pm,
                       int nnewton, double tol_newton,
                       double* const rho_new_arr,
                       double* const nu_new_ptr)
{
    int nsky = nskyx * nskyy;
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    double nu_pre = nu_init;
    double nu_new = nu_init;
    for(int iem = 0; iem < nem; iem ++){
        GetRhoNu_New(rho_pre_arr, nu_pre,
                     data_arr,
                     resp_norm_mat_arr,
                     bg_arr,
                     ndet, nskyx, nskyy,
                     mu,
                     npm, tol_pm,
                     nnewton, tol_newton,
                     rho_new_arr,
                     &nu_new);
        
        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre,
                                            rho_new_arr, nu_new, nsky);
        if (helldist < tol_em){
            printf("iem = %d, helldist = %e\n",
                   iem, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
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
    }
    delete [] rho_pre_arr;
    *nu_new_ptr = nu_new;
}

double GetHellingerDist(const double* const rho_arr, double nu, 
                        const double* const rho_new_arr, double nu_new,
                        int nsky)
{
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        double diff = sqrt(rho_arr[isky]) - sqrt(rho_new_arr[isky]);
        sum += diff * diff;
    }
    double diff = sqrt(nu) - sqrt(nu_new);
    sum += diff * diff;
    double ans = sqrt(sum);
    return (ans);
}


double GetFuncL(const double* const data_arr,
                const double* const bg_arr,
                const double* const rho_arr,
                double nu,
                const double* const resp_norm_mat_arr,
                int ndet, int nsky)
{
    // ans = - sum_v [ Y(v) log( sum_u t(v,u) rho_u + b(v) * nu ) ]
    double* tmp_arr = new double[ndet];
    GetDetArr(rho_arr, resp_norm_mat_arr,
              ndet, nsky, tmp_arr);
    daxpy_(ndet, nu, const_cast<double*>(bg_arr), 1, tmp_arr, 1);
    for(int idet = 0; idet < ndet; idet ++){
        tmp_arr[idet] = log(tmp_arr[idet]);
    }
    double ans = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1.0,
                              tmp_arr, 1.0);
    delete [] tmp_arr;
    return(ans);
}
