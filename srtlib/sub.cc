#include "sub.h"
#include "sub_pm.h"
#include "sub_newton.h"

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
    for(int ipm = 0; ipm < npm; ipm++){
        double* vval_arr = new double[nsky];

        double lip_const = 10.0;
        
        GetVvalArr(rho_pre_arr,
                   nskyx, nskyy,
                   mu, lip_const,
                   vval_arr);
        double wval = GetWval(nu_pre);

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

        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre,
                                            rho_new_arr, nu_new, nsky);
        if (helldist < tol_pm){
            printf("ipm = %d, helldist = %e\n",
                   ipm, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        nu_pre = nu_new;
        lambda = lambda_new;
    }
    *nu_new_ptr = nu_new;
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
