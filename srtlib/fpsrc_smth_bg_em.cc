#include "fpsrc_smth_bg_em.h"
#include "fpsrc_smth_bg_dc.h"
#include "fpsrc_smth_bg_statval.h"

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

// get mval_arr, nval_arr, pval
void GetMvalArrNvalArrPval(const double* const rho_arr,
                           const double* const nu_arr,
                           double phi,
                           const double* const data_arr,
                           const double* const bg_arr,
                           const double* const* const det_fpsrc_arr,
                           const double* const resp_norm_mat_arr,
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
    for(int isrc = 0; isrc < nsrc; isrc++){
        daxpy_(ndet, nu_arr[isrc],
               const_cast<double*>(det_fpsrc_arr[isrc]), 1,
               den_arr, 1);
    }
    double B_val = MibBlas::Sum(bg_arr, ndet);    
    daxpy_(ndet, phi/B_val, const_cast<double*>(bg_arr), 1, den_arr, 1);
    
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
    for(int isrc = 0; isrc < nsrc; isrc++){
        nval_arr[isrc] = ddot_(ndet, div_arr, 1,
                               const_cast<double*>(det_fpsrc_arr[isrc]), 1)
            * nu_arr[isrc];
    }
    double pval = ddot_(ndet, div_arr, 1, const_cast<double*>(bg_arr), 1)
        * phi / B_val;

    delete [] den_arr;
    delete [] div_arr;
    delete [] tmp_arr;
    *pval_ptr = pval;
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
    double B_val = MibBlas::Sum(bg_arr, ndet);
    int nph = MibBlas::Sum(data_arr, ndet);
    int nsky = nskyx * nskyy;
    double* rho_pre_arr = new double[nsky];
    double* nu_pre_arr = new double[nsrc];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_init_arr), 1, nu_pre_arr, 1);
    double phi_pre = phi_init;
    double phi_new = 0.0;
    for(int iem = 0; iem < nem; iem ++){
        printf("iem = %d\n", iem);

        double* mval_arr = new double[nsky];
        double* nval_arr = new double[nsrc];
        double pval = 0.0;
        GetMvalArrNvalArrPval(rho_pre_arr, nu_pre_arr, phi_pre,
                              data_arr, bg_arr, det_fpsrc_arr, resp_norm_mat_arr, 
                              ndet, nsky, nsrc,
                              mval_arr, nval_arr, &pval);

        //        for(int isky = 0; isky < nsky; isky ++){
        //    //if(rho_pre_arr[isky] <= 0.0){
        //    printf("rho_pre_arr[isky] = %e\n", rho_pre_arr[isky]);
        //    //}
        //}
        //for(int isrc = 0; isrc < nsrc; isrc ++){
        //    //if(nu_pre_arr[isrc] <= 0.0){
        //    printf("nu_pre_arr[isrc] = %e\n", nu_pre_arr[isrc]);
        //    //}
        //}
        GetRhoNuPhi_ByDC(rho_pre_arr, nu_pre_arr, phi_pre,
                         mval_arr, nval_arr, pval,
                         nph, B_val,
                         ndet, nskyx, nskyy, nsrc,
                         mu,
                         ndc, tol_dc,
                         npm, tol_pm,
                         nnewton, tol_newton,
                         rho_new_arr,
                         nu_new_arr,
                         &phi_new);
        delete [] mval_arr;
        delete [] nval_arr;
        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre_arr, phi_pre,
                                            rho_new_arr, nu_new_arr, phi_new,
                                            nsky, nsrc);
        if (helldist < tol_em){
            printf("iem = %d, helldist = %e\n",
                   iem, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nsrc, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
        phi_pre = phi_new;

        //double lval = 0.0;        
        if (iem % 100 == 0){
            //lval = GetFuncL(data_arr, bg_arr,
            //rho_new_arr, nu_new,
            //                resp_norm_mat_arr,
            //                ndet, nsky);
            //printf("iem = %d, helldist = %e, lval = %e\n",
            //       iem, helldist, lval);
        } else {
            printf("iem = %d, helldist = %e\n",
                   iem, helldist);
        }
    }
    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
    *phi_new_ptr = phi_new;
}
