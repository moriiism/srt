#include "fpsrc_smth_bg_em.h"
#include "fpsrc_smth_bg_dc.h"
#include "fpsrc_smth_bg_statval.h"
#include<unistd.h>

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

void RichlucyFpsrcSmthBg(FILE* const fp_log,
                         const double* const rho_init_arr,
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
        double* mval_arr = new double[nsky];
        double* nval_arr = new double[nsrc];
        double pval = 0.0;
        GetMvalArrNvalArrPval(rho_pre_arr, nu_pre_arr, phi_pre,
                              data_arr, bg_arr, det_fpsrc_arr,
                              resp_norm_mat_arr, 
                              ndet, nsky, nsrc,
                              mval_arr, nval_arr, &pval);
        double helldist_dc = 0.0;
        int flag_converge_dc = 0;
        GetRhoNuPhi_ByDC(fp_log,
                         rho_pre_arr, nu_pre_arr, phi_pre,
                         mval_arr, nval_arr, pval,
                         nph, B_val,
                         ndet, nskyx, nskyy, nsrc,
                         mu,
                         ndc, tol_dc,
                         npm, tol_pm,
                         nnewton, tol_newton,
                         rho_new_arr,
                         nu_new_arr,
                         &phi_new,
                         &helldist_dc,
                         &flag_converge_dc);
        delete [] mval_arr;
        delete [] nval_arr;
        double helldist  = GetHellingerDist(rho_pre_arr, nu_pre_arr, phi_pre,
                                            rho_new_arr, nu_new_arr, phi_new,
                                            nsky, nsrc);
        if (access( "/tmp/fpsrc_smth_bg_em_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/fpsrc_smth_bg_em_stop file is found, then stop.\n");
            break;
        }
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
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
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
        }
    }
    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
    *phi_new_ptr = phi_new;
}

// accerelation by SQUAREM
void RichlucyFpsrcSmthBg_Acc(FILE* const fp_log,
                             const double* const rho_init_arr,
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

    double* rho_0_arr  = new double[nsky];
    double* rho_1_arr  = new double[nsky];
    double* rho_2_arr  = new double[nsky];
    double* rho_dash_arr  = new double[nsky];
    double* r_rho_arr  = new double[nsky];
    double* r2_rho_arr  = new double[nsky];
    double* v_rho_arr  = new double[nsky];
    double* rho_0_new_arr  = new double[nsky];

    double* nu_0_arr  = new double[nsrc];
    double* nu_1_arr  = new double[nsrc];
    double* nu_2_arr  = new double[nsrc];
    double* nu_dash_arr  = new double[nsrc];
    double* r_nu_arr  = new double[nsrc];
    double* r2_nu_arr  = new double[nsrc];
    double* v_nu_arr  = new double[nsrc];
    double* nu_0_new_arr  = new double[nsrc];

    double phi_0 = 0.0;
    double phi_1 = 0.0;
    double phi_2 = 0.0;
    double phi_dash = 0.0;
    double r_phi = 0.0;
    double r2_phi = 0.0;
    double v_phi = 0.0;
    double phi_0_new = 0.0;

    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_0_arr, 1);
    dcopy_(nsrc, const_cast<double*>(nu_init_arr), 1, nu_0_arr, 1);
    phi_0 = phi_init;

    for(int iem = 0; iem < nem; iem ++){
        double* mval_arr = new double[nsky];
        double* nval_arr = new double[nsrc];
        double pval = 0.0;

        double helldist_dc1 = 0.0;
        int flag_converge_dc1 = 0;
        GetMvalArrNvalArrPval(rho_0_arr, nu_0_arr, phi_0,
                              data_arr, bg_arr, det_fpsrc_arr,
                              resp_norm_mat_arr,
                              ndet, nsky, nsrc,
                              mval_arr, nval_arr, &pval);
        GetRhoNuPhi_ByDC(fp_log,
                         rho_0_arr, nu_0_arr, phi_0,
                         mval_arr, nval_arr, pval,
                         nph, B_val,
                         ndet, nskyx, nskyy, nsrc,
                         mu,
                         ndc, tol_dc,
                         npm, tol_pm,
                         nnewton, tol_newton,
                         rho_1_arr,
                         nu_1_arr,
                         &phi_1,
                         &helldist_dc1,
                         &flag_converge_dc1);
        if (flag_converge_dc1 == 0){
            MiIolib::Printf2(fp_log,
                             "iem = %d: dc1: not converged: helldist_dc1 = %.2e\n",
                             iem,
                             helldist_dc1);
        }
        double helldist_dc2 = 0.0;
        int flag_converge_dc2 = 0;
        GetMvalArrNvalArrPval(rho_1_arr, nu_1_arr, phi_1,
                              data_arr, bg_arr, det_fpsrc_arr,
                              resp_norm_mat_arr,
                              ndet, nsky, nsrc,
                              mval_arr, nval_arr, &pval);
        GetRhoNuPhi_ByDC(fp_log,
                         rho_1_arr, nu_1_arr, phi_1,
                         mval_arr, nval_arr, pval,
                         nph, B_val,
                         ndet, nskyx, nskyy, nsrc,
                         mu,
                         ndc, tol_dc,
                         npm, tol_pm,
                         nnewton, tol_newton,
                         rho_2_arr,
                         nu_2_arr,
                         &phi_2,
                         &helldist_dc2,
                         &flag_converge_dc2);
        if (flag_converge_dc2 == 0){
            MiIolib::Printf2(fp_log,
                             "iem = %d: dc2: not converged: helldist_dc2 = %.2e\n",
                             iem,
                             helldist_dc2);
        }

        MibBlas::Sub(rho_1_arr, rho_0_arr, nsky, r_rho_arr);
        MibBlas::Sub(rho_2_arr, rho_1_arr, nsky, r2_rho_arr);
        MibBlas::Sub(r2_rho_arr, r_rho_arr, nsky, v_rho_arr);

        MibBlas::Sub(nu_1_arr, nu_0_arr, nsrc, r_nu_arr);
        MibBlas::Sub(nu_2_arr, nu_1_arr, nsrc, r2_nu_arr);
        MibBlas::Sub(r2_nu_arr, r_nu_arr, nsrc, v_nu_arr);

        r_phi = phi_1 - phi_0;
        r2_phi = phi_2 - phi_1;
        v_phi = r2_phi - r_phi;

        double r_norm2 = ddot_(nsky, r_rho_arr, 1, r_rho_arr, 1)
            + ddot_(nsrc, r_nu_arr, 1, r_nu_arr, 1)
            + r_phi * r_phi;
        double v_norm2 = ddot_(nsky, v_rho_arr, 1, v_rho_arr, 1)
            + ddot_(nsrc, v_nu_arr, 1, v_nu_arr, 1)
            + v_phi * v_phi;
        double alpha = -1.0 * sqrt(r_norm2 / v_norm2);


        int nk = 100;
        double eta = 0.8;
        int ifind_nonneg = 0;
        for (int ik = 0; ik < nk; ik ++){
            double alpha0 = alpha * pow(eta, ik);
            dcopy_(nsky, rho_0_arr, 1, rho_dash_arr, 1);
            dcopy_(nsrc, nu_0_arr, 1, nu_dash_arr, 1);
            phi_dash = phi_0;
            daxpy_(nsky, -2 * alpha0, r_rho_arr, 1, rho_dash_arr, 1);
            daxpy_(nsrc, -2 * alpha0, r_nu_arr, 1, nu_dash_arr, 1);
            phi_dash += -2 * alpha0 * r_phi;
            daxpy_(nsky, alpha0 * alpha0, v_rho_arr, 1, rho_dash_arr, 1);
            daxpy_(nsrc, alpha0 * alpha0, v_nu_arr, 1, nu_dash_arr, 1);
            phi_dash += alpha0 * alpha0 * v_phi;

            int nneg_tmp = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_dash_arr[isky] < 0.0){
                    nneg_tmp ++;
                }
            }
            for(int isrc = 0; isrc < nsrc; isrc ++){
                if(nu_dash_arr[isrc] < 0.0){
                    nneg_tmp ++;
                }
            }
            if(phi_dash < 0.0){
                nneg_tmp ++;
            }
            if (nneg_tmp > 0){
                continue;
            }

            double helldist_dc3 = 0.0;
            int flag_converge_dc3 = 0;
            GetMvalArrNvalArrPval(rho_dash_arr, nu_dash_arr, phi_dash,
                                  data_arr, bg_arr, det_fpsrc_arr,
                                  resp_norm_mat_arr,
                                  ndet, nsky, nsrc,
                                  mval_arr, nval_arr, &pval);
            GetRhoNuPhi_ByDC(fp_log,
                             rho_dash_arr, nu_dash_arr, phi_dash,
                             mval_arr, nval_arr, pval,
                             nph, B_val,
                             ndet, nskyx, nskyy, nsrc,
                             mu,
                             ndc, tol_dc,
                             npm, tol_pm,
                             nnewton, tol_newton,
                             rho_0_new_arr,
                             nu_0_new_arr,
                             &phi_0_new,
                             &helldist_dc3,
                             &flag_converge_dc3);
            if (flag_converge_dc3 == 0){
                MiIolib::Printf2(fp_log,
                                 "iem = %d: dc3: not converged: helldist_dc3 = %.2e\n",
                                 iem,
                                 helldist_dc3);
            }
            int nneg = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_0_new_arr[isky] < 0.0){
                    nneg ++;
                }
            }
            for(int isrc = 0; isrc < nsrc; isrc ++){
                if(nu_0_new_arr[isrc] < 0.0){
                    nneg ++;
                }
            }
            if(phi_0_new < 0.0){
                nneg ++;
            }
            if (nneg == 0){
                ifind_nonneg = 1;
                break;
            }
        }
        if(ifind_nonneg == 0){
            MiIolib::Printf2(fp_log, "warning: iem = %d, ifind_nonneg == 0\n",
                             iem);
        }
        
        delete [] mval_arr;
        delete [] nval_arr;

        double helldist  = GetHellingerDist(rho_0_arr, nu_0_arr, phi_0,
                                            rho_0_new_arr, nu_0_new_arr, phi_0_new,
                                            nsky, nsrc);
        if (access( "/tmp/fpsrc_smth_bg_em_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/fpsrc_smth_bg_em_stop file is found, then stop.\n");
            break;
        }
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
            break;
        }
        dcopy_(nsky, rho_0_new_arr, 1, rho_0_arr, 1);
        dcopy_(nsrc, nu_0_new_arr, 1, nu_0_arr, 1);
        phi_0 = phi_0_new;

        //double lval = 0.0;        
        if (iem % 100 == 0){
            //lval = GetFuncL(data_arr, bg_arr,
            //rho_new_arr, nu_new,
            //                resp_norm_mat_arr,
            //                ndet, nsky);
            //printf("iem = %d, helldist = %e, lval = %e\n",
            //       iem, helldist, lval);
        } else {
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
        }
    }

    dcopy_(nsky, rho_0_new_arr, 1, rho_new_arr, 1);
    dcopy_(nsrc, nu_0_new_arr, 1, nu_new_arr, 1);
    double phi_new = phi_0_new;

    delete [] rho_0_arr;
    delete [] rho_1_arr;
    delete [] rho_2_arr;
    delete [] rho_dash_arr;
    delete [] r_rho_arr;
    delete [] r2_rho_arr;
    delete [] v_rho_arr;
    delete [] rho_0_new_arr;

    delete [] nu_0_arr;
    delete [] nu_1_arr;
    delete [] nu_2_arr;
    delete [] nu_dash_arr;
    delete [] r_nu_arr;
    delete [] r2_nu_arr;
    delete [] v_nu_arr;
    delete [] nu_0_new_arr;
    
    *phi_new_ptr = phi_new;
}
