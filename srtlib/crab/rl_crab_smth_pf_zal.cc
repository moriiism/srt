#include "rl_crab.h"
#include "rl_crab_statval.h"
#include "rl_crab_smth_pf_zal.h"

void SrtlibRlCrabSmthPfZal::GetSkyNew(
    const double* const alpha_arr,
    const double* const beta_arr,
    const double* const mval_arr,
    double live_time_ratio_ave,
    int nsky, double mu,
    double* const sky_new_arr)
{
    for(int isky = 0; isky < nsky; isky++){    
        double bval = live_time_ratio_ave
            - mu * beta_arr[isky];
        double num = -1.0 * bval
            + sqrt(bval * bval
                   + 4.0 * mu * alpha_arr[isky]
                   * mval_arr[isky]);
        double den = 2.0 * mu * alpha_arr[isky];
        sky_new_arr[isky] = num / den;
    }
}

void SrtlibRlCrabSmthPfZal::GetSkyNewArr(
    const double* const sky_arr,
    const double* const mval_arr,
    double live_time_ratio_ave,    
    int nskyx, int nskyy, double mu,
    double* const sky_new_arr)
{
    int nsky = nskyx * nskyy;
    double* alpha_arr = new double[nsky];
    SrtlibSmthZal::GetDerivUAlphaArr(
        nskyx, nskyy, alpha_arr);
    double* beta_arr = new double[nsky];
    SrtlibSmthZal::GetDerivUBetaArr(
        sky_arr, nskyx, nskyy, beta_arr);
    SrtlibRlCrabSmthPfZal::GetSkyNew(
        alpha_arr, beta_arr, mval_arr, live_time_ratio_ave,
        nsky, mu, sky_new_arr);
    delete [] alpha_arr;
    delete [] beta_arr;
}

void SrtlibRlCrabSmthPfZal::GetFluxNewArr(
    const double* const nval_arr,
    const double* const flux_target_arr,
    const double* const phase_arr,
    const double* const live_time_ratio_arr,
    int nphase, double gamma,
    double* const flux_new_arr)
{
    for(int iphase = 0; iphase < nphase; iphase++){
        double bval =
            phase_arr[iphase] * live_time_ratio_arr[iphase]
            - 2.0 * gamma * flux_target_arr[iphase];
        double num = -1 * bval
            + sqrt(bval * bval
                   + 8.0 * gamma * nval_arr[iphase]);
        double den = 4.0 * gamma;
        flux_new_arr[iphase] = num / den;
    }
}

void SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
    const double* const sky_pre_arr,
    const double* const flux_pre_arr,
    const double* const* const data_arr,
    const double* const bg_arr,
    const double* const flux_target_arr,
    const double* const phase_arr,
    const double* const live_time_ratio_arr,
    const double* const det_0_arr,
    const double* const resp_norm_mat_arr,    
    int ndet, int nskyx, int nskyy, int nphase,
    double mu, double gamma,
    double* const sky_new_arr,
    double* const flux_new_arr)
{
    int nsky = nskyx * nskyy;
    double** den_arr = new double*[nphase];
    double** y_dash_arr = new double*[nphase];
    for(int iphase = 0; iphase < nphase; iphase++){
        den_arr[iphase] = new double[ndet];
        y_dash_arr[iphase] = new double[ndet];
    }
    double* mval_arr = new double[nsky];
    double* nval_arr = new double[nphase];
        
    SrtlibRlCrab::GetDenArr(sky_pre_arr,
                            flux_pre_arr,
                            det_0_arr,
                            bg_arr,
                            resp_norm_mat_arr,
                            ndet, nsky, nphase,
                            den_arr);
    SrtlibRlCrab::GetYDashArr(data_arr,
                              den_arr,
                              ndet, nphase,
                              y_dash_arr);
    SrtlibRlCrab::GetMvalArr(y_dash_arr,
                             resp_norm_mat_arr,
                             sky_pre_arr,
                             ndet, nsky, nphase,
                             mval_arr);
    SrtlibRlCrab::GetNvalArr(y_dash_arr,
                             flux_pre_arr,
                             det_0_arr,
                             ndet, nphase,
                             nval_arr);
    double live_time_ratio_ave = 0.0;
    live_time_ratio_ave = ddot_(
        nphase, const_cast<double*>(phase_arr), 1,
        const_cast<double*>(live_time_ratio_arr), 1);
    SrtlibRlCrabSmthPfZal::GetSkyNewArr(sky_pre_arr,
                                        mval_arr,
                                        live_time_ratio_ave,
                                        nskyx, nskyy, mu,
                                        sky_new_arr);
    SrtlibRlCrabSmthPfZal::GetFluxNewArr(nval_arr,
                                         flux_target_arr,
                                         phase_arr,
                                         live_time_ratio_arr,
                                         nphase, gamma,
                                         flux_new_arr);
    for(int iphase = 0; iphase < nphase; iphase++){
        delete [] den_arr[iphase];
        delete [] y_dash_arr[iphase];
    }
    delete [] den_arr;
    delete [] y_dash_arr;
    delete [] mval_arr;
    delete [] nval_arr;
}

void SrtlibRlCrabSmthPfZal::RichlucyCrabSmthPfZal(
    FILE* const fp_log,
    const double* const sky_init_arr,
    const double* const flux_init_arr,
    const double* const* const data_arr,
    const double* const bg_arr,
    const double* const flux_target_arr,
    const double* const phase_arr,
    const double* const live_time_ratio_arr,    
    const double* const det_0_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nskyx, int nskyy, int nphase,
    double mu, double gamma,
    string outdir,
    string outfile_head,
    int nem, double tol_em,
    double* const sky_new_arr,
    double* const flux_new_arr)
{
    int nsky = nskyx * nskyy;
    double* sky_pre_arr = new double[nsky];
    double* flux_pre_arr = new double[nphase];
    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1,
           sky_pre_arr, 1);
    dcopy_(nphase, const_cast<double*>(flux_init_arr), 1,
           flux_pre_arr, 1);
    for(int iem = 0; iem < nem; iem ++){
        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
            sky_pre_arr,
            flux_pre_arr,
            data_arr,
            bg_arr,
            flux_target_arr,
            phase_arr,
            live_time_ratio_arr,
            det_0_arr,
            resp_norm_mat_arr,    
            ndet, nskyx, nskyy, nphase,
            mu, gamma,
            sky_new_arr,
            flux_new_arr);
        double helldist
            = SrtlibCrabRlCrabStatval::GetHellingerDist(
            sky_pre_arr, flux_pre_arr,
            sky_new_arr, flux_new_arr,
            nsky, nphase);
        if (access( "/tmp/rl_crab_smth_pf_zal_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_crab_smth_pf_zal_stop file is found, "
                "then stop.\n");
            break;
        }
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n",
                             iem, helldist);
            break;
        }
        dcopy_(nsky, sky_new_arr, 1, sky_pre_arr, 1);
        dcopy_(nphase, flux_new_arr, 1, flux_pre_arr, 1);
        if(iem % 100 == 0){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n",
                             iem, helldist);
        }
    }
    delete [] sky_pre_arr;
    delete [] flux_pre_arr;
}

//// Accerelated by Zhou-Alexander-Lange,
//// which is H.Zhou, D.Alexander, K.Lange,
//// "A quasi-Newton acceleration for high-dimensional
//// optimization algorithms", Stat Comput (2011) 21, 261.
//// Case: q = 1
//void SrtlibRlCrabSmthPfZal::RichlucyCrabSmthPfZalQ1(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const flux_init_arr,
//    const double* const* const data_arr,
//    const double* const bg_arr,    
//    const double* const flux_target_arr,
//    const double* const phase_arr,
//    const double* const det_0_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nskyx, int nskyy, int nphase,
//    double mu, double gamma,
//    string outdir,
//    string outfile_head,
//    int nem, double tol_em,
//    double* const sky_new_arr,
//    double* const flux_new_arr)
//{
//    int nsky = nskyx * nskyy;
//    double* sky_0_arr  = new double[nsky];
//    double* sky_1_arr  = new double[nsky];
//    double* sky_2_arr  = new double[nsky];
//    double* u_sky_arr  = new double[nsky];
//    double* v_sky_arr  = new double[nsky];
//    double* diff_sky_arr  = new double[nsky];
//    double* sky_0_new_arr  = new double[nsky];
//
//    double* flux_0_arr  = new double[nphase];
//    double* flux_1_arr  = new double[nphase];
//    double* flux_2_arr  = new double[nphase];
//    double* u_flux_arr  = new double[nphase];
//    double* v_flux_arr  = new double[nphase];
//    double* diff_flux_arr  = new double[nphase];
//    double* flux_0_new_arr  = new double[nphase];
//    
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_0_arr, 1);
//    dcopy_(nphase, const_cast<double*>(flux_init_arr), 1, flux_0_arr, 1);
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_0_arr,
//            flux_0_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_1_arr,
//            flux_1_arr);
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_1_arr,
//            flux_1_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_2_arr,
//            flux_2_arr);
//
//        MibBlas::Sub(sky_1_arr, sky_0_arr, nsky, u_sky_arr);
//        MibBlas::Sub(flux_1_arr, flux_0_arr, nphase, u_flux_arr);        
//        MibBlas::Sub(sky_2_arr, sky_1_arr, nsky, v_sky_arr);
//        MibBlas::Sub(flux_2_arr, flux_1_arr, nphase, v_flux_arr);        
//        MibBlas::Sub(v_sky_arr, u_sky_arr, nsky, diff_sky_arr);
//        MibBlas::Sub(v_flux_arr, u_flux_arr, nphase, diff_flux_arr);
//
//        double num = ddot_(nsky, u_sky_arr, 1, u_sky_arr, 1)
//            + ddot_(nphase, u_flux_arr, 1, u_flux_arr, 1);
//        double den = ddot_(nsky, u_sky_arr, 1, diff_sky_arr, 1)
//            + ddot_(nphase, u_flux_arr, 1, diff_flux_arr, 1);
//        double cval = -1.0 * num / den;
//        MiIolib::Printf2(fp_log, "cval = %e\n", cval);
//
//        if (cval < 0.0){
//            // usual update
//            MiIolib::Printf2(fp_log, "cval < 0, then update by sky1 flux1\n");
//            dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//        } else {
//            int nk = 100;
//            double eta = 0.5;
//            int flag_find = 0;
//            for (int ik = 0; ik < nk; ik ++){
//                double cval0 = cval * pow(eta, ik);
//                dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//                dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//                dscal_(nsky, (1.0 - cval0), sky_0_new_arr, 1);
//                dscal_(nphase, (1.0 - cval0), flux_0_new_arr, 1);
//                daxpy_(nsky, cval0, sky_2_arr, 1, sky_0_new_arr, 1);
//                daxpy_(nphase, cval0, flux_2_arr, 1, flux_0_new_arr, 1);
//                
//                int nneg_tmp = 0;
//                for(int isky = 0; isky < nsky; isky ++){
//                    if(sky_0_new_arr[isky] < 0.0){
//                        nneg_tmp ++;
//                    }
//                }
//                for(int iphase = 0; iphase < nphase; iphase ++){
//                    if(flux_0_new_arr[iphase] < 0.0){
//                        nneg_tmp ++;
//                    }
//                }
//                if (nneg_tmp > 0){
//                    continue;
//                } else{
//                    MiIolib::Printf2(fp_log, "cval0 = %e, ik = %d\n",
//                                     cval0, ik);
//                    flag_find = 1;
//                    break;
//                }
//            }
//            if(flag_find == 0){
//                // usual update
//                MiIolib::Printf2(fp_log, "flag_find = 0, then update by sky2 flux2\n");
//                dcopy_(nsky, sky_2_arr, 1, sky_0_new_arr, 1);
//                dcopy_(nphase, flux_2_arr, 1, flux_0_new_arr, 1);
//            }
//        }
//        double helldist  = SrtlibCrabRlCrabStatval::GetHellingerDist(
//            sky_0_arr, flux_0_arr,
//            sky_0_new_arr, flux_0_new_arr,
//            nsky, nphase);
//        MiIolib::Printf2(
//            fp_log, "iem = %d, helldist = %e\n",
//            iem, helldist);
//        if (helldist < tol_em){
//            MiIolib::Printf2(
//                fp_log, "iem = %d, helldist = %e\n",
//                iem, helldist);
//            break;
//        }
//        dcopy_(nsky, sky_0_new_arr, 1, sky_0_arr, 1);
//        dcopy_(nphase, flux_0_new_arr, 1, flux_0_arr, 1);
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//    }
//    dcopy_(nsky, sky_0_new_arr, 1, sky_new_arr, 1);
//    dcopy_(nphase, flux_0_new_arr, 1, flux_new_arr, 1);
//
//    delete [] sky_0_arr;
//    delete [] sky_1_arr;
//    delete [] sky_2_arr;
//    delete [] u_sky_arr;
//    delete [] v_sky_arr;
//    delete [] diff_sky_arr;
//    delete [] sky_0_new_arr;
//
//    delete [] flux_0_arr;
//    delete [] flux_1_arr;
//    delete [] flux_2_arr;
//    delete [] u_flux_arr;
//    delete [] v_flux_arr;
//    delete [] diff_flux_arr;
//    delete [] flux_0_new_arr;
//}
//
//// Accerelated by Squarem S3
//void SrtlibRlCrabSmthPfZal::RichlucyCrabSmthPfSqS3(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const flux_init_arr,
//    const double* const* const data_arr,
//    const double* const bg_arr,    
//    const double* const flux_target_arr,
//    const double* const phase_arr,
//    const double* const det_0_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nskyx, int nskyy, int nphase,
//    double mu, double gamma,
//    string outdir,
//    string outfile_head,
//    int nem, double tol_em,
//    double* const sky_new_arr,
//    double* const flux_new_arr)
//{
//    int nsky = nskyx * nskyy;
//    double* sky_0_arr  = new double[nsky];
//    double* sky_1_arr  = new double[nsky];
//    double* sky_2_arr  = new double[nsky];
//    double* u_sky_arr  = new double[nsky];
//    double* v_sky_arr  = new double[nsky];
//    double* diff_sky_arr  = new double[nsky];
//    double* sky_0_new_arr  = new double[nsky];
//
//    double* flux_0_arr  = new double[nphase];
//    double* flux_1_arr  = new double[nphase];
//    double* flux_2_arr  = new double[nphase];
//    double* u_flux_arr  = new double[nphase];
//    double* v_flux_arr  = new double[nphase];
//    double* diff_flux_arr  = new double[nphase];
//    double* flux_0_new_arr  = new double[nphase];
//    
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_0_arr, 1);
//    dcopy_(nphase, const_cast<double*>(flux_init_arr), 1, flux_0_arr, 1);
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_0_arr,
//            flux_0_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_1_arr,
//            flux_1_arr);
//        SrtlibRlCrabSmthPfZal::GetSkyFluxNewArr(
//            sky_1_arr,
//            flux_1_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_0_arr,
//            resp_norm_mat_arr,    
//            ndet, nskyx, nskyy, nphase,
//            mu, gamma,
//            sky_2_arr,
//            flux_2_arr);
//
//        MibBlas::Sub(sky_1_arr, sky_0_arr, nsky, u_sky_arr);
//        MibBlas::Sub(flux_1_arr, flux_0_arr, nphase, u_flux_arr);        
//        MibBlas::Sub(sky_2_arr, sky_1_arr, nsky, v_sky_arr);
//        MibBlas::Sub(flux_2_arr, flux_1_arr, nphase, v_flux_arr);        
//        MibBlas::Sub(v_sky_arr, u_sky_arr, nsky, diff_sky_arr);
//        MibBlas::Sub(v_flux_arr, u_flux_arr, nphase, diff_flux_arr);
//
//        double num = ddot_(nsky, u_sky_arr, 1, u_sky_arr, 1)
//            + ddot_(nphase, u_flux_arr, 1, u_flux_arr, 1);
//        double den = ddot_(nsky, diff_sky_arr, 1, diff_sky_arr, 1)
//            + ddot_(nphase, diff_flux_arr, 1, diff_flux_arr, 1);
//        double sval = -1.0 * sqrt(num / den);
//        MiIolib::Printf2(fp_log, "sval = %e\n", sval);
//
//        int nk = 100;
//        double eta = 0.5;
//        int flag_find = 0;
//        for (int ik = 0; ik < nk; ik ++){
//            double sval0 = sval * pow(eta, ik);
//            dcopy_(nsky, sky_0_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_0_arr, 1, flux_0_new_arr, 1);
//            daxpy_(nsky, -2.0 * sval0, u_sky_arr, 1, sky_0_new_arr, 1);
//            daxpy_(nphase, -2.0 * sval0, u_flux_arr, 1, flux_0_new_arr, 1);
//            daxpy_(nsky, sval0 * sval0, diff_sky_arr, 1, sky_0_new_arr, 1);
//            daxpy_(nphase, sval0 * sval0, diff_flux_arr, 1, flux_0_new_arr, 1);
//
//            int nneg_tmp = 0;
//            for(int isky = 0; isky < nsky; isky ++){
//                if(sky_0_new_arr[isky] < 0.0){
//                    nneg_tmp ++;
//                }
//            }
//            for(int iphase = 0; iphase < nphase; iphase ++){
//                if(flux_0_new_arr[iphase] < 0.0){
//                    nneg_tmp ++;
//                }
//            }
//            if (nneg_tmp > 0){
//                continue;
//            } else{
//                MiIolib::Printf2(fp_log, "sval0 = %e, ik = %d\n",
//                                 sval0, ik);
//                flag_find = 1;
//                break;
//            }
//        }
//        if(flag_find == 0){
//            // usual update
//            MiIolib::Printf2(fp_log,
//                             "flag_find = 0, then update by sky2 flux2\n");
//            dcopy_(nsky, sky_2_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_2_arr, 1, flux_0_new_arr, 1);
//        }
//        double helldist  = SrtlibCrabRlCrabStatval::GetHellingerDist(
//            sky_0_arr, flux_0_arr,
//            sky_0_new_arr, flux_0_new_arr,
//            nsky, nphase);
//        MiIolib::Printf2(
//            fp_log, "iem = %d, helldist = %e\n",
//            iem, helldist);
//        if (helldist < tol_em){
//            MiIolib::Printf2(
//                fp_log, "iem = %d, helldist = %e\n",
//                iem, helldist);
//            break;
//        }
//        dcopy_(nsky, sky_0_new_arr, 1, sky_0_arr, 1);
//        dcopy_(nphase, flux_0_new_arr, 1, flux_0_arr, 1);
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//    }
//    dcopy_(nsky, sky_0_new_arr, 1, sky_new_arr, 1);
//    dcopy_(nphase, flux_0_new_arr, 1, flux_new_arr, 1);
//
//    delete [] sky_0_arr;
//    delete [] sky_1_arr;
//    delete [] sky_2_arr;
//    delete [] u_sky_arr;
//    delete [] v_sky_arr;
//    delete [] diff_sky_arr;
//    delete [] sky_0_new_arr;
//
//    delete [] flux_0_arr;
//    delete [] flux_1_arr;
//    delete [] flux_2_arr;
//    delete [] u_flux_arr;
//    delete [] v_flux_arr;
//    delete [] diff_flux_arr;
//    delete [] flux_0_new_arr;
//}

