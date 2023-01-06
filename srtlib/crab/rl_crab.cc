#include<unistd.h>
#include "mir_math.h"
#include "rl_crab.h"

void SrtlibRlCrab::GetDetArr(const double* const sky_arr,
                             const double* const resp_norm_mat_arr,
                             int ndet, int nsky,
                             double* const det_arr) // ndet
{
    // det_arr = R_mat %*% sky_arr
    char transa[2];
    strcpy(transa, "N");

    // y := alpha*A*x + beta*y
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_arr), ndet,
           const_cast<double*>(sky_arr), 1,
           0.0, det_arr, 1);
}

void SrtlibRlCrab::GetDenArr(const double* const sky_arr,
                             const double* const flux_arr,
                             const double* const det_0_arr,
                             const double* const bg_arr,
                             const double* const resp_norm_mat_arr,
                             int ndet, int nsky, int nphase,
                             double* const* const den_arr)
{
    double* det_arr = new double[ndet];
    SrtlibRlCrab::GetDetArr(sky_arr, resp_norm_mat_arr,
                            ndet, nsky, det_arr);
    for(int iphase = 0; iphase < nphase; iphase++){
        dcopy_(ndet, const_cast<double*>(bg_arr), 1, den_arr[iphase], 1);
        daxpy_(ndet, flux_arr[iphase], const_cast<double*>(det_0_arr), 1,
               den_arr[iphase], 1);
        daxpy_(ndet, 1.0, det_arr, 1, den_arr[iphase], 1);
    }
    delete [] det_arr;
}

void SrtlibRlCrab::GetYDashArr(const double* const* const data_arr,
                               const double* const* const den_arr,
                               int ndet, int nphase,
                               double* const* const y_dash_arr)
{
    for(int iphase = 0; iphase < nphase; iphase++){
        for(int idet = 0; idet < ndet; idet++){
            y_dash_arr[iphase][idet] = data_arr[iphase][idet]
                / den_arr[iphase][idet];
        }
    }
}

void SrtlibRlCrab::GetMvalArr(const double* const* const y_dash_arr,
                              const double* const resp_norm_mat_arr,
                              const double* const sky_arr,
                              int ndet, int nsky, int nphase,
                              double* const mval_arr)
{
    double* y_dash_sum_arr = new double[ndet];
    dcopy_(ndet, const_cast<double*>(y_dash_arr[0]), 1, y_dash_sum_arr, 1);
    for(int iphase = 1; iphase < nphase; iphase++){
        daxpy_(ndet, 1.0, const_cast<double*>(y_dash_arr[iphase]), 1,
               y_dash_sum_arr, 1);
    }

    double* coeff_arr = new double[nsky];
    char transa[2];
    strcpy(transa, "T");
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_arr), ndet,
           y_dash_sum_arr, 1,
           0.0, coeff_arr, 1);
    MibBlas::ElmWiseMul(nsky, 1.0,
                        coeff_arr, sky_arr,
                        mval_arr);
    delete [] y_dash_sum_arr;
    delete [] coeff_arr;
}

void SrtlibRlCrab::GetNvalArr(const double* const* const y_dash_arr,
                              const double* const flux_arr,
                              const double* const det_0_arr,
                              int ndet, int nphase,
                              double* const nval_arr)
{
    for(int iphase = 0; iphase < nphase; iphase++){
        nval_arr[iphase] = flux_arr[iphase]
            * ddot_(ndet, const_cast<double*>(y_dash_arr[iphase]), 1,
                    const_cast<double*>(det_0_arr), 1);
    }
}



//
//void SrtlibRlCrab::RichlucyCrab(FILE* const fp_log,
//                                const double* const sky_init_arr,
//                                const double* const flux_init_arr,
//                                const double* const* const data_arr,
//                                const double* const phase_arr,
//                                const double* const det_0_arr,
//                                const double* const resp_norm_mat_arr,
//                                int ndet, int nsky, int nphase,
//                                string outdir, string outfile_head,
//                                int nem, double tol_em,
//                                double* const sky_new_arr,
//                                double* const flux_new_arr)
//{
//    double* sky_pre_arr = new double[nsky];
//    double* flux_pre_arr = new double[nphase];
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_pre_arr, 1);
//    dcopy_(nphase, const_cast<double*>(flux_init_arr), 1, flux_pre_arr, 1);
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrab::GetSkyFluxNewArr(sky_pre_arr, flux_pre_arr,
//                                     data_arr, phase_arr, det_0_arr,
//                                     resp_norm_mat_arr,
//                                     ndet, nsky, nphase,
//                                     sky_new_arr, flux_new_arr);
//        double helldist  = SrtlibRlStatvalCrab::GetHellingerDist(
//            sky_pre_arr, flux_pre_arr,
//            sky_new_arr, flux_new_arr,
//            nsky, nphase);
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//        MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
//                         iem, helldist);
//        
//        if (helldist < tol_em){
//            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
//                             iem, helldist);
//            break;
//        }
//        dcopy_(nsky, const_cast<double*>(sky_new_arr), 1, sky_pre_arr, 1);
//        dcopy_(nphase, const_cast<double*>(flux_new_arr), 1, flux_pre_arr, 1);
//    }
//    delete [] sky_pre_arr;
//    delete [] flux_pre_arr;
//}
//
//
////
//// accerelated richardson lucy by Zhou-Alexander-Lange:
//// H.Zhou, D.Alexander, K.Lange,
//// "A quasi-Newton acceleration for high-dimensional
//// optimization algorithms", Stat Comput (2011) 21, 261.
//// Case: q = 1
//void SrtlibRlCrab::RichlucyCrabAccZALq1(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const flux_init_arr,
//    const double* const* const data_arr,
//    const double* const phase_arr,
//    const double* const det_0_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nsky, int nphase,
//    string outdir, string outfile_head,
//    int nem, double tol_em,
//    double* const sky_new_arr,
//    double* const flux_new_arr)
//{
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
//    
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrab::GetSkyFluxNewArr(sky_0_arr, flux_0_arr, data_arr,
//                                     phase_arr, det_0_arr,
//                                     resp_norm_mat_arr,
//                                     ndet, nsky, nphase,
//                                     sky_1_arr, flux_1_arr);
//        SrtlibRlCrab::GetSkyFluxNewArr(sky_1_arr, flux_1_arr, data_arr,
//                                     phase_arr, det_0_arr,
//                                     resp_norm_mat_arr,
//                                     ndet, nsky, nphase,
//                                     sky_2_arr, flux_2_arr);
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
//        printf("cval = %e\n", cval);
//
//        if (cval < 0.0){
//            // usual update
//            dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//        } else{
//            int nk = 10000;
//            double eta = 0.8;
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
//                    flag_find = 1;
//                    break;
//                }
//            }
//            if(flag_find == 0){
//                // usual update
//                dcopy_(nsky, sky_2_arr, 1, sky_0_new_arr, 1);
//                dcopy_(nphase, flux_2_arr, 1, flux_0_new_arr, 1);
//            }
//        }
//
//        double helldist  = SrtlibRlStatvalCrab::GetHellingerDist(
//            sky_0_arr, flux_0_arr,
//            sky_0_new_arr, flux_0_new_arr,
//            nsky, nphase);
//        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n", iem, helldist);
//        if (helldist < tol_em){
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
//       
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
//
////// accerelated richardson lucy (SQUAREM)
////void SrtlibRlCrab::RichlucyCrabAccSquarem(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const flux_init_arr,
//    const double* const* const data_arr,
//    const double* const phase_arr,
//    const double* const det_0_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nsky, int nphase,
//    string outdir, string outfile_head,
//    int nem, double tol_em,
//    double* const sky_new_arr,
//    double* const flux_new_arr)
//{
//    double* sky_0_arr  = new double[nsky];
//    double* sky_1_arr  = new double[nsky];
//    double* sky_2_arr  = new double[nsky];
//    double* sky_dash_arr  = new double[nsky];
//    double* r_sky_arr  = new double[nsky];
//    double* r2_sky_arr  = new double[nsky];
//    double* v_sky_arr  = new double[nsky];
//    double* sky_0_new_arr  = new double[nsky];
//    double* diff_sky_0_arr = new double[nsky];
//
//    double* flux_0_arr  = new double[nphase];
//    double* flux_1_arr  = new double[nphase];
//    double* flux_2_arr  = new double[nphase];
//    double* flux_dash_arr  = new double[nphase];
//    double* r_flux_arr  = new double[nphase];
//    double* r2_flux_arr  = new double[nphase];
//    double* v_flux_arr  = new double[nphase];
//    double* flux_0_new_arr  = new double[nphase];
//    double* diff_flux_0_arr = new double[nphase];
//    
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_0_arr, 1);
//    dcopy_(nphase, const_cast<double*>(flux_init_arr), 1, flux_0_arr, 1);
//    
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrab::GetSkyFluxNewArr(sky_0_arr, flux_0_arr, data_arr,
//                                     phase_arr, det_0_arr,
//                                     resp_norm_mat_arr,
//                                     ndet, nsky, nphase,
//                                     sky_1_arr, flux_1_arr);
//        SrtlibRlCrab::GetSkyFluxNewArr(sky_1_arr, flux_1_arr, data_arr,
//                                     phase_arr, det_0_arr,
//                                     resp_norm_mat_arr,
//                                     ndet, nsky, nphase,
//                                     sky_2_arr, flux_2_arr);
//        MibBlas::Sub(sky_1_arr, sky_0_arr, nsky, r_sky_arr);
//        MibBlas::Sub(flux_1_arr, flux_0_arr, nphase, r_flux_arr);
//        MibBlas::Sub(sky_2_arr, sky_1_arr, nsky, r2_sky_arr);
//        MibBlas::Sub(flux_2_arr, flux_1_arr, nphase, r2_flux_arr);
//        MibBlas::Sub(r2_sky_arr, r_sky_arr, nsky, v_sky_arr);
//        MibBlas::Sub(r2_flux_arr, r_flux_arr, nphase, v_flux_arr);        
//
//        double r_norm2 = ddot_(nsky, r_sky_arr, 1, r_sky_arr, 1)
//            + ddot_(nphase, r_flux_arr, 1, r_flux_arr, 1);
//        double v_norm2 = ddot_(nsky, v_sky_arr, 1, v_sky_arr, 1)
//            + ddot_(nphase, v_flux_arr, 1, v_flux_arr, 1);
//        double alpha = -1.0 * sqrt(r_norm2 / v_norm2);
//        // printf("alpha = %e\n", alpha);
//
//        // Zhou 2011 
//        double d_norm2 = ddot_(nsky, r_sky_arr, 1, v_sky_arr, 1)
//            + ddot_(nphase, r_flux_arr, 1, v_flux_arr, 1);
//        double beta = -1.0 * r_norm2 / d_norm2;
//        printf("beta = %e\n", beta);
//
//        if (beta < 0.0){
//            // usual update
//            dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//            dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//        } else{
//            int nk = 10000;
//            double eta = 0.8;
//            int flag_find = 0;
//            for (int ik = 0; ik < nk; ik ++){
//                double beta0 = beta * pow(eta, ik);
//                dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//                dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//                dscal_(nsky, (1.0 - beta0), sky_0_new_arr, 1);
//                dscal_(nphase, (1.0 - beta0), flux_0_new_arr, 1);
//                daxpy_(nsky, beta0, sky_2_arr, 1, sky_0_new_arr, 1);
//                daxpy_(nphase, beta0, flux_2_arr, 1, flux_0_new_arr, 1);
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
//                    flag_find = 1;
//                    break;
//                }
//            }
//            if(flag_find == 0){
//                // usual update
//                dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//                dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//            }
//        }
//
////        int nk = 10000;
////        double eta = 0.8;
////        //double negloglike_pre = 0.0;
////        //double epsilon = 1e-20;
////        //int ifind = 0;
////        int ifind_nonneg = 0;
////        for (int ik = 0; ik < nk; ik ++){
////            double alpha0 = alpha * pow(eta, ik);
////            dcopy_(nsky, sky_0_arr, 1, sky_dash_arr, 1);
////            dcopy_(nphase, flux_0_arr, 1, flux_dash_arr, 1);
////            daxpy_(nsky, -2 * alpha0, r_sky_arr, 1, sky_dash_arr, 1);
////            daxpy_(nphase, -2 * alpha0, r_flux_arr, 1, flux_dash_arr, 1);
////            daxpy_(nsky, alpha0 * alpha0, v_sky_arr, 1, sky_dash_arr, 1);
////            daxpy_(nphase, alpha0 * alpha0, v_flux_arr, 1, sky_dash_arr, 1);
////
////            int nneg_tmp = 0;
////            for(int isky = 0; isky < nsky; isky ++){
////                if(sky_dash_arr[isky] < 0.0){
////                    nneg_tmp ++;
////                }
////            }
////            for(int iphase = 0; iphase < nphase; iphase ++){
////                if(flux_dash_arr[iphase] < 0.0){
////                    nneg_tmp ++;
////                }
////            }
////            if (nneg_tmp > 0){
////                continue;
////            }
////            SrtlibRlCrab::GetSkyFluxNewArr(sky_dash_arr, flux_dash_arr, data_arr,
////                                         phase_arr, det_0_arr,
////                                         resp_norm_mat_arr,
////                                         ndet, nsky, nphase,
////                                         sky_0_new_arr, flux_0_new_arr);
////            int nneg = 0;
////            for(int isky = 0; isky < nsky; isky ++){
////                if(sky_0_new_arr[isky] < 0.0){
////                    nneg ++;
////                }
////            }
////            for(int iphase = 0; iphase < nphase; iphase ++){
////                if(flux_0_new_arr[iphase] < 0.0){
////                    nneg ++;
////                }
////            }
////            if (nneg == 0){
////                ifind_nonneg = 1;
////                break;
////                //                if (ifind == 0){
////                //                    negloglike_pre = GetNegLogLike(sky_0_new_arr,
////                //                                                   data_arr,
////                //                                                   resp_norm_mat_arr,
////                //                                                   ndet, nsky, epsilon);
////                //                    ifind = 1;
////                //                } else if (ifind == 1){
////                //                    double negloglike = GetNegLogLike(sky_0_new_arr,
////                //                                                      data_arr,
////                //                                                      resp_norm_mat_arr,
////                //                                                      ndet, nsky, epsilon);
////                //                    //MiIolib::Printf2(fp_log,
////                //                    //"ik = %d, alpha0 = %e, negloglike = %e\n",
////                //                    //ik, alpha0, negloglike);
////                //                    if(negloglike > negloglike_pre){
////                //                        break;
////                //                    }
////                //                    negloglike_pre = negloglike;
////                // }
////            }
////        }
//        
//        //if(ifind_nonneg == 0){
//        //    MiIolib::Printf2(fp_log, "warning: iem = %d, ifind_nonneg == 0\n", iem);
//        //// usual update
//        //    dcopy_(nsky, sky_1_arr, 1, sky_0_new_arr, 1);
//        //    dcopy_(nphase, flux_1_arr, 1, flux_0_new_arr, 1);
//        //}
//        
//        // double sum = MibBlas::Sum(sky_0_new_arr, nsky);
//        // printf("sum = %e\n", sum);
//        double helldist  = SrtlibRlStatvalCrab::GetHellingerDist(
//            sky_0_arr, flux_0_arr,
//            sky_0_new_arr, flux_0_new_arr,
//            nsky, nphase);
//        
//        //MibBlas::Sub(sky_0_new_arr, sky_0_arr, nsky, diff_sky_0_arr);
//        //double diff = sqrt(ddot_(nsky, diff_sky_0_arr, 1,
//        //                         diff_sky_0_arr, 1));
//        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n", iem, helldist);
//        if (helldist < tol_em){
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
//       
//    }
//    dcopy_(nsky, sky_0_new_arr, 1, sky_new_arr, 1);
//    dcopy_(nphase, flux_0_new_arr, 1, flux_new_arr, 1);
//
//    delete [] sky_0_arr;
//    delete [] sky_1_arr;
//    delete [] sky_2_arr;
//    delete [] sky_dash_arr;
//    delete [] r_sky_arr;
//    delete [] r2_sky_arr;
//    delete [] v_sky_arr;
//    delete [] sky_0_new_arr;
//    delete [] diff_sky_0_arr;
//    
//    delete [] flux_0_arr;
//    delete [] flux_1_arr;
//    delete [] flux_2_arr;
//    delete [] flux_dash_arr;
//    delete [] r_flux_arr;
//    delete [] r2_flux_arr;
//    delete [] v_flux_arr;
//    delete [] flux_0_new_arr;
//    delete [] diff_flux_0_arr;
//}
//
//
//
