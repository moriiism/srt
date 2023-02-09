#include<unistd.h>
#include "rl.h"
#include "rl_statval.h"

void SrtlibRl::GetDetArr(const double* const sky_arr,
                         const double* const resp_norm_mat_arr,
                         int ndet, int nsky,
                         double* const det_arr)
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

void SrtlibRl::GetDenArr(const double* const sky_arr,
                         const double* const bg_arr,
                         const double* const resp_norm_mat_arr,
                         int ndet, int nsky,
                         double* const den_arr)
{
    SrtlibRl::GetDetArr(sky_arr, resp_norm_mat_arr,
                        ndet, nsky, den_arr);
    daxpy_(ndet, 1.0, const_cast<double*>(bg_arr), 1,
           den_arr, 1);
}

void SrtlibRl::GetYDashArr(const double* const data_arr,
                           const double* const den_arr,
                           int ndet,
                           double* const y_dash_arr)
{
    for(int idet = 0; idet < ndet; idet++){
        y_dash_arr[idet] = data_arr[idet] / den_arr[idet];
    }
}

void SrtlibRl::GetSkyNewArr(const double* const sky_arr,
                            const double* const data_arr,
                            const double* const bg_arr,
                            const double* const resp_norm_mat_arr,
                            int ndet, int nsky,
                            double* const sky_new_arr)
{
    double* den_arr = new double[ndet];
    SrtlibRl::GetDenArr(sky_arr, bg_arr,
                        resp_norm_mat_arr,
                        ndet, nsky,
                        den_arr);
    double* y_dash_arr = new double[ndet];
    SrtlibRl::GetYDashArr(data_arr, den_arr,
                          ndet, y_dash_arr);
    double* tmp_arr = new double[nsky];
    char transa[2];
    strcpy(transa, "T");
    dgemv_(transa, ndet, nsky, 1.0,
           const_cast<double*>(resp_norm_mat_arr), ndet,
           y_dash_arr, 1,
           0.0, tmp_arr, 1);
    for(int isky = 0; isky < nsky; isky ++){
        sky_new_arr[isky] = tmp_arr[isky] * sky_arr[isky];
    }
    delete [] den_arr;
    delete [] y_dash_arr;
    delete [] tmp_arr;
}

void SrtlibRl::Richlucy(
    FILE* const fp_log,
    const double* const sky_init_arr,
    const double* const data_arr,
    const double* const bg_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nsky,
    string outdir, string outfile_head,
    int nem, double tol_em,
    double* const sky_new_arr)
{
    double* sky_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1,
           sky_pre_arr, 1);
    for(int iem = 0; iem < nem; iem ++){
        GetSkyNewArr(sky_pre_arr, data_arr, bg_arr,
                     resp_norm_mat_arr,
                     ndet, nsky, sky_new_arr);
        double helldist  = SrtlibRlStatval::GetHellingerDist(
            sky_pre_arr,
            sky_new_arr,
            nsky);
        if (access( "/tmp/rl_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_stop file is found, then stop.\n");
            break;
        }
        if(iem % 100 == 0){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
        }
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(sky_new_arr), 1,
               sky_pre_arr, 1);
    }
    delete [] sky_pre_arr;
}

// accerelated richardson lucy (SQUAREM)
void SrtlibRl::RichlucyAccSquarem(
    FILE* const fp_log,
    const double* const sky_init_arr,
    const double* const data_arr,
    const double* const bg_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nsky,
    string outdir, string outfile_head,
    int nem, double tol_em,
    double* const sky_new_arr)
{
    double* sky_0_arr  = new double[nsky];
    double* sky_1_arr  = new double[nsky];
    double* sky_2_arr  = new double[nsky];
    double* sky_dash_arr  = new double[nsky];
    double* r_arr  = new double[nsky];
    double* r2_arr  = new double[nsky];
    double* v_arr  = new double[nsky];
    double* sky_0_new_arr  = new double[nsky];
    double* diff_sky_0_arr = new double[nsky];
    
    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_0_arr, 1);    
    
    for(int iem = 0; iem < nem; iem ++){
        GetSkyNewArr(sky_0_arr, data_arr, bg_arr,
                     resp_norm_mat_arr,
                     ndet, nsky, sky_1_arr);
        GetSkyNewArr(sky_1_arr, data_arr, bg_arr,
                     resp_norm_mat_arr,
                     ndet, nsky, sky_2_arr);
        MibBlas::Sub(sky_1_arr, sky_0_arr, nsky, r_arr);
        MibBlas::Sub(sky_2_arr, sky_1_arr, nsky, r2_arr);
        MibBlas::Sub(r2_arr, r_arr, nsky, v_arr);

        double r_norm2 = ddot_(nsky, r_arr, 1, r_arr, 1);
        double v_norm2 = ddot_(nsky, v_arr, 1, v_arr, 1);
        double alpha = -1.0 * sqrt(r_norm2 / v_norm2);
        // printf("alpha = %e\n", alpha);

        int nk = 100;
        double eta = 0.8;
        //double negloglike_pre = 0.0;
        //double epsilon = 1e-20;
        //int ifind = 0;
        int ifind_nonneg = 0;
        for (int ik = 0; ik < nk; ik ++){
            double alpha0 = alpha * pow(eta, ik);
            dcopy_(nsky, sky_0_arr, 1, sky_dash_arr, 1);
            daxpy_(nsky, -2 * alpha0, r_arr, 1, sky_dash_arr, 1);
            daxpy_(nsky, alpha0 * alpha0, v_arr, 1, sky_dash_arr, 1);

            int nneg_tmp = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(sky_dash_arr[isky] < 0.0){
                    nneg_tmp ++;
                }
            }
            if (nneg_tmp > 0){
                continue;
            }
            GetSkyNewArr(sky_dash_arr, data_arr, bg_arr,
                         resp_norm_mat_arr,
                         ndet, nsky, sky_0_new_arr);
            int nneg = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(sky_0_new_arr[isky] < 0.0){
                    nneg ++;
                }
            }
            if (nneg == 0){
                ifind_nonneg = 1;
                break;
                //                if (ifind == 0){
                //                    negloglike_pre = GetNegLogLike(sky_0_new_arr,
                //                                                   data_arr,
                //                                                   resp_norm_mat_arr,
                //                                                   ndet, nsky, epsilon);
                //                    ifind = 1;
                //                } else if (ifind == 1){
                //                    double negloglike = GetNegLogLike(sky_0_new_arr,
                //                                                      data_arr,
                //                                                      resp_norm_mat_arr,
                //                                                      ndet, nsky, epsilon);
                //                    //MiIolib::Printf2(fp_log,
                //                    //"ik = %d, alpha0 = %e, negloglike = %e\n",
                //                    //ik, alpha0, negloglike);
                //                    if(negloglike > negloglike_pre){
                //                        break;
                //                    }
                //                    negloglike_pre = negloglike;
                // }
            }
        }
        if(ifind_nonneg == 0){
            MiIolib::Printf2(fp_log, "warning: iem = %d, ifind_nonneg == 0\n", iem);
        }
        
        // double sum = MibBlas::Sum(sky_0_new_arr, nsky);
        // printf("sum = %e\n", sum);
        double helldist  = SrtlibRlStatval::GetHellingerDist(
            sky_0_arr, sky_0_new_arr, nsky);
        
        //MibBlas::Sub(sky_0_new_arr, sky_0_arr, nsky, diff_sky_0_arr);
        //double diff = sqrt(ddot_(nsky, diff_sky_0_arr, 1,
        //                         diff_sky_0_arr, 1));
        if(iem % 100 == 0){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
        }
        if (helldist < tol_em){
            break;
        }
        dcopy_(nsky, sky_0_new_arr, 1, sky_0_arr, 1);
        if (access( "/tmp/rl_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_stop file is found, then stop.\n");
            break;
        }
       
    }
    dcopy_(nsky, sky_0_new_arr, 1, sky_new_arr, 1);

    delete [] sky_0_arr;
    delete [] sky_1_arr;
    delete [] sky_2_arr;
    delete [] sky_dash_arr;
    delete [] r_arr;
    delete [] r2_arr;
    delete [] v_arr;
    delete [] sky_0_new_arr;
    delete [] diff_sky_0_arr;
}

//// accerelated richardson lucy (Ikeda)
//void SrtlibRl::RichlucyAccIkeda(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const data_arr,
//    const double* const bg_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nsky,
//    string outdir, string outfile_head,
//    int nem, double tol_em,
//    int nph_data, int rand_seed,
//    double* const sky_new_arr)
//{
//    double* sky_arr = new double[nsky];
//    double* sky_dash_arr = new double[nsky];    
//    double* sky_rand_arr = new double[nsky];
//    double* sky_dash_bar_arr = new double[nsky];
//    double* diff_arr = new double[nsky];
//
//    int* evt_rand_arr = new int[nph_data];
//    double* data_rand_arr = new double[ndet];
//
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_arr, 1);
//    
//    for(int iem = 0; iem < nem; iem ++){
//        GetSkyNewArr(sky_arr, data_arr, bg_arr,
//                     resp_norm_mat_arr,
//                     ndet, nsky, sky_dash_arr);
//        SrtlibSim::GenRandomEvtFromProbDist(sky_dash_arr, nsky,
//                                            nph_data, rand_seed,
//                                            sky_rand_arr,
//                                            evt_rand_arr);
//        dscal_(nsky, 1.0 / nph_data, sky_rand_arr, 1);
//        GetDetArr(sky_rand_arr,
//                  resp_norm_mat_arr,
//                  ndet, nsky,
//                  data_rand_arr);
//        GetSkyNewArr(sky_arr, data_rand_arr, bg_arr,
//                     resp_norm_mat_arr,
//                     ndet, nsky, sky_dash_bar_arr);
//        MibBlas::Sub(sky_dash_arr, sky_dash_bar_arr, nsky, diff_arr);
//
//        
//        int nk = 100;
//        double eta = 0.8;
//        //double negloglike_pre = 0.0;
//        //double epsilon = 1e-20;
//        //int ifind = 0;
//        int ifind_nonneg = 0;
//        double alpha = 1.0;
//        for (int ik = 0; ik < nk; ik ++){
//            double alpha0 = alpha * pow(eta, ik);
//            dcopy_(nsky, sky_dash_arr, 1, sky_new_arr, 1);
//            daxpy_(nsky, alpha0, diff_arr, 1, sky_new_arr, 1);
//
//            int nneg = 0;
//            for(int isky = 0; isky < nsky; isky ++){
//                if(sky_new_arr[isky] < 0.0){
//                    nneg ++;
//                }
//            }
//            if (nneg == 0){
//                ifind_nonneg = 1;
//                break;
//                //                if (ifind == 0){
//                //                    negloglike_pre = GetNegLogLike(sky_0_new_arr,
//                //                                                   data_arr,
//                //                                                   resp_norm_mat_arr,
//                //                                                   ndet, nsky, epsilon);
//                //                    ifind = 1;
//                //                } else if (ifind == 1){
//                //                    double negloglike = GetNegLogLike(sky_0_new_arr,
//                //                                                      data_arr,
//                //                                                      resp_norm_mat_arr,
//                //                                                      ndet, nsky, epsilon);
//                //                    //MiIolib::Printf2(fp_log,
//                //                    //"ik = %d, alpha0 = %e, negloglike = %e\n",
//                //                    //ik, alpha0, negloglike);
//                //                    if(negloglike > negloglike_pre){
//                //                        break;
//                //                    }
//                //                    negloglike_pre = negloglike;
//                // }
//            }
//        }
//        if(ifind_nonneg == 0){
//            MiIolib::Printf2(fp_log, "warning: iem = %d, ifind_nonneg == 0\n", iem);
//        }
//        
//        // double sum = MibBlas::Sum(sky_0_new_arr, nsky);
//        // printf("sum = %e\n", sum);
//        double helldist  = SrtlibRlStatval::GetHellingerDist(sky_arr, sky_new_arr, nsky);
//        
//        //MibBlas::Sub(sky_0_new_arr, sky_0_arr, nsky, diff_sky_0_arr);
//        //double diff = sqrt(ddot_(nsky, diff_sky_0_arr, 1,
//        //                         diff_sky_0_arr, 1));
//        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n", iem, helldist);
//        if (helldist < tol_em){
//            break;
//        }
//        
//        dcopy_(nsky, sky_new_arr, 1, sky_arr, 1);
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//       
//    }
//
//    delete [] sky_arr;
//    delete [] sky_dash_arr;
//    delete [] sky_rand_arr;
//    delete [] sky_dash_bar_arr;
//    delete [] diff_arr;
//    delete [] evt_rand_arr;
//    delete [] data_rand_arr;
//}

void SrtlibRl::GetInvVec(const double* const vec_arr, int nelm,
                         double* const inv_arr)
{
    dcopy_(nelm, const_cast<double*>(vec_arr), 1, inv_arr, 1);
    double norm2 = ddot_(nelm, const_cast<double*>(vec_arr), 1,
                         const_cast<double*>(vec_arr), 1);
    dscal_(nelm, 1.0 / norm2, inv_arr, 1);
}


//// accerelated richardson lucy (Kuroda)
//void SrtlibRl::RichlucyAccKuroda(
//    FILE* const fp_log,
//    const double* const sky_init_arr,
//    const double* const data_arr,
//    const double* const bg_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nsky,
//    string outdir, string outfile_head,
//    int nem, double tol_em,
//    int k_restart, double delta_restart,
//    double* const sky_new_arr)
//{
//    if (delta_restart <= tol_em){
//        printf("delta_restart(=%e) <= tol_em(=%e), then abort.\n",
//               delta_restart, tol_em);
//        abort();
//    }
//    double* sky_0_arr  = new double[nsky];
//    double* sky_1_arr  = new double[nsky];
//    double* sky_2_arr  = new double[nsky];
//    double* sky_dot_old_arr = new double[nsky];
//    double* sky_dot_new_arr = new double[nsky];
//
//    dcopy_(nsky, const_cast<double*>(sky_init_arr), 1, sky_0_arr, 1);
//    GetSkyNewArr(sky_0_arr, data_arr, bg_arr,
//                 resp_norm_mat_arr,
//                 ndet, nsky, sky_1_arr);
//    dcopy_(nsky, sky_1_arr, 1, sky_dot_old_arr, 1);
//    
//    for(int iem = 0; iem < nem; iem ++){
//        GetSkyNewArr(sky_1_arr, data_arr, bg_arr,
//                     resp_norm_mat_arr,
//                     ndet, nsky, sky_2_arr);
//
//        // epsilon-acceleration step
//        double* delta_sky_0_arr = new double[nsky];
//        double* delta_sky_1_arr = new double[nsky];
//        MibBlas::Sub(sky_1_arr, sky_0_arr, nsky, delta_sky_0_arr);
//        MibBlas::Sub(sky_2_arr, sky_1_arr, nsky, delta_sky_1_arr);
//        double* delta_sky_0_inv_arr = new double[nsky];
//        double* delta_sky_1_inv_arr = new double[nsky];        
//        GetInvVec(delta_sky_0_arr, nsky, delta_sky_0_inv_arr);
//        GetInvVec(delta_sky_1_arr, nsky, delta_sky_1_inv_arr);
//
//        double* delta_delta_sky_01_inv_arr = new double[nsky];
//        MibBlas::Sub(delta_sky_1_inv_arr, delta_sky_0_inv_arr, nsky,
//                     delta_delta_sky_01_inv_arr);
//        double* delta_delta_sky_01_inv_inv_arr = new double[nsky];
//        GetInvVec(delta_delta_sky_01_inv_arr, nsky,
//                  delta_delta_sky_01_inv_inv_arr);
//        MibBlas::Add(sky_1_arr, delta_delta_sky_01_inv_inv_arr, nsky,
//                     sky_dot_new_arr);
//        delete [] delta_delta_sky_01_inv_arr;
//        delete [] delta_delta_sky_01_inv_inv_arr;
//        delete [] delta_sky_0_inv_arr;
//        delete [] delta_sky_1_inv_arr;
//        delete [] delta_sky_0_arr;
//        delete [] delta_sky_1_arr;
//
//        // re-starting procedure
//        double* delta_sky_dot_arr = new double[nsky];
//        MibBlas::Sub(sky_dot_new_arr, sky_dot_old_arr, nsky,
//                     delta_sky_dot_arr);
//        double diff2 = ddot_(nsky, delta_sky_dot_arr, 1,
//                             delta_sky_dot_arr, 1);
//        delete [] delta_sky_dot_arr;
//        MiIolib::Printf2(fp_log, "iem = %d, diff2 = %.2e\n",
//                         iem, diff2);
//        if (diff2 < delta_restart){
//            if (diff2 < tol_em){
//                MiIolib::Printf2(fp_log, "iem = %d, diff2 = %.2e\n",
//                                 iem, diff2);
//                break;
//            }
//            double* sky_tmp_arr = new double[nsky];
//            GetSkyNewArr(sky_dot_new_arr, data_arr, bg_arr,
//                         resp_norm_mat_arr,
//                         ndet, nsky, sky_tmp_arr);
//            double lval_tmp = -1.0 * SrtlibRlStatval::GetNegLogLike(sky_tmp_arr,
//                                                                    data_arr,
//                                                                    resp_norm_mat_arr,
//                                                                    ndet, nsky, 1.0e-10);
//            double lval_2 = -1.0 * SrtlibRlStatval::GetNegLogLike(sky_2_arr,
//                                                                  data_arr,
//                                                                  resp_norm_mat_arr,
//                                                                  ndet, nsky, 1.0e-10);
//            if (lval_tmp > lval_2){
//                dcopy_(nsky, sky_tmp_arr, 1, sky_2_arr, 1);
//                dcopy_(nsky, sky_dot_new_arr, 1, sky_1_arr, 1);
//                delta_restart *= pow(10.0, -1.0 * k_restart);
//                printf("delta_restart = %e\n", delta_restart);
//            }
//            delete [] sky_tmp_arr;
//        }
//        dcopy_(nsky, sky_dot_new_arr, 1, sky_dot_old_arr, 1);
//        dcopy_(nsky, sky_1_arr, 1, sky_0_arr, 1);
//        dcopy_(nsky, sky_2_arr, 1, sky_1_arr, 1);
//       
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//       
//    }
//    dcopy_(nsky, sky_dot_new_arr, 1, sky_new_arr, 1);
//    delete [] sky_0_arr;
//    delete [] sky_1_arr;
//    delete [] sky_2_arr;
//    delete [] sky_dot_old_arr;
//    delete [] sky_dot_new_arr;
//}
//

