#include<unistd.h>
#include "mir_math.h"
#include "rl_crab.h"
#include "rl_statval_crab.h"

void SrtlibRlCrab::GetRhoNuNewNumArr(
    const double* const rho_arr,
    const double* const nu_arr,
    const double* const* const data_arr,
    const double* const phase_arr,
    const double* const det_0_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nsky, int nphase,
    double* const rho_new_arr,
    double* const nu_new_arr)
{
    double* coeff_rho_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky++){
        coeff_rho_arr[isky] = 0.0;
    }
    double* den_arr = new double[ndet];
    double* div_arr = new double[ndet];    
    for(int iphase = 0; iphase < nphase; iphase++){
        SrtlibRl::GetDetArr(rho_arr, resp_norm_mat_arr,
                            ndet, nsky, den_arr);
        dscal_(ndet, phase_arr[iphase], den_arr, 1);
        daxpy_(ndet, nu_arr[iphase], const_cast<double*>(det_0_arr), 1,
               den_arr, 1);
        for(int idet = 0; idet < ndet; idet++){
            div_arr[idet] = data_arr[iphase][idet] / den_arr[idet];
        }

        char transa[1];
        strcpy(transa, "T");
        dgemv_(transa, ndet, nsky, phase_arr[iphase],
               const_cast<double*>(resp_norm_mat_arr), ndet,
               div_arr, 1,
               1.0, coeff_rho_arr, 1);

        // nu_new
        nu_new_arr[iphase] = ddot_(ndet, div_arr, 1,
                                   const_cast<double*>(det_0_arr), 1)
            * nu_arr[iphase];
    }
    delete [] den_arr;
    delete [] div_arr;
    
    // rho_new
    MibBlas::ElmWiseMul(nsky, 1.0,
                        coeff_rho_arr, rho_arr,
                        rho_new_arr);
    delete [] coeff_rho_arr;
    
}

void SrtlibRlCrab::GetRhoNuNewArr(const double* const rho_arr,
                                  const double* const nu_arr,
                                  const double* const* const data_arr,
                                  const double* const phase_arr,
                                  const double* const det_0_arr,
                                  const double* const resp_norm_mat_arr,
                                  int ndet, int nsky, int nphase,
                                  double* const rho_new_arr,
                                  double* const nu_new_arr)
{
    // denominator
    double nevt = 0;
    for(int iphase = 0; iphase < nphase; iphase++){
        nevt += MibBlas::Sum(data_arr[iphase], ndet);
    }
    // numerator
    SrtlibRlCrab::GetRhoNuNewNumArr(rho_arr,
                                    nu_arr,
                                    data_arr,
                                    phase_arr,
                                    det_0_arr,
                                    resp_norm_mat_arr,
                                    ndet, nsky, nphase,
                                    rho_new_arr,
                                    nu_new_arr);
    dscal_(nsky, 1.0/nevt, rho_new_arr, 1);
    dscal_(nphase, 1.0/nevt, nu_new_arr, 1);
}


void SrtlibRlCrab::RichlucyCrab(FILE* const fp_log,
                                const double* const rho_init_arr,
                                const double* const nu_init_arr,
                                const double* const* const data_arr,
                                const double* const phase_arr,
                                const double* const det_0_arr,
                                const double* const resp_norm_mat_arr,
                                int ndet, int nsky, int nphase,
                                string outdir, string outfile_head,
                                int nem, double tol_em,
                                double* const rho_new_arr,
                                double* const nu_new_arr)
{
    double* rho_pre_arr = new double[nsky];
    double* nu_pre_arr = new double[nphase];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    dcopy_(nphase, const_cast<double*>(nu_init_arr), 1, nu_pre_arr, 1);
    for(int iem = 0; iem < nem; iem ++){
        SrtlibRlCrab::GetRhoNuNewArr(rho_pre_arr, nu_pre_arr,
                                     data_arr, phase_arr, det_0_arr,
                                     resp_norm_mat_arr,
                                     ndet, nsky, nphase,
                                     rho_new_arr, nu_new_arr);
        double helldist  = SrtlibRlStatvalCrab::GetHellingerDist(
            rho_pre_arr, nu_pre_arr,
            rho_new_arr, nu_new_arr,
            nsky, nphase);
        if (access( "/tmp/rl_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_stop file is found, then stop.\n");
            break;
        }
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                         iem, helldist);
        
        if (helldist < tol_em){
            MiIolib::Printf2(fp_log, "iem = %d, helldist = %.2e\n",
                             iem, helldist);
            break;
        }
        dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, rho_pre_arr, 1);
        dcopy_(nphase, const_cast<double*>(nu_new_arr), 1, nu_pre_arr, 1);
    }
    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
}


//
// accerelated richardson lucy by Zhou-Alexander-Lange:
// H.Zhou, D.Alexander, K.Lange,
// "A quasi-Newton acceleration for high-dimensional
// optimization algorithms", Stat Comput (2011) 21, 261.
// Case: q = 1
void SrtlibRlCrab::RichlucyCrabAccZALq1(
    FILE* const fp_log,
    const double* const rho_init_arr,
    const double* const nu_init_arr,
    const double* const* const data_arr,
    const double* const phase_arr,
    const double* const det_0_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nsky, int nphase,
    string outdir, string outfile_head,
    int nem, double tol_em,
    double* const rho_new_arr,
    double* const nu_new_arr)
{
    double* rho_0_arr  = new double[nsky];
    double* rho_1_arr  = new double[nsky];
    double* rho_2_arr  = new double[nsky];
    double* u_rho_arr  = new double[nsky];
    double* v_rho_arr  = new double[nsky];
    double* diff_rho_arr  = new double[nsky];
    double* rho_0_new_arr  = new double[nsky];

    double* nu_0_arr  = new double[nphase];
    double* nu_1_arr  = new double[nphase];
    double* nu_2_arr  = new double[nphase];
    double* u_nu_arr  = new double[nphase];
    double* v_nu_arr  = new double[nphase];
    double* diff_nu_arr  = new double[nphase];
    double* nu_0_new_arr  = new double[nphase];
    
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_0_arr, 1);
    dcopy_(nphase, const_cast<double*>(nu_init_arr), 1, nu_0_arr, 1);
    
    for(int iem = 0; iem < nem; iem ++){
        SrtlibRlCrab::GetRhoNuNewArr(rho_0_arr, nu_0_arr, data_arr,
                                     phase_arr, det_0_arr,
                                     resp_norm_mat_arr,
                                     ndet, nsky, nphase,
                                     rho_1_arr, nu_1_arr);
        SrtlibRlCrab::GetRhoNuNewArr(rho_1_arr, nu_1_arr, data_arr,
                                     phase_arr, det_0_arr,
                                     resp_norm_mat_arr,
                                     ndet, nsky, nphase,
                                     rho_2_arr, nu_2_arr);
        MibBlas::Sub(rho_1_arr, rho_0_arr, nsky, u_rho_arr);
        MibBlas::Sub(nu_1_arr, nu_0_arr, nphase, u_nu_arr);
        MibBlas::Sub(rho_2_arr, rho_1_arr, nsky, v_rho_arr);
        MibBlas::Sub(nu_2_arr, nu_1_arr, nphase, v_nu_arr);
        MibBlas::Sub(v_rho_arr, u_rho_arr, nsky, diff_rho_arr);
        MibBlas::Sub(v_nu_arr, u_nu_arr, nphase, diff_nu_arr);

        double num = ddot_(nsky, u_rho_arr, 1, u_rho_arr, 1)
            + ddot_(nphase, u_nu_arr, 1, u_nu_arr, 1);
        double den = ddot_(nsky, u_rho_arr, 1, diff_rho_arr, 1)
            + ddot_(nphase, u_nu_arr, 1, diff_nu_arr, 1);
        double cval = -1.0 * num / den;
        printf("cval = %e\n", cval);

        if (cval < 0.0){
            // usual update
            dcopy_(nsky, rho_1_arr, 1, rho_0_new_arr, 1);
            dcopy_(nphase, nu_1_arr, 1, nu_0_new_arr, 1);
        } else{
            int nk = 10000;
            double eta = 0.8;
            int flag_find = 0;
            for (int ik = 0; ik < nk; ik ++){
                double cval0 = cval * pow(eta, ik);
                dcopy_(nsky, rho_1_arr, 1, rho_0_new_arr, 1);
                dcopy_(nphase, nu_1_arr, 1, nu_0_new_arr, 1);
                dscal_(nsky, (1.0 - cval0), rho_0_new_arr, 1);
                dscal_(nphase, (1.0 - cval0), nu_0_new_arr, 1);
                daxpy_(nsky, cval0, rho_2_arr, 1, rho_0_new_arr, 1);
                daxpy_(nphase, cval0, nu_2_arr, 1, nu_0_new_arr, 1);
                
                int nneg_tmp = 0;
                for(int isky = 0; isky < nsky; isky ++){
                    if(rho_0_new_arr[isky] < 0.0){
                        nneg_tmp ++;
                    }
                }
                for(int iphase = 0; iphase < nphase; iphase ++){
                    if(nu_0_new_arr[iphase] < 0.0){
                        nneg_tmp ++;
                    }
                }
                if (nneg_tmp > 0){
                    continue;
                } else{
                    flag_find = 1;
                    break;
                }
            }
            if(flag_find == 0){
                // usual update
                dcopy_(nsky, rho_2_arr, 1, rho_0_new_arr, 1);
                dcopy_(nphase, nu_2_arr, 1, nu_0_new_arr, 1);
            }
        }

        double helldist  = SrtlibRlStatvalCrab::GetHellingerDist(
            rho_0_arr, nu_0_arr,
            rho_0_new_arr, nu_0_new_arr,
            nsky, nphase);
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n", iem, helldist);
        if (helldist < tol_em){
            break;
        }
        dcopy_(nsky, rho_0_new_arr, 1, rho_0_arr, 1);
        dcopy_(nphase, nu_0_new_arr, 1, nu_0_arr, 1);
        if (access( "/tmp/rl_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_stop file is found, then stop.\n");
            break;
        }
       
    }
    dcopy_(nsky, rho_0_new_arr, 1, rho_new_arr, 1);
    dcopy_(nphase, nu_0_new_arr, 1, nu_new_arr, 1);

    delete [] rho_0_arr;
    delete [] rho_1_arr;
    delete [] rho_2_arr;
    delete [] u_rho_arr;
    delete [] v_rho_arr;
    delete [] diff_rho_arr;
    delete [] rho_0_new_arr;
    
    delete [] nu_0_arr;
    delete [] nu_1_arr;
    delete [] nu_2_arr;
    delete [] u_nu_arr;
    delete [] v_nu_arr;
    delete [] diff_nu_arr;
    delete [] nu_0_new_arr;
}


//// accerelated richardson lucy (SQUAREM)
//void SrtlibRlCrab::RichlucyCrabAccSquarem(
//    FILE* const fp_log,
//    const double* const rho_init_arr,
//    const double* const nu_init_arr,
//    const double* const* const data_arr,
//    const double* const phase_arr,
//    const double* const det_0_arr,
//    const double* const resp_norm_mat_arr,
//    int ndet, int nsky, int nphase,
//    string outdir, string outfile_head,
//    int nem, double tol_em,
//    double* const rho_new_arr,
//    double* const nu_new_arr)
//{
//    double* rho_0_arr  = new double[nsky];
//    double* rho_1_arr  = new double[nsky];
//    double* rho_2_arr  = new double[nsky];
//    double* rho_dash_arr  = new double[nsky];
//    double* r_rho_arr  = new double[nsky];
//    double* r2_rho_arr  = new double[nsky];
//    double* v_rho_arr  = new double[nsky];
//    double* rho_0_new_arr  = new double[nsky];
//    double* diff_rho_0_arr = new double[nsky];
//
//    double* nu_0_arr  = new double[nphase];
//    double* nu_1_arr  = new double[nphase];
//    double* nu_2_arr  = new double[nphase];
//    double* nu_dash_arr  = new double[nphase];
//    double* r_nu_arr  = new double[nphase];
//    double* r2_nu_arr  = new double[nphase];
//    double* v_nu_arr  = new double[nphase];
//    double* nu_0_new_arr  = new double[nphase];
//    double* diff_nu_0_arr = new double[nphase];
//    
//    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_0_arr, 1);
//    dcopy_(nphase, const_cast<double*>(nu_init_arr), 1, nu_0_arr, 1);
//    
//    for(int iem = 0; iem < nem; iem ++){
//        SrtlibRlCrab::GetRhoNuNewArr(rho_0_arr, nu_0_arr, data_arr,
//                                     phase_arr, det_0_arr,
//                                     resp_norm_mat_arr,
//                                     ndet, nsky, nphase,
//                                     rho_1_arr, nu_1_arr);
//        SrtlibRlCrab::GetRhoNuNewArr(rho_1_arr, nu_1_arr, data_arr,
//                                     phase_arr, det_0_arr,
//                                     resp_norm_mat_arr,
//                                     ndet, nsky, nphase,
//                                     rho_2_arr, nu_2_arr);
//        MibBlas::Sub(rho_1_arr, rho_0_arr, nsky, r_rho_arr);
//        MibBlas::Sub(nu_1_arr, nu_0_arr, nphase, r_nu_arr);
//        MibBlas::Sub(rho_2_arr, rho_1_arr, nsky, r2_rho_arr);
//        MibBlas::Sub(nu_2_arr, nu_1_arr, nphase, r2_nu_arr);
//        MibBlas::Sub(r2_rho_arr, r_rho_arr, nsky, v_rho_arr);
//        MibBlas::Sub(r2_nu_arr, r_nu_arr, nphase, v_nu_arr);        
//
//        double r_norm2 = ddot_(nsky, r_rho_arr, 1, r_rho_arr, 1)
//            + ddot_(nphase, r_nu_arr, 1, r_nu_arr, 1);
//        double v_norm2 = ddot_(nsky, v_rho_arr, 1, v_rho_arr, 1)
//            + ddot_(nphase, v_nu_arr, 1, v_nu_arr, 1);
//        double alpha = -1.0 * sqrt(r_norm2 / v_norm2);
//        // printf("alpha = %e\n", alpha);
//
//        // Zhou 2011 
//        double d_norm2 = ddot_(nsky, r_rho_arr, 1, v_rho_arr, 1)
//            + ddot_(nphase, r_nu_arr, 1, v_nu_arr, 1);
//        double beta = -1.0 * r_norm2 / d_norm2;
//        printf("beta = %e\n", beta);
//
//        if (beta < 0.0){
//            // usual update
//            dcopy_(nsky, rho_1_arr, 1, rho_0_new_arr, 1);
//            dcopy_(nphase, nu_1_arr, 1, nu_0_new_arr, 1);
//        } else{
//            int nk = 10000;
//            double eta = 0.8;
//            int flag_find = 0;
//            for (int ik = 0; ik < nk; ik ++){
//                double beta0 = beta * pow(eta, ik);
//                dcopy_(nsky, rho_1_arr, 1, rho_0_new_arr, 1);
//                dcopy_(nphase, nu_1_arr, 1, nu_0_new_arr, 1);
//                dscal_(nsky, (1.0 - beta0), rho_0_new_arr, 1);
//                dscal_(nphase, (1.0 - beta0), nu_0_new_arr, 1);
//                daxpy_(nsky, beta0, rho_2_arr, 1, rho_0_new_arr, 1);
//                daxpy_(nphase, beta0, nu_2_arr, 1, nu_0_new_arr, 1);
//                
//                int nneg_tmp = 0;
//                for(int isky = 0; isky < nsky; isky ++){
//                    if(rho_0_new_arr[isky] < 0.0){
//                        nneg_tmp ++;
//                    }
//                }
//                for(int iphase = 0; iphase < nphase; iphase ++){
//                    if(nu_0_new_arr[iphase] < 0.0){
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
//                dcopy_(nsky, rho_1_arr, 1, rho_0_new_arr, 1);
//                dcopy_(nphase, nu_1_arr, 1, nu_0_new_arr, 1);
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
////            dcopy_(nsky, rho_0_arr, 1, rho_dash_arr, 1);
////            dcopy_(nphase, nu_0_arr, 1, nu_dash_arr, 1);
////            daxpy_(nsky, -2 * alpha0, r_rho_arr, 1, rho_dash_arr, 1);
////            daxpy_(nphase, -2 * alpha0, r_nu_arr, 1, nu_dash_arr, 1);
////            daxpy_(nsky, alpha0 * alpha0, v_rho_arr, 1, rho_dash_arr, 1);
////            daxpy_(nphase, alpha0 * alpha0, v_nu_arr, 1, rho_dash_arr, 1);
////
////            int nneg_tmp = 0;
////            for(int isky = 0; isky < nsky; isky ++){
////                if(rho_dash_arr[isky] < 0.0){
////                    nneg_tmp ++;
////                }
////            }
////            for(int iphase = 0; iphase < nphase; iphase ++){
////                if(nu_dash_arr[iphase] < 0.0){
////                    nneg_tmp ++;
////                }
////            }
////            if (nneg_tmp > 0){
////                continue;
////            }
////            SrtlibRlCrab::GetRhoNuNewArr(rho_dash_arr, nu_dash_arr, data_arr,
////                                         phase_arr, det_0_arr,
////                                         resp_norm_mat_arr,
////                                         ndet, nsky, nphase,
////                                         rho_0_new_arr, nu_0_new_arr);
////            int nneg = 0;
////            for(int isky = 0; isky < nsky; isky ++){
////                if(rho_0_new_arr[isky] < 0.0){
////                    nneg ++;
////                }
////            }
////            for(int iphase = 0; iphase < nphase; iphase ++){
////                if(nu_0_new_arr[iphase] < 0.0){
////                    nneg ++;
////                }
////            }
////            if (nneg == 0){
////                ifind_nonneg = 1;
////                break;
////                //                if (ifind == 0){
////                //                    negloglike_pre = GetNegLogLike(rho_0_new_arr,
////                //                                                   data_arr,
////                //                                                   resp_norm_mat_arr,
////                //                                                   ndet, nsky, epsilon);
////                //                    ifind = 1;
////                //                } else if (ifind == 1){
////                //                    double negloglike = GetNegLogLike(rho_0_new_arr,
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
//        //    dcopy_(nsky, rho_1_arr, 1, rho_0_new_arr, 1);
//        //    dcopy_(nphase, nu_1_arr, 1, nu_0_new_arr, 1);
//        //}
//        
//        // double sum = MibBlas::Sum(rho_0_new_arr, nsky);
//        // printf("sum = %e\n", sum);
//        double helldist  = SrtlibRlStatvalCrab::GetHellingerDist(
//            rho_0_arr, nu_0_arr,
//            rho_0_new_arr, nu_0_new_arr,
//            nsky, nphase);
//        
//        //MibBlas::Sub(rho_0_new_arr, rho_0_arr, nsky, diff_rho_0_arr);
//        //double diff = sqrt(ddot_(nsky, diff_rho_0_arr, 1,
//        //                         diff_rho_0_arr, 1));
//        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n", iem, helldist);
//        if (helldist < tol_em){
//            break;
//        }
//        dcopy_(nsky, rho_0_new_arr, 1, rho_0_arr, 1);
//        dcopy_(nphase, nu_0_new_arr, 1, nu_0_arr, 1);
//        if (access( "/tmp/rl_stop", R_OK ) != -1){
//            MiIolib::Printf2(
//                fp_log,
//                "/tmp/rl_stop file is found, then stop.\n");
//            break;
//        }
//       
//    }
//    dcopy_(nsky, rho_0_new_arr, 1, rho_new_arr, 1);
//    dcopy_(nphase, nu_0_new_arr, 1, nu_new_arr, 1);
//
//    delete [] rho_0_arr;
//    delete [] rho_1_arr;
//    delete [] rho_2_arr;
//    delete [] rho_dash_arr;
//    delete [] r_rho_arr;
//    delete [] r2_rho_arr;
//    delete [] v_rho_arr;
//    delete [] rho_0_new_arr;
//    delete [] diff_rho_0_arr;
//    
//    delete [] nu_0_arr;
//    delete [] nu_1_arr;
//    delete [] nu_2_arr;
//    delete [] nu_dash_arr;
//    delete [] r_nu_arr;
//    delete [] r2_nu_arr;
//    delete [] v_nu_arr;
//    delete [] nu_0_new_arr;
//    delete [] diff_nu_0_arr;
//}
//
//
//
