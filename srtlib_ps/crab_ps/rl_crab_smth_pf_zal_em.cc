#include "rl_crab.h"
#include "rl_statval_crab.h"
#include "rl_crab_smth_pf_em.h"
#include "rl_crab_smth_pf_pm.h"

void SrtlibRlCrabSmthPfEm::RichlucyCrabSmthPf(
    FILE* const fp_log,
    const double* const rho_init_arr,
    const double* const nu_init_arr,
    const double* const* const data_arr,
    const double* const nu_target_arr,
    const double* const phase_arr,
    const double* const det_0_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nskyx, int nskyy, int nphase,
    double mu, double gamma,
    string outdir,
    string outfile_head,
    int nem, double tol_em,
    int npm, double tol_pm,
    int nnewton, double tol_newton,
    double* const rho_new_arr,
    double* const nu_new_arr)
{
    int nsky = nskyx * nskyy;
    double* rho_pre_arr = new double[nsky];
    double* nu_pre_arr = new double[nphase];
    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_pre_arr, 1);
    dcopy_(nphase, const_cast<double*>(nu_init_arr), 1, nu_pre_arr, 1);
    for(int iem = 0; iem < nem; iem ++){
        double* mval_arr = new double[nsky];
        double* nval_arr = new double[nphase];
        SrtlibRlCrab::GetRhoNuNewNumArr(rho_pre_arr, nu_pre_arr,
                                        data_arr,
                                        phase_arr, det_0_arr,
                                        resp_norm_mat_arr, 
                                        ndet, nsky, nphase,
                                        mval_arr, nval_arr);
        double helldist_pm = 0.0;
        int flag_converge_pm = 0;
        SrtlibRlCrabSmthPfPm::GetRhoNu_ByPm(
            fp_log, rho_pre_arr, nu_pre_arr,
            mval_arr, nval_arr, nu_target_arr,
            nskyx, nskyy, nphase, mu, gamma,
            npm, tol_pm, nnewton, tol_newton,
            rho_new_arr, nu_new_arr,
            &helldist_pm,
            &flag_converge_pm);
        if (flag_converge_pm == 0){
            MiIolib::Printf2(
                fp_log,
                "iem = %d: pm: not converged: helldist_pm = %.2e\n",
                iem,
                helldist_pm);
        }
        delete [] mval_arr;
        delete [] nval_arr;
        double helldist  = SrtlibRlStatvalCrab::GetHellingerDist(
            rho_pre_arr, nu_pre_arr,
            rho_new_arr, nu_new_arr,
            nsky, nphase);
        if (access( "/tmp/rl_crab_smth_pf_em_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_crab_smth_pf_em_stop file is found, then stop.\n");
            break;
        }
        if (helldist < tol_em){
            printf("iem = %d, helldist = %e\n",
                   iem, helldist);
            break;
        }
        dcopy_(nsky, rho_new_arr, 1, rho_pre_arr, 1);
        dcopy_(nphase, nu_new_arr, 1, nu_pre_arr, 1);
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e\n",
                         iem, helldist);
    }
    delete [] rho_pre_arr;
    delete [] nu_pre_arr;
}

// Accerelated by Zhou-Alexander-Lange and Nesterov.
// The former is H.Zhou, D.Alexander, K.Lange,
// "A quasi-Newton acceleration for high-dimensional
// optimization algorithms", Stat Comput (2011) 21, 261.
// Case: q = 1
void SrtlibRlCrabSmthPfEm::RichlucyCrabSmthPfAcc(
    FILE* const fp_log,
    const double* const rho_init_arr,
    const double* const nu_init_arr,
    const double* const* const data_arr,
    const double* const nu_target_arr,
    const double* const phase_arr,
    const double* const det_0_arr,
    const double* const resp_norm_mat_arr,
    int ndet, int nskyx, int nskyy, int nphase,
    double mu, double gamma,
    string outdir,
    string outfile_head,
    int nem, double tol_em,
    int npm, double tol_pm,
    int nnewton, double tol_newton,
    double* const rho_new_arr,
    double* const nu_new_arr)
{
    int nsky = nskyx * nskyy;

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
        double* mval_arr = new double[nsky];
        double* nval_arr = new double[nphase];
        double helldist_pm1 = 0.0;
        int flag_converge_pm1 = 0;
        SrtlibRlCrab::GetRhoNuNewNumArr(rho_0_arr, nu_0_arr,
                                        data_arr,
                                        phase_arr, det_0_arr,
                                        resp_norm_mat_arr, 
                                        ndet, nsky, nphase,
                                        mval_arr, nval_arr);        
        SrtlibRlCrabSmthPfPm::GetRhoNu_ByPm(
            fp_log, rho_0_arr, nu_0_arr,
            mval_arr, nval_arr, nu_target_arr,
            nskyx, nskyy, nphase, mu, gamma,
            npm, tol_pm, nnewton, tol_newton,
            rho_1_arr, nu_1_arr,
            &helldist_pm1,
            &flag_converge_pm1);
        if (flag_converge_pm1 == 0){
            MiIolib::Printf2(
                fp_log,
                "iem = %d: pm1: not converged: helldist_pm1 = %.2e\n",
                iem,
                helldist_pm1);
        }
        double helldist_pm2 = 0.0;
        int flag_converge_pm2 = 0;
        SrtlibRlCrab::GetRhoNuNewNumArr(rho_1_arr, nu_1_arr,
                                        data_arr,
                                        phase_arr, det_0_arr,
                                        resp_norm_mat_arr, 
                                        ndet, nsky, nphase,
                                        mval_arr, nval_arr);        
        SrtlibRlCrabSmthPfPm::GetRhoNu_ByPm(
            fp_log, rho_1_arr, nu_1_arr,
            mval_arr, nval_arr, nu_target_arr,
            nskyx, nskyy, nphase, mu, gamma,
            npm, tol_pm, nnewton, tol_newton,
            rho_2_arr, nu_2_arr,
            &helldist_pm2,
            &flag_converge_pm2);
        if (flag_converge_pm2 == 0){
            MiIolib::Printf2(
                fp_log,
                "iem = %d: pm2: not converged: helldist_pm2 = %.2e\n",
                iem,
                helldist_pm2);
        }

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
        MiIolib::Printf2(fp_log, "cval = %e\n", cval);

        if (cval < 0.0){
            // usual update
            dcopy_(nsky, rho_1_arr, 1, rho_0_new_arr, 1);
            dcopy_(nphase, nu_1_arr, 1, nu_0_new_arr, 1);
        } else {
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
                    MiIolib::Printf2(fp_log, "cval0 = %e\n", cval0);
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
        MiIolib::Printf2(
            fp_log, "iem = %d, helldist = %e\n",
            iem, helldist);
        if (helldist < tol_em){
            MiIolib::Printf2(
                fp_log, "iem = %d, helldist = %e\n",
                iem, helldist);
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

