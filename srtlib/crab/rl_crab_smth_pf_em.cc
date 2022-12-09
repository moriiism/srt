#include "rl.h"
#include "rl_statval.h"
#include "rl_crab_smth_pf_em.h"
#include "rl_crab_smth_pf_pm.h"

void SrtlibRlCrabSmthPfEm::RichlucyCrabSmthPf(
    FILE* const fp_log,
    const double* const rho_init_arr,
    const double* const nu_init_arr,
    const double* const* const data_arr,
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
            fp_log,
            rho_pre_arr, nu_pre_arr,
            mval_arr, nval_arr,
            nskyx, nskyy, nphase,
            mu, gamma,
            npm, tol_pm,
            nnewton, tol_newton,
            rho_new_arr,
            nu_new_arr,
            &helldist_pm,
            &flag_converge_pm);
        if (flag_converge_pm == 0){
            MiIolib::Printf2(fp_log,
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
        if (access( "/tmp/rl_bg2_smth_em_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_bg2_smth_em_stop file is found, then stop.\n");
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

// accerelation by SQUAREM
void SrtlibRlCrabSmthPfEm::RichlucyCrabSmthPf_Squarem(
    FILE* const fp_log,
    const double* const rho_init_arr,
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

    double* rho_0_arr  = new double[nsky];
    double* rho_1_arr  = new double[nsky];
    double* rho_2_arr  = new double[nsky];
    double* rho_dash_arr  = new double[nsky];
    double* r_rho_arr  = new double[nsky];
    double* r2_rho_arr  = new double[nsky];
    double* v_rho_arr  = new double[nsky];
    double* rho_0_new_arr  = new double[nsky];

    double nu_0 = 0.0;
    double nu_1 = 0.0;
    double nu_2 = 0.0;
    double nu_dash = 0.0;
    double r_nu = 0.0;
    double r2_nu = 0.0;
    double v_nu = 0.0;
    double nu_0_new = 0.0;

    dcopy_(nsky, const_cast<double*>(rho_init_arr), 1, rho_0_arr, 1);
    nu_0 = nu_init;
    for(int iem = 0; iem < nem; iem ++){
        double* mval_arr = new double[nsky];
        double nval = 0.0;

        double helldist_pm1 = 0.0;
        int flag_converge_pm1 = 0;
        GetMArrNval(rho_0_arr, nu_0,
                    data_arr, bg_arr,
                    resp_norm_mat_arr, 
                    ndet, nsky, mval_arr, &nval);
        SrtlibRlBg2SmthPm::GetRhoNu_ByPm_Nesterov(
            fp_log,
            rho_0_arr, nu_0,
            mval_arr, nval,
            nskyx, nskyy,
            mu,
            npm, tol_pm,
            nnewton, tol_newton,
            rho_1_arr,
            &nu_1,
            &helldist_pm1,
            &flag_converge_pm1);
        if (flag_converge_pm1 == 0){
            MiIolib::Printf2(fp_log,
                             "iem = %d: pm1: not converged: helldist_pm1 = %.2e\n",
                             iem,
                             helldist_pm1);
        }
        double helldist_pm2 = 0.0;
        int flag_converge_pm2 = 0;
        GetMArrNval(rho_1_arr, nu_1,
                    data_arr, bg_arr,
                    resp_norm_mat_arr, 
                    ndet, nsky, mval_arr, &nval);
        SrtlibRlBg2SmthPm::GetRhoNu_ByPm_Nesterov(
            fp_log,
            rho_1_arr, nu_1,
            mval_arr, nval,
            nskyx, nskyy,
            mu,
            npm, tol_pm,
            nnewton, tol_newton,
            rho_2_arr,
            &nu_2,
            &helldist_pm2,
            &flag_converge_pm2);
        if (flag_converge_pm2 == 0){
            MiIolib::Printf2(fp_log,
                             "iem = %d: pm2: not converged: helldist_pm2 = %.2e\n",
                             iem,
                             helldist_pm2);
        }

        MibBlas::Sub(rho_1_arr, rho_0_arr, nsky, r_rho_arr);
        MibBlas::Sub(rho_2_arr, rho_1_arr, nsky, r2_rho_arr);
        MibBlas::Sub(r2_rho_arr, r_rho_arr, nsky, v_rho_arr);

        r_nu = nu_1 - nu_0;
        r2_nu = nu_2 - nu_1;
        v_nu = r2_nu - r_nu;

        double r_norm2 = ddot_(nsky, r_rho_arr, 1, r_rho_arr, 1)
            + r_nu * r_nu;
        double v_norm2 = ddot_(nsky, v_rho_arr, 1, v_rho_arr, 1)
            + v_nu * v_nu;
        double alpha = -1.0 * sqrt(r_norm2 / v_norm2);

        int nk = 100;
        double eta = 0.8;
        int ifind_nonneg = 0;
        for (int ik = 0; ik < nk; ik ++){
            double alpha0 = alpha * pow(eta, ik);
            dcopy_(nsky, rho_0_arr, 1, rho_dash_arr, 1);
            nu_dash = nu_0;
            daxpy_(nsky, -2 * alpha0, r_rho_arr, 1, rho_dash_arr, 1);
            nu_dash += -2 * alpha0 * r_nu;
            daxpy_(nsky, alpha0 * alpha0, v_rho_arr, 1, rho_dash_arr, 1);
            nu_dash += alpha0 * alpha0 * v_nu;

            int nneg_tmp = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_dash_arr[isky] < 0.0){
                    nneg_tmp ++;
                }
            }
            if(nu_dash < 0.0){
                nneg_tmp ++;
            }
            if (nneg_tmp > 0){
                continue;
            }

            double helldist_pm3 = 0.0;
            int flag_converge_pm3 = 0;
            GetMArrNval(rho_dash_arr, nu_dash,
                        data_arr, bg_arr,
                        resp_norm_mat_arr, 
                        ndet, nsky, mval_arr, &nval);
            SrtlibRlBg2SmthPm::GetRhoNu_ByPm_Nesterov(
                fp_log,
                rho_dash_arr, nu_dash,
                mval_arr, nval,
                nskyx, nskyy,
                mu,
                npm, tol_pm,
                nnewton, tol_newton,
                rho_0_new_arr,
                &nu_0_new,
                &helldist_pm3,
                &flag_converge_pm3);
            if (flag_converge_pm3 == 0){
                MiIolib::Printf2(fp_log,
                                 "iem = %d: pm3: not converged: helldist_pm3 = %.2e\n",
                                 iem,
                                 helldist_pm3);
            }
            int nneg = 0;
            for(int isky = 0; isky < nsky; isky ++){
                if(rho_0_new_arr[isky] < 0.0){
                    nneg ++;
                }
            }
            if(nu_0_new < 0.0){
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
        double helldist  = SrtlibRlStatval::GetHellingerDist(
            rho_0_arr, nu_0,
            rho_0_new_arr, nu_0_new, nsky);
        if (access( "/tmp/rl_bg2_smth_em_stop", R_OK ) != -1){
            MiIolib::Printf2(
                fp_log,
                "/tmp/rl_bg2_smth_em_stop file is found, then stop.\n");
            break;
        }
        if (helldist < tol_em){
            printf("iem = %d, helldist = %e\n",
                   iem, helldist);
            break;
        }

        dcopy_(nsky, rho_0_new_arr, 1, rho_0_arr, 1);
        nu_0 = nu_0_new;
        MiIolib::Printf2(fp_log, "iem = %d, helldist = %e, nu_0 = %e\n",
                         iem, helldist, nu_0_new);
    }

    dcopy_(nsky, rho_0_new_arr, 1, rho_new_arr, 1);
    double nu_new = nu_0_new;

    delete [] rho_0_arr;
    delete [] rho_1_arr;
    delete [] rho_2_arr;
    delete [] rho_dash_arr;
    delete [] r_rho_arr;
    delete [] r2_rho_arr;
    delete [] v_rho_arr;
    delete [] rho_0_new_arr;
    
    *nu_new_ptr = nu_new;
}

