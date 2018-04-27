#include "sub_pmmn_case2.h"

void GetDetNpix(string respdir, int* npx_ptr, int* npy_ptr)
{
    char infile[kLineSize];
    sprintf(infile, "%s/gimage_000_000.img", respdir.c_str());
    int naxis = MifFits::GetNaxis(infile);
    printf("GetDetNpix: naxis = %d\n", naxis);
    if(2 != naxis){
        printf("GetDetNpix: bad naxis = %d\n", naxis);
        abort();
    }
    int* npix_arr = new int[naxis];
    for(int iaxis = 0; iaxis < naxis; iaxis ++){
        npix_arr[iaxis] = MifFits::GetAxisSize(infile, iaxis);
    }
    *npx_ptr = npix_arr[0];
    *npy_ptr = npix_arr[1];
    delete [] npix_arr;
}

void LoadModel(string respdir, double** mat_arr_ptr)
{
    printf("--- LoadModel --- \n");

//    row: detector
//    col: sky
    
    int npx = 0;
    int npy = 0;
    GetDetNpix(respdir, &npx, &npy);
    printf("LoadModel: npx = %d\n", npx);
    printf("LoadModel: npy = %d\n", npy);
    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, npx, npy);

    int nrow = npx * npy;
    int nibin = 60;
    int njbin = 60;
    int ncol = nibin * njbin;
    printf("LoadModel: nrow = %d, ncol = %d\n", nrow, ncol);
    double* mat_arr = new double [nrow * ncol];
    for(int jbin = 0; jbin < njbin; jbin ++){
        for(int ibin = 0; ibin < nibin; ibin ++){
            char infile[kLineSize];
            sprintf(infile, "%s/gimage_%3.3d_%3.3d.img", respdir.c_str(), ibin, jbin);

            int bitpix = 0;
            double* data_arr = NULL;
            MifFits::InFitsImageD(infile, img_info,
                                  &bitpix, &data_arr);
            double sum = 0.0;
            for(int ipx = 0; ipx < nrow; ipx ++){
                sum += data_arr[ipx];
            }
            
            int kbin = nibin * jbin + ibin;
            int imat = kbin * nrow;
            for(int ipx = 0; ipx < nrow; ipx ++){
                mat_arr[imat + ipx] = data_arr[ipx] / sum;
            }
            delete [] data_arr;
        }
    }
    delete img_info;
    
    *mat_arr_ptr = mat_arr;
}

void SolveByProxMapMN(const double* const rho_arr, int nph,
                      const double* const data_arr, const double* const resp_mat_arr,
                      double beta, double mu, double lconst,
                      double tol, int nstep,
                      string outdir, string outfile_head,
                      int ndet, int nskyx, int nskyy, double epsilon,
                      double* const out_arr)
{
    char outfile_moni[kLineSize];
    char outfile_timelog[kLineSize];
    sprintf(outfile_moni, "%s/%s_moni.dat",
            outdir.c_str(), outfile_head.c_str());
    sprintf(outfile_timelog, "%s/%s_time_logl.dat",
            outdir.c_str(), outfile_head.c_str());
    FILE* fp_moni = fopen(outfile_moni, "w");
    FILE* fp_timelog = fopen(outfile_timelog, "w");
    fprintf(fp_moni, "# istep, kldiv, logl, logl.inc, delta.logl, tdiff\n");
    fprintf(fp_timelog, "# tdiff logl.inc\n");
    
    double time_st = MiTime::GetTimeSec();
    double logl_init = GetFuncL(rho_arr, data_arr, resp_mat_arr, beta, mu,
                                ndet, nskyx, nskyy);
    double logl_pre = logl_init;
    fprintf(fp_moni, "0  0  %.10e  0.0  0.0  0.0\n", logl_init);
    
    double eta = 1.2;
    int nsky = nskyx * nskyy;
    double* rho_new_arr = new double[nsky];
    double* rho_pre_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_new_arr, 1);
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, rho_pre_arr, 1);
    for(int istep = 0; istep < nstep; istep ++){
        int ik = GetFindIk(rho_new_arr, data_arr, resp_mat_arr,
                           beta, mu, lconst, eta,
                           ndet, nskyx, nskyy, epsilon);
        double lconst_pre = lconst;
        lconst = pow(eta, ik) * lconst;
        printf("SolveByProxMap: ik = %d, L = %e, L.pre = %e\n", ik, lconst, lconst_pre);
        int flag_good = 0;
        GetProxMap(rho_new_arr, data_arr, resp_mat_arr,
                   beta, mu, lconst,
                   ndet, nskyx, nskyy, epsilon,
                   rho_new_arr, &flag_good);
        double kldiv = GetKLDiv(rho_pre_arr, rho_new_arr, resp_mat_arr, ndet, nsky);
        logl_pre    = GetFuncL(rho_pre_arr, data_arr, resp_mat_arr, beta, mu, ndet, nskyx, nskyy);
        double logl = GetFuncL(rho_new_arr, data_arr, resp_mat_arr, beta, mu, ndet, nskyx, nskyy);
        
        double delta_logl = logl - logl_pre;
        double logl_inc   = logl - logl_init;
        double time = MiTime::GetTimeSec();
        double tdiff = time - time_st;
        printf("%d  %e  %.10e  %e  %e  %e\n",
               istep, kldiv, logl, logl_inc, delta_logl, tdiff);
        fprintf(fp_moni, "%d  %e  %.10e  %e  %e  %e\n",
                istep, kldiv, logl, logl_inc, delta_logl, tdiff);
        fprintf(fp_timelog, "%e  %e\n", tdiff, logl_inc);

//        // save
//        if(istep % 100 == 0){
//            char outimg[kLinSize];
//            sprintf(outimg, "%s/outimg_%4.4d.fits", outdir.c_str(), istep);
//            array = array(rho.new.vec * Nph, dim=c(ncol, nrow))
//            writeFITSim(array, file=outimg)
//        }
        
        if( -1.0e-3 < delta_logl && delta_logl < 0.0){
            break;
        }
        if(kldiv < tol){
            break;
        }

        // for next
        dcopy_(nsky, rho_new_arr, 1, rho_pre_arr, 1);
        logl_pre = logl;

    }
    
}


int GetFindIk(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              double beta, double mu, double lconst, double eta,
              int ndet, int nskyx, int nskyy, double epsilon)
{
    int ik_max = 1000;
    int ik = 0;
    double lconst_org = lconst;
    int nsky = nskyx * nskyy;
    while(ik <= ik_max){
        lconst = pow(eta, ik) * lconst_org;
        double* pLy_arr = new double[nsky];
        int flag_good = 0;
        GetProxMap(rho_arr, data_arr, resp_mat_arr,
                   beta, mu, lconst,
                   ndet, nskyx, nskyy, epsilon,
                   pLy_arr, &flag_good);
        double qminusf = GetQMinusF(pLy_arr, rho_arr,
                                    beta, mu, lconst,
                                    nskyx, nskyy);
        printf("FindIk: (ik, L, qminusf) = (%d, %e, %e)\n", ik, lconst, qminusf);
        if(qminusf >= 0 && flag_good == 1){
            break;
        }
        ik = ik + 1;
    }
    return(ik);
}


void GetProxMap(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu, double lconst,
                int ndet, int nskyx, int nskyy, double epsilon,
                double* const out_arr, int* flag_good_ptr)
{
    // sigma.vec = FuncSigma(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)
    int nsky = nskyx * nskyy;
    double* sigma_arr = new double[nsky];
    GetFuncSigma(rho_arr, data_arr,
                 beta, mu, lconst, nskyx, nskyy,
                 sigma_arr);

    dcopy_(nsky, const_cast<double*>(rho_arr), 1, out_arr, 1);
    int nem = 1000;
    int flag_good = 1;
    double tau_pre = 1.0e-10;
    for(int iem = 0; iem < nem; iem ++){
        double lem = GetFuncLem(out_arr, data_arr, resp_mat_arr, sigma_arr,
                                ndet, nsky, lconst);
        double* mval_arr = new double[nsky];
        GetFuncM(out_arr, data_arr, resp_mat_arr, ndet, nsky, mval_arr);
        double* tau_thres_arr = new double[nsky];
        GetFuncTauThres(sigma_arr, mval_arr, nsky, lconst, epsilon, tau_thres_arr);
        double tau_thres_min = 0.0;
        double tau_thres_max = 0.0;
        GetMinMax(tau_thres_arr, nsky, &tau_thres_min, &tau_thres_max);

        int nnewton = 100;
        double tol_newton = 1.0e-3;
        double tau = 0.0;
        if(tau_pre > tau_thres_min){
            tau = tau_pre;
        } else {
            tau = tau_thres_max;
        }
        for(int inewton = 0; inewton < nnewton; inewton ++){
            tau = tau - GetFuncS(tau, sigma_arr, mval_arr, nsky, lconst, epsilon)
                / GetFuncDiffS(tau, sigma_arr, mval_arr, nsky, lconst, epsilon);
            if( fabs(GetFuncS(tau, sigma_arr, mval_arr, nsky, lconst, epsilon) ) < tol_newton){
                // printf("inewton = %d, tau = %e\n", inewton, tau)
                break;
            }
        }
        tau_pre = tau;

        double* rho_tmp_arr = new double[nsky];
        dcopy_(nsky, out_arr, 1, rho_tmp_arr, 1);
        GetFuncRho(tau, sigma_arr, mval_arr, nsky, lconst, epsilon, out_arr);
        // rho.new.vec = LineSearch(rho.tmp.vec, rho.new.vec, D.vec, R.mat, sigma.vec, L, epsilon)

        double lem_new = GetFuncLem(out_arr, data_arr, resp_mat_arr, sigma_arr,
                                    ndet, nsky, lconst);
        double diff_lem = lem_new - lem;
        double kldiv = GetKLDiv(rho_tmp_arr, out_arr, resp_mat_arr, ndet, nsky);
        printf("iem = %d, kldiv = %e, diff.lem = %e, lem = %e, lem.new = %e \n",
               iem, kldiv, diff_lem, lem, lem_new);
        if(diff_lem > 0){
            flag_good = 0;
            break;
        }
        if(kldiv < 0.0){
            flag_good = 0;
            break;
        }

        if(kldiv < 1.0e-10){
            break;
        }
    }

    *flag_good_ptr = flag_good;
}

//void GetLineSearch(const double* const xval_arr,
//                   const double* const xval_new_arr,
//                   const double* const data_arr,
//                   const double* const resp_mat_arr,
//                   const double* const sigma_arr,
//                   int ndet, int nsky, double lconst,
//                   double epsilon,
//                   double* const out_arr)
//{
//    double xval0     = xval_arr[0];
//    double xval0_new = xval_new_arr[0];
//    double* xval2_arr = new double[nsky - 1];
//    double* xval2_new_arr = new double[nsky - 1];
//    for(int isky = 0; isky < nsky - 1; isky ++){
//        xval2_arr[isky]     = xval_arr[isky + 1];
//        xval2_new_arr[isky] = xval_new_arr[isky + 1];
//        theta_arr[isky] = log( xval2_arr[isky] / (1 - xval0));
//        theta_new_arr[isky] = log( xval2_new_arr[isky] / (1 - xval0_new));
//    }
//    int nstep = 100;
//    double lem_init = GetFuncLem(xval_arr, data_arr, resp_mat_arr, sigma_arr,
//                                 ndet, nsky, lconst);
//    double lem_pre = lem_init;
//
//    double* xval_pre_arr = new double[nsky];
//    dcopy_(nsky, xval_arr, 1, xval_pre_arr, 1);
//
//    char outfile[kLineSize];
//    sprintf(outfile, "temp.dat"); 
//    FILE* fp_out = fopen(outfile, "w");
//    fprintf(fp_out, "skip sing\n");
//    fprintf(fp_out, "read\n");
//    
//    double eta = 3.0;
//    for(int istep = 0; istep < nstep; istep ++){
//        double factor = pow(eta, istep);
//        double lxval0_this = factor * (log(xval0_new) - log(xval0)) + log(xval0);
//        for(int isky = 0; isky < nsky - 1; isky ++){
//            xval2_this_arr = exp(lx1.this)
//
//        }
//
//        
//            theta.this.vec = factor * (theta.new.vec - theta.vec) + theta.vec
//        alpha = sum(exp(theta.this.vec))
//        x2.this.vec = (1 - x1.this) * exp(theta.this.vec) / alpha
//        x.this.vec = c(x1.this, x2.this.vec)
//
//        if(min(x.this.vec) < epsilon){
//            printf("min(x.this.vec) < epsilon: factor(istep) = %e (%d)\n", factor, istep)
//            if(istep != 1){
//                x.new.vec = x.pre.vec
//            }
//            break
//        }
//
//        lem = FuncLem(x.this.vec, D.vec, R.mat, sigma.vec, L)
//
//        fprintf(outfile, "%d  %e\n", istep, lem - lem.init)
//
//        if(lem.pre < lem){
//            ## printf("factor(istep) = %e (%d)\n", factor, istep)
//            if(istep != 1){
//                x.new.vec = x.pre.vec
//            }
//            break
//        }
//        lem.pre = lem
//        x.pre.vec = x.this.vec
//    }
//    return(x.new.vec)
//}
//

double GetKLDiv(const double* const rho_arr,
                const double* const rho_new_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky)
{
    char* transa = new char [1];
    strcpy(transa, "N");
    // q.vec = R.mat %*% y.vec
    double* q_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, q_arr, 1);
    // q.new.vec = R.mat %*% y.new.vec
    double* q_new_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_new_arr), 1,
           0.0, q_new_arr, 1);

    // q.vec = q.vec / sum(q.vec)
    // q.new.vec = q.new.vec / sum(q.new.vec)
    double sum_q = 0.0;
    double sum_q_new = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        sum_q += q_arr[idet];
        sum_q_new += q_new_arr[idet];
    }
    dscal_(ndet, 1.0/sum_q, q_arr, 1);
    dscal_(ndet, 1.0/sum_q_new, q_new_arr, 1);
    
    double ans = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        if(q_new_arr[idet] > 0.0){
            ans = ans + q_new_arr[idet] * log( q_new_arr[idet] / q_arr[idet] );
        }
    }
    return (ans);
}

//
//
//

void GetFuncM(const double* const rho_arr,
              const double* const data_arr,
              const double* const resp_mat_arr,
              int ndet, int nsky,
              double* const out_arr)
{
    char* transa = new char [1];
    strcpy(transa, "N");    
    // num.vec = R.mat %*% rho.vec
    double* det_arr = new double[ndet];
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);

    // ans.vec = t(R.mat) %*% (D.vec / num.vec) * rho.vec
    double* div_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        div_arr[idet] = data_arr[idet] / det_arr[idet];
    }
    strcpy(transa, "T");    
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(div_arr), 1,
           0.0, out_arr, 1);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = out_arr[isky] * rho_arr[isky];
    }
}

void GetFuncRho(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon,
                double* const out_arr)
{
    // termb.vec = sigma.vec + tau / L
    // ans.vec = ( termb.vec + sqrt( termb.vec * termb.vec + 4 * m.vec / L ) ) / 2.0
    // ans.vec = mapply(ThresEpsilon, ans.vec, epsilon)
    
    for(int isky = 0; isky < nsky; isky ++){
        double termb = sigma_arr[isky] + tau / lconst;
        out_arr[isky] = ( termb + sqrt( termb * termb + 4 * mval_arr[isky] / lconst ) ) / 2.0;
        if(out_arr[isky] < epsilon){
            out_arr[isky] = epsilon;
        }
    }
}


void GetFuncDiffRho(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double lconst, double epsilon,
                    double* const out_arr)
{
    double* tau_thres_arr = new double[nsky];
    GetFuncTauThres(sigma_arr, mval_arr, nsky, lconst, epsilon, tau_thres_arr);
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = 0.0;
        if(tau <= tau_thres_arr[isky]){
            out_arr[isky] = 0.0;
        } else {
            double termb = sigma_arr[isky] + tau / lconst;
            double root = sqrt( termb * termb + 4 * mval_arr[isky] / lconst );
            out_arr[isky] = (termb + root) / (2 * lconst * root);
        }
    }
}

void GetFuncTauThres(const double* const sigma_arr,
                     const double* const mval_arr,
                     int nsky, double lconst, double epsilon,
                     double* const out_arr)
{
    for(int isky = 0; isky < nsky; isky ++){
        out_arr[isky] = lconst * (epsilon - sigma_arr[isky]) - mval_arr[isky] / epsilon;
    }
}

double GetFuncS(double tau,
                const double* const sigma_arr,
                const double* const mval_arr,
                int nsky, double lconst, double epsilon)
{
    // ans = sum( FuncRho(tau, sigma.vec, m.vec, L, epsilon) ) - 1
    double* rho_arr = new double[nsky];
    GetFuncRho(tau, sigma_arr, mval_arr, nsky, lconst, epsilon, rho_arr);
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum += rho_arr[isky];
    }
    double ans = sum - 1.0;
    delete [] rho_arr;
    return(ans);
}

double GetFuncDiffS(double tau,
                    const double* const sigma_arr,
                    const double* const mval_arr,
                    int nsky, double lconst, double epsilon)
{
    // ans = sum(FuncDiffRhoVec(tau, sigma.vec, m.vec, L, epsilon))
    
    double* diff_rho_arr = new double[nsky];
    GetFuncDiffRho(tau, sigma_arr, mval_arr,
                   nsky, lconst, epsilon,
                   diff_rho_arr);
    double ans = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        ans += diff_rho_arr[isky];
    }
    return(ans);
}


//
//
//
 
double GetFuncLem(const double* const rho_arr,
                  const double* const data_arr,
                  const double* const resp_mat_arr,
                  const double* const sigma_arr,
                  int ndet, int nsky,
                  double lconst)
{
    double* diff_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, diff_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(sigma_arr), 1, diff_arr, 1);
    double term1 = lconst
        * ddot_(nsky, const_cast<double*>(diff_arr), 1, const_cast<double*>(diff_arr), 1) / 2.0;
    double term2 = GetFuncG(rho_arr, data_arr, resp_mat_arr, ndet, nsky);
    double ans = term1 + term2;
    return(ans);
}

    
double GetFuncL(const double* const rho_arr,
                const double* const data_arr,
                const double* const resp_mat_arr,
                double beta, double mu,
                int ndet, int nskyx, int nskyy)
{
    int nsky = nskyx * nskyy;
    double term1 = GetFuncF(rho_arr, beta, mu, nskyx, nskyy);
    double term2 = GetFuncG(rho_arr, data_arr, resp_mat_arr, ndet, nsky);
    double ans = term1 + term2;
    return(ans);
}


double GetFuncLSupp(const double* const rho_arr,
                    const int* const rho_supp_arr,
                    const double* const data_arr,
                    const double* const resp_mat_arr,
                    double beta, double mu,
                    int ndet, int nskyx, int nskyy)
{
    int nsky = nskyx * nskyy;
    double term1 = GetFuncFSupp(rho_arr, rho_supp_arr, beta, mu, nskyx, nskyy);
    double term2 = GetFuncG(rho_arr, data_arr, resp_mat_arr, ndet, nsky);
    double ans = term1 + term2;
    return(ans);
}

double GetQMinusF(const double* const rho_new_arr,
                  const double* const rho_arr,
                  double beta, double mu, double lconst,
                  int nskyx, int nskyy)
{
    //term1 =      FuncF(rho.vec, beta, mu, nrow, ncol)
    //term2 = -1 * FuncF(rho.new.vec, beta, mu, nrow, ncol)
    
    double term1 = GetFuncF(rho_arr, beta, mu, nskyx, nskyy);
    double term2 = -1 * GetFuncF(rho_new_arr, beta, mu, nskyx, nskyy);

    // term3 = sum( (rho.new.vec - rho.vec) * DiffF(rho.vec, beta, mu, nrow, ncol) )
    // term4 = L / 2.0 * sum( (rho.new.vec - rho.vec) * (rho.new.vec - rho.vec) )
    int nsky = nskyx * nskyy;
    double* diff_rho_arr = new double[nsky];
    dcopy_(nsky, const_cast<double*>(rho_new_arr), 1, diff_rho_arr, 1);
    daxpy_(nsky, -1.0, const_cast<double*>(rho_arr), 1, diff_rho_arr, 1);
    double* diff_f_arr = new double[nsky];
    GetDiffF(rho_arr, beta, mu, nskyx, nskyy, diff_f_arr);
    double term3 = ddot_(nsky, const_cast<double*>(diff_rho_arr), 1, const_cast<double*>(diff_f_arr), 1);
    double term4 = lconst *
        ddot_(nsky, const_cast<double*>(diff_rho_arr), 1, const_cast<double*>(diff_rho_arr), 1) / 2.0;
    double ans = term1 + term2 + term3 + term4;

    delete [] diff_rho_arr;
    delete [] diff_f_arr;
    return(ans);
}


void GetFuncSigma(const double* const rho_arr,
                  const double* const data_arr,
                  double beta, double mu,
                  double lconst, int nskyx, int nskyy,
                  double* out_arr)
{
//    term1 = rho.vec
//    term2 = -1.0 / L * DiffF(rho.vec, beta, mu, nrow, ncol)
//    ans.vec = term1 + term2

    int nsky = nskyx * nskyy;
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, out_arr, 1);
    double* difff_arr = new double[nsky];
    GetDiffF(rho_arr, beta, mu, nskyx, nskyy, difff_arr);
    daxpy_(nsky, -1.0 / lconst, const_cast<double*>(difff_arr), 1, out_arr, 1);
}

void GetFuncSigmaSupp(const double* const rho_arr,
                      const int* const rho_supp_arr,
                      const double* const data_arr,
                      double beta, double mu,
                      double lconst, int nskyx, int nskyy,
                      double* out_arr)
{
    int nsky = nskyx * nskyy;
    dcopy_(nsky, const_cast<double*>(rho_arr), 1, out_arr, 1);
    double* difff_arr = new double[nsky];
    GetDiffFSupp(rho_arr, rho_supp_arr, beta, mu, nskyx, nskyy, difff_arr);
    daxpy_(nsky, -1.0 / lconst, const_cast<double*>(difff_arr), 1, out_arr, 1);
}

double GetFuncF(const double* const rho_arr,
                double beta, double mu,
                int nskyx, int nskyy)
{
    //    term1 = (1.0 - beta) * sum(log(rho.vec))
    //    term2 = mu * TermV(rho.vec, nrow, ncol)
    //    ans = term1 + term2
    
    int nsky = nskyx * nskyy;
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum += log(rho_arr[isky]);
    }
    double ans = (1.0 - beta) * sum + mu * GetTermV(rho_arr, nskyx, nskyy);
    return(ans);
}

double GetFuncFSupp(const double* const rho_arr,
                    const int* const rho_supp_arr,
                    double beta, double mu,
                    int nskyx, int nskyy)
{
    // term1 = (1.0 - beta) * SumLogPlus(rho.vec, rho.supp.vec)
    // term2 = mu * TermV(rho.vec, nrow, ncol)
    // ans = term1 + term2

    int nsky = nskyx * nskyy;
    double ans
        = (1.0 - beta) * GetSumLogPlus(rho_arr, nsky, rho_supp_arr)
        + mu * GetTermV(rho_arr, nskyx, nskyy);
    return(ans);
}

void GetDiffF(const double* const rho_arr,
              double beta, double mu,
              int nskyx, int nskyy,
              double* const out_arr)
{
    // term1 = (1 - beta) / rho.vec
    int nsky = nskyx * nskyy;
    double* rho_inv_arr = new double[nsky];
    GetInvArr(rho_arr, nsky, rho_inv_arr);

    // term2 = mu * DiffTermV(rho.vec, nrow, ncol)
    GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);

    // ans = term1 + term2
    daxpy_(nsky, 1.0 - beta, rho_inv_arr, 1, out_arr, 1);
    delete [] rho_inv_arr;
}

void GetDiffFSupp(const double* const rho_arr,
                  const int* const rho_supp_arr,
                  double beta, double mu,
                  int nskyx, int nskyy,
                  double* const out_arr)
{
    // term1 = (1 - beta) / rho.vec
    int nsky = nskyx * nskyy;
    double* rho_inv_arr = new double[nsky];
    GetInvArrSupp(rho_arr, nsky, rho_supp_arr, rho_inv_arr);

    // term2 = mu * DiffTermV(rho.vec, nrow, ncol)
    GetDiffTermV(rho_arr, nskyx, nskyy, out_arr);
    dscal_(nsky, mu, out_arr, 1);

    // ans = term1 + term2
    daxpy_(nsky, 1.0 - beta, rho_inv_arr, 1, out_arr, 1);
    delete [] rho_inv_arr;
}


double GetFuncG(const double* const rho_arr, 
                const double* const data_arr,
                const double* const resp_mat_arr,
                int ndet, int nsky)
{
    //    num.vec = R.mat %*% (rho.vec)
    double* det_arr = new double[ndet];

    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rho_arr), 1,
           0.0, det_arr, 1);

    //ans = -1 * sum( D.vec * log( num.vec ) )
    for(int idet = 0; idet < ndet; idet ++){
        det_arr[idet] = log(det_arr[idet]);
    }
    double ans = -1.0 * ddot_(ndet, const_cast<double*>(data_arr), 1, const_cast<double*>(det_arr), 1);
    return(ans);
}


void GetDiffTermV(const double* const rho_arr, int nskyx, int nskyy,
                  double* const rho_diff_arr)
{
    // iskyx = 0, iskyy = 0
    // isky_plus_x, isky_plus_y
    {
        int iskyx = 0;
        int iskyy = 0;
        int isky = GetIbin(iskyx, iskyy, nskyx);
        int isky_plus_x = GetIbin(iskyx + 1, iskyy, nskyx);
        int isky_plus_y = GetIbin(iskyx, iskyy + 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_plus_x])
            + (rho_arr[isky] - rho_arr[isky_plus_y]);
    }

    // iskyx = 0, 1 <= iskyy <= nskyy - 2
    // isky_plus_x, isky_plus_y, isky_minus_y
    {
        int iskyx = 0;
        for(int iskyy = 1; iskyy < nskyy - 1; iskyy ++){
            int isky = GetIbin(iskyx, iskyy, nskyx);
            int isky_plus_x = GetIbin(iskyx + 1, iskyy, nskyx);
            int isky_plus_y = GetIbin(iskyx, iskyy + 1, nskyx);
            int isky_minus_y =  GetIbin(iskyx, iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }

    // iskyx = 0, iskyy = nskyy - 1
    // isky_plus_x, isky_minus_y
    {
        int iskyx = 0;
        int iskyy = nskyy - 1;
        int isky = GetIbin(iskyx, iskyy, nskyx);
        int isky_plus_x = GetIbin(iskyx + 1, iskyy, nskyx);
        int isky_minus_y = GetIbin(iskyx, iskyy - 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_plus_x])
            + (rho_arr[isky] - rho_arr[isky_minus_y]);
    }

    // 1 <= iskyx <= nskyx - 2, iskyy = 0
    // isky_plus_x, isky_minus_x, isky_plus_y
    {
        int iskyy = 0;
        for(int iskyx = 1; iskyx < nskyx - 1; iskyx ++){
            int isky = GetIbin(iskyx, iskyy, nskyx);
            int isky_plus_x  = GetIbin(iskyx + 1, iskyy    , nskyx);
            int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_plus_y  = GetIbin(iskyx    , iskyy + 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y]);
        }
    }

    // 1 <= iskyx <= nskyx - 2, 1 <= iskyy <= nskyy - 2
    // isky_plus_x, isky_minus_x, isky_plus_y, isky_minus_y
    for(int iskyx = 1; iskyx < nskyx - 1; iskyx ++){
        for(int iskyy = 1; iskyy < nskyy - 1; iskyy ++){
            int isky         = GetIbin(iskyx    , iskyy    , nskyx);
            int isky_plus_x  = GetIbin(iskyx + 1, iskyy    , nskyx);
            int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_plus_y  = GetIbin(iskyx    , iskyy + 1, nskyx);
            int isky_minus_y = GetIbin(iskyx    , iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }

    // 1 <= iskyx <= nskyx - 2, iskyy = nskyy - 1
    // isky_plus_x, isky_minus_x, isky_minus_y
    {
        int iskyy = nskyy - 1;
        for(int iskyx = 1; iskyx < nskyx - 1; iskyx ++){
            int isky          = GetIbin(iskyx    , iskyy    , nskyx);
            int isky_plus_x   = GetIbin(iskyx + 1, iskyy    , nskyx);
            int isky_minus_x  = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_minus_y  = GetIbin(iskyx    , iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_plus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }

    // iskyx = nskyx - 1, iskyy = 0
    // isky_minus_x, isky_plus_y
    {
        int iskyx = nskyx - 1;
        int iskyy = 0;
        int isky         = GetIbin(iskyx    , iskyy    , nskyx);
        int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
        int isky_plus_y  = GetIbin(iskyx    , iskyy + 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_minus_x])
            + (rho_arr[isky] - rho_arr[isky_plus_y]);
    }

    // iskyx = nskyx - 1, 1 <= iskyy <= nskyy - 2
    // isky_minus_x, isky_plus_y, isky_minus_y
    {
        int iskyx = nskyx - 1;
        for(int iskyy = 1; iskyy < nskyy - 1; iskyy ++){
            int isky          = GetIbin(iskyx    , iskyy    , nskyx);
            int isky_minus_x  = GetIbin(iskyx - 1, iskyy    , nskyx);
            int isky_plus_y   = GetIbin(iskyx    , iskyy + 1, nskyx);
            int isky_minus_y  = GetIbin(iskyx    , iskyy - 1, nskyx);
            rho_diff_arr[isky]
                = (rho_arr[isky] - rho_arr[isky_minus_x])
                + (rho_arr[isky] - rho_arr[isky_plus_y])
                + (rho_arr[isky] - rho_arr[isky_minus_y]);
        }
    }
    
    // iskyx = nskyx - 1, iskyy = nskyy - 1
    // isky_minus_x, isky_minus_y
    {
        int iskyx = nskyx - 1;
        int iskyy = nskyy - 1;
        int isky         = GetIbin(iskyx    , iskyy    , nskyx);
        int isky_minus_x = GetIbin(iskyx - 1, iskyy    , nskyx);
        int isky_minus_y = GetIbin(iskyx    , iskyy - 1, nskyx);
        rho_diff_arr[isky]
            = (rho_arr[isky] - rho_arr[isky_minus_x])
            + (rho_arr[isky] - rho_arr[isky_minus_y]);
    }
}

double GetTermV(const double* const rho_arr, int nskyx, int nskyy)
{
    double termv = 0.0;
    for(int iskyx = 0; iskyx < nskyx - 1; iskyx ++){
        for(int iskyy = 0; iskyy < nskyy - 1; iskyy ++){
            int isky = iskyx + iskyy * nskyx;
            int isky_plus_x = (iskyx + 1) + iskyy * nskyx;
            int isky_plus_y = iskyx + (iskyy + 1) * nskyx;
            double diff1 = rho_arr[isky] - rho_arr[isky_plus_x];
            double diff2 = rho_arr[isky] - rho_arr[isky_plus_y];
            termv += diff1 * diff1 + diff2 * diff2;
        }
    }
    for(int iskyx = 0; iskyx < nskyx - 1; iskyx ++){
        int isky = iskyx + (nskyy - 1) * nskyx;
        int isky_plus_x = (iskyx + 1) + (nskyy - 1) * nskyx;
        double diff1 = rho_arr[isky] - rho_arr[isky_plus_x];
        termv += diff1 * diff1;
    }
    for(int iskyy = 0; iskyy < nskyy - 1; iskyy ++){
        int isky = (nskyx - 1) + iskyy * nskyx;
        int isky_plus_y = (nskyx - 1) + (iskyy + 1) * nskyx;
        double diff2 = rho_arr[isky] - rho_arr[isky_plus_y];
        termv += diff2 * diff2;
    }
    return (termv);
}


double GetSumLogPlus(const double* const data_arr, int ndata,
                     const int* const supp_arr)
{
    double ans = 0.0;
    for(int idata = 0; idata < ndata; idata ++){
        if(supp_arr[idata] > 0){
            ans = ans + log(data_arr[idata]);
        }
    }
    return(ans);
}

void GetInvArr(const double* const data_arr, int ndata,
               double* const out_arr)
{
    for(int idata = 0; idata < ndata; idata ++){
        out_arr[idata] = 1.0 / data_arr[idata];
    }
}

void GetInvArrSupp(const double* const data_arr, int ndata,
                   const int* const supp_arr, 
                   double* const out_arr)
{
    for(int idata = 0; idata < ndata; idata ++){
        if(supp_arr[idata] > 0){
            out_arr[idata] = 1.0 / data_arr[idata];
        } else {
            out_arr[idata] = 0.0;
        }
    }
}
    
void GetSuppArr(const double* const data_arr, int ndata,
                double epsilon,
                int* const out_arr)
{
    for(int idata = 0; idata < ndata; idata ++){
        out_arr[idata] = 0;
        if(data_arr[idata] > epsilon){
            out_arr[idata] = 1;
        }
    }
}

int GetNZeroArr(const double* const data_arr, int ndata,
                double epsilon)
{
    int nzero = 0;
    for(int idata = 0; idata < ndata; idata ++){
        if(fabs(data_arr[idata]) <= epsilon){
            nzero = nzero + 1;
        }
    }
    return(nzero);
}

int GetIbin(int ibinx, int ibiny, int nbinx)
{
    int ibin = ibinx + ibiny * nbinx;
    return(ibin);
}

void GetMinMax(const double* const data_arr, int ndata,
               double* min_ptr, double* max_ptr)
{
    double min = data_arr[0];
    double max = data_arr[0];
    for(int idata = 0; idata < ndata; idata ++){
        if(min < data_arr[idata]){
            min = data_arr[idata];
        }

        if(max > data_arr[idata]){
            max = data_arr[idata];
        }
    }
    *min_ptr = min;
    *max_ptr = max;
}
