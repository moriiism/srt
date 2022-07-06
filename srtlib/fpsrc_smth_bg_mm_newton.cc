#include "mir_math.h"
#include "fpsrc_smth_bg_mm_newton.h"

void FpsrcSmthBgMmNewton::GetRhoArr_FromLambda_Mm(double lambda,
                                                  double lip_const,
                                                  const double* const vval_arr,
                                                  const double* const mval_arr,
                                                  int nsky,
                                                  double* const rho_arr)
{
    for(int isky = 0; isky < nsky; isky++){
        double bval = lambda - lip_const * vval_arr[isky];
        double num = sqrt(bval * bval + 4 * lip_const * mval_arr[isky]) - bval;
        double den = 2 * lip_const;
        rho_arr[isky] = num / den;
    }
}

void FpsrcSmthBgMmNewton::GetNuArr_FromLambda_Mm(double lambda,
                                                 double lip_const,
                                                 const double* const wval_arr,
                                                 const double* const nval_arr,
                                                 int nsrc,
                                                 double* const nu_arr)
{
    for(int isrc = 0; isrc < nsrc; isrc++){
        double bval = lambda - lip_const * wval_arr[isrc];
        double num = sqrt(bval * bval + 4 * lip_const * nval_arr[isrc]) - bval;
        double den = 2 * lip_const;
        nu_arr[isrc] = num / den;
    }
}

double FpsrcSmthBgMmNewton::GetPhi_FromLambda_Mm(double lambda,
                                                 double nph,
                                                 double lip_const,
                                                 double zval,
                                                 double phi_val,
                                                 double pval)
{
    double bval = lambda - lip_const * zval + nph / phi_val;
    double num = sqrt(bval * bval + 4 * lip_const * pval) - bval;
    double den = 2 * lip_const;
    double phi = num / den;
    return phi;
}

// get derivative of rho from lambda
void FpsrcSmthBgMmNewton::GetDerivRhoArr_FromLambda_Mm(double lambda,
                                                       double lip_const,
                                                       const double* const vval_arr,
                                                       const double* const mval_arr,
                                                       int nsky,
                                                       double* const deriv_rho_arr)
{
    double* rho_arr = new double[nsky];
    GetRhoArr_FromLambda_Mm(lambda,
                            lip_const,
                            vval_arr,
                            mval_arr,
                            nsky,
                            rho_arr);
    for(int isky = 0; isky < nsky; isky++){
        double bval = lambda - lip_const * vval_arr[isky];
        double den = sqrt(bval * bval + 4 * lip_const * mval_arr[isky]);
        double num = -1 * rho_arr[isky];
        deriv_rho_arr[isky] = num / den;

        if(1 == isnan(deriv_rho_arr[isky])){
            printf("GetDerivRhoArr_FromLambda_Mm: vval_arr[%d] = %e\n",
                   isky, vval_arr[isky]);
            printf("GetDerivRhoArr_FromLambda_Mm: lip_const = %e\n",
                   lip_const);
            printf("GetDerivRhoArr_FromLambda_Mm: lambda = %e\n",
                   lambda);
            printf("GetDerivRhoArr_FromLambda_Mm: mval_arr[%d] = %e\n",
                   isky, mval_arr[isky]);
            printf("GetDerivRhoArr_FromLambda_Mm: rho_arr[%d] = %e\n",
                   isky, rho_arr[isky]);
            printf("GetDerivRhoArr_FromLambda_Mm: bval = %e\n",
                   bval);
        }
    }
    
    delete [] rho_arr;
}


// get derivative of nu from lambda
void FpsrcSmthBgMmNewton::GetDerivNuArr_FromLambda_Mm(double lambda,
                                                      double lip_const,
                                                      const double* const wval_arr,
                                                      const double* const nval_arr,
                                                      int nsrc,
                                                      double* const deriv_nu_arr)
{
    double* nu_arr = new double[nsrc];
    GetNuArr_FromLambda_Mm(lambda,
                           lip_const,
                           wval_arr,
                           nval_arr,
                           nsrc,
                           nu_arr);
    for(int isrc = 0; isrc < nsrc; isrc++){
        double bval = lambda - lip_const * wval_arr[isrc];
        double den = sqrt(bval * bval + 4 * lip_const * nval_arr[isrc]);
        double num = -1 * nu_arr[isrc];
        deriv_nu_arr[isrc] = num / den;
    }
    delete [] nu_arr;
}

// get derivative of phi from lambda
double FpsrcSmthBgMmNewton::GetDerivPhi_FromLambda_Mm(double lambda,
                                                      double nph,
                                                      double lip_const,
                                                      double zval,
                                                      double phi_val,
                                                      double pval)
{
    double phi = GetPhi_FromLambda_Mm(lambda,
                                      nph,
                                      lip_const,
                                      zval,
                                      phi_val,
                                      pval);
    double bval = lambda - lip_const * zval + nph / phi_val;
    double den = sqrt(bval * bval + 4 * lip_const * pval);
    double num = -1 * phi;
    double deriv_phi = num / den;
    return deriv_phi;
}

double FpsrcSmthBgMmNewton::GetSval_FromLambda_Mm(double lambda,
                                                  double lip_const,
                                                  const double* const vval_arr,
                                                  const double* const wval_arr,
                                                  double zval,
                                                  const double* const mval_arr,
                                                  const double* const nval_arr,
                                                  double pval,
                                                  double phi_val,
                                                  int nsky, int nsrc, int nph)
{
    double* rho_arr = new double[nsky];
    GetRhoArr_FromLambda_Mm(lambda,
                            lip_const,
                            vval_arr,
                            mval_arr,
                            nsky,
                            rho_arr);
    double* nu_arr = new double[nsrc];
    GetNuArr_FromLambda_Mm(lambda,
                           lip_const,
                           wval_arr,
                           nval_arr,
                           nsrc,
                           nu_arr);
    double phi = GetPhi_FromLambda_Mm(lambda,
                                      nph,
                                      lip_const,
                                      zval,
                                      phi_val,
                                      pval);
    double sval = phi - 1.0;
    for(int isky = 0; isky < nsky; isky ++){
        sval += rho_arr[isky];
    }
    for(int isrc = 0; isrc < nsrc; isrc ++){
        sval += nu_arr[isrc];
    }
    delete [] rho_arr;
    delete [] nu_arr;
    return sval;
}

double FpsrcSmthBgMmNewton::GetDerivSval_FromLambda_Mm(
    double lambda,
    double lip_const,
    const double* const vval_arr,
    const double* const wval_arr,
    double zval,
    const double* const mval_arr,
    const double* const nval_arr,
    double pval,
    double phi_val,
    int nsky, int nsrc, double nph)
{
    double* deriv_rho_arr = new double[nsky];
    GetDerivRhoArr_FromLambda_Mm(lambda,
                                 lip_const,
                                 vval_arr,
                                 mval_arr,
                                 nsky,
                                 deriv_rho_arr);
    double* deriv_nu_arr = new double[nsrc];
    GetDerivNuArr_FromLambda_Mm(lambda,
                                lip_const,
                                wval_arr,
                                nval_arr,
                                nsrc,
                                deriv_nu_arr);
    double deriv_phi = GetDerivPhi_FromLambda_Mm(lambda,
                                                 nph,
                                                 lip_const,
                                                 zval,
                                                 phi_val,
                                                 pval);
    double deriv_sval = deriv_phi;
    for(int isky = 0; isky < nsky; isky ++){
        deriv_sval += deriv_rho_arr[isky];
    }
    for(int isrc = 0; isrc < nsrc; isrc ++){
        deriv_sval += deriv_nu_arr[isrc];
    }
    if(1 == isnan(deriv_sval)){
        printf("            : GetDerivSval_FromLambda_Mm: deriv_sval = %e\n",
               deriv_sval);
        double max_deriv_rho = MirMath::GetMax(nsky, deriv_rho_arr);
        double min_deriv_rho = MirMath::GetMin(nsky, deriv_rho_arr);
        double max_deriv_nu = MirMath::GetMax(nsrc, deriv_nu_arr);
        double min_deriv_nu = MirMath::GetMin(nsrc, deriv_nu_arr);
        printf("            : GetDerivSval_FromLambda_Mm: max_deriv_rho = %e\n",
               max_deriv_rho);
        printf("            : GetDerivSval_FromLambda_Mm: min_deriv_rho = %e\n",
               min_deriv_rho);
        printf("            : GetDerivSval_FromLambda_Mm: max_deriv_nu = %e\n",
               max_deriv_nu);
        printf("            : GetDerivSval_FromLambda_Mm: min_deriv_nu = %e\n",
               min_deriv_nu);
        printf("            : GetDerivSval_FromLambda_Mm: deriv_phi = %e\n",
               deriv_phi);
        for(int isky = 0; isky < nsky; isky ++){
            if(1 == isnan(deriv_rho_arr[isky])){
                printf("            : GetDerivSval_FromLambda_Mm: deriv_rho_arr[%d] = %e\n",
                       isky, deriv_rho_arr[isky]);
            }
        }
        for(int isrc = 0; isrc < nsrc; isrc ++){
            if(1 == isnan(deriv_nu_arr[isrc])){
                printf("            : GetDerivSval_FromLambda_Mm: deriv_nu_arr[%d] = %e\n",
                       isrc, deriv_nu_arr[isrc]);
            }
        }
        abort();
    }

    delete [] deriv_rho_arr;
    delete [] deriv_nu_arr;
    
    return deriv_sval;
}

double FpsrcSmthBgMmNewton::GetLambdaUpdate_ByNewton_Mm(
    double lambda,
    double lip_const,
    const double* const vval_arr,
    const double* const wval_arr,
    double zval,
    const double* const mval_arr,
    const double* const nval_arr,
    double pval,
    double phi_val,
    int nsky, int nsrc, int nph)
{
    double sval = GetSval_FromLambda_Mm(lambda,
                                        lip_const,
                                        vval_arr, wval_arr, zval,
                                        mval_arr, nval_arr, pval,
                                        phi_val,
                                        nsky, nsrc, nph);
    double deriv_sval = GetDerivSval_FromLambda_Mm(lambda,
                                                   lip_const,
                                                   vval_arr, wval_arr, zval,
                                                   mval_arr, nval_arr, pval,
                                                   phi_val,
                                                   nsky, nsrc, nph);
    double lambda_new = lambda - sval / deriv_sval;

    if(1 == isnan(lambda_new)){
        printf("            : GetLambdaUpdate_ByNewton_Mm: sval = %e\n",
               sval);
        printf("            : GetLambdaUpdate_ByNewton_Mm: deriv_sval = %e\n",
               deriv_sval);
        printf("            : GetLambdaUpdate_ByNewton_Mm: lambda = %e\n",
               lambda);
    }
    return lambda_new;
}

double FpsrcSmthBgMmNewton::GetLambda_ByNewton_Mm(
    FILE* const fp_log,
    double lambda_init,
    const double* const vval_arr,
    const double* const wval_arr,
    double zval,
    const double* const mval_arr,
    const double* const nval_arr,
    double pval,
    double phi_val,
    int nsky, int nsrc, int nph,
    double lip_const,
    int nnewton, double tol_newton)
{
    double lambda = lambda_init;
    double lambda_new = lambda_init;
    int flag_converge = 0;
    double sval = 0.0;
    for(int inewton = 0; inewton < nnewton; inewton ++){
        lambda_new = GetLambdaUpdate_ByNewton_Mm(lambda,
                                                 lip_const,
                                                 vval_arr, wval_arr, zval,
                                                 mval_arr, nval_arr, pval,
                                                 phi_val,
                                                 nsky, nsrc, nph);
        sval = GetSval_FromLambda_Mm(lambda_new,
                                     lip_const,
                                     vval_arr, wval_arr, zval,
                                     mval_arr, nval_arr, pval,
                                     phi_val,
                                     nsky, nsrc, nph);
        if(fabs(sval) < tol_newton){
            flag_converge = 1;
            break;
        }
        lambda = lambda_new;
    }
    if(flag_converge != 1){
        MiIolib::Printf2(fp_log,
                         "      newton: not converge: sval = %e\n", sval);
        MiIolib::Printf2(fp_log,
                         "            : lambda_new = %e\n", lambda_new);
        MiIolib::Printf2(fp_log,
                         "            : lambda_init = %e\n", lambda_init);
        MiIolib::Printf2(fp_log,
                         "            : lip_const = %e\n", lip_const);        
        double max_vval = MirMath::GetMax(nsky, vval_arr);
        double min_vval = MirMath::GetMin(nsky, vval_arr);
        double max_wval = MirMath::GetMax(nsrc, wval_arr);
        double min_wval = MirMath::GetMin(nsrc, wval_arr);
        double max_mval = MirMath::GetMax(nsky, mval_arr);
        double min_mval = MirMath::GetMin(nsky, mval_arr);
        double max_nval = MirMath::GetMax(nsrc, nval_arr);
        double min_nval = MirMath::GetMin(nsrc, nval_arr);
        MiIolib::Printf2(fp_log,
                         "            : max_vval = %e\n", max_vval);
        MiIolib::Printf2(fp_log,
                         "            : min_vval = %e\n", min_vval);
        MiIolib::Printf2(fp_log,
                         "            : max_wval = %e\n", max_wval);
        MiIolib::Printf2(fp_log,
                         "            : min_wval = %e\n", min_wval);
        MiIolib::Printf2(fp_log,
                         "            : zval = %e\n", zval);
        MiIolib::Printf2(fp_log,
                         "            : max_mval = %e\n", max_mval);
        MiIolib::Printf2(fp_log,
                         "            : min_mval = %e\n", min_mval);        
        MiIolib::Printf2(fp_log,
                         "            : max_nval = %e\n", max_nval);
        MiIolib::Printf2(fp_log,
                         "            : min_nval = %e\n", min_nval);
        MiIolib::Printf2(fp_log,
                         "            : pval = %e\n", pval);
        MiIolib::Printf2(fp_log,
                         "            : phi_val = %e\n", phi_val);
        MiIolib::Printf2(fp_log,
                         "      abort.\n");
        abort();
    }
    
    return lambda_new;
}


void FpsrcSmthBgMmNewton::GetRhoArrNuArrPhi_ByNewton_Mm(
    FILE* const fp_log,
    const double* const vval_arr,
    const double* const wval_arr,
    double zval,
    const double* const mval_arr,
    const double* const nval_arr,
    double pval,
    double phi_val,
    int nsky, int nsrc, int nph,
    double lip_const,
    int nnewton, double tol_newton,
    double lambda,
    double* const rho_new_arr,
    double* const nu_new_arr,
    double* const phi_new_ptr,
    double* const lambda_new_ptr)
{
    double lambda_new = GetLambda_ByNewton_Mm(fp_log,
                                              lambda,
                                              vval_arr, wval_arr, zval,
                                              mval_arr, nval_arr, pval,
                                              phi_val,
                                              nsky, nsrc, nph,
                                              lip_const,
                                              nnewton, tol_newton);
    GetRhoArr_FromLambda_Mm(lambda_new,
                            lip_const,
                            vval_arr,
                            mval_arr,
                            nsky,
                            rho_new_arr);
    GetNuArr_FromLambda_Mm(lambda_new,
                           lip_const,
                           wval_arr,
                           nval_arr,
                           nsrc,
                           nu_new_arr);
    double phi_new = GetPhi_FromLambda_Mm(lambda_new,
                                          nph,
                                          lip_const,
                                          zval,
                                          phi_val,
                                          pval);
    *phi_new_ptr = phi_new;
    *lambda_new_ptr = lambda_new;
}
