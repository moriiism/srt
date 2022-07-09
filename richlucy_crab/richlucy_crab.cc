// Image reconstruction by Richardson-Lucy method
// with one point source fixed position
// and with smoothness regularization

#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_richlucy_crab.h"
#include "rl.h"
#include "rl_bg2_smth_em.h"
#include "TRandom3.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValRichlucyCrab* argval = new ArgValRichlucyCrab;
    argval->Init(argc, argv);
    argval->Print(stdout);

    char logfile[kLineSize];
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        system(cmd);
    }
    sprintf(logfile, "%s/%s_%s.log",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str(),
            argval->GetProgname().c_str());
    FILE* fp_log = fopen(logfile, "w");
    argval->Print(fp_log);
    MiIolib::Printf2(fp_log, "-----------------------------\n");    

    int nskyx = argval->GetNskyx();
    int nskyy = argval->GetNskyy();
    int ndetx = argval->GetNdetx();
    int ndety = argval->GetNdety();
    int nsky = nskyx * nskyy;
    int ndet = ndetx * ndety;    
    
    // load data
    MifImgInfo* img_info_data = new MifImgInfo;
    img_info_data->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_data = 0;
    double* data_arr = NULL;
    MifFits::InFitsImageD(argval->GetDatafile(), img_info_data,
                          &bitpix_data, &data_arr);
    int nph_data = MirMath::GetSum(ndet, data_arr);
    MiIolib::Printf2(fp_log, "N photon = %d\n", nph_data);    

    // load sky image of fixed source with normalized flux
    MifImgInfo* img_info_fixed_src_norm = new MifImgInfo;
    img_info_fixed_src_norm->InitSetImg(1, 1, nskyx, nskyy);
    int bitpix_fixed_src_norm = 0;
    double* sky_fixed_src_norm_arr = NULL;
    MifFits::InFitsImageD(argval->GetFixedSrcNormFile(),
                          img_info_fixed_src_norm,
                          &bitpix_fixed_src_norm,
                          &sky_fixed_src_norm_arr);
    int nph_fixed_src_norm = MirMath::GetSum(nsky, sky_fixed_src_norm_arr);
    MiIolib::Printf2(fp_log, "N photon fixed source with normalized flux = %d\n",
                     nph_fixed_src_norm);
    
    // load response file
    int naxis0 = MifFits::GetAxisSize(argval->GetRespFile(), 0);
    int naxis1 = MifFits::GetAxisSize(argval->GetRespFile(), 1);
    if ((naxis0 != ndet) || (naxis1 != nsky)){
        MiIolib::Printf2(fp_log, "Error: response file size error.\n");
        abort();
    }
    double* resp_mat_arr = NULL;
    int bitpix_resp = 0;
    MifImgInfo* img_info_resp = new MifImgInfo;
    img_info_resp->InitSetImg(1, 1, ndet, nsky);
    MifFits::InFitsImageD(argval->GetRespFile(), img_info_resp,
                          &bitpix_resp, &resp_mat_arr);

    // load efficiency file
    double* eff_mat_arr = NULL;
    int bitpix_eff = 0;
    MifImgInfo* img_info_eff = new MifImgInfo;
    img_info_eff->InitSetImg(1, 1, nskyx, nskyy);
    MifFits::InFitsImageD(argval->GetEffFile(), img_info_eff,
                          &bitpix_eff, &eff_mat_arr);

    // normalize response file
    double* resp_norm_mat_arr = new double [ndet * nsky];
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;
            for(int idet = 0; idet < ndet; idet ++){
                resp_norm_mat_arr[imat + idet]
                    = resp_mat_arr[imat + idet] / eff_mat_arr[isky];
            }
        }
    }
    // check
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;            
            double resp_norm_sum = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                resp_norm_sum += resp_norm_mat_arr[imat + idet];
            }
            if ( fabs(resp_norm_sum - 1.0) > 1.0e-10){
                // printf("warning: resp_norm_sum = %e\n", resp_norm_sum);
            }
        }
    }

    // det image of fixed source with normalized flux
    double* det_fixed_src_norm_arr = new double[ndet];
    SrtlibRl::GetDetArr(sky_fixed_src_norm_arr,
                        resp_norm_mat_arr,
                        ndet, nsky,
                        det_fixed_src_norm_arr);
    
    // sky image to be reconstructed
    double* rho_init_arr = new double[nsky];
    double nu_init = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        rho_init_arr[isky] = 0.5 / nsky;
        nu_init = 0.5;
    }
    if("none" != argval->GetSkyfile()){
        printf("not implimented, yet.");
        abort();
        // load
        // reference sky image 
        //double* sky_ref_arr = new double[nsky];
        //MifImgInfo* img_info_sky = new MifImgInfo;
        //img_info_sky->InitSetImg(1, 1, nskyx, nskyy);
        //int bitpix_sky = 0;
        //MifFits::InFitsImageD(argval->GetSkyfile(), img_info_sky,
        //                      &bitpix_sky, &sky_ref_arr);
        //for(int isky = 0; isky < nsky; isky ++){
        //    sky_ref_arr[isky] *= eff_mat_arr[isky];
        //}
        //double N_val = MirMath::GetSum(nsky, sky_ref_arr);
        //double N_B = N_val + B_val;
        //for(int isky = 0; isky < nsky; isky ++){
        //    rho_init_arr[isky] = sky_ref_arr[isky] / N_B;
        // }
        //nu_init = B_val / N_B;
        //delete [] sky_ref_arr;
        //delete img_info_sky;
    }

    double* rho_new_arr = new double[nsky];
    double nu = 0.0;
    SrtlibRlBg2SmthEm::RichlucyBg2Smth_Acc(
        fp_log,
        rho_init_arr, nu_init,
        data_arr, det_fixed_src_norm_arr,
        resp_norm_mat_arr,
        ndet, nskyx, nskyy, argval->GetMu(),
        argval->GetOutdir(),
        argval->GetOutfileHead(),
        argval->GetNem(), argval->GetTolEm(),
        argval->GetNpm(), argval->GetTolPm(),
        argval->GetNnewton(), argval->GetTolNewton(),
        rho_new_arr, &nu);

    double N_B = MibBlas::Sum(data_arr, ndet);
    MiIolib::Printf2(fp_log, "N_B = %e\n", N_B);
    double B_val = nu * N_B;
    MiIolib::Printf2(fp_log, "B_val = %e\n", B_val);

    double sum_rho_new = MirMath::GetSum(nsky, rho_new_arr);
    MiIolib::Printf2(fp_log, "sum_rho_new = %e\n", sum_rho_new);
    
    // output reconstructed sky image: lambda for nebula
    double* sky_new_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_new_arr[isky] = rho_new_arr[isky] * N_B;
    }
    double sum_sky_new = MirMath::GetSum(nsky, sky_new_arr);
    MiIolib::Printf2(fp_log, "sum_sky_new = %e\n", sum_sky_new);

    // nebula + pulsar: sky_new_arr + B_val * sky_fixed_src_norm_arr
    double* sky_total_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_total_arr[isky] = sky_new_arr[isky]
            + B_val * sky_fixed_src_norm_arr[isky];
    }
    
    // div by eff_arr
    for(int isky = 0; isky < nsky; isky ++){
        sky_new_arr[isky] /= eff_mat_arr[isky];
        sky_total_arr[isky] /= eff_mat_arr[isky];
    }

    long naxes[2];
    naxes[0] = nskyx;
    naxes[1] = nskyy;
    int bitpix_out = -32;    
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "nebula", 2,
                           bitpix_out,
                           naxes, sky_new_arr);

    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "pulsar+nebula", 2,
                           bitpix_out,
                           naxes, sky_total_arr);

    double time_ed = MiTime::GetTimeSec();
    MiIolib::Printf2(fp_log, "duration = %e sec.\n", time_ed - time_st);

    fclose(fp_log);

    return status_prog;
}
