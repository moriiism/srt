// Image reconstruction by Richardson-Lucy method for
// diffuse component with smoothness prior,
// with fixed point sources, and with non-X-ray background.

#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_richlucy_fpsrc_smth_bg.h"
#include "fpsrc.h"
#include "TRandom3.h"
#include "fpsrc_smth_bg_mm_em.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;

    double time_st = MiTime::GetTimeSec();
    
    ArgValRichlucyFpsrcSmthBg* argval = new ArgValRichlucyFpsrcSmthBg;
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

    // b_v^{(k)}
    double** det_fpsrc_arr = NULL;
    int* xpos_src_arr = NULL;
    int* ypos_src_arr = NULL;
    int nsrc = 0;
    GenFixedPointSrcDetImg(fp_log,
                           argval->GetFixedSrcList(),
                           resp_norm_mat_arr,
                           nskyx, nskyy, ndet,
                           &nsrc,
                           &xpos_src_arr,
                           &ypos_src_arr,
                           &det_fpsrc_arr);
    
    // sky image to be reconstructed
    double* rho_init_arr = new double[nsky];
    double* nu_init_arr = new double[nsrc];
    double phi_init = 1.0/3.0;
    for(int isky = 0; isky < nsky; isky ++){
        rho_init_arr[isky] = (1.0/3.0)/nsky;
    }
    for(int isrc = 0; isrc < nsrc; isrc ++){
        nu_init_arr[isrc] = (1.0/3.0)/nsrc;
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

    // load bg data
    int bitpix_bg = 0;
    double* bg_arr = NULL;
    MifImgInfo* img_info_bg = new MifImgInfo;
    img_info_bg->InitSetImg(1, 1, ndetx, ndety);
    MifFits::InFitsImageD(argval->GetBgfile(), img_info_bg,
                          &bitpix_bg, &bg_arr);
    int nph_bg = MirMath::GetSum(ndet, bg_arr);
    MiIolib::Printf2(fp_log, "N bg = %d\n", nph_bg);

    
    double* rho_new_arr = new double[nsky];
    double* nu_new_arr = new double[nsrc];
    double phi_new = 0.0;
//    RichlucyFpsrcSmthBgMM(fp_log,
//                          rho_init_arr, nu_init_arr, phi_init,
//                          data_arr, bg_arr, det_fpsrc_arr,
//                          resp_norm_mat_arr,
//                          ndet, nskyx, nskyy, nsrc,
//                          argval->GetMu(),
//                          argval->GetOutdir(),
//                          argval->GetOutfileHead(),
//                          argval->GetNem(), argval->GetTolEm(),
//                          argval->GetNpm(), argval->GetTolPm(),
//                          argval->GetNnewton(), argval->GetTolNewton(),
//                          rho_new_arr, nu_new_arr, &phi_new);
    RichlucyFpsrcSmthBgMM_Acc(fp_log,
                              rho_init_arr, nu_init_arr, phi_init,
                              data_arr, bg_arr, det_fpsrc_arr,
                              resp_norm_mat_arr,
                              ndet, nskyx, nskyy, nsrc,
                              argval->GetMu(),
                              argval->GetOutdir(),
                              argval->GetOutfileHead(),
                              argval->GetNem(), argval->GetTolEm(),
                              argval->GetNpm(), argval->GetTolPm(),
                              argval->GetNnewton(), argval->GetTolNewton(),
                              rho_new_arr, nu_new_arr, &phi_new);
    
    double B_val = MibBlas::Sum(bg_arr, ndet);
    MiIolib::Printf2(fp_log, "B_val = %e\n", B_val);
    double N_prime = B_val / phi_new;
    MiIolib::Printf2(fp_log, "N_prime = %e\n", N_prime); 
    
    // output reconstructed sky image: lambda for diffuse
    double* sky_new_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_new_arr[isky] = rho_new_arr[isky] * N_prime;
    }
    double sum_sky_new = MirMath::GetSum(nsky, sky_new_arr);
    MiIolib::Printf2(fp_log, "sum_sky_new = %e\n", sum_sky_new);

    double* flux_arr = new double[nsrc];
    for(int isrc = 0; isrc < nsrc; isrc++){
        flux_arr[isrc] = nu_new_arr[isrc] * N_prime;
    }
    
    // diffuse + fpsrc: sky_new_arr
    // + sum_k flux_k * 
    double* sky_total_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_total_arr[isky] = sky_new_arr[isky];
    }

    // debug
    for(int isrc = 0; isrc < nsrc; isrc ++){
        int isky_pos = ypos_src_arr[isrc] * nskyx + xpos_src_arr[isrc];
        sky_total_arr[isky_pos] += flux_arr[isrc];
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
                           "diffuse", 2,
                           bitpix_out,
                           naxes, sky_new_arr);

    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "diffuse+fpsrc", 2,
                           bitpix_out,
                           naxes, sky_total_arr);

    double time_ed = MiTime::GetTimeSec();
    MiIolib::Printf2(fp_log, "duration = %e sec.\n", time_ed - time_st);
    
    fclose(fp_log);
    
    return status_prog;
}
