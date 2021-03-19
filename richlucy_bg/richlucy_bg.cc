#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_richlucy_bg.h"
#include "sub_bg.h"
#include "TRandom3.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValRichlucyBg* argval = new ArgValRichlucyBg;
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
    MiIolib::Printf2(fp_log, "-----------------------------\n");
    argval->Print(fp_log);


    int nskyx = argval->GetNskyx();
    int nskyy = argval->GetNskyy();
    int ndetx = argval->GetNdetx();
    int ndety = argval->GetNdety();
    int nsky = nskyx * nskyy;
    int ndet = ndetx * ndety;    

    
    // load data
    MifImgInfo* img_info_data = new MifImgInfo;
    img_info_data->InitSetImg(1, 1, ndetx, ndety);
    int bitpix = 0;
    double* data_arr = NULL;
    MifFits::InFitsImageD(argval->GetDatafile(), img_info_data,
                          &bitpix, &data_arr);
    int nph = MirMath::GetSum(ndet, data_arr);
    printf("N photon = %d\n", nph);


    // sky image
    double* rho_arr = new double[nsky];
    double* rho_new_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        rho_arr[isky] = 1.0 / nsky;
        rho_new_arr[isky] = 1.0 / nsky;
    }
    if("none" != argval->GetSkyfile()){
        // load
        double* rho_ref_arr = new double[nsky];
        MifImgInfo* img_info_sky = new MifImgInfo;
        img_info_sky->InitSetImg(1, 1, nskyx, nskyy);
        int bitpix_sky = 0;
        MifFits::InFitsImageD(argval->GetSkyfile(), img_info_sky,
                              &bitpix_sky, &rho_ref_arr);
        double nph_ref = MirMath::GetSum(nsky, rho_ref_arr);
        for(int isky = 0; isky < nsky; isky ++){
            rho_arr[isky] = (rho_arr[isky] + rho_ref_arr[isky] / nph_ref ) / 2.0;
            // rho_arr[isky] = rho_ref_arr[isky] / nph_ref;
            rho_new_arr[isky] = rho_arr[isky];
        }
        delete [] rho_ref_arr;
        delete img_info_sky;
    }


    // load response file
    int naxis = 2;
    int* naxes_arr = new int[naxis];
    for(int iaxis = 0; iaxis < naxis; iaxis ++){
        naxes_arr[iaxis] = MifFits::GetAxisSize(argval->GetRespFile(), iaxis);
    }
    if ((naxes_arr[0] != ndet) || (naxes_arr[1] != nsky)){
        abort();
    }
    double* resp_mat_arr = NULL;
    int bitpix_resp = 0;
    MifImgInfo* img_info_resp = new MifImgInfo;
    img_info_resp->InitSetImg(1, 1, ndet, nsky);
    MifFits::InFitsImageD(argval->GetRespFile(), img_info_resp,
                          &bitpix_resp, &resp_mat_arr);

    // load efficiency file
    for(int iaxis = 0; iaxis < naxis; iaxis ++){
        naxes_arr[iaxis] = MifFits::GetAxisSize(argval->GetEffFile(), iaxis);
    }
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

    // load bg model
    bitpix = 0;
    double* bg_arr = NULL;
    MifImgInfo* img_info_bg = new MifImgInfo;
    img_info_bg->InitSetImg(1, 1, ndetx, ndety);
    MifFits::InFitsImageD(argval->GetBgfile(), img_info_bg,
                          &bitpix, &bg_arr);
    int nph_bg = MirMath::GetSum(ndet, bg_arr);
    printf("N bg = %d\n", nph_bg);


    double N_B = 0.0;
    bitpix = -32;
    RichlucyBg(rho_arr,
               data_arr, resp_norm_mat_arr, bg_arr,
               argval->GetNloopMain(),
               argval->GetNloopEm(),
               argval->GetNloopNewton(),
               argval->GetOutdir(), argval->GetOutfileHead(),
               ndet, nskyx, nskyy,
               argval->GetTolMain(),
               argval->GetTolEm(),               
               argval->GetTolNewton(),
               rho_new_arr, &N_B);

    //printf("nph = %d\n", nph);
    //double sum_image = MirMath::GetSum(nsky, rho_arr);
    //printf("sum_image = %e\n", sum_image);
 
    // output reconstructed sky image
    for(int isky = 0; isky < nsky; isky ++){
        rho_new_arr[isky] *= N_B;
    }
    // div by eff_arr
    for(int isky = 0; isky < nsky; isky ++){
        rho_new_arr[isky] /= eff_mat_arr[isky];
    }    
    naxis = 2;
    long* naxes = new long[naxis];
    naxes[0] = nskyx;
    naxes[1] = nskyy;
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "rec", 2,
                           bitpix,
                           naxes, rho_new_arr);

    double time_ed = MiTime::GetTimeSec();
    printf("duration = %e sec.\n", time_ed - time_st);
    double sum_image = MirMath::GetSum(nsky, rho_new_arr);
    printf("sum_image = %e\n", sum_image);
    
    return status_prog;
}
