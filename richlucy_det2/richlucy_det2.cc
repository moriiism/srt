// Image reconstruction by Richardson-Lucy method
// for HXI1 + HXI2

#include "mi_time.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "arg_richlucy_det2.h"
#include "rl.h"
#include "rl_det2.h"
#include "srtmathlib.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValRichlucyDet2* argval = new ArgValRichlucyDet2;
    argval->Init(argc, argv);
    argval->Print(stdout);

    char logfile[kLineSize];
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        int ret = system(cmd);
        (void) ret;
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
    MifImgInfo* img_info_data1 = new MifImgInfo;
    img_info_data1->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_data1 = 0;
    double* data1_arr = NULL;
    MifFits::InFitsImageD(argval->GetDatafile1(), img_info_data1,
                          &bitpix_data1, &data1_arr);
    int nph_data1 = SrtMathlib::GetSum(ndet, data1_arr);
    MiIolib::Printf2(fp_log, "N photon (1) = %d\n", nph_data1);

    MifImgInfo* img_info_data2 = new MifImgInfo;
    img_info_data2->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_data2 = 0;
    double* data2_arr = NULL;
    MifFits::InFitsImageD(argval->GetDatafile2(), img_info_data2,
                          &bitpix_data2, &data2_arr);
    int nph_data2 = SrtMathlib::GetSum(ndet, data2_arr);
    MiIolib::Printf2(fp_log, "N photon (2) = %d\n", nph_data2);

    
    // load bg model data
    double* bg1_arr = NULL;
    if(argval->GetBgFile1() == "none"){
        int ndet = ndetx * ndety;
        bg1_arr = new double[ndet];
        for(int idet = 0; idet < ndet; idet++){
            bg1_arr[idet] = 0.0;
        }
    } else {
        MifImgInfo* img_info_bg1 = new MifImgInfo;
        img_info_bg1->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_bg1 = 0;
        MifFits::InFitsImageD(argval->GetBgFile1(), img_info_bg1,
                              &bitpix_bg1, &bg1_arr);
        int nph_bg1 = SrtMathlib::GetSum(ndet, bg1_arr);
        MiIolib::Printf2(fp_log, "N bg (1) = %d\n", nph_bg1);
        delete img_info_bg1;
    }

    double* bg2_arr = NULL;
    if(argval->GetBgFile2() == "none"){
        int ndet = ndetx * ndety;
        bg2_arr = new double[ndet];
        for(int idet = 0; idet < ndet; idet++){
            bg2_arr[idet] = 0.0;
        }
    } else {
        MifImgInfo* img_info_bg2 = new MifImgInfo;
        img_info_bg2->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_bg2 = 0;
        MifFits::InFitsImageD(argval->GetBgFile2(), img_info_bg2,
                              &bitpix_bg2, &bg2_arr);
        int nph_bg2 = SrtMathlib::GetSum(ndet, bg2_arr);
        MiIolib::Printf2(fp_log, "N bg (2) = %d\n", nph_bg2);
        delete img_info_bg2;
    }

    // load response file
    int naxis0 = MifFits::GetAxisSize(argval->GetRespNormFile1(), 0);
    int naxis1 = MifFits::GetAxisSize(argval->GetRespNormFile1(), 1);
    if ((naxis0 != ndet) || (naxis1 != nsky)){
        MiIolib::Printf2(fp_log,
                         "Error: response file size error.\n");
        abort();
    }

    double* resp_norm_mat_det1_arr = NULL;
    int bitpix_resp_norm_det1 = 0;
    MifImgInfo* img_info_resp_norm_det1 = new MifImgInfo;
    img_info_resp_norm_det1->InitSetImg(1, 1, ndet, nsky);
    MifFits::InFitsImageD(argval->GetRespNormFile1(),
                          img_info_resp_norm_det1,
                          &bitpix_resp_norm_det1,
                          &resp_norm_mat_det1_arr);
    double* resp_norm_mat_det2_arr = NULL;
    int bitpix_resp_norm_det2 = 0;
    MifImgInfo* img_info_resp_norm_det2 = new MifImgInfo;
    img_info_resp_norm_det2->InitSetImg(1, 1, ndet, nsky);
    MifFits::InFitsImageD(argval->GetRespNormFile2(),
                          img_info_resp_norm_det2,
                          &bitpix_resp_norm_det2,
                          &resp_norm_mat_det2_arr);
    
    // load efficiency file
    double* eff_mat_arr = NULL;
    int bitpix_eff = 0;
    MifImgInfo* img_info_eff = new MifImgInfo;
    img_info_eff->InitSetImg(1, 1, nskyx, nskyy);
    MifFits::InFitsImageD(argval->GetEffFile(), img_info_eff,
                          &bitpix_eff, &eff_mat_arr);

    // check
    //    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
    //        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
    //            int isky = nskyx * iskyy + iskyx;
    //            int imat = isky * ndet;            
    //            double resp_norm_sum = 0.0;
    //            for(int idet = 0; idet < ndet; idet ++){
    //                resp_norm_sum += resp_norm_mat_arr[imat + idet];
    //            }
    //            if ( fabs(resp_norm_sum - 1.0) > 1.0e-10){
    //                // printf("warning: resp_norm_sum = %e\n",
    //                // resp_norm_sum);
    //            }
    //        }
    //    }

    // sky image to be reconstructed
    double* sky_init_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_init_arr[isky] = float(nph_data1 + nph_data2) / nsky;
    }
    printf("nph_data1 = %d, nph_data2 = %d, nsky = %d\n",
           nph_data1, nph_data2, nsky);
    
    double* sky_new_arr = new double[nsky];
    if (argval->GetAccMethod() == "none"){
        SrtlibRlDet2::Richlucy(
            fp_log,
            sky_init_arr,
            data1_arr, data2_arr,
            bg1_arr, bg2_arr,
            resp_norm_mat_det1_arr,
            resp_norm_mat_det2_arr,
            ndet, nsky,
            argval->GetOutdir(), argval->GetOutfileHead(),
            argval->GetNem(), argval->GetTolEm(),
            sky_new_arr);
//    } else if (argval->GetAccMethod() == "kuroda"){
//        int k_restart = 2;
//        double delta_restart = 1.0;
//        SrtlibRl::RichlucyAccKuroda(
//            fp_log,
//            sky_init_arr,
//            data_arr, bg_arr,
//            resp_norm_mat_arr,
//            ndet, nsky,
//            argval->GetOutdir(), argval->GetOutfileHead(),
//            argval->GetNem(), argval->GetTolEm(),
//            k_restart, delta_restart,
//            sky_new_arr);
//    } else if (argval->GetAccMethod() == "squarem"){
//        SrtlibRl::RichlucyAccSquarem(
//            fp_log,
//            sky_init_arr,
//            data_arr, bg_arr,
//            resp_norm_mat_arr,
//            ndet, nsky,
//            argval->GetOutdir(), argval->GetOutfileHead(),
//            argval->GetNem(), argval->GetTolEm(),
//            sky_new_arr);
//    } else if (argval->GetAccMethod() == "ikeda"){
//        int rand_seed = 1;
//        SrtlibRl::RichlucyAccIkeda(
//            fp_log,
//            sky_init_arr,
//            data_arr, bg_arr,
//            resp_norm_mat_arr,
//            ndet, nsky,
//            argval->GetOutdir(), argval->GetOutfileHead(),
//            argval->GetNem(), argval->GetTolEm(),
//            nph_data, rand_seed,
//            sky_new_arr);
    } else {
        printf("bad acc_method\n");
        abort();
    }
    
    double sum_sky_new = SrtMathlib::GetSum(nsky, sky_new_arr);
    MiIolib::Printf2(fp_log, "sum_sky_new = %e\n", sum_sky_new);

    // div by eff_arr
    for(int isky = 0; isky < nsky; isky ++){
        sky_new_arr[isky] /= eff_mat_arr[isky];
    }
    
    long naxes[2];
    naxes[0] = nskyx;
    naxes[1] = nskyy;
    int bitpix_out = -32;
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "rec", 2,
                           bitpix_out,
                           naxes, sky_new_arr);

    double time_ed = MiTime::GetTimeSec();
    MiIolib::Printf2(fp_log,
                     "duration = %e sec.\n", time_ed - time_st);

    fclose(fp_log);
    
    return status_prog;
}
