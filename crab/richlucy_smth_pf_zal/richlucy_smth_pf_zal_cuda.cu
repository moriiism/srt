// Simultaneous image reconstruction by Richardson-Lucy method
// for multiple pulse phase data of Crab pulsar with Crab nebula
// under regularizations of smoothness for Crab nebula and
// pulse flux of Crab pulsar. Non X-ray background is also considered.

#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "srtmathlib.h"
#include "arg_richlucy_smth_pf_zal.h"
#include "rl_crab.h"
#include "rl_crab_cuda.h"
#include "rl_crab_smth_pf_zal_cuda.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValRichlucySmthPfZal* argval = new ArgValRichlucySmthPfZal;
    argval->Init(argc, argv);
    argval->Print(stdout);

    char logfile[kLineSize];
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        int ret = system(cmd);
        (void) ret;
    }
    sprintf(logfile, "%s/%s.log",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    FILE* fp_log = fopen(logfile, "w");
    argval->Print(fp_log);
    MiIolib::Printf2(fp_log, "-----------------------------\n");    

    int nskyx = argval->GetNskyx();
    int nskyy = argval->GetNskyy();
    int ndetx = argval->GetNdetx();
    int ndety = argval->GetNdety();
    int nsky = nskyx * nskyy;
    int ndet = ndetx * ndety;    

    // load data list
    long nphase_long = 0;
    string* line_data_list_arr = NULL;
    MiIolib::GenReadFileSkipComment(argval->GetDataList(),
                                    &line_data_list_arr,
                                    &nphase_long);
    int nphase = nphase_long;
    printf("nphase = %d\n", nphase);
    string* data_list_arr = new string[nphase];
    string* data_vl_list_arr = new string[nphase];    
    string* phase_tag_arr = new string[nphase];
    double* phase_arr = new double[nphase];
    double* live_time_ratio_arr = new double[nphase];
    double* flux_target_arr = new double[nphase];
    for(int iphase = 0; iphase < nphase; iphase ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_data_list_arr[iphase],
                        &nsplit, &split_arr);
        if(nsplit != 6){
            printf("error: bad nsplit(=%d)\n", nsplit);
            abort();
        }
        data_list_arr[iphase] = split_arr[0];
        data_vl_list_arr[iphase] = split_arr[1];      
        phase_tag_arr[iphase] = split_arr[2];
        phase_arr[iphase] = atof(split_arr[3].c_str());
        live_time_ratio_arr[iphase] = atof(split_arr[4].c_str());
        flux_target_arr[iphase] = atof(split_arr[5].c_str());
        MiStr::DelSplit(split_arr);
    }
    MiIolib::DelReadFile(line_data_list_arr);

    for(int iphase = 0; iphase < nphase; iphase ++){
        printf("phase_arr[%d] = %e\n",
               iphase, phase_arr[iphase]);
    }
    for(int iphase = 0; iphase < nphase; iphase ++){
        printf("flux_target_arr[%d] = %e\n",
               iphase, flux_target_arr[iphase]);
    }
    
    // load image data
    double** data_arr = new double*[nphase];
    int* nph_data_arr = new int[nphase];
    int nph_data = 0;
    for(int iphase = 0; iphase < nphase; iphase ++){
        MifImgInfo* img_info_data = new MifImgInfo;
        img_info_data->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_data = 0;
        MifFits::InFitsImageD(data_list_arr[iphase], img_info_data,
                              &bitpix_data, &data_arr[iphase]);
        nph_data_arr[iphase] = SrtMathlib::GetSum(
            ndet, data_arr[iphase]);
        MiIolib::Printf2(fp_log, "N photon = %d\n",
                         nph_data_arr[iphase]);
        nph_data += nph_data_arr[iphase];
    }
    MiIolib::Printf2(fp_log, "N photon = %d\n", nph_data);

    // load bg model data
    double* bg_arr = NULL;
    if(argval->GetBgFile() == "none"){
        int ndet = ndetx * ndety;
        bg_arr = new double[ndet];
        for(int idet = 0; idet < ndet; idet++){
            bg_arr[idet] = 0.0;
        }
    } else {
        MifImgInfo* img_info_bg = new MifImgInfo;
        img_info_bg->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_bg = 0;
        MifFits::InFitsImageD(argval->GetBgFile(), img_info_bg,
                              &bitpix_bg, &bg_arr);
        int nph_bg = SrtMathlib::GetSum(ndet, bg_arr);
        MiIolib::Printf2(fp_log, "N bg = %d\n", nph_bg);
        delete img_info_bg;
    }
    
    // load sky image of fixed source with normalized flux
    MifImgInfo* img_info_fixed_src_norm = new MifImgInfo;
    img_info_fixed_src_norm->InitSetImg(1, 1, nskyx, nskyy);
    int bitpix_fixed_src_norm = 0;
    double* sky_fixed_src_norm_arr = NULL;
    MifFits::InFitsImageD(argval->GetFixedSrcNormFile(),
                          img_info_fixed_src_norm,
                          &bitpix_fixed_src_norm,
                          &sky_fixed_src_norm_arr);
    int nph_fixed_src_norm = SrtMathlib::GetSum(
        nsky, sky_fixed_src_norm_arr);
    MiIolib::Printf2(
        fp_log,
        "N photon fixed source with normalized flux = %d\n",
        nph_fixed_src_norm);

    
    // load response file
    int naxis0 = MifFits::GetAxisSize(argval->GetRespNormFile(), 0);
    int naxis1 = MifFits::GetAxisSize(argval->GetRespNormFile(), 1);
    if ((naxis0 != ndet) || (naxis1 != nsky)){
        MiIolib::Printf2(
            fp_log,
            "Error: normalized response file size error.\n");
        abort();
    }
    double* resp_norm_mat_arr = NULL;
    int bitpix_resp = 0;
    MifImgInfo* img_info_resp = new MifImgInfo;
    img_info_resp->InitSetImg(1, 1, ndet, nsky);
    MifFits::InFitsImageD(argval->GetRespNormFile(),
                          img_info_resp,
                          &bitpix_resp,
                          &resp_norm_mat_arr);

    // load efficiency file
    double* eff_mat_arr = NULL;
    int bitpix_eff = 0;
    MifImgInfo* img_info_eff = new MifImgInfo;
    img_info_eff->InitSetImg(1, 1, nskyx, nskyy);
    MifFits::InFitsImageD(argval->GetEffFile(), img_info_eff,
                          &bitpix_eff, &eff_mat_arr);

    // check response file
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;            
            double resp_norm_sum = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                resp_norm_sum += resp_norm_mat_arr[imat + idet];
            }
            if ( fabs(resp_norm_sum - 1.0) > 1.0e-10){
                //printf("warning: resp_norm_sum = %e\n",
                //       resp_norm_sum);
            }
        }
    }

    // det image of fixed source with normalized flux
    double* det_fixed_src_norm_arr = new double[ndet];
    SrtlibRlCrab::GetDetArr(sky_fixed_src_norm_arr,
                            resp_norm_mat_arr,
                            ndet, nsky,
                            det_fixed_src_norm_arr);
    
    // nebula sky image to be reconstructed
    double* sky_init_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_init_arr[isky] = nph_data / (2.0 * nsky);
    }
    // flux to be reconstructed
    double* flux_init_arr = new double[nphase];
    for(int iphase = 0; iphase < nphase; iphase ++){
        flux_init_arr[iphase] = nph_data / (2.0 * nphase);
    }
    double* sky_new_arr = new double[nsky];
    double* flux_new_arr = new double[nphase];
    if (argval->GetAccMethod() == "none"){
        SrtlibRlCrabSmthPfZalCuda::RichlucyCrabSmthPfZal(
            fp_log,
            sky_init_arr,
            flux_init_arr,
            data_arr,
            bg_arr,
            flux_target_arr,
            phase_arr,
            live_time_ratio_arr,
            det_fixed_src_norm_arr,
            resp_norm_mat_arr,
            ndet, nskyx, nskyy, nphase,
            argval->GetMu(), argval->GetGamma(),
            argval->GetOutdir(),
            argval->GetOutfileHead(),
            argval->GetNem(), argval->GetTolEm(),
            sky_new_arr, flux_new_arr);
//    } else if (argval->GetAccMethod() == "zalq1"){
//        SrtlibRlCrabSmthPfZalCuda::RichlucyCrabSmthPfZalQ1(
//            fp_log,
//            sky_init_arr,
//            flux_init_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_fixed_src_norm_arr,
//            resp_norm_mat_arr,
//            ndet, nskyx, nskyy, nphase,
//            argval->GetMu(), argval->GetGamma(),
//            argval->GetOutdir(),
//            argval->GetOutfileHead(),
//            argval->GetNem(), argval->GetTolEm(),
//            sky_new_arr, flux_new_arr);
//    } else if (argval->GetAccMethod() == "sqs3"){
//        SrtlibRlCrabSmthPfZalCuda::RichlucyCrabSmthPfSqS3(
//            fp_log,
//            sky_init_arr,
//            flux_init_arr,
//            data_arr,
//            bg_arr,
//            flux_target_arr,
//            phase_arr,
//            det_fixed_src_norm_arr,
//            resp_norm_mat_arr,
//            ndet, nskyx, nskyy, nphase,
//            argval->GetMu(), argval->GetGamma(),
//            argval->GetOutdir(),
//            argval->GetOutfileHead(),
//            argval->GetNem(), argval->GetTolEm(),
//            sky_new_arr, flux_new_arr);
    } else {
        printf("bad acc_method\n");
        abort();
    }
    
    // output reconstructed sky image of nebula
    double sum_sky = SrtMathlib::GetSum(nsky, sky_new_arr);
    MiIolib::Printf2(fp_log, "sum_sky = %e\n", sum_sky);

    // output reconstructed flux
    for(int iphase = 0; iphase < nphase; iphase ++){
        MiIolib::Printf2(fp_log, "flux [%d] = %e\n",
                         iphase, flux_new_arr[iphase]);
    }
    // pulse profile
    char qdp_file[kLineSize];
    sprintf(qdp_file, "%s/%s_pulseprof.qdp",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    FILE* fp_qdp = NULL;
    fp_qdp = fopen(qdp_file, "w");
    fprintf(fp_qdp, "skip sing\n");
    fprintf(fp_qdp, "! flux_new__arr\n");
    for(int iphase = 0; iphase < nphase; iphase ++){
        fprintf(fp_qdp, "%d  %e\n",
                iphase, flux_new_arr[iphase]);
    }
    fprintf(fp_qdp, "\n");
    fprintf(fp_qdp, "no\n");
    fprintf(fp_qdp, "\n");
    fprintf(fp_qdp, "! flux_target_arr\n");
    for(int iphase = 0; iphase < nphase; iphase ++){
        fprintf(fp_qdp, "%d  %e\n",
                iphase, flux_target_arr[iphase]);
    }
    fprintf(fp_qdp, "\n");
    fprintf(fp_qdp, "la file\n");
    fprintf(fp_qdp, "time off\n");    
    fprintf(fp_qdp, "lw 5\n");
    fprintf(fp_qdp, "csize 1.2\n");
    fprintf(fp_qdp, "la rot \n");
    fprintf(fp_qdp, "loc 0.05 0.05 0.95 0.95\n");
    fprintf(fp_qdp, "la pos y 3.0\n");
    fprintf(fp_qdp, "\n");
    fclose(fp_qdp);

    // pulse phase resolved reconstructed
    // sky image of nebula + pulsar
    double** sky_pulse_arr = new double*[nphase];
    for(int iphase = 0; iphase < nphase; iphase ++){
        sky_pulse_arr[iphase] = new double[nsky];
        for(int isky = 0; isky < nsky; isky ++){
            sky_pulse_arr[iphase][isky] = sky_new_arr[isky]
                + flux_new_arr[iphase]
                * sky_fixed_src_norm_arr[isky];
        }
    }

    // if(nfold > 1): cross validation
    // otherwise:     no cross validation
    if(argval->GetNfoldCv() > 1){ 
        // evaluate resultant images with validation images
        // load validation image data
        double** data_vl_arr = new double*[nphase];
        int* nph_data_vl_arr = new int[nphase];
        int nph_data_vl = 0;
        for(int iphase = 0; iphase < nphase; iphase ++){
            MifImgInfo* img_info_data_vl = new MifImgInfo;
            img_info_data_vl->InitSetImg(1, 1, ndetx, ndety);
            int bitpix_data_vl = 0;
            MifFits::InFitsImageD(data_vl_list_arr[iphase],
                                  img_info_data_vl,
                                  &bitpix_data_vl,
                                  &data_vl_arr[iphase]);
            nph_data_vl_arr[iphase] = MirMath::GetSum(
                ndet, data_vl_arr[iphase]);
            MiIolib::Printf2(fp_log, "N photon (vl) = %d\n",
                             nph_data_vl_arr[iphase]);
            nph_data_vl += nph_data_vl_arr[iphase];
        }
        MiIolib::Printf2(fp_log, "N photon (vl) = %d\n", nph_data_vl);

        double num_rmse = 0.0;
        double den_rmse = 0.0;
        //   det images of resultant images
        double** det_pulse_arr = new double*[nphase];
        for(int iphase = 0; iphase < nphase; iphase ++){
            det_pulse_arr[iphase] = new double[ndet];
            SrtlibRlCrab::GetDetArr(sky_pulse_arr[iphase],
                                    resp_norm_mat_arr,
                                    ndet, nsky,
                                    det_pulse_arr[iphase]);
            // add non X-ray background
            daxpy_(ndet, 1.0, bg_arr, 1, det_pulse_arr[iphase], 1);
            // multiply phase ratio
            dscal_(ndet, phase_arr[iphase],
                   det_pulse_arr[iphase], 1);
            // scale det_pulse_arr by 1.0/(nfold - 1) for
            // evalution between validation data
            dscal_(ndet, 1.0 / (argval->GetNfoldCv() - 1),
                   det_pulse_arr[iphase], 1);
            double rmse = SrtMathlib::GetRootMeanSquareError(
                det_pulse_arr[iphase], data_vl_arr[iphase], ndet);
            printf("iphase = %d, rmse = %e\n", iphase, rmse);
            num_rmse += rmse * rmse * ndet;
            den_rmse += ndet;
        }
        double rmse_tot = sqrt(num_rmse / den_rmse);
        printf("rmse_tot = %e\n", rmse_tot);

        char outfile[kLineSize];
        sprintf(outfile, "%s/%s_rmse.txt",
                argval->GetOutdir().c_str(),
                argval->GetOutfileHead().c_str());
        FILE* fp_out = fopen(outfile, "w");
        fprintf(fp_out, "%e\n", rmse_tot);
        fclose(fp_out);
    } else {
        printf("nfold <= 1, then no cross validation.\n");
    }
    
    // div by eff_arr
    for(int isky = 0; isky < nsky; isky ++){
        sky_new_arr[isky] /= eff_mat_arr[isky];
    }
    for(int iphase = 0; iphase < nphase; iphase ++){
        for(int isky = 0; isky < nsky; isky ++){
            sky_pulse_arr[iphase][isky] /= eff_mat_arr[isky];
        }
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
    for(int iphase = 0; iphase < nphase; iphase ++){
        char tag[kLineSize];
        sprintf(tag, "rec_%2.2d", iphase);
        MifFits::OutFitsImageD(argval->GetOutdir(),
                               argval->GetOutfileHead(),
                               tag, 2,
                               bitpix_out,
                               naxes, sky_pulse_arr[iphase]);
    }

    double time_ed = MiTime::GetTimeSec();
    MiIolib::Printf2(fp_log, "duration = %e sec.\n",
                     time_ed - time_st);
    fclose(fp_log);
    return status_prog;
}
