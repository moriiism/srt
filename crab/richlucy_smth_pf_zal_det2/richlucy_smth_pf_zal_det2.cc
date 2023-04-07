// Simultaneous image reconstruction by Richardson-Lucy method
// for multiple pulse phase data of Crab pulsar with Crab nebula
// under regularizations of smoothness for Crab nebula and
// pulse flux of Crab pulsar.
// Non X-ray background is also considered.
// HXI1 and HXI2 are used at once.

#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "srtmathlib.h"
#include "arg_richlucy_smth_pf_zal.h"
#include "rl_crab.h"
#include "rl_crab_smth_pf_zal.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValRichlucySmthPfZalDet2* argval = new ArgValRichlucySmthPfZalDet2;
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

    // load data list (HXI1 + HXI2)
    long nphase_long = 0;
    string* line_data_list_arr = NULL;
    MiIolib::GenReadFileSkipComment(argval->GetDataList(),
                                    &line_data_list_arr,
                                    &nphase_long);
    int nphase = nphase_long;
    printf("nphase = %d\n", nphase);
    string* data_list_det1_arr = new string[nphase];
    string* data_vl_list_det1_arr = new string[nphase];
    string* data_list_det2_arr = new string[nphase];
    string* data_vl_list_det2_arr = new string[nphase];
    string* phase_tag_arr = new string[nphase];
    double* phase_arr = new double[nphase];
    double* live_time_ratio_det1_arr = new double[nphase];
    double* live_time_ratio_det2_arr = new double[nphase];
    double* flux_target_arr = new double[nphase];
    double* flux_target_err_arr = new double[nphase];    
    for(int iphase = 0; iphase < nphase; iphase ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_data_list_arr[iphase],
                        &nsplit, &split_arr);
        if(nsplit != 10){
            printf("error: bad nsplit(=%d)\n", nsplit);
            abort();
        }
        data_list_det1_arr[iphase] = split_arr[0];
        data_vl_list_det1_arr[iphase] = split_arr[1];
        data_list_det2_arr[iphase] = split_arr[2];
        data_vl_list_det2_arr[iphase] = split_arr[3]; 
        phase_tag_arr[iphase] = split_arr[4];
        phase_arr[iphase] = atof(split_arr[5].c_str());
        live_time_ratio_det1_arr[iphase] = atof(split_arr[6].c_str());
        live_time_ratio_det2_arr[iphase] = atof(split_arr[7].c_str());        
        flux_target_arr[iphase] = atof(split_arr[8].c_str());
        flux_target_err_arr[iphase] = atof(split_arr[9].c_str());
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
    for(int iphase = 0; iphase < nphase; iphase ++){
        printf("flux_target_err_arr[%d] = %e\n",
               iphase, flux_target_err_arr[iphase]);
    }
    

    // load image data (HXI1, HXI2)
    double** data1_arr = new double*[nphase];
    double** data2_arr = new double*[nphase];
    int* nph_data1_arr = new int[nphase];
    int* nph_data2_arr = new int[nphase];
    int nph_data = 0;
    for(int iphase = 0; iphase < nphase; iphase ++){
        MifImgInfo* img_info_data = new MifImgInfo;
        img_info_data->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_data = 0;
        MifFits::InFitsImageD(data_list_det1_arr[iphase], img_info_data,
                              &bitpix_data, &data1_arr[iphase]);
        MifFits::InFitsImageD(data_list_det2_arr[iphase], img_info_data,
                              &bitpix_data, &data2_arr[iphase]);
        nph_data1_arr[iphase] = SrtMathlib::GetSum(
            ndet, data1_arr[iphase]);
        nph_data2_arr[iphase] = SrtMathlib::GetSum(
            ndet, data2_arr[iphase]);        
        MiIolib::Printf2(fp_log, "N photon (det1) = %d\n",
                         nph_data1_arr[iphase]);
        MiIolib::Printf2(fp_log, "N photon (det2) = %d\n",
                         nph_data2_arr[iphase]);        
        nph_data += nph_data1_arr[iphase];
        nph_data += nph_data2_arr[iphase];
    }
    MiIolib::Printf2(fp_log, "N photon = %d\n", nph_data);

    // load bg model data (HXI1)
    double* bg1_arr = NULL;
    if(argval->GetBgFile1() == "none"){
        int ndet = ndetx * ndety;
        bg1_arr = new double[ndet];
        for(int idet = 0; idet < ndet; idet++){
            bg1_arr[idet] = 0.0;
        }
    } else {
        MifImgInfo* img_info_bg = new MifImgInfo;
        img_info_bg->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_bg = 0;
        MifFits::InFitsImageD(argval->GetBgFile1(), img_info_bg,
                              &bitpix_bg, &bg1_arr);
        int nph_bg1 = SrtMathlib::GetSum(ndet, bg1_arr);
        MiIolib::Printf2(fp_log, "N bg (det1) = %d\n", nph_bg1);
        delete img_info_bg;
    }

    // load bg model data (HXI2)
    double* bg2_arr = NULL;
    if(argval->GetBgFile2() == "none"){
        int ndet = ndetx * ndety;
        bg2_arr = new double[ndet];
        for(int idet = 0; idet < ndet; idet++){
            bg2_arr[idet] = 0.0;
        }
    } else {
        MifImgInfo* img_info_bg = new MifImgInfo;
        img_info_bg->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_bg = 0;
        MifFits::InFitsImageD(argval->GetBgFile2(), img_info_bg,
                              &bitpix_bg, &bg2_arr);
        int nph_bg2 = SrtMathlib::GetSum(ndet, bg2_arr);
        MiIolib::Printf2(fp_log, "N bg (det2) = %d\n", nph_bg2);
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

    
    // load normalized response file
    int naxis0 = MifFits::GetAxisSize(argval->GetRespNormFile1(), 0);
    int naxis1 = MifFits::GetAxisSize(argval->GetRespNormFile1(), 1);
    if ((naxis0 != ndet) || (naxis1 != nsky)){
        MiIolib::Printf2(
            fp_log,
            "Error: normalized response file size error.\n");
        abort();
    }
    double* resp_norm_mat_det1_arr = NULL;
    double* resp_norm_mat_det2_arr = NULL;    
    int bitpix_resp = 0;
    MifImgInfo* img_info_resp = new MifImgInfo;
    img_info_resp->InitSetImg(1, 1, ndet, nsky);
    MifFits::InFitsImageD(argval->GetRespNormFile1(),
                          img_info_resp,
                          &bitpix_resp,
                          &resp_norm_mat_det1_arr);
    MifFits::InFitsImageD(argval->GetRespNormFile2(),
                          img_info_resp,
                          &bitpix_resp,
                          &resp_norm_mat_det2_arr);
    
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
                resp_norm_sum += resp_norm_mat_det1_arr[imat + idet];
                resp_norm_sum += resp_norm_mat_det2_arr[imat + idet];
            }
            if ( fabs(resp_norm_sum - 1.0) > 1.0e-10){
                //printf("warning: resp_norm_sum = %e\n",
                //       resp_norm_sum);
            }
        }
    }

    // det image of fixed source with normalized flux (HXI1)
    double* det_fixed_src_norm_det1_arr = new double[ndet];
    SrtlibRlCrab::GetDetArr(sky_fixed_src_norm_arr,
                            resp_norm_mat_det1_arr,
                            ndet, nsky,
                            det_fixed_src_norm_det1_arr);

    // det image of fixed source with normalized flux (HXI2)
    double* det_fixed_src_norm_det2_arr = new double[ndet];
    SrtlibRlCrab::GetDetArr(sky_fixed_src_norm_arr,
                            resp_norm_mat_det2_arr,
                            ndet, nsky,
                            det_fixed_src_norm_det2_arr);
    
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
        SrtlibRlCrabSmthPfZalDet2::RichlucyCrabSmthPfZalDet2(
            fp_log,
            sky_init_arr,
            flux_init_arr,
            data1_arr,
            data2_arr,
            bg1_arr,
            bg2_arr,
            flux_target_arr,
            phase_arr,
            live_time_ratio_det1_arr,
            live_time_ratio_det2_arr,
            det_fixed_src_norm_det1_arr,
            det_fixed_src_norm_det2_arr,
            resp_norm_mat_det1_arr,
            resp_norm_mat_det2_arr,
            ndet, nskyx, nskyy, nphase,
            argval->GetMu(), argval->GetGamma(),
            argval->GetOutdir(),
            argval->GetOutfileHead(),
            argval->GetNem(), argval->GetTolEm(),
            sky_new_arr, flux_new_arr);
//    } else if (argval->GetAccMethod() == "zalq1"){
//        SrtlibRlCrabSmthPfZal::RichlucyCrabSmthPfZalQ1(
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
//        SrtlibRlCrabSmthPfZal::RichlucyCrabSmthPfSqS3(
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
    fprintf(fp_qdp, "read serr 2\n");    
    fprintf(fp_qdp, "! flux_new_arr\n");
    for(int iphase = 0; iphase < nphase; iphase ++){
        fprintf(fp_qdp, "%d  %e  0.0\n",
                iphase, flux_new_arr[iphase]);
    }
    fprintf(fp_qdp, "\n");
    fprintf(fp_qdp, "no\n");
    fprintf(fp_qdp, "\n");
    fprintf(fp_qdp, "! flux_target_arr\n");
    for(int iphase = 0; iphase < nphase; iphase ++){
        fprintf(fp_qdp, "%d  %e  %e\n",
                iphase, flux_target_arr[iphase],
                flux_target_err_arr[iphase]);
    }
    fprintf(fp_qdp, "\n");
    fprintf(fp_qdp, "la file\n");
    fprintf(fp_qdp, "time off\n");    
    fprintf(fp_qdp, "lw 7\n");
    fprintf(fp_qdp, "csize 1.2\n");
    fprintf(fp_qdp, "la rot \n");
    fprintf(fp_qdp, "loc 0.05 0.05 0.95 0.95\n");
    fprintf(fp_qdp, "la pos y 3.0\n");
    fprintf(fp_qdp, "ma 6 on\n");
    fprintf(fp_qdp, "line on\n");    
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

    // pulse phase integrated reconstructed
    // sky image of nebula + pulsar
    double* sky_integ_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_integ_arr[isky] = 0.0;
        for(int iphase = 0; iphase < nphase; iphase ++){
            sky_integ_arr[isky] += phase_arr[iphase] *
                sky_pulse_arr[iphase][isky];
        }
    }


    ////----------------------
    
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
            nph_data_vl_arr[iphase] = SrtMathlib::GetSum(
                ndet, data_vl_arr[iphase]);
            MiIolib::Printf2(fp_log, "N photon (vl) = %d\n",
                             nph_data_vl_arr[iphase]);
            nph_data_vl += nph_data_vl_arr[iphase];
        }
        MiIolib::Printf2(fp_log, "N photon (vl) = %d\n",
                         nph_data_vl);

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
            // multiply phase_ratio * live_time_ratio
            dscal_(ndet, phase_arr[iphase]
                   * live_time_ratio_arr[iphase],
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
    for(int isky = 0; isky < nsky; isky ++){
        sky_integ_arr[isky] /= eff_mat_arr[isky];
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
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "rec_integ", 2,
                           bitpix_out,
                           naxes, sky_integ_arr);
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
