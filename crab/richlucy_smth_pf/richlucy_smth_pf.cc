// Simultaneous image reconstruction by Richardson-Lucy method
// for multiple pulse phase data of Crab pulsar with Crab nebula
// under regularizations of smoothness for Crab nebula and
// pulse flux of Crab pulsar.

#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_richlucy_smth_pf.h"
#include "rl_crab.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValRichlucySmthPf* argval = new ArgValRichlucySmthPf;
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

    // load data list
    long nphase_long = 0;
    string* line_data_list_arr = NULL;
    MiIolib::GenReadFileSkipComment(argval->GetDataList(),
                                    &line_data_list_arr,
                                    &nphase_long);
    int nphase = nphase_long;
    printf("nphase = %d\n", nphase);
    string* data_list_arr = new string[nphase];
    string* phase_tag_arr = new string[nphase];
    double* phase_arr = new double[nphase];
    for(int iphase = 0; iphase < nphase; iphase ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_data_list_arr[iphase], &nsplit, &split_arr);
        data_list_arr[iphase] = split_arr[0];
        phase_tag_arr[iphase] = split_arr[1];
        phase_arr[iphase] = atof(split_arr[2].c_str());
        MiStr::DelSplit(split_arr);
    }
    MiIolib::DelReadFile(line_data_list_arr);

    for(int iphase = 0; iphase < nphase; iphase ++){
        printf("phase_arr[%d] = %e\n", iphase, phase_arr[iphase]);
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
        nph_data_arr[iphase] = MirMath::GetSum(ndet, data_arr[iphase]);
        MiIolib::Printf2(fp_log, "N photon = %d\n", nph_data_arr[iphase]);
        nph_data += nph_data_arr[iphase];
    }
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

    // load nu_0_file
    long nphase_long_nu_0 = 0;
    string* line_nu_0_arr = NULL;
    MiIolib::GenReadFileSkipComment(argval->GetNu0File(),
                                    &line_nu_0_arr,
                                    &nphase_long_nu_0);
    // here!
    





    
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
    for(int isky = 0; isky < nsky; isky ++){
        rho_init_arr[isky] = 1.0 / (nsky + nphase);
    }
    // nu to be reconstructed
    double* nu_init_arr = new double[nphase];
    for(int iphase = 0; iphase < nphase; iphase ++){
        nu_init_arr[iphase] = 1.0 / (nsky + nphase);
    }
    double* rho_new_arr = new double[nsky];
    double* nu_new_arr = new double[nphase];
    if (argval->GetAccMethod() == "none"){

        // SrtlibRlBg2SmthEm::RichlucyBg2Smth_Acc(
        SrtlibRlCrabSmthPf::RichlucyCrabSmthPf(
            fp_log,
            rho_init_arr,
            nu_init_arr,
            data_arr,
            nu_0_arr,
            phase_arr,
            det_fixed_src_norm_arr,
            resp_norm_mat_arr,
            ndet, nsky, nphase,
            argval->GetOutdir(),
            argval->GetOutfileHead(),
            argval->GetNem(), argval->GetTolEm(),
            argval->GetNpm(), argval->GetTolPm(),
            argval->GetNnewton(), argval->GetTolNewton(),
            rho_new_arr, nu_new_arr);
    } else if (argval->GetAccMethod() == "squarem"){
        SrtlibRlCrab::RichlucyCrabAccSquarem(
            fp_log,
            rho_init_arr,
            nu_init_arr,
            data_arr,
            phase_arr,
            det_fixed_src_norm_arr,
            resp_norm_mat_arr,
            ndet, nsky, nphase,
            argval->GetOutdir(),
            argval->GetOutfileHead(),
            argval->GetNloop(),
            argval->GetTol(),
            rho_new_arr, nu_new_arr);
    } else {
        printf("bad acc_method\n");
        abort();
    }
    
    // output reconstructed sky image
    double* sky_new_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_new_arr[isky] = rho_new_arr[isky] * nph_data;
    }
    double sum_sky_new = MirMath::GetSum(nsky, sky_new_arr);
    MiIolib::Printf2(fp_log, "sum_sky_new = %e\n", sum_sky_new);

    // output reconstructed flux
    double* flux_pulsar_arr = new double[nphase];
    for(int iphase = 0; iphase < nphase; iphase ++){
        flux_pulsar_arr[iphase] = nu_new_arr[iphase] * nph_data
            / phase_arr[iphase];
        MiIolib::Printf2(fp_log, "flux_pulsar[%d] = %e\n",
                         iphase, flux_pulsar_arr[iphase]);
    }

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
    MiIolib::Printf2(fp_log, "duration = %e sec.\n", time_ed - time_st);

    fclose(fp_log);
    
    return status_prog;
}
