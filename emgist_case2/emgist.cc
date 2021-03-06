#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_emgist.h"
#include "sub_emgist.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValEmgist* argval = new ArgValEmgist;
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

    // load response 
    int nskyx = 60;
    int nskyy = 60;
    int ndetx = 0;
    int ndety = 0;
    double* resp_mat_arr = NULL;
    LoadResp(argval->GetRespdir(), nskyx, nskyy,
             &resp_mat_arr, &ndetx, &ndety);
    int nsky = nskyx * nskyy;
    int ndet = ndetx * ndety;    

    // load data
    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, ndetx, ndety);
    int bitpix = 0;
    double* data_arr = NULL;
    MifFits::InFitsImageD(argval->GetDatafile(), img_info,
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

    bitpix = -32;
    SolveByEM(rho_arr, nph,
              data_arr, resp_mat_arr,
              argval->GetBeta(), argval->GetMu(), argval->GetLconst(),
              argval->GetNem(), argval->GetTolEm(), argval->GetNpm(), argval->GetTolPm(),
              argval->GetFlagLineSearch(),
              argval->GetOutdir(), argval->GetOutfileHead(),
              ndet, nskyx, nskyy, argval->GetEpsilon(),
              bitpix,
              rho_new_arr);

    // output reconstructed sky image
    for(int isky = 0; isky < nsky; isky ++){
        rho_new_arr[isky] *= nph;
    }
    int naxis = 2;
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
    
    return status_prog;
}
