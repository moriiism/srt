//
// mkobs_cv.cc
//
// make observation detector images as
// training and validation detector images for cross-validation
//
// detector image --> randomized event list --> nfold images
//

#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkobs_cv.h"
#include "sim.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMkobsCv* argval = new ArgValMkobsCv;
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
    
    // get ndetx, ndety
    int ndetx = MifFits::GetAxisSize(argval->GetOrgfile(), 0);
    int ndety = MifFits::GetAxisSize(argval->GetOrgfile(), 1);
    int ndet = ndetx * ndety;

    // load orgfile
    double* org_arr = NULL;
    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, ndetx, ndety);
    int bitpix = 0;
    MifFits::InFitsImageD(argval->GetOrgfile(), img_info,
                          &bitpix, &org_arr);
    delete img_info;
    double nph = MirMath::GetSum(ndet, org_arr);
    printf("nph = %e\n", nph);
    int nph_int = ceil(nph);

    // norm orgfile
    double* org_norm_arr = new double [ndet];
    for(int idet = 0; idet < ndet; idet ++){
        org_norm_arr[idet] = org_arr[idet] / nph;
    }
    delete [] org_arr;

    double* obs_bin_arr = new double[ndet];
    int*    obs_evt_arr = new int[nph_int];
    GenRandomEvtFromProbDist(org_norm_arr, ndet,
                             nph_int,
                             argval->GetRandSeed(),
                             obs_bin_arr,
                             obs_evt_arr);
    {
        long naxes[2];
        naxes[0] = ndetx;
        naxes[1] = ndety;
        char tag[kLineSize];
        sprintf(tag, "obs_%4.4d",
                argval->GetRandSeed());
        int bitpix = -64;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                               tag, 2,
                               bitpix,
                               naxes, obs_bin_arr);
    }

    double** obs_tr_arr = new double*[argval->GetNfold()];
    double** obs_vl_arr = new double*[argval->GetNfold()];
    for(int ifold = 0; ifold < argval->GetNfold(); ifold ++){
        obs_tr_arr[ifold] = new double[ndet];
        obs_vl_arr[ifold] = new double[ndet];
    }
    GenCVImageByPartition(obs_evt_arr, nph_int,
                          argval->GetNfold(),
                          argval->GetRandSeed(),
                          ndet,
                          obs_tr_arr,
                          obs_vl_arr);
    for(int ifold = 0; ifold < argval->GetNfold(); ifold ++){
        long naxes[2];
        naxes[0] = ndetx;
        naxes[1] = ndety;
        char tag[kLineSize];
        sprintf(tag, "obs_%4.4d_%2.2dfold%2.2d_tr",
                argval->GetRandSeed(),
                argval->GetNfold(), ifold);
        int bitpix = -64;
        MifFits::OutFitsImageD(argval->GetOutdir(),
                               argval->GetOutfileHead(),
                               tag, 2,
                               bitpix,
                               naxes, obs_tr_arr[ifold]);
        sprintf(tag, "obs_%4.4d_%2.2dfold%2.2d_vl",
                argval->GetRandSeed(),
                argval->GetNfold(), ifold);
        MifFits::OutFitsImageD(argval->GetOutdir(),
                               argval->GetOutfileHead(),
                               tag, 2,
                               bitpix,
                               naxes, obs_vl_arr[ifold]);
    }

    // check
    for(int ifold = 0; ifold < argval->GetNfold(); ifold ++){
        int nevt_out_tr = 0;
        int nevt_out_vl = 0;
        for(int idet = 0; idet < ndet; idet ++){
            nevt_out_tr += obs_tr_arr[ifold][idet];
            nevt_out_vl += obs_vl_arr[ifold][idet];
        }
        printf("obs_tr_arr[%d], obs_vl_arr[%d] = (%d, %d)\n",
               ifold, ifold, nevt_out_tr, nevt_out_vl);
    }
        
    for(int ifold = 0; ifold < argval->GetNfold(); ifold ++){
        delete [] obs_tr_arr[ifold];
        delete [] obs_vl_arr[ifold];
    }
    delete [] obs_tr_arr;
    delete [] obs_vl_arr;
    
    return status_prog;
}
