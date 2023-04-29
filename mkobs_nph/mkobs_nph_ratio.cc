//
// mkobs_nph_ratio.cc
//
// make observation detector images with total photons (N_tot),
// where the ratio of the number of ouput photons (N_tot) to
// that of original photons (N_org) is set to input ratio.

#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkobs_nph_ratio.h"
#include "sim.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMkobsNphRatio* argval = new ArgValMkobsNphRatio;
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

    int nphoton = (int) (nph_int * argval->GetRatio());
    printf("output nphoton = %d\n", nphoton);
    double* obs_bin_arr = new double[ndet];
    int*    obs_evt_arr = new int[nph_int];
    SrtlibSim::GenRandomEvtFromProbDist(org_norm_arr, ndet,
                                        nphoton,
                                        argval->GetRandSeed(),
                                        obs_bin_arr,
                                        obs_evt_arr);
    long naxes[2];
    naxes[0] = ndetx;
    naxes[1] = ndety;
    char tag[kLineSize];
    sprintf(tag, "%4.4d", argval->GetRandSeed());
    bitpix = -64;
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           tag, 2,
                           bitpix,
                           naxes, obs_bin_arr);
    delete [] org_norm_arr;
    delete [] obs_bin_arr;
    delete [] obs_evt_arr;
    
    return status_prog;
}
