#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkresp_det2.h"
#include "load_resp_det2.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    ArgValMkrespDet2* argval = new ArgValMkrespDet2;
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
    int nphoton_input = argval->GetNphotonInput();
    int nskyx = argval->GetNskyx();
    int nskyy = argval->GetNskyy();
    int ndetx = 0;
    int ndety = 0;
    double* resp_norm_mat_det1_arr = NULL;
    double* resp_norm_mat_det2_arr = NULL;
    double* eff_mat_arr = NULL;
    LoadRespDet2(argval->GetRespdir1(),
                 argval->GetRespdir2(),
                 nskyx, nskyy,
                 nphoton_input,
                 &resp_norm_mat_det1_arr,
                 &resp_norm_mat_det2_arr,
                 &eff_mat_arr,
                 &ndetx, &ndety);

    // output
    int nsky = nskyx * nskyy;
    int ndet = ndetx * ndety;

    int bitpix = -32;
    int naxis = 2;
    long* naxes = new long[naxis];
    naxes[0] = ndet;
    naxes[1] = nsky;
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "resp_norm_det1", naxis, bitpix,
                           naxes, resp_norm_mat_det1_arr);

    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "resp_norm_det2", naxis, bitpix,
                           naxes, resp_norm_mat_det2_arr);
    
    naxes[0] = nskyx;
    naxes[1] = nskyy;
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "eff", naxis, bitpix,
                           naxes, eff_mat_arr);
    return status_prog;
}
