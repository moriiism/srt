#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_simobs.h"
#include "sub.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValSimobs* argval = new ArgValSimobs;
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

    // load infile
    double* in_arr = NULL;
    MifImgInfo* img_info_in = new MifImgInfo;
    img_info_in->InitSetImg(1, 1, nskyx, nskyy);
    int bitpix_in = 0;
    MifFits::InFitsImageD(argval->GetInfile(), img_info_in,
                          &bitpix_in, &in_arr);
    delete img_info_in;
    double nph_in = MirMath::GetSum(nsky, in_arr);
    printf("nph_in = %e\n", nph_in);

    // norm infile
    double* in_norm_arr = new double [nsky];
    for(int isky = 0; isky < nsky; isky ++){
        in_norm_arr[isky] = in_arr[isky] / nph_in;
    }
    delete [] in_arr;

    // sky image
    double* sky_rand_arr = new double[nsky];
    GenRandomEvtFromProbDist(in_norm_arr, nsky,
                             argval->GetNevt(), argval->GetRandSeedSky(),
                             sky_rand_arr);
    {
        int naxis = 2;
        long* naxes = new long[naxis];
        naxes[0] = nskyx;
        naxes[1] = nskyy;
        char tag[kLineSize];
        sprintf(tag, "sky_%4.4d", argval->GetRandSeedSky());
        int bitpix = -64;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                               tag, 2,
                               bitpix,
                               naxes, sky_rand_arr);
        delete [] naxes;
    }
    
    // det_arr = R_mat %*% sky_rand_arr
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(sky_rand_arr), 1,
           0.0, det_arr, 1);

    double sum_det = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        sum_det +=det_arr[idet];
    }
    printf("sum_det = %e\n", sum_det);
    double* det_norm_arr = new double [ndet];
    for(int idet = 0; idet < ndet; idet ++){
        det_norm_arr[idet] = det_arr[idet] / sum_det;
    }


    for(int iobs = 0; iobs < argval->GetNobs(); iobs++){

        printf("simobs iobs = %d\n", iobs);

        
        int rand_seed_det = argval->GetRandSeedDet() + iobs;
        double* obs_arr = new double[ndet];
        GenRandomEvtFromProbDist(det_norm_arr, ndet,
                                 argval->GetNevt(), rand_seed_det,
                                 obs_arr);
        int naxis = 2;
        long* naxes = new long[naxis];
        naxes[0] = ndetx;
        naxes[1] = ndety;
        char tag[kLineSize];
        sprintf(tag, "sky_%4.4d_obs_%4.4d",
                argval->GetRandSeedSky(), rand_seed_det);

        int bitpix = -64;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                               tag, 2,
                               bitpix,
                               naxes, obs_arr);
        delete [] naxes;
        delete [] obs_arr;
    }

    
    return status_prog;
}
