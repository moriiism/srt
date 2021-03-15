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
    int*    sky_evt_arr = new int[argval->GetNevt()];
    GenRandomEvtFromProbDist(in_norm_arr, nsky,
                             argval->GetNevt(), argval->GetRandSeedSky(),
                             sky_rand_arr,
                             sky_evt_arr);
    delete [] sky_evt_arr;
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
        sum_det += det_arr[idet];
    }
    printf("sum_det = %e\n", sum_det);
    double* det_norm_arr = new double [ndet];
    for(int idet = 0; idet < ndet; idet ++){
        det_norm_arr[idet] = det_arr[idet] / sum_det;
    }

    double* obs_bin_arr = new double[ndet];
    int*    obs_evt_arr = new int[argval->GetNevt()];
    GenRandomEvtFromProbDist(det_norm_arr, ndet,
                             argval->GetNevt(), argval->GetRandSeedDet(),
                             obs_bin_arr,
                             obs_evt_arr);
    {
        int naxis = 2;
        long* naxes = new long[naxis];
        naxes[0] = ndetx;
        naxes[1] = ndety;
        char tag[kLineSize];
        sprintf(tag, "sky_%4.4d_obs_%4.4d",
                argval->GetRandSeedSky(),
                argval->GetRandSeedDet());
        int bitpix = -64;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                               tag, 2,
                               bitpix,
                               naxes, obs_bin_arr);
        delete [] naxes;
    }

    for(int ipart = 0; ipart < argval->GetNpartition(); ipart ++){
        double** obs_tr_arr = new double*[argval->GetNfold()];
        double** obs_vl_arr = new double*[argval->GetNfold()];
        for(int ifold = 0; ifold < argval->GetNfold(); ifold ++){
            obs_tr_arr[ifold] = new double[ndet];
            obs_vl_arr[ifold] = new double[ndet];
        }
        int rand_seed_part = argval->GetRandSeedPartition() + ipart;
        GenCVImageByPartition(obs_evt_arr, argval->GetNevt(),
                              argval->GetNfold(), rand_seed_part,
                              ndet,
                              obs_tr_arr,
                              obs_vl_arr);
        for(int ifold = 0; ifold < argval->GetNfold(); ifold ++){
            int naxis = 2;
            long* naxes = new long[naxis];
            naxes[0] = ndetx;
            naxes[1] = ndety;
            char tag[kLineSize];
            sprintf(tag, "sky_%4.4d_obs_%4.4d_part%2.2d_%2.2dfold%2.2d_tr",
                    argval->GetRandSeedSky(), argval->GetRandSeedDet(), rand_seed_part,
                    argval->GetNfold(), ifold);
            int bitpix = -64;
            MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                                   tag, 2,
                                   bitpix,
                                   naxes, obs_tr_arr[ifold]);
            sprintf(tag, "sky_%4.4d_obs_%4.4d_part%2.2d_%2.2dfold%2.2d_vl",
                    argval->GetRandSeedSky(), argval->GetRandSeedDet(), rand_seed_part,
                    argval->GetNfold(), ifold);
            MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                                   tag, 2,
                                   bitpix,
                                   naxes, obs_vl_arr[ifold]);
            delete [] naxes;
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
    }
    
    return status_prog;
}
