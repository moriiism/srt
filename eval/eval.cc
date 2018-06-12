#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_eval.h"
#include "sub.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValEval* argval = new ArgValEval;
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

    // load recfile
    double* rec_arr = NULL;
    MifImgInfo* img_info_rec = new MifImgInfo;
    img_info_rec->InitSetImg(1, 1, nskyx, nskyy);
    int bitpix_rec = 0;
    MifFits::InFitsImageD(argval->GetRecfile(), img_info_rec,
                          &bitpix_rec, &rec_arr);
    delete img_info_rec;
    double nph_rec = MirMath::GetSum(nsky, rec_arr);
    printf("nph_rec = %e\n", nph_rec);

    double* rec_norm_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky++){
        rec_norm_arr[isky] = rec_arr[isky] / nph_rec;
    }

    // load valfile
    double* val_arr = NULL;
    MifImgInfo* img_info_val = new MifImgInfo;
    img_info_val->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_val = 0;
    MifFits::InFitsImageD(argval->GetValfile(), img_info_val,
                          &bitpix_val, &val_arr);
    double nph_val = MirMath::GetSum(ndet, val_arr);
    printf("nph_val = %e\n", nph_val);

    double* val_norm_arr = new double[ndet];
    for(int idet = 0; idet < ndet; idet++){
        val_norm_arr[idet] = val_arr[idet] / nph_val;
    }
    
    // det_arr = R_mat %*% rec_norm_arr
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rec_norm_arr), 1,
           0.0, det_arr, 1);

    double heldist = GetHellingerDist(det_arr, val_norm_arr, ndet);
    printf("heldist = %e\n", heldist);

    printf("nph_rec = %e\n", nph_rec);
    printf("nph_val = %e\n", nph_val);

    int naxis = 2;
    long* naxes = new long[naxis];
    naxes[0] = nskyx;
    naxes[1] = nskyy;
    char tag[kLineSize];
    sprintf(tag, "val");

    int bitpix = -64;
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                           tag, 2,
                           bitpix,
                           naxes, val_norm_arr);
    
    sprintf(tag, "det");
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                           tag, 2,
                           bitpix,
                           naxes, det_arr);
    
    delete [] naxes;

    char outfile[kLineSize];
    sprintf(outfile, "%s/%s_heldist.txt",
            argval->GetOutdir().c_str(), argval->GetOutfileHead().c_str());
    FILE* fp_out = fopen(outfile, "w");
    fprintf(fp_out, "%e\n", heldist);
    fclose(fp_out);
    
    return status_prog;
}
