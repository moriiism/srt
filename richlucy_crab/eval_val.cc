#include "mir_math.h"
#include "mib_blas.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_eval_val.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

double GetRootMeanSquareError(const double* const det_arr, 
                              const double* const val_arr,
                              int ndet);

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValEvalVal* argval = new ArgValEvalVal;
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

    // load recfile
    int nskyx = MifFits::GetAxisSize(argval->GetRecfile(), 0);
    int nskyy = MifFits::GetAxisSize(argval->GetRecfile(), 1);
    int nsky = nskyx * nskyy;
    double* rec_arr = NULL;
    MifImgInfo* img_info_rec = new MifImgInfo;
    img_info_rec->InitSetImg(1, 1, nskyx, nskyy);
    int bitpix_rec = 0;
    MifFits::InFitsImageD(argval->GetRecfile(), img_info_rec,
                          &bitpix_rec, &rec_arr);
    delete img_info_rec;
    double nph_rec = MirMath::GetSum(nsky, rec_arr);
    printf("nph_rec = %e\n", nph_rec);

    // load valfile
    int ndetx = MifFits::GetAxisSize(argval->GetValfile(), 0);
    int ndety = MifFits::GetAxisSize(argval->GetValfile(), 1);
    int ndet = ndetx * ndety;
    double* val_arr = NULL;
    MifImgInfo* img_info_val = new MifImgInfo;
    img_info_val->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_val = 0;
    MifFits::InFitsImageD(argval->GetValfile(), img_info_val,
                          &bitpix_val, &val_arr);
    double nph_val = MirMath::GetSum(ndet, val_arr);
    printf("nph_val = %e\n", nph_val);

    // load response 
    double* resp_mat_arr = NULL;
    int bitpix_resp = 0;
    MifImgInfo* img_info_resp = new MifImgInfo;
    img_info_resp->InitSetImg(1, 1, ndet, nsky);
    MifFits::InFitsImageD(argval->GetRespFile(), img_info_resp,
                          &bitpix_resp, &resp_mat_arr);

    // det_arr = R_mat %*% rec_arr
    double* det_arr = new double[ndet];
    char* transa = new char [1];
    strcpy(transa, "N");
    dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
           const_cast<double*>(rec_arr), 1,
           0.0, det_arr, 1);

    // scale det_arr
    int nfold = argval->GetNfold();
    dscal_(ndet, 1.0/(nfold - 1), det_arr, 1);

    double rmse = GetRootMeanSquareError(det_arr, val_arr, ndet);
    printf("rmse = %e\n", rmse);

    long naxes[2];
    naxes[0] = ndetx;
    naxes[1] = ndety;
    char tag[kLineSize];
    sprintf(tag, "val");
    int bitpix = -64;
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                           tag, 2,
                           bitpix,
                           naxes, val_arr);
    sprintf(tag, "det");
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                           tag, 2,
                           bitpix,
                           naxes, det_arr);
    char outfile[kLineSize];
    sprintf(outfile, "%s/%s_rmse.txt",
            argval->GetOutdir().c_str(), argval->GetOutfileHead().c_str());
    FILE* fp_out = fopen(outfile, "w");
    fprintf(fp_out, "%e\n", rmse);
    fclose(fp_out);
    
    return status_prog;
}

double GetRootMeanSquareError(const double* const det_arr,
                              const double* const val_arr,
                              int ndet)
{
    double sum = 0.0;
    for(int idet = 0; idet < ndet; idet ++){
        double diff = det_arr[idet] - val_arr[idet];
        sum += diff * diff;
    }
    double ave = sum / ndet;
    double ans = sqrt(ave);
    return ans;
}



