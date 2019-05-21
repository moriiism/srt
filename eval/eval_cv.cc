#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_eval_cv.h"
#include "sub.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValEvalCv* argval = new ArgValEvalCv;
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

    // cvlist
    string* lines_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetCvlist(),
                                    &lines_arr, &nline);

    double* heldist_arr = new double[nline];
    
    for(int iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);
        printf("nsplit = %d\n", nsplit);
        printf("%s\n", lines_arr[iline].c_str());

        string recfile = split_arr[0];
        string valfile = split_arr[1];

        printf("recfile = %s\n", recfile.c_str());
        printf("valfile = %s\n", valfile.c_str());
        
        // load recfile
        double* rec_arr = NULL;
        MifImgInfo* img_info_rec = new MifImgInfo;
        img_info_rec->InitSetImg(1, 1, nskyx, nskyy);
        int bitpix_rec = 0;
        MifFits::InFitsImageD(recfile, img_info_rec,
                              &bitpix_rec, &rec_arr);
        double nph_rec = MirMath::GetSum(nsky, rec_arr);
        printf("nph_rec = %e\n", nph_rec);

       
        double* rec_norm_arr = new double[nsky];
        if(fabs(nph_rec) > 1.0e-20){
            for(int isky = 0; isky < nsky; isky++){
                rec_norm_arr[isky] = rec_arr[isky] / nph_rec;
            }
        } else {
            // adhoc
            for(int isky = 0; isky < nsky; isky++){
                rec_norm_arr[isky] = 1./ nsky; 
            }
        }
        delete img_info_rec;
        delete [] rec_arr;
        
        // load valfile
        double* val_arr = NULL;
        MifImgInfo* img_info_val = new MifImgInfo;
        img_info_val->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_val = 0;
        MifFits::InFitsImageD(valfile, img_info_val,
                              &bitpix_val, &val_arr);
        double nph_val = MirMath::GetSum(ndet, val_arr);
        printf("nph_val = %e\n", nph_val);

        double* val_norm_arr = new double[ndet];
        for(int idet = 0; idet < ndet; idet++){
            val_norm_arr[idet] = val_arr[idet] / nph_val;
        }
        delete img_info_val;
        delete [] val_arr;

        // det_arr = R_mat %*% rec_norm_arr
        double* det_arr = new double[ndet];
        char* transa = new char [1];
        strcpy(transa, "N");
        dgemv_(transa, ndet, nsky, 1.0, const_cast<double*>(resp_mat_arr), ndet,
               const_cast<double*>(rec_norm_arr), 1,
               0.0, det_arr, 1);
        delete [] transa;

        heldist_arr[iline] = GetHellingerDist(det_arr, val_norm_arr, ndet);
        printf("heldist_arr[%d] = %e\n", iline, heldist_arr[iline]);

        delete [] split_arr;
        delete [] rec_norm_arr;
        delete [] val_norm_arr;
        delete [] det_arr;
    }
    delete [] resp_mat_arr;
    delete [] lines_arr;

    double heldist_ave = 0.0;
    for(int iline = 0; iline < nline; iline ++){
        printf("heldist_arr[%d] = %e\n", iline, heldist_arr[iline]);
        heldist_ave += heldist_arr[iline];
    }
    heldist_ave /= nline;

    double heldist_stddev = MirMath::GetStddev(nline, heldist_arr);
    printf("heldist_ave = %e\n", heldist_ave);
    printf("heldist_stddev = %e\n", heldist_stddev);
    
    char outfile[kLineSize];
    sprintf(outfile, "%s/%s_heldist.txt",
            argval->GetOutdir().c_str(), argval->GetOutfileHead().c_str());
    FILE* fp_out = fopen(outfile, "w");
    fprintf(fp_out, "%e  %e\n", heldist_ave, heldist_stddev);
    fclose(fp_out);
    
    return status_prog;
}
