#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_pmmn_case2.h"
#include "sub_pmmn_case2.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValPmmnCase2* argval = new ArgValPmmnCase2;
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

    double* resp_mat_arr = NULL;
    LoadModel(argval->GetRespdir(), &resp_mat_arr);



    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, 60, 60);
    int bitpix = 0;
    double* data_arr = NULL;
    MifFits::InFitsImageD(argval->GetDatafile(), img_info,
                          &bitpix, &data_arr);
    int nph = MirMath::GetSum(3600, data_arr);
    printf("N photon = %d\n", nph);

    



    

//    printf("--- img_info ---\n");
//    MifImgInfo* img_info = new MifImgInfo;
//    img_info->Load(argval->GetSubcube());
//    img_info->PrintInfo();
//    printf("=== img_info ===\n");
//    printf("sparsity ratio = ncard / npixel = %e\n",
//           (double) argval->GetNcard() / (double)(img_info->GetNpixelImg() * img_info->GetNframe()));
//
//    int bitpix = 0;
//    double* X_mat = NULL;
//    int status = MifFits::InFitsCubeD(argval->GetInFile(),
//                                      img_info, &bitpix, &X_mat);
//    printf("MifFits::InFitsCubeF: status = %d\n", status);
//
//    double* S_mat = NULL;
//    double* L_mat = NULL;
//    double* G_mat = NULL;
//    double* L_U_mat = NULL;
//    double* L_S_vec = NULL;
//    double* L_V_mat = NULL;
//    GenGodec(X_mat,
//             img_info->GetNpixelImg(),
//             img_info->GetNframe(),
//             argval->GetNrank(), argval->GetNcard(),
//             argval->GetQpar(), argval->GetEpsilon(),
//             argval->GetNiter(), argval->GetSeed(),
//             &S_mat, &L_mat, &G_mat,
//             &L_U_mat, &L_S_vec, &L_V_mat);
//
//    double* LS_mat = new double [img_info->GetNpixelTotal()];
//    for(long ielem = 0; ielem < img_info->GetNpixelTotal(); ielem ++){
//        LS_mat[ielem] = L_mat[ielem] + S_mat[ielem];
//    }
//    MifFits::OutFitsCubeD(argval->GetOutdir(),
//                          argval->GetOutfileHead(),
//                          "S", img_info->GetNaxis(),
//                          bitpix,
//                          img_info->GetNaxesArr(), S_mat);
//    MifFits::OutFitsCubeD(argval->GetOutdir(),
//                          argval->GetOutfileHead(),
//                          "L", img_info->GetNaxis(),
//                          bitpix,
//                          img_info->GetNaxesArr(), L_mat);
//    MifFits::OutFitsCubeD(argval->GetOutdir(),
//                          argval->GetOutfileHead(),
//                          "G", img_info->GetNaxis(),
//                          bitpix,
//                          img_info->GetNaxesArr(), G_mat);
//    MifFits::OutFitsCubeD(argval->GetOutdir(),
//                          argval->GetOutfileHead(),
//                          "O", img_info->GetNaxis(),
//                          bitpix,
//                          img_info->GetNaxesArr(), X_mat);
//    MifFits::OutFitsCubeD(argval->GetOutdir(),
//                          argval->GetOutfileHead(),
//                          "LS", img_info->GetNaxis(),
//                          bitpix,
//                          img_info->GetNaxesArr(), LS_mat);
//
//    // svd image
//    long* naxes_L_U_mat = new long [3];
//    naxes_L_U_mat[0] = img_info->GetNaxesArrElm(0);
//    naxes_L_U_mat[1] = img_info->GetNaxesArrElm(1);
//    naxes_L_U_mat[2] = argval->GetNrank();
//    MifFits::OutFitsCubeD(argval->GetOutdir(),
//                          argval->GetOutfileHead(),
//                          "L_U", img_info->GetNaxis(),
//                          bitpix,
//                          naxes_L_U_mat, L_U_mat);
//    delete [] naxes_L_U_mat;
//
//    // svd time series
//    long* naxes_L_V_mat = new long [2];
//    naxes_L_V_mat[0] = img_info->GetNaxesArrElm(2);
//    naxes_L_V_mat[1] = argval->GetNrank();
//    MifFits::OutFitsImageD(argval->GetOutdir(),
//                           argval->GetOutfileHead(),
//                           "L_V", 2,
//                           bitpix,
//                           naxes_L_V_mat, L_V_mat);
//    delete [] naxes_L_V_mat;
//
//    // svd singular values
//    char outfile_L_S_vec[kLineSize];
//    sprintf(outfile_L_S_vec, "%s/%s_L_S.txt",
//            argval->GetOutdir().c_str(),
//            argval->GetOutfileHead().c_str());
//    FILE* fp = fopen(outfile_L_S_vec, "w");
//    for(int irank = 0; irank < argval->GetNrank(); irank ++){
//        fprintf(fp, "%-4d %e\n", irank + 1, L_S_vec[irank]);
//    }
//    fclose(fp);
//
//    
//    fclose(fp_log);
//    delete argval;
//    delete img_info;
//    delete [] X_mat;
//    delete [] S_mat;
//    delete [] L_mat;
//    delete [] G_mat;
//    delete [] L_U_mat;
//
//    double time_ed = MiTime::GetTimeSec();
//    printf("duration = %e sec.\n", time_ed - time_st);
    
    return status_prog;
}
